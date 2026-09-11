#!/usr/bin/env python3
"""Audit the public candidate and optionally refresh its code inventory."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
CODE_SUFFIXES = {".R": "R", ".r": "R", ".py": "Python", ".sh": "Shell"}
MAX_PUBLIC_BYTES = 5 * 1024 * 1024
BANNED_PARTS = {".venv", ".Rlib", "__pycache__", ".pytest_cache", "node_modules"}
TEXT_SKIP_SUFFIXES = {".png", ".jpg", ".jpeg", ".pdf", ".rds", ".rdata", ".h5", ".h5ad"}
PATTERNS = {
    "macOS user path": re.compile(r"/Users/[A-Za-z0-9._-]+/"),
    "Linux user path": re.compile(r"/home/[A-Za-z0-9._-]+/"),
    "Windows user path": re.compile(r"[A-Za-z]:\\\\Users\\\\[^\\\\]+\\\\", re.I),
    "private key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    "AWS access key": re.compile(r"\bAKIA[0-9A-Z]{16}\b"),
    "literal secret assignment": re.compile(
        r"(?i)(?:api[_-]?key|access[_-]?token|password|passwd)\s*[:=]\s*['\"][^'\"\n]{8,}['\"]"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def public_files() -> list[Path]:
    command = ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"]
    result = subprocess.run(command, cwd=ROOT, check=True, capture_output=True)
    names = [name for name in result.stdout.decode().split("\0") if name]
    return sorted((ROOT / name for name in names if (ROOT / name).is_file()), key=lambda p: str(p.relative_to(ROOT)))


def tier(relative: Path) -> str:
    parts = relative.parts
    if parts and parts[0] == "revision":
        return "revision"
    if parts and parts[0] == "scripts":
        return "historical"
    return "infrastructure"


def scan(files: list[Path]) -> list[str]:
    errors: list[str] = []
    for path in files:
        relative = path.relative_to(ROOT)
        if any(part in BANNED_PARTS for part in relative.parts):
            errors.append(f"banned runtime/cache path: {relative}")
        size = path.stat().st_size
        if size > MAX_PUBLIC_BYTES:
            errors.append(f"public file exceeds {MAX_PUBLIC_BYTES} bytes: {relative} ({size})")
        if path.suffix.lower() in TEXT_SKIP_SUFFIXES:
            continue
        payload = path.read_bytes()
        if b"\0" in payload:
            continue
        text = payload.decode("utf-8", errors="replace")
        for label, pattern in PATTERNS.items():
            if pattern.search(text):
                errors.append(f"{label}: {relative}")
    return errors


def write_inventory(files: list[Path]) -> int:
    rows = []
    for path in files:
        language = CODE_SUFFIXES.get(path.suffix)
        if not language:
            continue
        relative = path.relative_to(ROOT)
        rows.append(
            {
                "path": relative.as_posix(),
                "language": language,
                "tier": tier(relative),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    inventory = ROOT / "docs" / "CODE_INVENTORY.csv"
    inventory.parent.mkdir(parents=True, exist_ok=True)
    with inventory.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "language", "tier", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)
    checksum_file = ROOT / "checksums" / "code_sha256.txt"
    checksum_file.parent.mkdir(parents=True, exist_ok=True)
    checksum_file.write_text(
        "".join(f"{row['sha256']}  {row['path']}\n" for row in rows), encoding="utf-8"
    )
    return len(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-only", action="store_true", help="Do not refresh inventory/checksum files")
    args = parser.parse_args()
    files = public_files()
    errors = scan(files)
    code_count = sum(path.suffix in CODE_SUFFIXES for path in files)
    if not args.check_only and not errors:
        code_count = write_inventory(files)
    summary = {
        "status": "PASS" if not errors else "FAIL",
        "public_candidate_files": len(files),
        "code_files": code_count,
        "inventory_refreshed": not args.check_only and not errors,
        "errors": errors,
    }
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
