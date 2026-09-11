"""Build the validation atlas from the author's public count matrix.

All published cells are retained for visualization; QC flags are exported for
sensitivity analyses. Authors' cluster labels define malignancy independently
of ROS/OSARS. Expression scores must use uncorrected RNA, never Harmony values.
"""
from pathlib import Path
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ.setdefault(key,'4')
import gzip, hashlib, json, importlib.metadata
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
import scanpy as sc
import harmonypy

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/processed'
OUT.mkdir(parents=True,exist_ok=True)
SEED=20260905
RAW=Path(os.environ.get('HCC_OS_GSE149614_COUNTS', ROOT/'data/raw/GSE149614_HCC.scRNAseq.S71915.count.fresh.txt.gz'))
METADATA=Path(os.environ.get('HCC_OS_GSE149614_METADATA_CSV', OUT/'GSE149614_metadata_annotated.csv'))
sc.settings.n_jobs=4
np.random.seed(SEED)

def load_counts():
    if not RAW.is_file():
        raise FileNotFoundError(f'GSE149614 count matrix not found: {RAW}. Set HCC_OS_GSE149614_COUNTS.')
    if not METADATA.is_file():
        raise FileNotFoundError(f'Annotated metadata not found: {METADATA}. Set HCC_OS_GSE149614_METADATA_CSV.')
    cached=OUT/'GSE149614_counts.h5ad'
    if cached.exists():
        return ad.read_h5ad(cached)
    blocks=[];indices=[];genes=[];indptr=[0]
    with gzip.open(RAW,'rt') as f:
        cells=next(f).rstrip('\n\r').split('\t')
        for k,line in enumerate(f):
            gene,values=line.rstrip('\n\r').split('\t',1)
            x=np.fromstring(values,sep='\t',dtype=np.int32)
            if len(x)!=len(cells): raise ValueError(f'Malformed count row {k}: {gene}')
            if np.any(x<0): raise ValueError('Negative count')
            idx=np.flatnonzero(x).astype(np.int32)
            indices.append(idx);blocks.append(x[idx]);indptr.append(indptr[-1]+len(idx));genes.append(gene)
            if k%3000==0: print('Read genes',k,flush=True)
    mat=sparse.csr_matrix((np.concatenate(blocks),np.concatenate(indices),np.array(indptr,dtype=np.int64)),shape=(len(genes),len(cells))).T.tocsr()
    meta=pd.read_csv(METADATA).set_index('Cell')
    assert len(meta)==len(cells) and set(meta.index)==set(cells)
    assert len(set(genes))==len(genes), 'Duplicated genes require explicit resolution'
    a=ad.AnnData(mat,obs=meta.loc[cells].copy(),var=pd.DataFrame(index=pd.Index(genes,name='gene')))
    a.var['mt']=a.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(a,qc_vars=['mt'],percent_top=None,log1p=False,inplace=True)
    a.obs['technical_qc_pass']=(a.obs.n_genes_by_counts>=200)&(a.obs.n_genes_by_counts<=6000)&(a.obs.pct_counts_mt<20)
    a.obs.to_csv(OUT/'GSE149614_metadata_with_QC.csv')
    a.write_h5ad(cached,compression='gzip')
    audit={'source_id':RAW.name,'raw_sha256':hashlib.sha256(RAW.read_bytes()).hexdigest(),'genes':a.n_vars,'cells':a.n_obs,'seed':SEED,'count_sum':int(a.X.sum()),'nnz':int(a.X.nnz),'versions':{p:importlib.metadata.version(p) for p in ['scanpy','anndata','numpy','scipy','harmonypy']}}
    (ROOT/'audit/atlas_input_manifest.json').write_text(json.dumps(audit,indent=2))
    return a

def embedding(a,label):
    if (OUT/f'{label}_embedding.csv').exists(): return
    x=a.copy()
    sc.pp.normalize_total(x,target_sum=1e4)
    sc.pp.log1p(x)
    sc.pp.highly_variable_genes(x,n_top_genes=2500,flavor='seurat',batch_key='sample')
    pd.DataFrame({'gene':x.var_names,'highly_variable':x.var.highly_variable}).to_csv(OUT/f'{label}_hvg.csv',index=False)
    x=x[:,x.var.highly_variable].copy()
    sc.pp.scale(x,max_value=10)
    sc.tl.pca(x,n_comps=30,svd_solver='arpack',random_state=SEED)
    hm=harmonypy.run_harmony(x.obsm['X_pca'],x.obs,'sample',random_state=SEED,max_iter_harmony=20,verbose=True)
    z=np.asarray(hm.Z_corr)
    if z.shape!=x.obsm['X_pca'].shape: z=z.T
    assert z.shape==x.obsm['X_pca'].shape
    x.obsm['X_harmony']=z
    sc.pp.neighbors(x,n_neighbors=20,n_pcs=30,use_rep='X_harmony',random_state=SEED)
    sc.tl.umap(x,min_dist=.35,random_state=SEED)
    table=x.obs.copy()
    table['UMAP1']=x.obsm['X_umap'][:,0];table['UMAP2']=x.obsm['X_umap'][:,1]
    table.to_csv(OUT/f'{label}_embedding.csv')
    np.savez_compressed(OUT/f'{label}_pca_harmony.npz',cells=np.asarray(x.obs_names,dtype=str),pca=x.obsm['X_pca'],harmony=z)
    print('Saved embedding',label,len(table),flush=True)

if __name__=='__main__':
    a=load_counts()
    print('Validated count matrix',a.shape,flush=True)
    embedding(a,'GSE149614_all_cells')
    embedding(a[a.obs.primary_analysis.astype(bool)].copy(),'GSE149614_primary_malignant')
    print('ATLAS_COMPLETE',flush=True)
