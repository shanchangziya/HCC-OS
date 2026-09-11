args <- commandArgs(trailingOnly=TRUE)
src <- normalizePath(args[1], mustWork=TRUE)
out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
e <- new.env(parent=globalenv())
load(file.path(src,"Resnet.Rdata"), envir=e)
stopifnot(all(c("resnet2","rp") %in% ls(e)))
stopifnot(ncol(e$resnet2)==2048L)
write.csv(data.frame(patient_id=rownames(e$resnet2), e$resnet2,
                     check.names=FALSE), file.path(out,"original_patient_features.csv"), row.names=FALSE)
write.csv(e$rp, file.path(out,"original_molecular_clinical.csv"),row.names=FALSE)
p <- new.env(parent=globalenv())
load(file.path(src,"pathdat.Rdata"), envir=p)
stopifnot("dat_all" %in% ls(p))
write.csv(data.frame(original_row_id=rownames(p$dat_all),p$dat_all,
                     check.names=FALSE),file.path(out,"original_pathology_split.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"export_R_session.txt"))
