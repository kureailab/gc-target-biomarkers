#task: load GC bulk mRNA-seq salmon matrices files and perform Wilcoxon rank-sum test between tumor vs. non-tumor (When non-tumor is not avail, combine those with GTEx by ComBat)

graphics.off()
closeAllConnections()
rm(list=ls())

proj_d="../.."
setwd(file.path(proj_d,"code","de_genes"))
source(file.path(proj_d,'code','lib_project.r'))

fpath_dt <- get_current_script_fpath()
proj_info = get_proj_info(task_d=fpath_dt$parentd,debug2=0)

# main()
# >==========
args <- data.table(reuse=0,
                   exptag="cohort23",
                   assay='gene',
                   ncpu=8,
                   debug=0,
                   outd=file.path(proj_info$results,fpath_dt$fbase0))

message(str(args))
create_dir_if_not_exist(args$outd)
assay.2 = args$assay
# <--------


jmessage('load cohorts ...')
# >==========
# select cohorts where both tumor and normal samples are available
deg.lst = deg_test_blkRNAseq_pairedTvN(assay=assay.2)

# register it into the local resource
db_entry = "blk_gc_deg_by_paired_tn"
rds_fn.1 = file.path(args$outd,sprintf("%s.rds",db_entry))
update_db_resource(query_name=db_entry,
                   rds_fpath = rds_fn.1, 
                   file.obj = deg.lst)
# <-----------


jmessage('including tumor only and compare them w/ gtex (reuse the prev source code); gene annotation')
# > ===========
deg_obj = deg_test_blkRNAseq_TvGTEx(assay=assay.2)

db_entry = "blk_gc_deg_tumor.only_vs_gtex"

rds_fn.1 = file.path(args$outd,sprintf("%s.rds",db_entry))

update_db_resource(query_name=db_entry,
                   rds_fpath = rds_fn.1, 
                   file.obj = deg_obj)
# <-------------
