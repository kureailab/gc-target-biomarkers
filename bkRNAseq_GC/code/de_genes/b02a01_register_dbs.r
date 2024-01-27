#task: register salmon gene/tx level assay/slots to project resource table

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
									 exptag="register_dbs",
									 rds_pats='../../../results/de_genes/b02a02_update_cohorts/*.rds',
									 gene_annot_d = '~/projects/refdb/gene_annotation_for_anitbody',
									 ncpu=1,
									 debug=0)

message(str(args))
# <-----------


jmessage('registering GC mRNA-bulk s2m files ...')
# >==========
load_db_table()

rds_fpaths = Sys.glob(args$rds_pats)
names(rds_fpaths) = sapply(rds_fpaths,function(rds_fpath) {
  decompose_fname(rds_fpath)$fbase0
})

imap(rds_fpaths, function(rds_fpath,db_entry) {
  message(db_entry)
  rds_fpath = R.utils::getAbsolutePath(rds_fpath)
  update_db_resource(query_name=db_entry, rds_fpath=rds_fpath)
  NA
})
# <-----------

jmessage('register gene annotation and filtering info ...')
# >=========
proj_db.dt <- load_db_table()

db_quries = c('cosmic_gi_ann',
              'gtex_ntpm',
              'hpa_stomach_ann',
              'omnipath',
              'coota7_expr_prof')

names(db_quries)= db_quries

imap(db_quries,function(db_query,dmy){
  
  rds_fpath.1 = file.path(args$gene_annot_d,sprintf('%s.rds',db_query)) %>%
    R.utils::getAbsolutePath()
  stopifnot(file.exists(rds_fpath.1))
  has_entry = proj_db.dt %>% 
    dplyr::filter(name == `db_query`)
  
  if (nrow(has_entry)==0) {
    update_db_resource(query_name=db_query,rds_fpath = rds_fpath.1)
  }
})
# <---------
