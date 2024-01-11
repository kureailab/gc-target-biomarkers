#keywords: annotated the DEG/marker table by gene description, cosmic, hpa, omnipath(druggable/surface/endogenous), coota(7 cell types expressed in critical organs), and GTEx to secure the safety and reduce off-targets
# ref: motivated by https://www.nature.com/articles/s41587-023-01684-0, 

graphics.off()
closeAllConnections()
rm(list=ls())

proj_d=".."
setwd(file.path(proj_d,"code"))
source(file.path(proj_d,'code','lib_project.r'))

fpath_dt <- get_current_script_fpath()
proj_info = get_proj_info(debug2=0)

args <- data.table(exptag=get_pipeline_step(fpath_dt$fpath),
                   reuse=1,
                   deg_pvco = 0.001,
                   coota_min_pct = 0.05,
                   coota_nz.expr = 0.25,
                   gtex_ntpm_co = 25,
                   ncpu=4,
                   gene_annot_d = '~/projects/refdb/gene_annotation_for_anitbody',
                   deg_rds.fpath = '../results/a02b01_deg_epithelial/epi_fndmkr_tm.l_p0.001.rds',
                   outd=file.path(proj_info$results,fpath_dt$fbase0))

message(str(args))

# create output directory
args$outd=create_dir_if_not_exist(args$outd) %>%
  R.utils::getAbsolutePath()

# register annotation resource if not done prevly
proj_db.dt <- load_db_table()

db_quries = c('cosmic_gi_ann','gtex_ntpm','hpa_stomach_ann','omnipath','coota7_expr_prof')
names(db_quries)= db_quries

imap(db_quries,function(db_query,dmy){
  
  rds_fpath.1 = file.path(args$gene_annot_d,sprintf('%s.rds',db_query))
  stopifnot(file.exists(rds_fpath.1))
  has_entry = proj_db.dt %>% 
    dplyr::filter(name == `db_query`)
  
  if (nrow(has_entry)==0) {
    update_db_resource(query_name=db_query,rds_fpath = rds_fpath.1)
  }
})
# <-----


jmessage('combine all the DEG table together and stratify by each cluster (celltype) ...')
# ======>
stopifnot(file.exists(args$deg_rds.fpath))

deg_freq.dt = readRDS(args$deg_rds.fpath) %>%
  imap(function(deg.dt,cohort) {
    deg.dt %>%
      dplyr::filter((avg_log2FC>0. & avg_expr_nz>0.8) & p_val<args$deg_pvco) %>%
      dplyr::group_by(gene) %>%
      dplyr::top_n(n=1,wt=mlog10pv) # %>%
  }) %>% rbindlist() %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(n_cnt=n(),
                   observed_in=unique(paste0(cohort,collapse=",")),
                   m.mlog10pv = mean(mlog10pv),
                   max.mlog10pv = max(mlog10pv),
                   sum.mlog10pv = sum(mlog10pv),
                   m.avg_expr=mean(avg_expr),
                   m.avg_expr_nz=mean(avg_expr_nz),
                   m.avg_log2FC=mean(avg_log2FC),
                   m.pct.1=mean(pct.1),
                   m.pct.2=mean(pct.2),
                   m.p_val_adj=mean(p_val_adj)) %>%
  dplyr::filter(n_cnt>=2) %>%
  dplyr::arrange(-n_cnt, -max.mlog10pv, -m.avg_expr)
  
# <------


# >======
# annoate gene description
deg_freq_annot.l = annotate_deg_w_potential_target(list(epi=deg_freq.dt),
                                                   gtex_ntpm_co = args$gtex_ntpm_co)

xslx_fpath.2 = file.path(args$outd,'epithelial_deg_tier0.xlsx')

#tier0: gtex_brain_GI<co & either surface or stomach_cancer annotated from HPA
#tier1: gtex_brain_GI<co & either surface and druggable (more strigent)
imap(deg_freq_annot.l,function(deg_freq_annot,cluster) {
  deg_freq_annot %>%
    dplyr::filter(grepl(sprintf('gtex_brain_gi<%d',args$gtex_ntpm_co),priority)) %>%
    dplyr::filter(grepl('(surface|stomach_cancer)',priority))
}) %>%
  dt_to_xlsx.sheet(xlsx_fpath=xslx_fpath.2)
#<------