#task: to register the preprocessed seurat objec file from each GC scRNA-seq cohort into the local placeholder

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
                   outd=file.path(proj_info$results,fpath_dt$fbase0))

message(str(args))

# create output directory
args$outd=create_dir_if_not_exist(args$outd) %>%
  R.utils::getAbsolutePath()

db_entry='cellid_iseu_fpaths'
rds_fpath.2 = file.path(args$outd,sprintf('%s.rds',db_entry))
if (args$reuse==1 & file.exists(rds_fpath.2)) {
  seu_fpaths = readRDS(rds_fpath.2)
} else {
  seu_fpaths = list()
  #locate each integrated seurat objec file path
  
  seu_fpaths[['pzhang19']]='~/projects/2019/04_P.Zhang_early.GC_scRNA/results/sctomach_repr/dt_cellid.rds'
  
  seu_fpaths[['asathe20']]='~/projects/2020/2020.sathe_gc_scell/results/a02a02/seu_ro_cid_by_author.rds'
  
  seu_fpaths[['mkwon21']]='~/projects/2021/M.Kwon21_GC.msiH.ICI/results/d01a01_cellid_stomach/mkwon_seu_sct.rds'
  
  seu_fpaths[['hjeong21']]='~/projects/2021/HY.Jeong21_GC.CCL2.inflammatory/results/b01a01_cca_integration/seu_cca.ams.rds'
  
  seu_fpaths[['jkim22']]='~/projects/2022/20220127.GC_intratumor.jhKim/results/b02b02_cellid/harmony_seu_cellid.rds'
  
  seu_fpaths[['vkumar22']]='~/projects/2022/2022.vkumar_GC_scr/results/a01c01_cca_integration/v.kumar22_rpca_cellid_ams.rds'
  
  seu_fpaths[['ksun22']]='~/projects/2022/08_kSun_GC_scr/results/a02a01_iseu_from_author/seu_cellid.rds'
  
  seu_fpaths[['hjiang22']]='~/projects/2022/H.Jiang22_GC.hetero.mets/results/b01a02_cca_cellid/seu_rpca.ams.rds'
  
  seu_fpaths[['rwang23']]='~/projects/2023/09_rWang_GSE234129_GC_scr/results/a02a01_cellid_cca/iseu_cca.rds'
  
  seu_fpaths[['healthy20']] = '~/projects.scratch/2023/08_A.Gottschlich_AML/results/b02a02_build_coota/HCL/han20.AdultStomach.rds'
  
  do.call(file.exists, seu_fpaths) %>%
    stopifnot(all())
  
  update_db_resource(query_name=db_entry,
                     file.obj = seu_fpaths,
                     rds_fpath = rds_fpath.2)
}

# >==========

jmessage('assigning tumor(like) for an expriment group and normal(like) for a control group for a downstream comparive analysis ...')
db_entry = "comp_orig.ident_map"
rds_fpath.2 =  file.path(args$outd,sprintf('%s.rds',db_entry))
if (args$reuse==1 & file.exists(rds_fpath.2)) {
  cnt_l = readRDS(rds_fpath.2)
} else {
  cohort = "pzhang19"
  seu=readRDS(seu_fpaths[[cohort]])
  
  cmeta = get_cmeta_from_seu(seu) %>%
    dplyr::mutate(comp_group=gsub('\\d+','',orig.ident))
  
  celltype_n.list = list()
  comp_group.list = list()
  
  celltype_n.list[[cohort]] = cmeta %>%
    dplyr::count(celltype) %>%
    dplyr::mutate(cohort=cohort)
  
  comp_group.list[[cohort]] = cmeta %>% #IMx & EGC vs. CAG & NAG
    dplyr::count(orig.ident, comp_group)%>%
    dplyr::mutate(cohort=cohort)
  
  # ----
  cohort = 'asathe20'
  seu=readRDS(seu_fpaths[[cohort]])
  cmeta = get_cmeta_from_seu(seu) %>%
    dplyr::mutate(comp_group=tstrsplit(orig.ident,'_')[[2]])
  
  celltype_n.list[[cohort]] = cmeta %>%
    dplyr::count(celltype) %>%
    dplyr::mutate(cohort=cohort)
  
  comp_group.list[[cohort]] = cmeta %>%
    dplyr::count(orig.ident, comp_group) %>%
    dplyr::mutate(cohort=cohort)
  
  # ---
  cohort = 'hjeong21'
  seu = readRDS(seu_fpaths[[cohort]])
  
  cmeta = get_cmeta_from_seu(seu)
  
  cmeta[,comp_group:=tstrsplit(orig.ident,'\\.')[[2]]]
  
  celltype_n.list[[cohort]]=cmeta %>%
    dplyr::count(celltype) %>%
    dplyr::mutate(cohort=cohort)
  
  comp_group.list[[cohort]] =cmeta %>%
    dplyr::count(orig.ident, comp_group) %>%
    dplyr::mutate(cohort=cohort)
  
  # ----
  cohort = 'jkim22'
  seu = readRDS(seu_fpaths[[cohort]])
  
  cmeta = get_cmeta_from_seu(seu) %>%
    dplyr::mutate(comp_group=lesion)
  
  celltype_n.list[[cohort]]=cmeta %>%
    dplyr::count(celltype) %>%
    dplyr::mutate(cohort=cohort)
  
  comp_group.list[[cohort]] =cmeta %>%
    dplyr::count(orig.ident, comp_group) %>%
    dplyr::mutate(cohort=cohort)
  
  # ----
  cohort = 'vkumar22'
  seu = readRDS(seu_fpaths[[cohort]])
  
  cmeta = get_cmeta_from_seu(seu)
  cmeta[,comp_group:=lesion]
  
  celltype_n.list[[cohort]]=cmeta %>%
    dplyr::count(celltype) %>%
    dplyr::mutate(cohort=cohort)
  
  comp_group.list[[cohort]] = cmeta %>%
    dplyr::count(orig.ident, comp_group) %>%
    dplyr::mutate(cohort=cohort)
  
  # ----
  cohort = 'ksun22'
  seu = readRDS(seu_fpaths[[cohort]])
  
  cmeta = get_cmeta_from_seu(seu)
  
  cmeta[,comp_group:=tissue]
  
  celltype_n.list[[cohort]]=cmeta %>%
    dplyr::count(celltype) %>%
    dplyr::mutate(cohort=cohort)
  
  comp_group.list[[cohort]] = cmeta %>%
    dplyr::count(orig.ident, comp_group) %>%
    dplyr::mutate(cohort=cohort)
  
  # ----
  cohort = 'hjiang22'
  seu = readRDS(seu_fpaths[[cohort]])
  
  cmeta = get_cmeta_from_seu(seu)
  # PT1, PT2, NT1, PT3
  cmeta[,comp_group:="Else"]
  cmeta[orig.ident %in% c('PT1','PT2','PT3'),comp_group:="Primary_Tumor"]
  cmeta[orig.ident %in% c('NT1'),comp_group:="Non_Tumor"]
  
  celltype_n.list[[cohort]]=cmeta %>%
    dplyr::count(celltype) %>%
    dplyr::mutate(cohort=cohort)
  
  comp_group.list[[cohort]] = cmeta %>%
    dplyr::count(orig.ident, comp_group) %>%
    dplyr::mutate(cohort=cohort)
  
  # ---
  # load healthy donor
  cohort = 'healthy20'
  seu = readRDS(seu_fpaths[[cohort]])
  cmeta = get_cmeta_from_seu(seu)
  setnames(cmeta,'cell_type','celltype')
  cmeta[,comp_group:='normal']
  
  celltype_n.list[[cohort]]=cmeta %>%
    dplyr::count(celltype) %>%
    dplyr::mutate(cohort=cohort)
  
  comp_group.list[[cohort]] = cmeta %>%
    dplyr::count(orig.ident, comp_group) %>%
    dplyr::mutate(cohort=cohort)
  
  # ---
  cohort = 'mkwon21' # this is all tumor; use annot_2nd_half; use healthy donor samples for ctrl
  seu = readRDS(seu_fpaths[[cohort]])
  if (F) {
    unique(seu$annot_3rd)
    cmeta = get_cmeta_from_seu(seu)
    cmeta[,celltype:=annot_2nd_half]
    cmeta[annot_3rd=='SMC',celltype:='SMC']
    cmeta[celltype=="Fibroblast_SMC" & annot_3rd!='SMC',celltype:='fibroblasts']
    
    seu[['celltype']] = seu$celltype
    saveRDS(seu, file = seu_fpaths[[cohort]])
  }
  
  cmeta = get_cmeta_from_seu(seu)
  cmeta = cmeta[TimePoint=="Pre",]
  cmeta[,comp_group:='Tumor']
  
  celltype_n.list[[cohort]]=cmeta %>%
    dplyr::count(celltype) %>%
    dplyr::mutate(cohort=cohort)
  
  comp_group.list[[cohort]] = cmeta %>%
    dplyr::count(orig.ident, comp_group) %>%
    dplyr::mutate(cohort=cohort)
  
  # ----
  cnt_l = list()
  
  cnt_l[['celltype']]=celltype_n.list %>%
    rbindlist()
  
  cnt_l[['comp_group']]=comp_group.list %>%
    rbindlist()
  
  update_db_resource(query_name=db_entry,
                     file.obj = cnt_l,
                     rds_fpath = rds_fpath.2)
}