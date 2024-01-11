#task: to extract epithelial cells; perform subclustering in each cohort; identify cancerous epithelial population based on prev known tumor/cancer tx properties (HH, glycolysis, hipoxia, and etc.) and tumor-cell dominant clusters

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
                   ncpu=4,
                   test_schedule_xlsx='../results/a01b01_locate_gc_cellid_seus/gc_scr_celltype_lesion_cnt.xlsx', #mapping to shared celltype lineage
                   outd=file.path(proj_info$results,fpath_dt$fbase0))

message(str(args))

# create output directory
args$outd=create_dir_if_not_exist(args$outd) %>%
  R.utils::getAbsolutePath()

# >======
# load known cancerous gene sets
gset.lst = load_gc_cancer_gene_sets()
# <------

jmessage('select tumor (-like) samples ...')
# ======>
db_entry = "comp_orig.ident_map"
cnt_l=load_obj_from_db(query_name = db_entry)

db_entry='cellid_iseu_fpaths'
rds_fpaths = load_obj_from_db(query_name = db_entry)

# determine common higher level cell lineages
ctype.group_tests = args$test_schedule_xlsx %>%
  read.xlsx(sheet=1) %>%
  as.data.table() %>%
  split(by="cohort")

# define tumor vs. non-tumor cell group to compare with
sample.group_tests = args$test_schedule_xlsx %>%
  read.xlsx(sheet=2) %>%
  as.data.table() %>%
  split(by="cohort")

sample_groups = split(cnt_l$comp_group,by="cohort")

# extracting epithelials and perform a subclustering w/o the batch effect removal
# >====
iseus = imap(sample_groups[names(sample_groups)!='healthy20'],function(sample_group, cohort2) {
  
  message(sprintf('loading seu [%s] ...', cohort2))
  
  # subset cells to analyze
  iseu = extract_subset_seu_to_comp(cohort2,
                                    rds_fpaths,
                                    sample_groups,
                                    sample.group_tests, 
                                    ctype.group_tests,
                                    cohort1 = 'healthy20',
                                    is_epithelial=1)
  
  # perform a subclustering on RNA assay
  iseu = sc_cluster(iseu,
                    scale_to_regress=c('nCount_RNA'),
                    reuse_sct=0,
                    default_assay = 'RNA',
                    fc.resolution = 0.5)
  
  # to compute ams on the cancerous gene sets
  iseu = AddModuleScore_Assay(iseu, 
                              gset_lst = gset.lst, 
                              max_iter = 20,
                              cell_group_by = 'seurat_clusters',
                              new_assay_name = 'CancerProbes')
  iseu
})

db_entry = 'epitheial_subclusters'
rds_fn.3=file.path(args$outd,sprintf('%s.rds',db_entry))
update_db_resource(query_name=db_entry, rds_fpath=rds_fn.3,file.obj=iseus)
# <----


jmessage('in each cohort epithelial lineage, quantify known cancer tx property and cell contribution frequency at each cluster ...')
# >=====
cluster_prof.l = imap(iseus,function(iseu,cohort.2) {
  message(cohort.2)
  DefaultAssay(iseu) = 'CancerProbes'
  cluster_prof = list()
  
  # ---
  # to generate average known cancer tx property
  ams.l=AverageExpression(iseu)
  
  cluster_prof[['ams']]=ams.l$CancerProbes
  cluster_prof[['ams.sum']]=colSums(ams.l$CancerProbes)
  
  ams.m=ams.l$CancerProbes %>%
    t() %>%
    scale() %>% #scale across clusters in each feature
    t()
  
  cluster_prof[['ams_scaled']]=ams.m
  
  cluster_prof[['ams_scaled.sum']]=colSums(ams.m)
  
  # ---
  # to generate cell population frequency
  cmeta = get_cmeta_from_seu(iseu)
  
  #comp = c('non-tumor','tumor')
  cell_c.m = cmeta[,.N,by=c('seurat_clusters','comp')] %>%
    acast(comp ~ seurat_clusters, value.var = 'N')
  
  cell_c.m[is.na(cell_c.m)] = 0
  cluster_prof[['cell_cnt']] = cell_c.m
  
  cell_f.m=sweep(cell_c.m,1,rowSums(cell_c.m),'/')
  cluster_prof[['cell_freq']] = cell_f.m
  
  # compute a fold change w.r.t tumor/normal frequency
  cluster_prof[['cell_freq_l2fc']] = apply(cell_f.m, 2, function(c2) {
    fc_adj = mean(c2)*1e-5
    log2((c2[['expr']]+fc_adj)/(c2[['ctrl']]+fc_adj))
  })
  
  cluster_prof
})

db_entry = 'epitheial_subclusters_prof'
rds_fn.3=file.path(args$outd,sprintf('%s.rds',db_entry))
update_db_resource(query_name=db_entry, rds_fpath=rds_fn.3,file.obj=cluster_prof.l)
# <------