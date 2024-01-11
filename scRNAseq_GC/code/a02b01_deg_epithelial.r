#task: to find marker genes in each cancerous epithelial clusters that are dominated by tumor-cells AND represented by a high accumulated known cancer tx property

graphics.off()
closeAllConnections()
rm(list=ls())

proj_d=".."
setwd(file.path(proj_d,"code"))
source(file.path(proj_d,'code','lib_project.r'))

fpath_dt <- get_current_script_fpath()
proj_info = get_proj_info(debug2=0)

args <- data.table(exptag=get_pipeline_step(fpath_dt$fpath),
                   deg_pvco = 0.001,
                   reuse=1,
                   ncpu=4,
                   dge_method = 'LR',
                   outd=file.path(proj_info$results,fpath_dt$fbase0))

message(str(args))

# create output directory
args$outd=create_dir_if_not_exist(args$outd) %>%
  R.utils::getAbsolutePath()

jmessage('loading epithelial seurat objects and cluster profile computed in the prev step ...')
# ======>
db_entry = "epitheial_subclusters"
iseu.l = load_obj_from_db(query_name = db_entry)

db_entry = "epitheial_subclusters_prof"
cluster_prof.l=load_obj_from_db(query_name = db_entry)
# <------


jmessage('perform fndmkr ...')
# >======
db_entry = 'epi_fndmkr.l'
rds_fn.1=file.path(args$outd, sprintf('%s.rds',db_entry))
if (args$reuse==1 & file.exists(rds_fn.1)) {
  fndmkr.l = readRDS(rds_fn.1)
} else {
  fndmkr.l = imap(iseu.l,function(iseu, cohort.2) {
    FindAllMarkers(iseu, test.use = 'LR') %>%
      as.data.table(keep.rownames=TRUE)
  })
  update_db_resource(query_name=db_entry, rds_fpath = rds_fn.1, file.obj = fndmkr.l)
}
# <------

jmessage('in each cohort, select tumor-dominated clusters, select the ones among those such that an accumulated cancer tx score > 0 among the subclusters')
# >======
db_entry = 'epi_fndmkr_tm.l'
rds_fn.2=file.path(args$outd, sprintf('%s.rds',db_entry))
if (args$reuse==1 & file.exists(rds_fn.2)) {
  fndmkr_tm.l = readRDS(rds_fn.2)
} else {
  fndmkr_tm.l = imap(fndmkr.l,function(fndmkr, cohort.2) {
    cluster_prof=cluster_prof.l[[cohort.2]]
    cprof.dt = data.table(seurat_clusters=names(cluster_prof$cell_freq_l2fc),
                          cfreq_l2fc = cluster_prof$cell_freq_l2fc,
                          ams_scaled.sum=cluster_prof$ams_scaled.sum) %>%
      dplyr::filter(cfreq_l2fc>0 & ams_scaled.sum>0)
    
    fndmkr %>%
      dplyr::filter(p_val < args$deg_pvco) %>%
      dplyr::filter(cluster %in% cprof.dt$seurat_clusters) %>%
      dplyr::mutate(cluster = droplevels(cluster))
  })
  update_db_resource(query_name=db_entry, rds_fpath = rds_fn.2, file.obj = fndmkr_tm.l)
}
# <------  


jmessage('taking average gene expression on each cell group ...')
# >=======
expr_prof.lst = imap(fndmkr_tm.l, function(deg.dt, cohort2) {
  gc()
  message(cohort2)
  iseu = iseu.l[[cohort2]]
  
  iseu = subset(iseu, seurat_clusters %in% unique(deg.dt$cluster))
  
  deg.dts = split(deg.dt, by = "cluster")
  
  DefaultAssay(iseu) = 'RNA'
  
  expr_summary.dt = split_seurats(seu=iseu, split_by = 'seurat_clusters') %>%
    imap(function(by_cluster,scluster){
      message(sprintf('averaging gexpr/pct on [%s] ...',scluster))
      mean_expr.dt = gexpr_prof_by_cgroup(seu = by_cluster,
                                          features = deg.dts[[scluster]]$gene,
                                          debug2=0,
                                          group.by='comp') %>%
        dplyr::mutate(cluster=scluster)
      mean_expr.dt
    }) %>% rbindlist()
  expr_summary.dt
})
# <-------

jmessage('appending the average gene expression of the tumor group into the DEG/marker table ...')
# >========
fndmkr_tm.l = add_avgexpr_to_sdeg(fndmkr_tm.l,
                                  expr_prof.lst, 
                                  deg_pvco=args$deg_pvco, 
                                  debug2=0)

db_entry = sprintf('epi_fndmkr_tm.l_p%g',args$deg_pvco)
rds_fn.4 = file.path(args$outd,sprintf('%s.rds',db_entry))
saveRDS(fndmkr_tm.l, file=rds_fn.4)
# <--------