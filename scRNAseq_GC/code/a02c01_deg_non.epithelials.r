#task: to perform DEG between tumor vs. non-tumor on each dominant cell-type lineage (except epithelial)

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
                   deg_pvco=0.001,
                   ncpu=4,
                   dge_method = 'LR',
                   default_assay = 'RNA',
                   test_schedule_xlsx='../results/a01b01_locate_gc_cellid_seus/gc_scr_celltype_lesion_cnt.xlsx',
                   outd=file.path(proj_info$results,fpath_dt$fbase0))

message(str(args))

# create output directory
args$outd=create_dir_if_not_exist(args$outd) %>%
  R.utils::getAbsolutePath()

# db_entry='cellid_iseu_fpaths'

jmessage('select tumor (-like) samples and perform the copykat ...')
# ======>
db_entry = "comp_orig.ident_map"
cnt_l=load_obj_from_db(query_name = db_entry)

db_entry='cellid_iseu_fpaths'
rds_fpaths = load_obj_from_db(query_name = db_entry)

ctype.group_tests = args$test_schedule_xlsx %>%
  read.xlsx(sheet=1) %>%
  as.data.table() %>%
  split(by="cohort")

sample.group_tests = args$test_schedule_xlsx %>%
  read.xlsx(sheet=2) %>%
  as.data.table() %>%
  split(by="cohort")

sample_groups = split(cnt_l$comp_group,by="cohort")

# ---
# deg.l =  imap(sample_groups[names(sample_groups)=='jkim22'], 
db_entry = 'non.epithelial_deg'
rds_fn.2 = file.path(args$outd,sprintf('%s.rds',db_entry))

if (args$reuse==1 & file.exists(rds_fn.3)) {
  deg.l = readRDS(rds_fn.2)
} else {
  deg.l = imap(sample_groups[names(sample_groups)!='healthy20'],function(sample_group, cohort2) {
    # browser()
    message(sprintf('loading seu [%s] ...', cohort2))
    
    iseu = extract_subset_seu_to_comp(cohort2,
                                      rds_fpaths,
                                      sample_groups,
                                      sample.group_tests, 
                                      ctype.group_tests,
                                      cohort1 = 'healthy20')
    # ---
    deg.dt = diff_from_integration(iseu,
                                   cluster_col="celltype_test",
                                   min_cell_cnts=3,
                                   method = args$dge_method,
                                   ncpu=4,
                                   debug=0)
    deg.dt
  })
  
  update_db_resource(query_name=db_entry, rds_fpath = rds_fn.2, file.obj=deg.l) 
}
#<------

gc()

db_entry = 'non.epithelial_deg.gexp'
rds_fn.3 = file.path(args$outd,sprintf('%s.rds',db_entry))
if (args$reuse==1 & file.exists(rds_fn.3)) {
  expr_prof.lst = readRDS(rds_fn.3)
} else {
  expr_prof.lst = imap(deg.l, function(deg.dt,cohort2) {
    message(cohort2)
    iseu = extract_subset_seu_to_comp(cohort2,
                                      rds_fpaths, 
                                      sample_groups,
                                      sample.group_tests, 
                                      ctype.group_tests,
                                      cohort1 = 'healthy20')
    gc()
    
    deg.dt = deg.dt[p_val<args$deg_pvco,]
    deg.dts = split(deg.dt, by = "cluster_id")
    
    DefaultAssay(iseu) = 'RNA'
    
    expr_summary.dt = split_seurats(seu=iseu,split_by = 'celltype_test') %>%
      imap(function(by_ctype,ctype){
        mean_expr.dt = gexpr_prof_by_cgroup(seu = by_ctype,
                                            features = deg.dts[[ctype]]$gene,
                                            debug2=0,
                                            group.by='comp') %>%
          dplyr::mutate(cluster=ctype)
        mean_expr.dt
      }) %>% rbindlist()
    expr_summary.dt
  })
  
  update_db_resource(query_name=db_entry,
                     rds_fpath = rds_fn.3, 
                     file.obj=expr_prof.lst)
}

gc()

# >========
deg.l = add_avgexpr_to_sdeg(deg.l, expr_prof.lst, deg_pvco=args$deg_pvco, debug2=0)
db_entry = sprintf('non.epithelial_deg_p%g',args$deg_pvco)
rds_fn.4 = file.path(args$outd,sprintf('%s.rds',db_entry))
saveRDS(deg.l, file=rds_fn.4)
# <--------