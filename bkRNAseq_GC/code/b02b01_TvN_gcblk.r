#task: load GC bulk mRNA-seq salmon matrices files and perform Wilcoxon rank-sum test between tumor vs. non-tumor (When non-tumor is not avail, combine those with GTEx by ComBat)

graphics.off()
closeAllConnections()
rm(list=ls())

proj_d=".."
setwd(file.path(proj_d,"code"))
source(file.path(proj_d,'code','lib_project.r'))

fpath_dt <- get_current_script_fpath()
proj_info = get_proj_info(debug2=0)

# main()
# >==========
args <- data.table(reuse=0,
									 exptag="cohort23",
									 ncpu=8,
									 debug=0,
									 outd=file.path(proj_info$results,fpath_dt$fbase0))

message(str(args))
create_dir_if_not_exist(args$outd)
# <--------


jmessage('load cohorts ...')
# >==========
# select cohorts where both tumor and normal samples are available
bgc_lookup = load_obj_from_db(query_name="lookup.cohort")

query_db = bgc_lookup$db %>%
	dplyr::filter(grepl('T',sample_key) & grepl('N',sample_key))

row_index = 1:nrow(query_db)
names(row_index) = query_db$cohort

plan2(ncpu=length(row_index))
deg.lst = future_imap(row_index, function(r2,cohort) {
	message(cohort)

	s2m.lst = load_obj_from_db(query_name=query_db$db_entry[r2][[1]])
	
	# sync sample meta to perform wilcox-based DEG
	s2m.l = sync_s2m_for_comparison(cohort=cohort,
																	s2m = s2m.lst, 
																	comp_meta="tumor_normal", 
																	celltype="PanCK+", #applies only spark23
																	debug2=0)
	
	# browser()
	stopifnot(identical(dim(s2m.l$fxs.m)[2],
											dim(s2m.l$smeta)[1]))
	
	if (cohort != "spark23") { # Nanostring GeoMx is normalized by vst which is already log scaled.
		s2m.l$fxs.m=s2m.l$fxs.m %>%
			log1p()
	}
	
	cnames=colnames(s2m.l$smeta)
	resp_column="tumor_normal"
	
	deg_obj = prep_wilcox_test(fxs.mtx = s2m.l$fxs.m,
														 smeta.dt = s2m.l$smeta,
														 comp_col = resp_column,
														 ctrl_var = 'Normal',
														 expr_var = 'Tumor',
														 test_method = "wilcox",
														 debug2 = 0,
														 comp_name = query_db$cohort[r2])
	
	deg_obj$comp_var.orig = ''
	deg_obj = run_wilcoxon_test(deg_obj=deg_obj,
															logExpr=1)
	deg_obj$deg.dt
	
})
plan2()

# register it into the local resource
db_entry = "blk_gc_deg_by_paired_tn"
rds_fn.1 = file.path(args$outd,sprintf("%s.rds",db_entry))
update_db_resource(query_name=db_entry,rds_fpath = rds_fn.1, file.obj = deg.lst)
# <-----------


jmessage('including tumor only and compare them w/ gtex (reuse the prev source code); gene annotation')
# > ===========
bgc_lookup = load_obj_from_db(query_name="lookup.cohort")

# look for tumor only cohorts
query_db = bgc_lookup$db %>%
	dplyr::filter((grepl('T',sample_key) & !grepl('N',sample_key))|cohort=="gtex")

row_index = 1:nrow(query_db)
names(row_index) = query_db$cohort
assay="gene"
slot="abundance"
feat_name="gene"

s2m.lst = load_obj_from_db(query_name = "gtex_s2m")
feats=rownames(s2m.lst[[assay]]$txi[[slot]])

plan2(ncpu=length(row_index))

#combine the tumor-only cohorts and generate a union gene expression matrix
fxs.m = future_imap(row_index, function(r2, cohort) {
	s2m.lst = load_obj_from_db(query_name = query_db$db_entry[r2])
	reshape2::melt(s2m.lst[[assay]]$txi[[slot]], varnames=c(feat_name,'sample'),value.name = 'expr')
}) %>% rbindlist() %>%
	dplyr::filter(gene %in% feats) %>%
	acast(gene ~ sample, value.var = 'expr')

smeta.dt = future_imap(row_index, function(r2, cohort) {
	s2m.lst = load_obj_from_db(query_name = query_db$db_entry[r2])
	smeta=s2m.lst$smeta[,list(sample)]
	smeta$cohort=cohort
	smeta
}) %>% rbindlist()
smeta.dt[,tumor_normal:="Tumor"]
smeta.dt[cohort=="gtex",tumor_normal:="Normal"]
plan2()
stopifnot(!any(duplicated(smeta.dt$sample)))

#batch effect removal
dim(fxs.m)
dim(smeta.dt)

# 
fxs.m[is.na(fxs.m)] = 0
fxs.m=fxs.m[(rowSums(fxs.m)>0.),]
fxs.m=fxs.m[,colSums(fxs.m)>0]

snames = intersect(colnames(fxs.m),smeta.dt$sample)
fxs.m = fxs.m[,snames]
smeta.dt = smeta.dt[match(snames,sample),]
stopifnot(identical(smeta.dt$sample, colnames(fxs.m)))

# use combat to align the tumor-only cohorts to GTEx
fxs_log1p_combat.m = sva::ComBat(dat=log1p(fxs.m),
										batch=smeta.dt$cohort,
										ref.batch = 'gtex')


resp_column="tumor_normal"
deg_obj = prep_wilcox_test(fxs.mtx = fxs_log1p_combat.m,
													 smeta.dt = smeta.dt,
													 comp_col = resp_column,
													 ctrl_var = 'Normal',
													 expr_var = 'Tumor',
													 test_method = "wilcox",
													 debug2 = 0,
													 comp_name = "tumor.only_vs_gtex")

deg_obj$comp_var.orig = ''
deg_obj = run_wilcoxon_test(deg_obj=deg_obj,
														logExpr=1)

db_entry = "blk_gc_deg_tumor.only_vs_gtex"

rds_fn.1 = file.path(args$outd,sprintf("%s.rds",db_entry))

update_db_resource(query_name=db_entry,
									 rds_fpath = rds_fn.1, 
									 file.obj = deg_obj)
# <-------------