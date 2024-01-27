#task: collapse DEG tables (Tumor vs. Normal) and annotate them by FDA (potentially) druggable and surface proteins, filter out any gene expressed in healthy critical organs (>5% & non-zero mean expr>0.25)

graphics.off()
closeAllConnections()
rm(list=ls())

proj_d="../.."
setwd(file.path(proj_d,"code","de_transcripts"))
source(file.path(proj_d,'code','lib_project.r'))

fpath_dt <- get_current_script_fpath()
proj_info = get_proj_info(task_d=fpath_dt$parentd,debug2=0)

# main()
# >==========
args <- data.table(reuse=0,
									 exptag="cohort23",
									 topk_nd = 20,
									 expr_co=0.5,
									 min_cohort_n=2,
									 coota_min_pct = 0.02,
									 coota_nz.expr = 0.25,
									 gtex_ntpm_co=25,
									 ncpu=8,
									 debug=0,
									 outd=file.path(proj_info$results,fpath_dt$fbase0))

message(str(args))
create_dir_if_not_exist(args$outd)
# <--------


jmessage('load TvN deg table and combine all together ...')
# >==========
# select cohorts where both tumor and normal samples are available
load_db_table()
bgc_lookup = load_obj_from_db(query_name="lookup.cohort")

deg.lst = list()
deg.lst[['paired']] = load_obj_from_db(query_name="blk_gc_deg_by_paired_tn.transcript") %>%
	imap(function(deg.dt,cohort) {
		deg.dt %>%
			dplyr::mutate(cohort=cohort)
	}) %>% rbindlist()

deg.lst[['tumor_only']] = load_obj_from_db(query_name="blk_gc_deg_tumor.only_vs_gtex.transcript")$deg.dt %>%
	dplyr::mutate(cohort="tumor_vs_gtex")

deg.dt=load_obj_from_db(query_name='chida_s2m')$transcript$tx2gene %>%
  dplyr::mutate(gene_tx=tstrsplit(Name,'\\|')[[5]]) %>%
  rename(gene='GENESYMBOL') %>%
  dplyr::select(TXNAME,gene_tx,gene,EffectiveLength) %>%
  merge(rbindlist(deg.lst),by.x='TXNAME',by.y='feat',all.y=TRUE) %>%
  rename(feat='TXNAME')

rm(deg.lst)
gc()
# <-----------


jmessage('select frequent significant DEG across the cohorts ...')
# >=============
# deg.dt[lfc_group=="expr",hist(exprMean)]

deg_freq.dt = deg.dt %>%
	dplyr::mutate(mlog10pv = -log10(P.Value)) %>%
	dplyr::filter((lfc_group=="expr" & exprMean>args$expr_co) & P.Value<0.05) %>%
	dplyr::group_by(feat) %>%
	dplyr::summarize(n_cnt=n(),
	                 gene = unique(gene),
	                 gene_tx = unique(gene_tx),
									 observed_in=unique(paste0(cohort,collapse=",")),
									 m.P.Value=mean(P.Value),
									 sum.mlog10pv = sum(mlog10pv),
									 m.mlog10pv = mean(mlog10pv),
									 m.ctrlMean=mean(ctrlMean),
									 m.exprMean=mean(exprMean),
									 max.mlog10pv = max(mlog10pv),
									 m.meanExpr_diff=mean(meanExpr_diff)) %>%
	dplyr::filter(n_cnt>=args$min_cohort_n) %>%
	dplyr::arrange(-n_cnt, -sum.mlog10pv)
# <------------


# >=============

deg_freq_annot.l = annotate_deg_w_potential_target(list(gc_bulk_target=deg_freq.dt),
                                                   gtex_ntpm_co = args$gtex_ntpm_co)

deg_freqs = list()

deg_freq_annot = deg_freq_annot.l$gc_bulk_target %>%
  dplyr::filter(grepl(sprintf('gtex_brain_gi<%d',args$gtex_ntpm_co),priority)) %>%
  dplyr::filter(grepl('surface',priority))

deg_freqs[['druggable']] = deg_freq_annot %>%
  dplyr::filter(grepl('druggable',priority)) %>%
  dplyr::arrange(-max.mlog10pv) %>%
  dplyr::mutate(ranking_idx = row_number())

deg_freqs[['no.drug_info']] = deg_freq_annot %>%
  dplyr::filter(!grepl('druggable',priority)) %>%
  dplyr::arrange(-max.mlog10pv) %>%
  dplyr::slice(1:args$topk_nd) %>%
  dplyr::mutate(ranking_idx = row_number())

deg_annot_head.dt=add_deg_annot_header()
deg_freqs[['col_desc']] = deg_annot_head.dt
xslx_fpath.2 = file.path(args$outd,'TvN_gcblk_annot_tier1.xlsx')
dt_to_xlsx.sheet(deg_freqs,xlsx_fpath=xslx_fpath.2)
# <------------