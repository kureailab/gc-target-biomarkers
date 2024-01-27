stopifnot(exists('proj_d'))

suppressWarnings(suppressPackageStartupMessages({
  library(openxlsx)
  library(randomcoloR)
  library(data.table)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(VennDiagram)
  library(ggpubr)
  library(parallel)
  library(purrr)
  library(rlang)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(Biobase)
  library(preprocessCore)
  library(Seurat)
  library(sva)
  library(tidyverse)
  library(dplyr)
  
  source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
  source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))
  # source(file.path(Sys.getenv('R_UTIL'),'lib_rnaseq.r'))
  source(file.path(Sys.getenv('R_UTIL'),'lib_clustering.r'))
  source(file.path(Sys.getenv('R_UTIL'),'lib_seurat3.r'))
  source(file.path(Sys.getenv('R_UTIL'),'local_resource.r'))
}))

get_proj_info=function(task_d=NA,debug2=0) {
  if (debug2==1){browser()}
  sub_dirs = c('data','results','code')
  names(sub_dirs) = sub_dirs
  proj_dirs = imap(sub_dirs,function(sub_d,dmy) {
    if (debug2==1){browser()}
    if (is.na(task_d)) {
      target_d=file.path(proj_d,sub_d) %>%
        R.utils::getAbsolutePath()
      create_dir_if_not_exist(target_d)
    } else {
      target_d=file.path(proj_d,sub_d,task_d) %>%
        R.utils::getAbsolutePath()
      
      if (sub_d=='results') {
        create_dir_if_not_exist(target_d)
      }
    }

  })
  
  a=data.table(data=proj_dirs$data,
               results=proj_dirs$results,
               code=proj_dirs$code)
  
  a$file_receive_date="20220516"
  a[,resource_db_path:=R.utils::getAbsolutePath(file.path(results,'resource_db.tsv'))]
  a$project_alias="bkRNAseq_GC"
  a$project_name="bkRNAseq_GC"
  a$email="changjin@kurea.ai"
  a$pmid=""
  a
}

get_gc_bulk_meta_to_compare <- function() {
  list(stad=c('sample','race','gender','paper_Lauren.Class','paper_Molecular.Subtype','tissue_or_organ_of_origin','vital_status','paper_Total.Mutation.Rate','ajcc_pathologic_stage','paper_MSI.status','overall_survival','sgroup'), #sgroup (non.tumor)
       yonsei_precl_mcarray = c('sample','OS.months','Status','Stage','Adjuvant.Chemotherapy','Molecular Subgroup','Age','Gender','Lauren Type','Lymphovascular Invasion','Perineural Invasion','Tumor Location'), #all
       yonsei_precl_rnaseq = c('sample','source_type','tissue_source1','seq_source','TNMstage2','recur','Lauren'), #source_type (Normal)
       hslat20 = c('sample','Age','Gender','Lauren type','Death','overall_survival','Recurrence','Molecular Subtype','Total Mutation Rate','TP53 mutation','PIK3CA mutation','KRAS mutation','ARIDA mutation','RHOA mutation','Tumor purity','Tumor location','CDH1 somatic mutation','Clinical.Stage','H.pylori.Infection','lesion','age_group'), #lesion (normal)
       spgc = c('sample','stained','tumor_status','MSI','response'), #tumor_status (Non-Tumor)
       S.Kim = c('sample','msi_type','ebv_in_situ','num_snv','tcga','Mesenchymal_subtype','immsig','response'), #all
       Kor_St.Mary=c('sample','Age','Gender','Disease.Status','MSI.status','response','HER2.status'), #all
       Kor_Yonsei=c('sample','response','MSI.status'), # all
       K.Chida=c('sample','OS_event','OS_time','Age','sex','TMB','response'), # all
       C.Park=c('sample','Age','sgroup','acquisition','sex'), # sgroup (healthy,post_hp_eradication,hp_associated_gastritis)
       D.Mun=c('sample','histology','tissue','msi','ebv'), # all
       M.Kwon=c('sample','timepoint','response','tumor_normal','HER2')) # all
}

split_gc_5ici <- function(gexp_dt) {
  gc_cohort='immgc_target'
  smeta.dt = load_gc_bulks(cohort_name=gc_cohort,smeta_only = 1)$smeta
  
  tpm.lst = merge(x=gexp_dt,y=smeta.dt[,list(sample,author)],by="sample") %>%
    as.data.table() %>%
    dplyr::mutate(author=gsub('Kor:','Kor_',author)) %>%
    split(by="author")
  
  pmeta_lst = get_gc_cohort_pmeta()
  
  cohort_map = list(samsung="S.Kim",stmary="Kor_St.Mary",yonsei="Kor_Yonsei",chida="K.Chida",cpark="C.Park",dmun="D.Mun",jlee="M.Kwon",shared="shared")
  cohort_names = cohort_map[match(names(pmeta_lst),names(cohort_map))] %>%
    unlist()
  names(cohort_names)=NULL
  names(pmeta_lst) = cohort_names
  
  stopifnot(all(names(tpm.lst) %in% cohort_names))
  
  pmeta_lst$shared=NULL
  
  pmeta_lst = imap(pmeta_lst,function(pmeta.dt,cohort){
    pmeta.dt$is_normal=0
    if (cohort=="C.Park") {
      pmeta.dt[sgroup %in% c('healthy','post_hp_eradication','hp_associated_gastritis'),is_normal:=1]
    } 
    pmeta.dt
  })
  
  list(tpm=tpm.lst,
       smeta=pmeta_lst)
}

cvd_gene_lst <- function(gene_lst, assay_genes) {
  imap(gene_lst,function(genes,gs_name) {
    genes[genes %in% assay_genes]
  })
}

get_resp_col <- function(resp="both") {
  if (resp=="NR") {
    col_in_hex="#04BCC4"
  } else if (resp=="R") {
    col_in_hex="#F4786E"
  } else {
    col_in_hex="#000000"
  }
  col_in_hex
}

color_columns <- function(dt2,col2="rev_gset_name",debug2=0) {
  if (debug2==1){browser()}
  dt2[resp=="both",colored1:=sprintf("<span style='color: %s'>%s</span>",get_resp_col(resp="both"),get(col2))]
  dt2[resp=="NR",colored1:=sprintf("<span style='color: %s'>%s</span>",get_resp_col(resp="NR"),get(col2))]
  dt2[resp=="R",colored1:=sprintf("<span style='color: %s'>%s</span>",get_resp_col(resp="R"),get(col2))]
  dt2
}

fisher.z <- function (r1,r2,n1,n2,debug2=0) {
  if (debug2==1) {browser()}
  ((0.5*log((1+r1)/(1-r1)))-(0.5*log((1+r2)/(1-r2))))/((1/(n1-3))+(1/(n2-3)))^0.5

}

sync_s2m_for_comparison <- function(cohort, s2m, assay="gene", slot="abundance", comp_meta="tumor_normal", celltype=NULL, debug2=0) {
  if (debug2==1){browser()}
  
  message(sprintf('cohort[%s], comp_meta[%s]',cohort, comp_meta))
  
  if (cohort=="spark23") {
    if (comp_meta=="tumor_normal") {
      stopifnot(!invalid(celltype))
      s2m$smeta[,tumor_normal:=tumor_status]
      s2m$smeta[tumor_normal=="Non-Tumor",tumor_normal:="Normal"]
      s2m$smeta[tumor_normal=="Tumor",tumor_normal:="Tumor"]
      smeta = s2m$smeta[tumor_normal %in% c('Normal','Tumor'),]
      smeta = smeta[stained==celltype,]
    }
  } else if (cohort=="stad") {
    if (comp_meta=="tumor_normal") {
      s2m$smeta[,tumor_normal:=sgroup]
      s2m$smeta[tumor_normal=="non.tumor",tumor_normal:="Normal"]
      s2m$smeta[tumor_normal=="tumor",tumor_normal:="Tumor"]
      smeta = s2m$smeta[tumor_normal %in% c('Normal','Tumor'),]
      # fxs.m = s2m$txi[[assay]][,smeta$sample]
    }
  } else if (cohort=="yonsei_pdx_ogd") {
    if (comp_meta=="tumor_normal") {
      s2m$smeta[,tumor_normal:=source_type]
      s2m$smeta[tumor_normal=="Normal",tumor_normal:="Normal"]
      s2m$smeta[tumor_normal=="Tumor",tumor_normal:="Tumor"]
      smeta = s2m$smeta[tumor_normal %in% c('Normal','Tumor'),]
    }
  } else if (cohort=="jcheong22") {
    if (comp_meta=="tumor_normal") {
      s2m$smeta$tumor_normal="Tumor"
      smeta = s2m$smeta
    }
  } else if (cohort=="hslat20") { #hslat20=c('lesion','tumor','normal')
    if (comp_meta=="tumor_normal") {
      s2m$smeta[,tumor_normal:=lesion]
      s2m$smeta[tumor_normal=="normal",tumor_normal:="Normal"]
      s2m$smeta[tumor_normal=="tumor",tumor_normal:="Tumor"]
      smeta = s2m$smeta[tumor_normal %in% c('Normal','Tumor'),]
    }
  } else if (cohort=="cooi09") {#cooi09=c('lesion','Tumor',NA)# need to add lesion="Tumor"?
    if (comp_meta=="tumor_normal") {
      s2m$smeta$tumor_normal="Tumor"
      smeta = s2m$smeta
    }
  } else if (cohort=="gtex") { #gtex=c('tumor_or_normal',NA,'normal')
    if (comp_meta=="tumor_normal") {
      s2m$smeta$tumor_normal="Normal"
      smeta = s2m$smeta
    }
  } else if (cohort=="samsung") { #samsung=c('lesion','Tumor',NA)# need to add lesion="Tumor"?
    if (comp_meta=="tumor_normal") {
      s2m$smeta$tumor_normal="Tumor"
      smeta = s2m$smeta
    }
  } else if (cohort=="stmary") { #stmary=c('lesion','Tumor',NA)# sample ending w/ T: tumor, w/ N: normal
    if (comp_meta=="tumor_normal") {
      s2m$smeta[,tumor_normal:="Tumor"]
      s2m$smeta[grepl('N$',sample),tumor_normal:="Normal"]
      smeta=s2m$smeta
    }
  } else if (cohort=="yonsei") {
    if (comp_meta=="tumor_normal") {
      s2m$smeta[,tumor_normal:="Tumor"]
      smeta=s2m$smeta
    }
  } else if (cohort=="chida") {
    if (comp_meta=="tumor_normal") {
      s2m$smeta[,tumor_normal:="Tumor"]
      smeta=s2m$smeta
    }
  } else if (cohort=="dmun") {#dmun = c('tissue2','T','N')
    if (comp_meta=="tumor_normal") {
      s2m$smeta[,tumor_normal:=tissue2]
      s2m$smeta[tumor_normal=="T",tumor_normal:="Tumor"]
      s2m$smeta[tumor_normal=="N",tumor_normal:="Normal"]
      smeta=s2m$smeta
    }
  } else if (cohort=="cpark.surgery") { #cpark =  surgery: gastric_cancer vs. adj_severe_gastritis, biopsy: post_hp_eradication+hp_associated_gastritis+gastric_cancer vs. healthy;normal
    if (comp_meta=="tumor_normal") {
      s2m$smeta[,tumor_normal:="Tumor"]
      s2m$smeta[group=="Severe_gastritis",tumor_normal:="Normal"]
      smeta=s2m$smeta
    }
  } else if (cohort=="cpark.biopsy") {
    if (comp_meta=="tumor_normal") {
      s2m$smeta[,tumor_normal:="Tumor"]
      s2m$smeta[group=="Healthy",tumor_normal:="Normal"]
      smeta=s2m$smeta
    }
  } else if (cohort=="mkwon") { #mkwon = c('tumor_normal','tumor','normal')
    if (comp_meta=="tumor_normal") {
      s2m$smeta[tumor_normal=="tumor",tumor_normal:="Tumor"]
      s2m$smeta[tumor_normal=="normal",tumor_normal:="Normal"]
      smeta=s2m$smeta
    }
  } else if (cohort=="rkim") { #rkim = c('tumor_normal','tumor',NA)
    if (comp_meta=="tumor_normal") {
      s2m$smeta[tumor_normal=="tumor",tumor_normal:="Tumor"]
      s2m$smeta[tumor_normal=="normal",tumor_normal:="Normal"]
      smeta=s2m$smeta
    }
  }
  
  fxs.m = s2m[[assay]]$txi[[slot]][,smeta$sample]
  list(smeta=smeta,
       fxs.m=fxs.m)
}


combine_stmary_yonsei <- function(slot="abundance", lesion=c("Tumor")) {
  s2m.lst = list()
  s2m.lst[["stmary"]]=load_obj_from_db(query_name="stmary_s2m")
  s2m.lst[["yonsei"]]=load_obj_from_db(query_name="yonsei_s2m")
  
  s2m_sy.lst = list()
  s2m_sy.lst$txi[[slot]] = cbind(s2m.lst$stmary$gene$txi[[slot]],
                                 s2m.lst$yonsei$gene$txi[[slot]])
  
  smeta=s2m.lst$stmary$smeta[,list(sample, response, MSI.status, patient_id, Disease.Status)]
  smeta$cohort='stmary'
  smeta[,tumor_normal:="Tumor"]
  smeta[Disease.Status=="Normal",tumor_normal:="Normal"]
  smeta$Disease.Status=NULL
  
  smeta2=s2m.lst$yonsei$smeta[,list(sample, response, MSI.status, patient_id)]
  smeta2$cohort="yonsei"
  smeta2$tumor_normal = "Tumor"
  
  s2m_sy.lst$smeta = rbind(smeta,smeta2)
  
  s2m_sy.lst$smeta = s2m_sy.lst$smeta[tumor_normal %in% lesion,]
  s2m_sy.lst$txi[[slot]] = s2m_sy.lst$txi[[slot]][,s2m_sy.lst$smeta$sample]
  s2m_sy.lst
}


prep_ici_s2m_lst <- function() {
  bgc_lookup = load_obj_from_db(query_name="lookup.cohort")
  
  s2m_sy.lst = combine_stmary_yonsei()
  
  query_db = bgc_lookup$db %>%
    dplyr::filter(cohort %in% c('samsung','chida','mkwon'))
  
  db_entries=query_db$db_entry
  names(db_entries) = query_db$cohort
  
  db_entries = c(db_entries,preclin='stmary_yonsei')
  
  # prep only pre-treatment and GC samples
  # plan2(ncpu=length(db_entries))
  
  s2m_processed.lst = imap(db_entries, function(db_entry, cohort) {
    message(cohort)
    if (cohort=="preclin") {
      smeta = s2m_sy.lst$smeta %>%
        rename(msi_type='MSI.status')
      smeta[msi_type %in% c('MSI-H','MSI'),msi_type:="MSI"]
      smeta[msi_type=="N/A",msi_type:="Unk"]
      
      mod = model.matrix(~as.factor(response), data=smeta)
      
      fxs.m = sva::ComBat(dat = s2m_sy.lst$txi$abundance,
                          batch = smeta$cohort,
                          mod = mod)
      
    } else {
      s2m.l = load_obj_from_db(query_name=db_entry)
      fxs.m = s2m.l$gene$txi$abundance
      smeta = s2m.l$smeta
      
      if (cohort=="chida") {
        smeta = smeta %>%
          dplyr::filter(cancer_type %in% c('gastric')) %>% 
          dplyr::mutate(msi_type='MSI')
      } else if (cohort=="mkwon") {
        smeta = smeta %>%
          dplyr::filter(timepoint=="B" & tumor_normal=="tumor") %>% 
          dplyr::mutate(msi_type='MSI')
      }
    }
    fxs.m = fxs.m[,smeta$sample]
    
    list(smeta=smeta,
         fxs.m=fxs.m)
  })
  s2m_processed.lst
}


ici4_cohort_names <- function() {
  author.map=data.table(author=c('samsung','mkwon','preclin','chida'),
                        cohort_desc=c('S.Kim18','M.Kwon21','PreClinical','K.Chida22'))
}


annotate_deg_w_potential_target <- function(deg_freq.l,
                                            order_by1='n_cnt',
                                            order_by2='max.mlog10pv',
                                            coota_min_pct=0.05,
                                            coota_nz.expr=0.25,
                                            gtex_ntpm_co=25,
                                            debug2=0) {
  if (debug2==1) {browser()}
  
  stopifnot('gene' %in% colnames(deg_freq.l[[1]]))
  
  genes = rbindlist(deg_freq.l) %>%
    dplyr::pull(gene) %>%
    unique() %>%
    sort()
  
  gene_desc.dt = data.table(gene = genes)
  gene_desc.dt$gene_desc = get_gene_desc(gene_desc.dt$gene)
  
  coota7 = load_obj_from_db(query_name='coota7_expr_prof')
  
  expressed_coota = coota7 %>%
    dplyr::filter(pct.1 >= coota_min_pct &
                    avg_expr_nz>=coota_nz.expr) %>%
    dplyr::pull(gene) %>%
    unique() %>%
    sort()
  
  # ---
  # load stomach HPA
  hpa.ann = load_obj_from_db(query_name="hpa_stomach_ann")
  hpa.ann = load_stomach_hpa(hpa.ann)
  
  # ---
  # load Cosmic
  cosmic.ann = load_obj_from_db(query_name="cosmic_gi_ann")
  cosmic.ann = load_cosmic(cosmic.ann)
  
  # ---
  omnipath = load_obj_from_db(query_name="omnipath")
  omnipath.ann = load_omnipath(omnipath)
  
  # ---
  gtex_ntpm <- load_obj_from_db(query_name="gtex_ntpm")
  gtex_ntpm.ann <- load_gtex_ntpm(gtex_ntpm)
  
  # ---
  deg_freq.l = imap(deg_freq.l,function(deg_freq.dt, cohort2){
    deg_freq.dt = deg_freq.dt %>%
      dplyr::filter(!(gene %in% expressed_coota)) %>%
      merge(gene_desc.dt,by="gene",all.x=T) %>%
      merge(omnipath.ann,by="gene",all.x=T) %>%
      merge(hpa.ann, by="gene", all.x=T) %>%
      merge(cosmic.ann, by="gene", all.x=T) %>%
      merge(gtex_ntpm.ann,by="gene",all.x=T) %>%
      dplyr::arrange(-get(order_by1),-get(order_by2)) %>%
      as.data.table()
    
    # <-------------
    
    deg_freq.dt[,priority:=""]
    deg_freq.dt[!is.na(hpa_pathology.prog.stomach.cancer),priority:="stomach_cancer"]
    deg_freq.dt[!is.na(surface_protein),priority:=sprintf('%s;surface',priority)]
    deg_freq.dt[!is.na(druggable),priority:=sprintf('%s;druggable',priority)]
    deg_freq.dt[is.na(hpa_rna.brain.regional.specificity),hpa_rna.brain.regional.specificity:="Unk"]
    
    deg_freq.dt[brain<gtex_ntpm_co & stomach<gtex_ntpm_co & colon<gtex_ntpm_co & esophagus<gtex_ntpm_co, 
                priority:=sprintf('%s;gtex_brain_gi<%d',priority,gtex_ntpm_co)]
    deg_freq.dt[,priority:=gsub("^;","",priority)]

    deg_freq.dt
  })
  deg_freq.l
}
# <----

deg_test_blkRNAseq_TvGTEx <- function(assay='gene') {
  # look for tumor only cohorts
  bgc_lookup = load_obj_from_db(query_name="lookup.cohort")
  query_db = bgc_lookup$db %>%
    dplyr::filter((grepl('T',sample_key) & !grepl('N',sample_key))|cohort=="gtex")
  
  row_index = 1:nrow(query_db)
  names(row_index) = query_db$cohort
  
  slot="abundance"
  s2m.lst = load_obj_from_db(query_name = "gtex_s2m")
  feats=rownames(s2m.lst[[assay]]$txi[[slot]])
  
  plan2(ncpu=length(row_index))
  #combine the tumor-only cohorts and generate a union gene expression matrix
  fxs.m = future_imap(row_index, function(r2, cohort) {
    s2m.lst = load_obj_from_db(query_name = query_db$db_entry[r2])
    ret = NA
    if (assay %in% names(s2m.lst)) {
      ret = reshape2::melt(s2m.lst[[assay]]$txi[[slot]],
                           varnames=c('feat','sample'),
                           value.name = 'expr')
    }
    ret
  }) %>% `[`(!is.na(.)) %>% 
    rbindlist() %>%
    dplyr::filter(feat %in% feats) %>%
    acast(feat ~ sample, value.var = 'expr')
  
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
  deg_obj
}


deg_test_blkRNAseq_pairedTvN <- function(assay='gene') {
  bgc_lookup = load_obj_from_db(query_name="lookup.cohort")
  
  query_db = bgc_lookup$db %>%
    dplyr::filter(grepl('T',sample_key) & grepl('N',sample_key))
  
  row_index = 1:nrow(query_db)
  names(row_index) = query_db$cohort
  
  plan2(ncpu=length(row_index))
  deg.lst = future_imap(row_index, function(r2,cohort) { #future_
    message(cohort)
    
    s2m.lst = load_obj_from_db(query_name=query_db$db_entry[r2][[1]])
    
    ret = NA
    if (assay %in% names(s2m.lst)) {
      message(sprintf('found the assay [%s] from %s ...',assay,cohort))
      # sync sample meta to perform wilcox-based DEG
      s2m.l = sync_s2m_for_comparison(cohort=cohort,
                                      s2m = s2m.lst,
                                      assay = assay,
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
      ret = deg_obj$deg.dt
    }
    ret
  }) %>% `[`(!is.na(.))
  plan2()
  deg.lst
}

add_deg_annot_header <- function() {
  col_head_str = "gene: gene symbol,
feat:ensembl transcript id (isoform),
n_cnt:# of cohorts with P.Value<0.05,
gene_tx:gene transcript id (isoform),
observed_in: cohorts with P.Value<0.05,
m.P.Value:mean(P.Value),
sum.mlog10pv:sum(-log10(P.Value)),
m.mlog10pv:mean(-log10(P.Value)),
m.ctrlMean:mean(mRNA expression level in normal-like samples in log(TPM)),
m.exprMean:mean(mRNA expression level in tumor-like samples in log(TPM)),
max.mlog10pv:max(P.Value),
m.meanExpr_diff:mean(exprMean-ctrlMean),
gene_desc:gene description,
druggable:omnipath - drug information,
surface_protein:omnpath - surface/endogeneous proteins,
hpa_rna.specific.nTPM:HPA - nTPM specfic to certain organ,
hpa_rna.brain.regional.specificity:HPA - brain regional specificity,
hpa_rna.cancer.specific.fpkm:HPA - FPKM specific to cancer organ,
hpa_pathology.prog.stomach.cancer:HPA - pathology prognosis in stomach cancer,
cosmic_gene.ids:COSMIC,
cosmic_hallmarks:COSMIC,
aorta:GTEx - TPM,
blood:GTEx - TPM,
brain:GTEx - TPM,
breast:GTEx - TPM,
colon:GTEx - TPM,
coronary.artery:GTEx - TPM,
EBV.transformed.lymphocyte:GTEx - TPM,
esophagus:GTEx - TPM,
greater.omentum:GTEx - TPM,
heart:GTEx - TPM,
kidney:GTEx - TPM,
liver:GTEx - TPM,
lower.leg.skin:GTEx - TPM,
lung:GTEx - TPM,
minor.salivary.gland:GTEx - TPM,
pancreas:GTEx - TPM,
prostate.gland:GTEx - TPM,
skeletal.muscle.tissue:GTEx - TPM,
spinal.chord:GTEx - TPM,
spleen:GTEx - TPM,
stomach:GTEx - TPM,
subcutaneous.adipose.tissue:GTEx - TPM,
suprapubic.skin:GTEx - TPM,
testis:GTEx - TPM,
thyroid.gland:GTEx - TPM,
tibial.artery:GTEx - TPM,
tibial.nerve:GTEx - TPM,
transformed.skin.fibroblast:GTEx - TPM,
urinary.bladder:GTEx - TPM,
vagina:GTEx - TPM,
priority:gtex_brain_gi < n (GTEx.TPM@brain & colon & stomach<n); surface (omnipath); druggable (omnipath),
ranking_idx:sorted by max.mlog10pv w/ descreasing"
  
  heads=comma_string_to_list(col_head_str, sep = ',')
  head_desc.dt = lapply(heads,function(head){
    data.table(column=tstrsplit(head,':')[[1]],
               desc=tstrsplit(head,':')[[2]])
  }) %>% rbindlist()
  head_desc.dt
}