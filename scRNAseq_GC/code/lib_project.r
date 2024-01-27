stopifnot(exists('proj_d'))

library(data.table)
library(openxlsx)
library(randomcoloR)
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
# library(VennDiagram)
library(ggpubr)
library(parallel)
library(purrr)
# library(RUVSeq)
# library(DESeq2)
library(patchwork)
library(ComplexHeatmap)
library(furrr)
# library(survival)
# library(survminer)
# library(maftools)
library(R.utils)
library(spam)
library(fields)
library(reticulate)
library(openxlsx)
library(gtools)

source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))
# source(file.path(Sys.getenv('R_UTIL'),'lib_rnaseq.r'))
source(file.path(Sys.getenv('R_UTIL'),'lib_seurat3.r'))
# source(file.path(Sys.getenv('R_UTIL'),'lib_clustering.r'))
source(file.path(Sys.getenv('R_UTIL'),'local_resource.r'))
source(file.path(Sys.getenv('R_UTIL'),'lib_workflow.r'))

project_name = "scRNAseq_GC"

get_proj_info=function(task_d=NA,debug2=0) {
	if (debug2==1){browser()}
	if (!is.na(task_d)) {
		if (task_d=="code") {
			task_d=NA
		}
	}
	
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
	
	a$file_receive_date="20231220"
	a[,resource_db_path:=R.utils::getAbsolutePath(file.path(results,'resource_db.tsv'))]
	a$owner="kureai"
	a$project_name=project_name
	a$email="changjin@kure.ai"
	a$pmid="pmid_kureai"
	a
}

extract_subset_seu_to_comp <- function(cohort2, 
                                       rds_fpaths, 
                                       sample_groups, 
                                       sample.group_tests, 
                                       ctype.group_tests,
                                       cohort1 = 'healthy20',
                                       is_epithelial=0) {
  
  sample_group=sample_groups[[cohort2]]
  
  sample.group_test.2 = sample.group_tests[[cohort2]] %>%
    merge(sample_group, by="comp_group")
  
  iseu = rds_fpaths[[cohort2]] %>%
    readRDS() %>%
    subset(orig.ident %in% sample.group_test.2$orig.ident)
  
  ctype.group_test.2 = ctype.group_tests[[cohort2]]
  
  cmeta = get_cmeta_from_seu(iseu) %>%
    merge(ctype.group_test.2[,list(celltype,celltype_2nd)],
          by='celltype') %>%
    merge(sample.group_test.2[,list(orig.ident,test_group)],
          by="orig.ident") %>%
    dplyr::filter(!is.na(celltype_2nd))
  
  iseu = subset(iseu, cells = cmeta$cbc)
  default_assay='RNA'
  
  # check if this cohort has ctrl group. Otherwise, append a ctrl group from ref
  test_n = sample.group_test.2 %>%
    dplyr::count(test_group)
  
  stopifnot('expr' %in% test_n$test_group)
  
  if (!('ctrl' %in% test_n$test_group)) {
    sample_group.c = sample_groups[[cohort1]]
    sample.group_test.c = sample.group_tests[[cohort1]] %>%
      merge(sample_group.c, by="comp_group")
    
    cseu = rds_fpaths[[cohort1]] %>%
      readRDS() %>%
      subset(orig.ident %in% sample.group_test.c$orig.ident)
    
    ctype.group_test.c = ctype.group_tests[[cohort1]]
    
    cmeta.c = get_cmeta_from_seu(cseu) %>%
      merge(y=.,x=ctype.group_test.c[,list(celltype,celltype_2nd)],
            by.y='cell_type',by.x='celltype') %>%
      merge(sample.group_test.c[,list(orig.ident,test_group)],
            by="orig.ident") %>%
      dplyr::filter(!is.na(celltype_2nd))
    
    cseu = subset(cseu, cells = cmeta.c$cbc)
    
    iseu = merge_seurats(list(iseu,cseu),add_sample_idx=F)
    
    cmeta = rbind(cmeta[,list(cbc,celltype_2nd,test_group)],
                  cmeta.c[,list(cbc,celltype_2nd,test_group)])
    
  }
  
  if (is_epithelial==1){
    cmeta = cmeta[celltype_2nd=='epithelial',]
  } else {
    cmeta = cmeta[celltype_2nd!='epithelial',]
  }
  
  iseu = subset(iseu, cells = cmeta$cbc) %>%
    NormalizeData(assay='RNA')
  
  iseu[['celltype_test']]=cmeta[match(colnames(iseu),cbc),celltype_2nd]
  iseu[['comp']]=cmeta[match(colnames(iseu),cbc),test_group]
  iseu
}

load_gc_cancer_gene_sets <- function(){
  cancer_hmk_w_stat3=comma_string_to_list('
HALLMARK_HYPOXIA,
KEGG_PATHWAYS_IN_CANCER,
WP_HIF1A_AND_PPARG_REGULATION_OF_GLYCOLYSIS,
HALLMARK_HEDGEHOG_SIGNALING')
  
  msigdb_hs.cancer_hallmark = get_msigdbr_hs(sp2="hs",ucache_dir="~/temp") %>%
    dplyr::filter(gs_name %in% cancer_hmk_w_stat3) %>%
    as.data.table()
  
  gset.lst = msigdb_hs.cancer_hallmark %>%
    split(by="gs_name") %>%
    imap(.,function(by_gs,gs) {
      by_gs$gene_symbol
    })
  
  gset.lst[['AUNG_GASTRIC_CANCER']] = 'ALDH7A1,alpha4GnT,APIN,AQR,ATE1,ATPIF1,BIRC5,BRD4,C4orf9,CBFA2T3,CCT3,CYP2W1,DEFA5,DEFA6,DKK4,ETS2,FLJ10036,FLJ36666,FXYD3,GITA,GPP34R,GW112,HORMAD1,HOXA10,IFRD1,IL16,JUN,KIF4A,LMO6,MAPK13,MGC20806,MIA,MLL4,MMP10,MYBL2,NEK9,NIPSNAP3B,PEGASUS,PPARBP,PRKAG1,REG4,RPL8,SEC31L2,SFRS9,SH3BGRL2,STAT2,SULT1C1,TAPBP,TD-60,THBS3,TMLHE,TPT1,TRAG3,TYRO3' %>%
    comma_string_to_list()
  
  gset.lst[['pzhang19']] = c('CEACAM5','CEACAM6','OLFM4','EPHB2','SOX9')
  gset.lst[['jkim21']] = c('EPCAM','CDH17','CDH4')
  gset.lst
}

add_avgexpr_to_sdeg <- function(deg.l, avg.l, deg_pvco=0.001, debug2=0) {
  
  stopifnot(identical(names(deg.l),names(avg.l)))
  stopifnot('cluster' %in% colnames(deg.l[[1]]))
  stopifnot('cluster' %in% colnames(avg.l[[1]]))
  
  deg_mg.l = imap(deg.l, function(deg.dt, sname) {
    if (debug2==1){browser()}
    deg_by_cg = deg.dt %>%
      dplyr::filter(p_val<deg_pvco) %>%
      split(by='cluster')
    
    avg_by_cg =  avg.l[[sname]] %>%
      split(by='cluster')
    
    imap(deg_by_cg,function(deg.dt,cgroup) {
      if (debug2==1){browser()}
      
      #add expr_prof
      p_val.nz = deg.dt[p_val>0.,min(p_val)]
      
      deg.dt[p_val==0.,p_val:=p_val.nz*0.99]
      deg.dt[,mlog10pv:= -log10(p_val)]
      deg.dt[, gexpr := NULL]
      
      #append an average expression of the expr group only
      avg_gexp.dt.2 = avg_by_cg[[cgroup]] %>%
        dplyr::filter(cell_group=='expr') %>%
        dplyr::filter(gene %in% deg.dt$gene) %>%
        dplyr::select(-cell_group, -pct.1, -cluster) %>%
        merge(deg.dt, by='gene', all.y=T) %>%
        dplyr::mutate(cohort = sname)
      
      gc()
      
      avg_gexp.dt.2
    })%>% rbindlist()
  })
  deg_mg.l
}


# >====
annotate_deg_w_potential_target <- function(deg_freq.l,
                                            coota_min_pct=0.05,
                                            coota_nz.expr=0.25,
                                            gtex_ntpm_co=25) {
  
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
      dplyr::arrange(-n_cnt, -sum.mlog10pv) %>%
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
    deg_freq.dt[order(-n_cnt, -max.mlog10pv)]
    deg_freq.dt
  })
  deg_freq.l
}
# <----
