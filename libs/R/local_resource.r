library(data.table)
library(purrr)
library(gtools)
library(clusterProfiler)
source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))
source(file.path(Sys.getenv('R_UTIL'),'lib_apps.R'))


get_msigdbr_hs <- function(sp2="hs",ucache_dir="~/temp",reuse=1) {
	
	rds_fn = file.path(ucache_dir,sprintf("msigdbr_%s.rds",sp2))
	if (reuse==1 & file.exists(rds_fn)) {
		message(sprintf("reusing %s",rds_fn))
		msigdbr_hs_m_t2g=readRDS(rds_fn)
	} else {
		if (!file.exists(ucache_dir)){
			create_dir_if_not_exist(ucache_dir)
		}
		if (sp2=="hs") {
			msigdbr_hs_m_t2g = msigdbr(species = "Homo sapiens")
		} else {
			msigdbr_hs_m_t2g = msigdbr(species = "Mus musculus")
		}
		saveRDS(msigdbr_hs_m_t2g,file=rds_fn)
	}
	msigdbr_hs_m_t2g
}
