#!/usr/bin/env Rscript

#Authors: Chanjing Hong
#Date Created: 1/2/2018

library(data.table)
library(stringr)
library(tools)
library(future)
library(openxlsx)
library(gtools)

jmessage <- function(msg,task="message_name") {
	message(sprintf("[%s] %s|%s ...",Sys.time(), msg, task))
}

print_table <- function(u_title,u_fontSize,ms3a) {
	
	title = textGrob(u_title,gp=gpar(fontsize=12))
	padding = unit(5,'mm')
	table = gtable_add_rows(tableGrob(ms3a,rows=NULL,theme=ttheme_default(base_size=u_fontSize)),
													heights = grobHeight(title) + padding,
													pos=0)
	
	table = gtable_add_grob(table,
													title,
													1,1,1,ncol(table),
													clip = 'off')
	grid.newpage()
	grid.draw(table)
}

file_filter_rows <- function(fn,colidx1,key) {
	fn_head = sprintf('%s.head',fn)
	system(sprintf('head -n1 %s > %s',fn,fn_head))
	fn_tmp = sprintf('%s.tmp',fn)
	system(sprintf('awk \'$%d == \"%s\" {print $0}\' %s > %s',colidx1,key,fn,fn_tmp))
	system(sprintf('cat %s >> %s',fn_tmp,fn_head))
	unlink(fn_tmp)
	fn_head
}

get_max_val_in_list <- function(my_list,absolute=F){
  M <- length(my_list)
  maxVals <- rep(-1e7,M)
  for (i in 1:M){
    maxVals[i] <- max(my_list[[i]],na.rm = T)
  }
  return(max(maxVals))
}


get_min_val_in_list <- function(my_list,absolute=F){
  M <- length(my_list)
  minVals <- rep(1e7,M)
  for (i in 1:M){
    minVals[i] <- min(my_list[[i]],na.rm = T)
  }
  return(min(minVals))
}

comma1k <- function(my_integers){
  
  my_integer_dc <- data.class(my_integers[1])
  
  if (my_integer_dc=='numeric'){
    char_integers_with_comma <- format(my_integers, big.mark=",", scientific=FALSE)
  } else if (my_integer_dc == 'character') {
    char_integers_with_comma <- format(as.numeric(my_integers), big.mark=",", scientific=FALSE)
  } else {
    stop(sprintf('cannot support the data type [%s]',my_integer_dc))
  }
  return(char_integers_with_comma)
}

wo_bgnd_ggplot <- function() {
	wo_bgnd <- theme_bw() + 
		theme(panel.border = element_blank(), panel.grid.major = element_blank(),
											 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	
	return(wo_bgnd)
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
# 
# b37_to_hg19_bam <- function(bam){
# 	chr_bam <- sprintf('%s.chr.bam',bam)
# 	cmd <- sprintf("samtools view -H %s | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | grep -vP '(SN:GL)|(SN:hs)' | samtools reheader - %s > %s",bam,bam,chr_bam)
# 	system(cmd)
# 	return(chr_bam)
# }

decompose_fname <- function(filepath2,debug2=0) {
  if (debug2==1){browser()}
  fbase <- basename(filepath2)
  if (endsWith(fbase,'.gz')) {
    fbase <- sub(".gz","",fbase)
  }
  
  #finally strip off the original file extension
  fbase0 <- tools::file_path_sans_ext(fbase)
  
  parent_dir=sapply(filepath2,function(fpath) {
  	dirsegs <- tstrsplit(dirname(fpath),.Platform$file.sep)
  	tail(dirsegs, n=1)[[1]]
  })
  
  fname <- data.table(dir=dirname(filepath2),
  										parentd=parent_dir,
                      fbase0=file_path_sans_ext(fbase),
                      ext=file_ext(fbase))
  
  return(fname)
}

get_file_base <- function(filepath2) {
  # message(filepath2)
  fbase <- basename(filepath2)
  if (endsWith(fbase,'.gz')) {
    fbase <- sub(".gz","",fbase)
  }
  
  #finally strip off the original file extension
  fbase0 <- tools::file_path_sans_ext(fbase)
  return(fbase0)
}

tag_file <- function(filepath2,tag=NA,new_ext=NA){
  
  message(sprintf('in:%s',filepath2))
  fname <- decompose_fname(filepath2)
  if (invalid(tag) & invalid(new_ext)) {
  	stop('assign a value to either tag or new_ext!')
  }
  
  if (invalid(new_ext)) {
    tagged_fname <- sprintf(file.path(fname$dir,sprintf('%s_%s.%s',fname$fbase0,tag,fname$ext)))
  } else {
  	if (invalid(tag)) {
  		tagged_fname <- sprintf(file.path(fname$dir,sprintf('%s.%s',fname$fbase0,new_ext)))
  	} else {
  		tagged_fname <- sprintf(file.path(fname$dir,sprintf('%s_%s.%s',fname$fbase0,tag,new_ext)))
  	}
  }
  
  message(sprintf('in:%s',tagged_fname))
  return(tagged_fname)
}

page_to_print <- function(range_str) {
  
  range2 <- strsplit(range_str,',')[[1]]
  range2a <- strsplit(range2,'-')
  range0 <- list()
  for (i in length(range2a)) {
    range3 <- range2a[[i]]
    range4 <- as.numeric(unlist(range3))
    if (length(range4)==2) {
      range0[[i]] <- range4[1]:range4[2]
    } else {
      range0[[i]] <- range4
    }
  }
  page_to_print <- unique(unlist(range0))
  return(page_to_print)
}

annotate_reads <- function(read_dt,is_paired) {
  read_dt$sample_label1 <- sapply(read_dt$r1,function(read_fpath) {
    fname <- decompose_fname(read_fpath)
    return(substr(fname$fbase0, 1, nchar(fname$fbase0)-3))
  })
  if (is_paired==1) {
    read_dt$sample_label2 <- sapply(read_dt$r2,function(read_fpath) {
      fname <- decompose_fname(read_fpath)
      return(substr(fname$fbase0, 1, nchar(fname$fbase0)-3))
    })
  }
  
  read_dt$dir1 <- sapply(read_dt$r1,function(read_fpath) {
    fname <- decompose_fname(read_fpath)
    return(fname$dir)
  })
  
  if (is_paired==1) {
    read_dt$dir2 <- sapply(read_dt$r2,function(read_fpath) {
      fname <- decompose_fname(read_fpath)
      return(fname$dir)
    })
  }
  
  if (is_paired==1) {
    is_valid <- all(read_dt[,sample_label1==sample_label2]) & all(read_dt[,dir1==dir2])
    if (is_valid) {
      read_dt$dir <- read_dt$dir1
      read_dt$sample_label <- read_dt$sample_label1
    } else {
      read_dt$dir <- NA
      read_dt$sample_label <- NA
    }
    read_dt[,sample_label1:=NULL]
    read_dt[,sample_label2:=NULL]
    read_dt[,dir1:=NULL]
    read_dt[,dir2:=NULL]
  } else {
    is_valid <- TRUE
    read_dt$dir <- read_dt$dir1
    read_dt$sample_label <- read_dt$sample_label1
  }
  
  return(list(is_valid=is_valid,read_dt=read_dt))
}

get_reads_from_prefix <- function(dir_pref,is_paired=1,suffix='fastq.gz') {
  message(sprintf('%s/*%s',dir_pref,suffix))
  ret <- sort(Sys.glob(sprintf('%s*/*%s',dir_pref,suffix)))
  reads <- list()
  if (is_paired==1) {
    reads[[1]] <- ret[endsWith(ret,'R1.fastq.gz')]
    reads[[2]] <- ret[endsWith(ret,'R2.fastq.gz')]
    
    read_dt <- data.table(r1=reads[[1]],r2=reads[[2]])
    ret <- annotate_reads(read_dt,is_paired)
    stopifnot(ret$is_valid)
    read_dt <- ret$read_dt
  } else {
    read_dt <- data.table(r1=ret)
    ret <- annotate_reads(read_dt,is_paired)
    read_dt <- ret$read_dt
  }
  return(read_dt)
}


run_system <- function(cmd,syslog_fn=NA,intern=TRUE,verbose=0,debug2=0) {
  if (debug2==1){browser()}
	if (verbose>0) {message(cmd)}
	
  sys_out <- system(cmd,intern=intern)
  
  if (intern & 'Error' %in% sys_out) {
    stop('error occurs while runing %s\nRefer %s',cmd,sys_out)
  }
  if (!invalid(syslog_fn)) {
    fwrite(sys_out,file=syslog_fn)
  }
  return(sys_out)
}


comma_string_to_list <- function(comma_string,sort_flag="",sep=",",keep_unique=1,debug2=0) {
	if (debug2==1) {browser()}
	items = str_trim(tstrsplit(comma_string,sep))
	if (keep_unique==1) {
		items <- unique(items)
	}
	
	if (sort_flag=="ascending") {
		items <- sort(items,decreasing = FALSE)
	} else if (sort_flag=="descending") {
		items <- sort(items,decreasing = TRUE)
	}
	if (debug2==1) {browser()}
	return(items)
}

try2 <- function(code, silent = FALSE) {
  tryCatch(code, error = function(c) {
    msg <- conditionMessage(c)
    if (!silent) message(c)
    invisible(structure(msg, class = "try-error"))
  })
}

concat_txt_files <- function(text_files,out_fn) {
  cmd <- sprintf("cat %s | gzip -fc > %s",paste0(text_files,sep=" ",collapse=""),out_fn)
  system(cmd)
}

copy_to_local <- function(remote_fpath,local_cached) {
  
	if (!file.exists(local_cached)) {dir.create(local_cached,showWarnings = F,recursive = T)}
  local_cpy <- file.path(local_cached,basename(remote_fpath))
  cmd <- sprintf("rsync -ravP %s %s/",remote_fpath,local_cached)
  system(cmd)
  message(sprintf("copied [%s] to [%s]",remote_fpath,local_cpy))
	
  return(local_cpy)
}

split_vec <- function(my_list,chunk_len=10) {
	split(my_list,ceiling(seq_along(my_list) / chunk_len))
}


chunk2 <- function(my_list,num_batches=3) {
	split(my_list, rep_len(1:num_batches, length(my_list)))
}

chunk_by_unitlen <- function(list2split,unit_len=10) {
	A <- length(list2split)
	if (A <= unit_len) {
		mbatches <- list(list2split)
	} else {
		num_batches <- ceiling(A/unit_len)
		b <- lapply(1:num_batches, function(j) {
			rep(j,unit_len)
		})
		f <- unlist(b)
		f <- f[1:A]
		mbatches <- split(list2split, f)
	}
	message(sprintf("total [%d] batches w/ [%d] elements!",length(mbatches),unit_len))
	mbatches
}

dt_to_xlsx.sheet <- function(dt.lst, xlsx_fpath = "./dts_in_wks.xlsx") {
  wb <- createWorkbook("wbk")
  for (name1 in names(dt.lst)) {
    sheet_name=name1
    if (nchar(name1)>30){
      sheet_name=strtrim(name1,30)
    }
    message(name1)
    addWorksheet(wb,sheetName = sheet_name)
    freezePane(wb, sheet=sheet_name,firstRow = TRUE, firstCol = TRUE)
    writeData(wb,sheet_name,dt.lst[[name1]],withFilter=TRUE)
  }
  saveWorkbook(wb, file=xlsx_fpath, overwrite = TRUE)
}

dt_to_xlsx <- function(input_dt,col1,xlsx_fpath,debug2=0) {
	if(debug2==1){browser()}
	dts1 <- split(input_dt,by=col1)
	dt_to_xlsx.sheet(dts1, xlsx_fpath=xlsx_fpath)
	if(debug2==1){browser()}
}

matdf_to_scina_marker_fmt <- function(matdf_tsv,out_csv) {
	dt <- fread(matdf_tsv)
	N <- dim(dt)[2]
	genes <- dt$marker
	dt$marker <- NULL
	M <- max(colSums(dt))
	
	markers2 <- apply(dt,2,function(mycol) {
		genes2 <- genes[mycol>0]
		M2 <- M - length(genes2)
		return(c(genes2,rep(NA,M2)))
	})
	
	fwrite(as.data.table(markers2),file=out_csv,sep=',')
}


fread_between_two_lines <- function(txt_file,from_word,to_word=NA,debug2=0) {
	if (debug2==1){browser()}
	dt2 <- NA
	if (file.exists(txt_file)) {
		#check if 'from_word' exists
		cmd<-sprintf("grep '%s' %s",from_word,txt_file)
		exit_status <- run_system(cmd,intern = FALSE)
		if (exit_status==0) {
			cmd <- sprintf("sed -n '/%s/=' %s",from_word,txt_file)
			line_str <- run_system(cmd)
			message(line_str)
			line_from <- as.numeric(line_str)
			
			if (!invalid(to_word)) {
				cmd<-sprintf("grep '%s' %s",to_word,txt_file)
				exit_status <- run_system(cmd,intern = FALSE)
				if (exit_status==1) {
					to_word <- NA
				}
			}
			
			if (invalid(to_word)) {
				dt <- fread(file=txt_file,skip=line_from-1)
			} else {
				cmd <- sprintf("sed -n '/%s/=' %s",to_word,txt_file)
				
				line_str <- run_system(cmd)
				message(line_str)
				line_to <- as.numeric(line_str) - 1
				# browser()
				dt <- fread(file=txt_file,skip=line_from-1,nrows=(line_to-line_from))
			}
			dt2<-as.data.table(dt)
		}
	}
	return(dt2)
}

rbindlist2 <- function(res_stat) {
	if (all(is.na(res_stat)|invalid(res_stat))) {
		res_stats <- NA
	} else {
		res_stats <- rbindlist(res_stat[!is.null(res_stat) & !is.na(res_stat)])
	}
	return(res_stats)
}

cut_nice <- function(x, lower = 0, upper, by,
										 sep = "-", above.char = "+",debug2=0) {
	
	if(debug2==1){browser()}
	labs <- c(paste(sprintf("%3.4f",seq(lower, upper - by, by = by)),
									sprintf("%3.4f",seq(lower + by, upper, by = by)),
									sep = sep),
						paste(sprintf("%3.3f",upper), above.char, sep = ""))
	
	cut(x, breaks = c(seq(lower, upper, by = by), Inf),
			right = TRUE, labels = labs)
}

get_current_script_fpath <- function() {
	
	if (rstudioapi::isAvailable()) {
		script_fpath <- rstudioapi::getSourceEditorContext()$path
	} else {
		args = commandArgs()
		script_fpath = args[substr(args,1,7) == '--file=']
		script_fpath <- substr(script_fpath, 8, nchar(script_fpath))
	}
	
	script_fpath <- normalizePath(script_fpath)
	
	message(script_fpath)
	
	fpath_loc <- data.table(fpath=script_fpath)
	fpath_loc <- cbind(fpath_loc,decompose_fname(script_fpath))
	fpath_loc[,step_name:=tstrsplit(fbase0,"_")[[1]]]
	# fpath_loc[,parent_dir:=basename(dir)]
	return(fpath_loc)
}

create_outdirectory <- function(args,script_loc) {
	outd<-file.path(args$baseoutd,script_loc$fbase0)
	if (!file.exists(outd)){dir.create(outd,showWarnings = F,recursive = T)}
	message("=====================")
	message(sprintf("output directory [%s]",outd))
	message("=====================")
	return(outd)
}

check_fpath <- function(fpath) {
	stopifnot(file.exists(fpath))
}

create_dir_if_not_exist <- function(newd,active=T) {
	if (!file.exists(newd)){
		if (active==T) {
			message(sprintf("creating %s ...",newd))
			dir.create(newd,showWarnings = F,recursive = T)
		}
	}
	message(newd)
	newd
}


get_outd <- function(outd0,script_fpath,create=1,debug2=0) {
	if (debug2==1){browser()}
	decfn <- decompose_fname(script_fpath)

	outd<-file.path(outd0,decfn$fbase0)
	
	if (create==1 & !file.exists(outd)) {dir.create(outd,showWarnings = F,recursive = T)}
	message(sprintf("output directory[%s]",outd))
	outd
}

get_wkd <- function(projd,script_fpath,create=1) {
	
	decfn <- decompose_fname(script_fpath)
	
	wkd<-file.path(projd,decfn$fbase0)
	
	if (create==1 & !file.exists(wkd)) {dir.create(wkd,showWarnings = F,recursive = T)}
	message(sprintf("output directory[%s]",wkd))
	wkd
}

get_pipeline_step <- function(script_fpath) {
	decfn <- decompose_fname(script_fpath)
	step_name <- tstrsplit(decfn$fbase0,"_")[[1]]
	message(sprintf("step_name[%s]",step_name))
	step_name
}

get_out_fpath <- function(outd,fname) {
	stopifnot(dir.exists(outd))
	out_fpath = file.path(outd,fname)
	message(sprintf("chekcing %s",out_fpath))
	out_fpath
}

get_proj_out0 <- function(proj_name) {
	file.path(Sys.getenv('HOME'),'projects',proj_name)
}

read_rd_or_rds <- function(rds_fpath) {
	message(sprintf("loading %s",rds_fpath))
	if (endsWith(rds_fpath,'.rd')) {
		ret <- get(load(rds_fpath))
	} else if (endsWith(rds_fpath,'.rds')) {
		ret <- readRDS(rds_fpath)
	} else {
		stop(sprintf("cannot support the file format %s",rds_fpath))
	}
	ret
}

check_prog_in_path <- function(bin_fpath) {
	cmd_str=sprintf("which %s",bin_fpath)
	ret <- run_system(cmd_str,intern = FALSE)
	bin_in_path = TRUE
	if (ret==1){
		bin_in_path=FALSE
	}
	bin_in_path
}

replace_env_variables <- function(str_w_senv,debug2=0) {
	if (debug2==1) {browser()}
	str_wo_senv <- tryCatch (
		{
			senv=str_extract_all(str_w_senv,"(?<=\\$\\{)\\S+(?=\\})",simplify=TRUE)[[1]]
			str_wo_senv=gsub(sprintf("\\$\\{%s\\}",senv),Sys.getenv(`senv`),str_w_senv)
			message(sprintf("replacing [%s] by [%s]",str_w_senv,str_wo_senv))
			str_wo_senv
		},
		error=function(cond) {
			return(NA)
		},
		finally={
		}
	)
	return(str_wo_senv)
}

list2_to_list1 <- function(list_of_dt){
	list_names = names(list_of_dt)
	fields = names(list_of_dt[[1]])
	dt.lst = lapply(fields,function(field) {
		message(field)
		list1=imap(list_of_dt,function(entry,list_name){
			entry[[field]]
		})
		names(list1) = list_names
		list1
	})
	names(dt.lst) = fields
	dt.lst
}

has_value <- function(item) {
	has_flag = TRUE
	if (isEmpty(item)) {
		has_flag = FALSE
	} else if (invalid(item)) {
		has_flag = FALSE
	} else if (invalid(item)) {
		has_flag = FALSE
	}
	has_flag
}

##################
topk_by_group <- function(my.dt,ugroup_by,field1,field2,field3=NA,topk=3,debug2=0) {
	
	if (debug2==1){browser()}
	cnames = colnames(my.dt)
	stopifnot(ugroup_by %in% cnames)
	stopifnot(field1 %in% cnames)
	stopifnot(field2 %in% cnames)
	
	topk_deg.dt = imap(split(my.dt,by=`ugroup_by`),function(by_gr1,gr1) {
	  message(gr1)
		N=nrow(by_gr1)
		topkj = topk
		if (topkj>N) {topkj = N}
		if (debug2==1){browser()}
		# if (gr1==2){browser()}
		if (invalid(field3)){
		  by_gr1=by_gr1[order(get(field1),get(field2),decreasing = TRUE),][1:topkj,]
		} else {
		  by_gr1=by_gr1[order(get(field1),get(field2),get(field3),decreasing = TRUE),][1:topkj,]
		}
		by_gr1[,ranking:=1:topkj]
	}) %>% rbindlist()
	topk_deg.dt[!is.na(get(ugroup_by)),]
}

sync_two_mtx <- function(mtx1,mtx2,debug2=0){
	if (debug2==1){browser()}
	shared_rnames=intersect(rownames(mtx1),rownames(mtx2))
	shared_cnames=intersect(colnames(mtx1),colnames(mtx2))
	
	mtx1=mtx1[shared_rnames,shared_cnames]
	mtx2=mtx2[shared_rnames,shared_cnames]
	
	list(mtx1=mtx1,mtx2=mtx2)
}

###############
approx_nz_pval <- function(deg.dt,pval_col="p_val",nz_pct=0.9) {
	deg.dt$old.p_val = rlang::duplicate(deg.dt[[pval_col]])
	deg.dt$pval_nz = deg.dt[[pval_col]]
	
	p_val_mz = deg.dt[get(pval_col)>0,min(get(pval_col))] * nz_pct
	deg.dt[pval_nz==0.,pval_nz:=p_val_mz]
	deg.dt[[pval_col]] = deg.dt$pval_nz
	deg.dt$mlog10pv = -log10(deg.dt$pval_nz)
	deg.dt$pval_nz = NULL
	deg.dt
}

approx_inf_logfc <- function(deg.dt,logfc="logFC") {
	deg.dt$old.logfc = rlang::duplicate(deg.dt[[logfc]])
	deg.dt$logfc1234 = deg.dt[[logfc]]
	logFC.max = 1.1 * deg.dt[!is.infinite(logfc1234) & logfc1234>0, max(logfc1234)]
	logFC.min = 1.1 * deg.dt[!is.infinite(logfc1234) & logfc1234<0, min(logfc1234)]
	
	deg.dt[is.infinite(logfc1234)&logfc1234>0,logfc1234:=logFC.max]
	deg.dt[is.infinite(logfc1234)&logfc1234<0,logfc1234:=logFC.min]
	deg.dt[[logfc]] = deg.dt$logfc1234
	deg.dt$alogfc = abs(deg.dt$logfc1234)
	deg.dt$logfc1234= NULL
	deg.dt
}
################
# wrap a string of column from input data table (obsolete)
wrap_long_names <- function(dt2,col2wrap) {
	L2 = round(dt2[,max(nchar(get(col2wrap)))]/2)
	stri_sub(dt2[[col2wrap]],L2,(L2-1)) = "\n"
	dt2
}

wrap_n <- function(x,width=35) {
	y <- gsub("_"," ",x)
	# y <- sub("\\s+$", "", gsub('(.{35})', '\\1 ', x))
	stringr::str_wrap(y,width=width)
}

# read matrix
# ##############
fread_as_matrix <- function(tsv_fn,rowname_col=1) {
	# browser()
	message(tsv_fn)
	stopifnot(file.exists(tsv_fn))
	my.m = fread(tsv_fn)
	if (is.character(rowname_col)) {
		stopifnot(rowname_col %in% colnames(my.m))
	}
	row_names = my.m[[rowname_col]]
	my.m[[rowname_col]] = NULL
	my.m = as.matrix(my.m)
	rownames(my.m) = row_names
	my.m
}

load_smk_echo_out <- function(smk_echo_fn) {
	stopifnot(file.exists(smk_echo_fn))
	out=fread(smk_echo_fn) %>% 
		colnames() %>%
		unique()
	
	data.table(out1=out,
						 sample=basename(out))
}

plan2 <- function(ncpu=1,nG=8) {
	if (!inherits(plan(), "sequential")) plan(sequential)
	if (ncpu>1){
		options(future.globals.maxSize = nG*1e3*1024^2)
		plan(multisession,workers=ncpu,gc=TRUE)
	}
}

# Function for computing Jaccard Similarity
jaccard_similarity <- function(A, B) {
	intersection = length(intersect(A, B))
	union = length(A) + length(B) - intersection
	return (intersection/union)
}


# shorten long string in a certain column in a given data.table
shorten_long_str_column <- function(dt2,col_orig,col_rev,trim_len=75,debug2=0) {
	
	if (debug2==1){browser()}
	
	if (nrow(dt2)>0) {
		stopifnot(col_orig %in% colnames(dt2))
		if ('col_long_tmp' %in% colnames(dt2)) {dt2$col_long_tmp = NULL}
		if (col_rev %in% colnames(dt2)) {dt2[[col_rev]] = NULL}
		
		dt2[,col_long_tmp:=get(col_orig)]
		
		trim_uniq = dt2[nchar(col_long_tmp)>trim_len,] %>%
			nrow() %>%
			1:.
		
		dt2[nchar(col_long_tmp)>trim_len, col_long_tmp:=sprintf("%s.%d",substr(col_long_tmp,1,trim_len),trim_uniq)]
		
		setnames(dt2,'col_long_tmp',col_rev)
	}
	dt2
}

get_comp_schs <- function(entries_to_comp,min_cnt=2,debug2=0) {
	if (debug2==1){browser()}
	stat.1=table(entries_to_comp)
	members = names(stat.1[stat.1>=min_cnt])
	comp.m = t(combn(members,2))
	lapply(1:dim(comp.m)[1],function(r2){
		c(comp.m[r2,])
	})
}

merge_fxs_with_smeta <- function(fxs.m, smeta.dt, feat="gene", sample_col="sample", value_name="data", debug2=0) {
	if (debug2==1){browser()}
	stopifnot(sample_col %in% colnames(smeta.dt))
	stopifnot(!(value_name %in% colnames(smeta.dt)))
	
	melt(fxs.m, varnames=c(feat,sample_col),value.name=value_name) %>%
		merge(smeta.dt,by=sample_col)
}

clip_fxs_matrix <- function(fxs.m,qlo=0.1,qhi=0.9,debug2=0) {
  if (debug2==1){browser()}
	apply(fxs.m, 1, function(r2) {
	  qco = quantile(r2, c(qlo, qhi))
	  r2[r2<=qco[[1]]]=qco[[1]]
	  r2[r2>=qco[[2]]]=qco[[2]]
	  r2
	}) %>%
	  t()
	
}

clip_and_scale<-function(fxs.m,qlo=0.05,qhi=0.95,scale.center=TRUE,debug2=0) {
  
  fxs.m.cs=clip_fxs_matrix(fxs.m,qlo=qlo,qhi=qhi,debug2=debug2) %>%
    t() %>%
    scale(center = scale.center) %>%
    t()
  
  #get rid of any NaN
  j <- rowSums(is.nan(fxs.m.cs))==0
  fxs.m.cs[j,]
}

reduce_fxs_matrix_clip <- function(fxs_log.m,min_frac_expr=0.2,var_method="var_mean_ratio",qlo=.05,qhi=.95,n_feats=2000,incl_feats=NA,corr_cutoff=0.5,ncpu=4,debug2=0) {
	
	if (debug2==1){browser()}
	
	fxs_log.m.bkp = rlang::duplicate(fxs_log.m)
	
	fxs_log.m = fxs_log.m[apply(fxs_log.m, 1, function(x) sum(x > 0)/length(x) > min_frac_expr), , drop = FALSE]
	N = dim(fxs_log.m)[1]
	
	if (n_feats>N) {
		n_feats = N
	}
	
	n_hvgs = n_feats
	if (corr_cutoff>0){
		n_hvgs = n_feats * 1.25
		if (n_hvgs>N) {
			n_hvgs=N
		}
	}
	
	if (var_method=="mad"){
		mads = apply(fxs_log.m,1,mad)
		ind = rev(order(mads))[1:n_hvgs]
	} else {
		ind = order(apply(fxs_log.m, 1, function(x) {
			q = quantile(x, c(qlo, qhi))
			x = x[x > q[1] & x < q[2]]
			var(x)/mean(x)
		}), decreasing = TRUE)[1:n_hvgs]
	}
	
	fxs_log.m = fxs_log.m[ind, , drop = FALSE]
	
	if (corr_cutoff>0.){
		dt = WGCNA::cor(t(fxs_log.m), nThreads = ncpu)
		diag(dt) = 0
		dt[abs(dt) < corr_cutoff] = 0
		dt[dt < 0] = -1
		dt[dt > 0] = 1
		
		# consider two variables to select top-k genes
		ind = order(colSums(abs(dt)),decreasing=T)[1:n_feats]
		rm(dt)
		fxs_log.m = fxs_log.m[ind, ,drop = FALSE]
	}
	
	if (!invalid(incl_feats)) {
		fxs_log.m = fxs_log.m.bkp[unique(c(incl_feats,rownames(fxs_log.m))),]
	}
	
	fxs_log.m = apply(fxs_log.m, 1, function(x) {
		q = quantile(x, c(qlo, qhi))
		x[x<q[1]] = q[1]
		x[x>q[2]] = q[2]
		x
	}) %>% t()
	fxs_log.m
}

impute_NA_fxs_mat <- function(fxs.m) {
  fxs.m = apply(fxs.m,1,function(r) {
    mean.r=mean(r[!is.na(r)])
    r[is.na(r)] = mean.r
    r
  }) %>% t()
  fxs.m
}
