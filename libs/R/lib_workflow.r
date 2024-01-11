#!/usr/bin/env Rscript

#Authors: Chanjing Hong
#Date Created: 2/1/2019

library(data.table)
library(stringr)
library(yaml)
# library(configr) #RcppTOML <mayo_setup>
library(httr)
library(gtools)
source(file.path(Sys.getenv('R_UTIL'),'lib_util.R'))

get_prog_info <- function(method,is_paired=1,reuse=0,debug=0,action='script') {
  if (debug==1) {browser()}

  #prog_config_fname <- get_prog_config_filename()
	rwkf.ini <- Sys.getenv('PROG_CONF')
	
	message(sprintf('prog_config[%s]',rwkf.ini))
	stopifnot(file.exists(rwkf.ini))

  prog <- list()
  prog$method <- method
  prog$exptag <- method
  prog$user_email <- ''
  prog$bin_path <- eval.config(file = rwkf.ini, config = method, value = "bin_path")
  prog$option1 <- eval.config(file = rwkf.ini, config = method, value = "default_opt")
  prog$option2 <- eval.config(file = rwkf.ini, config = method, value = "config_path_opt")
  prog$reqcmd <- eval.config(file = rwkf.ini, config = method, value = "reqcmd")
  prog$preamble <- eval.config(file = rwkf.ini, config = method, value = "preamble")
  prog$reuse <- reuse
  prog$debug <- debug
  prog$action <- action
  prog$is_paired <- is_paired
  return(prog)
}

set_user_option <- function(cmethod,args){
	
	opts <- names(args)
	
	message(opts)
	
	stopifnot('exptag' %in% opts)
	stopifnot('user_email' %in% opts)
	stopifnot('outd' %in% opts)

	if (!('reuse' %in% opts)) {args$reuse=0}
	if (!('is_paired' %in% opts)) {args$is_paired=1}
	if (!('ncpu' %in% opts)) {args$ncpu=2}
	if (!('hours' %in% opts)) {args$hours=24}
	if (!('reuse' %in% opts)) {args$reuse=0}
	if (!('memG' %in% opts)) {args$memG=48}
	if (!('run_opt' %in% opts)) {args$run_opt=""}
	
	cmethod$jobname = make.names(args$exptag)
	cmethod$user_email = args$user_email
	cmethod$outd = args$outd
	cmethod$reuse=  args$reuse
	cmethod$is_paired = args$is_paired
	cmethod$ncpu = args$ncpu
	cmethod$hours = args$hours
	cmethod$reuse = args$reuse
	cmethod$memG = args$memG
	cmethod$run_opt = args$run_opt
	
	if ('jobq' %in% names(args)) {
		cmethod$jobq = args$jobq
	} else {
		cmethod$jobq = ""
	}
	
	if ('nscripts' %in% names(args)) {
	  cmethod$nscripts = args$nscripts
	} else {
	  cmethod$nscripts = 1
	}
	
	return(cmethod)
}

locate_prog_config <- function(config.ini=NA,debug2=0) {
  
	if (debug2==1){browser()}
  if (invalid(config.ini)) {
  	#mount_prefix <- get_mount_dir()
    if (Sys.getenv("PROG_CONF")=="") {
  	  conf_file <- file.path(Sys.getenv("PPLINE"),'configs','local',sprintf('hwanglab_prog_conf_%s.yml',Sys.getenv("HOSTNAME")))
    } else {
      conf_file <- Sys.getenv("PROG_CONF")
    }
    if (length(conf_file)==0) {
      stop('check if [%s] exists and is valid',conf_file)
    }
  } else {
    conf_file <- config.ini
  }
  
  if (!file.exists(conf_file)) {
    stop(sprintf('check if [%s] exists !',conf_file))
  }
  return(conf_file)
}

load_prog_config <- function(prog_conf_file=NA,debug2=0) {
	
	prog_conf_file <- locate_prog_config(debug2=debug2)
	message(sprintf('prog_config[%s]',prog_conf_file))
	prog_yaml <- read_yaml(prog_conf_file)

	return(prog_yaml)
}

prep_update_worksheet <- function(in_wksheet_file,out_wksheet_file) {
  
  jmessage('prepare updated worksheet template ...')
  cmd <- sprintf("grep -vP '^#' %s > %s",in_wksheet_file,out_wksheet_file)
  message(cmd)
  system(cmd)
}

replace_env_vars <- function(str1) {
  # browser()
  str2 <- gsub(sprintf("\\$\\{%s\\}","APP"),Sys.getenv('APP'),str1)
  str2 <- gsub(sprintf("\\$\\{%s\\}","LABSHARED"),Sys.getenv('LABSHARED'),str2)
  str2 <- gsub(sprintf("\\$\\{%s\\}","TMP"),Sys.getenv('TMP'),str2)
  
  return(str2)
}

get_prog_info_yaml <- function(method=NA,debug2=0) {
	if (debug2==1) {browser()}
	prog_yaml <- load_prog_config(debug2=debug2)
	
	# stopifnot(method %in% names(prog_yaml))
	if ((method %in% names(prog_yaml))) {
		cmethod <- prog_yaml[[method]]
		cmethod$required_input_columns <- comma_string_to_list(replace_env_vars(cmethod$required_input_columns))
		cmethod$optional_input_columns <- comma_string_to_list(replace_env_vars(cmethod$optional_input_columns))
		
		if (!invalid(cmethod[['bin_path']])) {
			if (!startsWith(cmethod$bin_path,'/')) {
				subd_componets <- tstrsplit(replace_env_vars(cmethod$bin_path),'/')
				if (length(subd_componets)>1) {
					cmethod$bin_path <- replace_env_vars(cmethod$bin_path)
				}
			}
		}
		
		if (!invalid(cmethod[['resource']])) {
			for (entry in names(cmethod$resource)) {
				if (!invalid(cmethod$resource[[entry]])){
					if (endsWith(entry,'_path') & !startsWith(cmethod$resource[[entry]],'/')) {
						subd_componets <- tstrsplit(replace_env_vars(cmethod$resource[[entry]]),'/')
						if (length(subd_componets)>1) {
							cmethod$resource[[entry]] <- replace_env_vars(cmethod$resource[[entry]])
						}
					}
				}
			}
		}
	} else {
		cmethod = list()
		cmethod$bin_path = method
		cmethod$required_input_columns = NA
		cmethod$optional_input_columns = NA
		cmethod$resource = NA
	}
	
	return(cmethod)
}

load_method_in_wksheet <- function(worksheet_file) {
  # browser()
  cmd_str <- sprintf("grep -P '^#' %s | cut -f2,3,4,5",worksheet_file)
  message(cmd_str)
  sdt <- fread(cmd = cmd_str, header=F, col.names = c('step', 'item','value','condition'))
  return(sdt)
}


load_method_in_wksheet_old <- function(worksheet_file,stepi) {
  # browser()
  cmd_str <- sprintf("grep -P '^#' %s | cut -f2,3,4,5 | grep -w '^%d'",worksheet_file,stepi)
  message(cmd_str)
  sdt <- fread(cmd = cmd_str, header=F, col.names = c('step', 'item','value','condition'))
  return(sdt)
}

load_samples_in_wksheet <- function(worksheet_file) {
  # browser()
  cmd_str <- sprintf("grep -vP '^#' %s",worksheet_file)
  message(cmd_str)
  sdt <- fread(cmd = cmd_str, header=T)
  return(sdt)
}

check_consistency <- function(cmethod,worksheet) {
	stopifnot(all(cmethod$required_input_columns %in% colnames(worksheet)))
}

check_consistency_in_wksheet <- function(mdt,sdt) {
  
  # browser()
  
  mdt_req <- mdt[condition=='required' & is.na(value),]
  
  if (dim(mdt_req)[1]>0) {
    stop('empty value in a required field in the worksheet file [method section]')
  }
  
  sample_cols <- colnames(sdt)
  mdt_sreq <- mdt[condition=='s.required',]
  
  for (i in 1:dim(mdt_sreq)[1]) {
    
    sreq_row <- mdt_sreq[i,]
    vals <- strsplit(sreq_row$value,';')[[1]]
    if (length(vals) > 0) {
      
      for (j in 1:length(vals)) {
        if (vals[j] %in% sample_cols) {
          message(sprintf('column[%s] verified',vals[j]))
        } else {
          stop(sprintf('check if column[%s] exist in the worksheet file',
                       vals[j]))
        }
      }
    } else {
      stop(sprintf('empty value at the row[%s] in the worksheet file:method',
                   sreq_row$item))
    }
  }
  
}

launch_pipeline <- function(runinst,sample,inputs,output_dir,debug=0,reuse=0) {
  # browser()
  if (debug==0) {
    message(output_dir)
    if (!dir.exists(output_dir)){
      dir.create(output_dir,recursive=TRUE)
    }
  } else {
    # browser()
    x<-1
  }
  
  if (endsWith(runinst$prog_path,'.py')){
    cmd <- sprintf("python %s",runinst$prog_path)
  } else if (endsWith(runinst$prog_path,'.r')) {
    cmd <- sprintf("Rscript %s",runinst$prog_path)
  } else {
    cmd <- runinst$prog_path
  }
  
  outputs <- list()
  
  if (!invalid(runinst$common_opt)) {
    cmd <- sprintf("%s %s",cmd,runinst$common_opt)
  }
  
  message(sprintf('generating a commandline[%s]',runinst$rstep$method))
  
  if (runinst$rstep$method %in% c('pathoqc','pathoscope2','salmon')) {
    I <- length(inputs)
    cmd <- sprintf("%s -1 %s",cmd,inputs[[1]])
    if (I==2) {
      cmd <- sprintf("%s -2 %s",cmd,inputs[[2]])
    }
    cmd <- sprintf("%s -o %s",cmd,output_dir)
   
  } else if (runinst$rstep$method == 'fastqc') {
    cmd <- sprintf("%s %s",cmd,paste(inputs,collapse = ' '))
    cmd <- sprintf("%s -o %s",cmd,output_dir)
  } else if (runinst$rstep$method == 'cutadapt') {
    
    cmd <- sprintf("%s -o %s",cmd,file.path(output_dir,basename(inputs[[1]])))
    if (length(inputs)==2) {
      cmd <- sprintf("%s -p %s",cmd,file.path(output_dir,basename(inputs[[2]])))
    }
    cmd <- sprintf("%s %s",cmd,paste(inputs,collapse = ' '))
    
  } else if (runinst$rstep$method == 'cellranger') {
    run_dir <- dirname(inputs[[1]])
    cmd <- sprintf("%s --id=%s",cmd,sample)
    cmd <- sprintf("%s --run=%s",cmd,run_dir)
    cmd <- sprintf("%s --samplesheet=%s",cmd,inputs[[1]])
    cmd <- sprintf("%s --output-dir %s",cmd,output_dir)
    
  } else {
    stop(sprintf('method [%s] is not supported yet',runinst$rstep$method))
  }
  
  if (!invalid(runinst$add_opt)) {
    cmd <- sprintf("%s %s",cmd,runinst$add_opt)
  }
  
  outfpat <- runinst$rstep$mdt[item=='output_file_suffix',value]
  fpats <- strsplit(outfpat,';')[[1]]
  
  for (j in 1:length(fpats)) {
    outfpat_n <- fpats[j]
    if (startsWith(outfpat_n,'*')) { #check if this is for a suffix type
      outputs[[j]] <- file.path(output_dir,outfpat_n)
      message(outputs[[j]])
    } else {
      outputs[[j]] <- outfpat
    }
  }
  
  log_file <- file.path(output_dir,sprintf("%s.log",sample))
  if (F) {
    cmd <- sprintf("%s 2>&1 | tee %s",cmd,log_file)
  }
  
  message(cmd)
  if (debug==0) {
    out_fpath <- file.path(output_dir,outputs[[1]])
    message(sprintf('expected output file[%s]',out_fpath))
    if (reuse==1 & file.exists(out_fpath)) {
      message('reuse prev result ...')
    } else {
      #message(sprintf('************check <%s>**************',out_fpath))
      system(cmd) #debug
    }
  }
  return(outputs)
}

gen_cmd <- function(cmethod,cmd,out1) {
  if (cmethod$action=="run") {
    if (cmethod$reuse==1) {
      if (file.exists(out1)) {
        message(sprintf('The expected output file [%s] already exists and skip running',out1))
      } else {
        message(cmd)
        system(cmd)
      }
      cmd <- ""
      status <- 'successful'
    }
  } else {
    status <- cmethod$action
    if (file.exists(out1)) {
      cmd <- ""
      status <- 'successful'
    }
  }
  return(list(cmd=cmd,status=status))
}


run_module <- function(cmethod,args) {
  # browser()
  if (cmethod$method == "cutadapt") {
    cmethod <- cutadapt2(cmethod)
  } else if (cmethod$method == "salmon") {
    cmethod <- salmon2(cmethod)
  } else {
    message('not impleted yet...')
  }
  return(cmethod)
  
}

prep_pipeline_shortreads <- function(method,args,action) {
  
	if (args$debug==1){browser()}
  cmethod <- get_prog_info(method,
                           args$is_paired,
                           reuse=args$reuse,
                           debug=args$debug,
                           action=action)

  cmethod <- set_user_runoption(args$runopt,cmethod)
  
  cmethod$input <- get_reads_from_prefix(args$in_pat,
                                         args$is_paired,
                                         args$suffix)
  
  cmethod$input$cmd <- ""
  if (args$debug==1){browser()}
  return(cmethod)
}

salmon_quant <- function(args,wks,action,debug2=args$debug) { #prev_cmethod,method,action
	if (debug2==1){browser()}
	
	#$SALMON quant -i $SALMON_IDX -l IU -p 20 -1 $SALMON_IN/$mMDSC_1_p1 -2 $SALMON_IN/$mMDSC_1_p2 --numBootstraps 100 -o $SALMON_OUT/mMDSC_1
	#
	message('salmon_quant ...')
	cmethod <- get_prog_info_yaml('salmon')
	cmethod <- set_user_option(cmethod,args)
	cmethod$action <- action
	
	refidx <- cmethod$config_path_opt
	if (args$refidx!="") {refidx <- args$refidx}
	
	runopt <- cmethod$default_opt
	if (args$runopt!="") {runopt <- args$runopt}
	
	M <- dim(wks)[1]
	
	ret <- lapply(1:M,function(i) {
		# if (debug2==1){browser()}
		cmd <- cmethod$bin_path
		cmd <- sprintf("%s quant",cmd)
		
		cmd <- sprintf("%s -i %s",cmd,refidx)
		
		if (runopt!="") {
			cmd <- sprintf("%s %s",cmd,runopt)
		}
		
		cmd <- sprintf("%s -p %d",cmd,cmethod$ncpu)
		
		r1fpath <- wks$read1[i]
		r2fpath <- wks$read2[i]
		
		cmd <- sprintf("%s -1 %s",cmd,r1fpath)
		cmd <- sprintf("%s -2 %s",cmd,r2fpath)
		
		if (("outd" %in% colnames(wks))) {
			outd <- normalizePath(wks$outd[i])
		} else {
			stopifnot(args$outd!="")
			outd <- normalizePath(args$outd)
		}
		
		if (!file.exists(outd)) {
			message(sprintf("creating [%s]...",outd))
			dir.create(outd,recursive = T)
		}
		
		quant_outd <- file.path(outd,'salmon_quant',wks$sample[i])
		if (!file.exists(quant_outd)) {
			dir.create(quant_outd,recursive = T)
		}
		
		cmd <- sprintf("%s -o %s",cmd,quant_outd)
		
		out1 <- file.path(quant_outd,"quant.sf")
		
		if (cmethod$reuse==1 & file.exists(out1)) {
			cmd <- sprintf("## reuse [%s]",out1)
		}
		
		if (action == "run" & !startsWith(cmd,"## reuse")) {
			run_system(cmd)
		}
		# browser()
		return(list(sample=wks$sample[i],read1=r1fpath,read2=r2fpath,cmd=cmd,out1=out1))
	})
	
	cmethod$wks <- data.table(sample=sapply(1:M,function(m) {ret[[m]]$sample}),
														cmd=sapply(1:M,function(m) {ret[[m]]$cmd}),
														read1=sapply(1:M,function(m) {ret[[m]]$read1}),
														read2=sapply(1:M,function(m) {ret[[m]]$read2}),
														out1=sapply(1:M,function(m) {ret[[m]]$out1}))
	
	if (debug2==1){browser()}
	message('Done.')
	
	return(cmethod)
}


kallisto <- function(args,wks,kallisto_opt_str="",action="script",debug2=0) {
	# kallisto index -i hs_b38_v105 /fs/ess/PCCF0022/refdb/kallisto/Homo_sapiens.GRCh38.cdna.all.fa.gz
	
	if (debug2==1){browser()}
	message("kallisto ...")
	cmethod <- get_prog_info_yaml('kallisto')
	cmethod <- set_user_option(cmethod,args)
	cmethod$action <- action

	stopifnot(args$outd!="")
	outd <- normalizePath(args$outd)
	create_dir_if_not_exist(outd)
	
	refidx <- cmethod$resource[["index_path"]]
	if (args$refidx!="") {refidx <- args$refidx}
	
	M <- dim(wks)[1]
	ret <- lapply(1:M,function(i) {
		if (debug2==1){browser()}
		
		cmd <- cmethod$bin_path
		
		cmd <- sprintf("%s quant -i %s",cmd,refidx)
		
		outdj <- file.path(outd,wks$sample[i])
		if (!file.exists(outdj)) {
			dir.create(outdj,recursive = T)
		}
		cmd <- sprintf("%s -o %s",cmd,outdj)
		
		cmd <- sprintf("%s -b 100",cmd) #bootstrap-samples
		cmd <- sprintf("%s -t %d",cmd,args$ncpu)

		if (args$gen_bam==1) {
			cmd <- sprintf("%s --pseudobam",cmd)
			cmd <- sprintf("%s --gtf %s", cmd, cmethod$resource[['gtf_path']])
			cmd <- sprintf("%s --chromosomes %s",cmd,cmethod$resource[['chrom_dict']])
		}
		
		cmd <- sprintf("%s %s %s",cmd,wks$read1[i],wks$read2[i])
		return(list(sample=wks$sample[i],read1=wks$read1[i],read2=wks$read2[i],cmd=cmd,outputd=outdj))
	})
	
	cmethod$wks <- data.table(sample=sapply(1:M,function(m) {ret[[m]]$sample}),
														cmd=sapply(1:M,function(m) {ret[[m]]$cmd}),
														read1=sapply(1:M,function(m) {ret[[m]]$read1}),
														read2=sapply(1:M,function(m) {ret[[m]]$read2}),
														outputd=sapply(1:M,function(m) {ret[[m]]$outputd}))
	
	return(cmethod)
	
	message("Done")
}


parse_run_option <- function(arg_runopt) {
  runopts <- strsplit(arg_runopt,';')[[1]]
  runopt_list <- lapply(runopts,function(runopt){
    items <- strsplit(runopt,':')[[1]]
    method <- items[[1]]
    user_opt <- trimws(items[[2]])
    return(list(method,user_opt))
    })
  
  runopt_dt <- rbindlist(runopt_list)
  colnames(runopt_dt) <- c('method','option1')
  return(runopt_dt)
}

set_user_runoption <- function(user_runopt_str,cmethod) {
	if (user_runopt_str!="") {
	  user_opt_dt <- parse_run_option(user_runopt_str)
	  option1 <- user_opt_dt[method==cmethod$method,option1]
	  if (!is_empty(option1)) {
	    cmethod$option1 <- option1
	  }
	}
  return(cmethod)
}

get_job_scheduler_type <- function() {
	
	ret_success0 = system("which sbatch",ignore.stderr = TRUE)
	job_scheduler = "NotAvailable"
	if (ret_success0==0) {
		job_scheduler = "slurm"
	} else {
		ret_success0 = system("which qsub",ignore.stderr = TRUE)
		if (ret_success0==0) {
			job_scheduler = "torque"
		}
	}
	message(sprintf("job scheduler type[%s]",job_scheduler))
	job_scheduler
}

assign_batch_run <- function(tn4_dt,num_scripts=8,batch_run_prefix="br") {
	tn4_dt$sidx <- 1:nrow(tn4_dt)
	run_batches <- chunk_by_unitlen(1:nrow(tn4_dt),unit_len = ceiling(nrow(tn4_dt)/num_scripts))
	names(run_batches) <- sprintf("%s.%s",batch_run_prefix,names(run_batches))
	run_schedule = rbindlist(lapply(names(run_batches),function(brun){
		data.table(brun,sidx=run_batches[[brun]])
	}))
	tn4_dt$brun <- run_schedule[match(tn4_dt$sidx,run_schedule$sidx),brun]
	tn4_dt$sidx <- NULL
	tn4_dt
}

get_db_fpath <- function(proj_info.lst=get_proj_info(),debug2=0) {
	if (debug2==1){browser()}
	stopifnot('resource_db_path' %in% names(proj_info.lst))
	res_fn = proj_info.lst$resource_db_path
	create_dir_if_not_exist(dirname(file.path(res_fn)))
	res_fn
}

load_db_table <- function(debug2=0) {
	if (debug2==1){browser()}
	message('use load_obj_from_db() to load obj file ...')
	res_fn = get_db_fpath(debug2=debug2)
	res.dt = NA
	if (file.exists(res_fn)) {
		res.dt = fread(res_fn)
	}
	res.dt
}

prep_resource_entry <- function(idx, name, rds_fpath, tag="initial", comment="default") {
	
	stopifnot(dir.exists(dirname(rds_fpath)))
	
	data.table(idx=idx,
						 proj_name = project_name,
						 name=name,
						 tag=tag,
						 fpath=rds_fpath,
						 timestamp=Sys.time(),
						 is_active=1,
						 comment=comment)
}

register_db_table <- function(entry.dt, file.obj=NULL, res.dt=NULL, update=1,debug2=0) {
	if (debug2==1){browser()}
	
	if (update==1 | !file.exists(entry.dt$fpath)) {
		if (!invalid(file.obj)) {
			message(sprintf('saving rds file [%s] ...',entry.dt$fpath))
			saveRDS(file.obj, file = entry.dt$fpath)
		}
	}
	
	if (invalid(res.dt)) {
		res.dt = rlang::duplicate(entry.dt)
	} else {
		if (entry.dt$idx==0) {
			entry.dt$idx = max(res.dt$idx) + 1
		} else {
			res.dt = res.dt[idx != entry.dt$idx,]
		}
		res.dt = rbind(res.dt, entry.dt)
	}
	res.dt
}

load_obj_from_db <- function(query_name, debug2=0) {
	if (debug2==1){browser()}
	target_idx = get_db_idx(query_name, debug2=debug2)
	file.obj = list()
	if (!isEmpty(target_idx)) {
		file.obj = load_db_table() %>%
			dplyr::filter(idx==target_idx) %>%
			dplyr::slice(1) %>%
			dplyr::pull(fpath) %>%
			readRDS()
	}
	file.obj
}

get_db_idx <- function(query_name, allow_multiple=0,debug2=0) {
	if (debug2==1){browser()}
	db.dt = load_db_table()
	
	if (!invalid(query_name)) {
		target_idx = db.dt[name == query_name, idx]
		if (!isEmpty(target_idx)) {
		  if (allow_multiple==0) {
			  stopifnot(length(target_idx)<2)
		  }
		}
	} else {
		stop('query_name (+query_tag) should be given!')
	}
	target_idx
}

drop_db_resource <- function(query_name, rm_file=0,  debug2=0) {
	if (debug2==1){browser()}
	db_fpath = get_db_fpath()
	if (file.exists(db_fpath)) {
		res.dt = load_db_table()
		target_idx = get_db_idx(query_name, allow_multiple=1,debug2=0)
		if (!isEmpty(target_idx)) {
			if (rm_file==1) {
				message('deleting the actual file is not implemented yet ...')
				# res.dt[idx %in% target_idx, unlink(fpath)]
			}
			res.dt = res.dt[!(idx %in% target_idx),]
			fwrite(res.dt, file=db_fpath, sep="\t")
		}
	} else {
		stop('drop_db_resource() requires db_table()!')
	}
}

update_db_resource <- function(query_name, rds_fpath, file.obj=NULL, query_tag=NA, comment="desc", update=1, debug2=0) {
	if (debug2==1){browser()}
	db_fpath = get_db_fpath()
	resource_entry = prep_resource_entry(idx=0,name=query_name,rds_fpath=rds_fpath,tag=query_tag,comment=comment)
	if (file.exists(db_fpath)) {
		res.dt = load_db_table()
		target_idx = get_db_idx(query_name, debug2=debug2)
		if (isEmpty(target_idx)) {
			#add new entry
			res.dt = register_db_table(resource_entry, file.obj, res.dt, update=update, debug2=debug2)
			fwrite(res.dt, file=db_fpath, sep="\t")
		} else if (update==1) {
			#update
			resource_entry$idx = target_idx
			res.dt = register_db_table(resource_entry, file.obj, res.dt, update=update, debug2=debug2)
			fwrite(res.dt, file=db_fpath, sep="\t")
		} else {
			message('already exists ...')
		}
	} else {
		resource_entry$idx=1
		res.dt = register_db_table(resource_entry, file.obj,update=update, debug2=debug2)
		fwrite(res.dt, file=db_fpath, sep="\t")
	}
	res.dt
}
