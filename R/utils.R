#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Utility Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(biomaRt)
library(rBLAST)
library(stringr)
library(seqinr)
library(D3GB)
library(cowplot)
library(TmCalculator)
library(stringi)
library(Hmisc)
library(pbmcapply)
library(pbapply)
library(Biostrings)
library(intervals)
library(IRanges)
library(future)
library(future.apply)
library(reshape2)

break_str <- function(sequence){
  return(toupper(strsplit(sequence,"")[[1]]))
}

Tm_calc <- function(probe, ambiguous = FALSE, comSeq = NULL, shift = 0, nn_table = "DNA_NN4", 
                tmm_table = "DNA_TMM1", imm_table = "DNA_IMM1",de_table = "DNA_DE1", dnac1 = 5,
                dnac2 = 5, selfcomp = FALSE, Na = 330, K = 0, Tris = 0, Mg = 0, dNTPs = 0, saltcorr = 5,Tm, DMSO = 0, fmd = 20, DMSOfactor = 0.75, 
                fmdfactor = 0.65, fmdmethod = "concentration", ptGC = NULL){
  tm <- Tm_NN(probe, ambiguous, comSeq, shift, nn_table, tmm_table, imm_table,de_table, dnac1, dnac2, selfcomp, Na, K, Tris, Mg, dNTPs, saltcorr)
  tm_salt <- chem_correction(tm, DMSO, fmd, DMSOfactor, fmdfactor, fmdmethod, ptGC)
  if(tm_salt < 0){
    warning("Negative Tm resets to zero \n", call. = FALSE)
    tm_salt <- 0
  }
  return(tm_salt)
}

gc_content <- function(sequence){
  strings <- break_str(sequence)
  return(length(which(strings %in% c("G","C")))/length(strings))
}

find_consecut <- function(sequence, times = 5, base = NULL){
  if(is.null(base)){
    a_n <- strrep("A",times)
    t_n <- strrep("T",times)
    g_n <- strrep("G",times)
    c_n <- strrep("C",times)
    return(str_count(sequence, pattern = paste(a_n,t_n,g_n,c_n,sep = "|")))
  } else {
    base_n <- strrep(base, times)
    return(str_count(sequence, pattern = base_n))
  }
}

get_sequences <- function(transcript_id,mart,genetable,ensembl_table,species = "human",seqtype="cdna"){
  seq_all <- data.frame(NULL)
    if(species == "mouse"){
      seq_all <- biomaRt::getSequence(id = transcript_id, type="ensembl_transcript_id", mart = mart, seqType = seqtype,verbose = F)
      ensembl.names <- genetable[match(seq_all[["ensembl_transcript_id"]], genetable[["transcript_id"]]),"gene_id"]
      gene.names <- ensembl_table[match(ensembl.names, ensembl_table[["ensembl_gene_id"]]),"mgi_symbol"]
    } else {
      seq_all <- biomaRt::getSequence(id = transcript_id, type="ensembl_transcript_id", mart = mart, seqType = seqtype,verbose = F)
      ensembl.names <- genetable[match(seq_all[["ensembl_transcript_id"]], genetable[["transcript_id"]]),"gene_id"]
      gene.names <- ensembl_table[match(ensembl.names, ensembl_table[["ensembl_gene_id"]]),"hgnc_symbol"]
    }
  if(nrow(seq_all) == 0){
    if(nrow(seq_all) == 0) warning("No ENSEMBL sequence returned for ", transcript_id)
    return(NULL)
  } else {
    seq_meta <- biomaRt::getBM(attributes=c('chromosome_name', 'transcript_start','transcript_end','exon_chrom_start','exon_chrom_end','strand'),
                               filters = 'ensembl_transcript_id',
                               values = transcript_id,
                               mart = mart)
    if(seq_meta[["strand"]][1] == 1){
      seq_meta <- seq_meta[order(seq_meta[["exon_chrom_start"]],decreasing = FALSE),]
    } else {
      seq_meta <- seq_meta[order(seq_meta[["exon_chrom_start"]],decreasing = TRUE),]
    }
    seq_all[["ensembl_gene_id"]] <- ensembl.names
    seq_all[["GeneName"]] <- gene.names
    seq_all[["Chromosome"]] <- paste("chr", seq_meta[["chromosome_name"]][1], sep = "")
    seq_all[["TranscriptStart"]] <- seq_meta[["transcript_start"]][1]
    seq_all[["TranscriptEnd"]] <- seq_meta[["transcript_end"]][1]
    seq_all[["ExonStarts"]] <- paste(seq_meta[["exon_chrom_start"]], collapse = ",")
    seq_all[["ExonEnds"]] <- paste(seq_meta[["exon_chrom_end"]], collapse = ",")
    seq_all[["Strand"]] <- seq_meta[["strand"]][1]
    seq_all <- seq_all[which(seq_all[,1] != "Sequence unavailable"),]
    if(nrow(seq_all) == 0){
      warning("No ENSEMBL sequence returned for ", transcript_id)
      return(NULL)
    } else {
      return(seq_all)
    }
  }
}

prob_extract <- function(sequence,probe_len = 25, consecutive_repeats = 5){
  probe_len <- probe_len-1
  if(str_length(sequence) < probe_len){
    warning("The probe length is longer than the sequences!\n", call. = FALSE)
    return(data.frame(NULL))
  } else {
    if(str_length(sequence) == probe_len){
      start_seq <- 1
      stop_seq <- start_seq + probe_len
    } else {
      start_seq <- 1:(str_length(sequence)-probe_len)
      stop_seq <- start_seq + probe_len
    }
    gc.list <- lapply(as.list(1:length(start_seq)), function(i){
      seq.temp <- substr(sequence,start_seq[i],stop_seq[i])
      gc.temp <- gc_content(seq.temp)
      consecut.temp <- find_consecut(seq.temp,times = consecutive_repeats)
      return(c(Sequence=seq.temp,GC=gc.temp,Consecutive=consecut.temp))
    })
    gc.data <- data.frame(do.call(rbind,gc.list),stringsAsFactors = F)
    gc.data[,2] <- as.numeric(gc.data[,2])
    gc.data[,3] <- as.numeric(gc.data[,3])
    gc.data[["Start"]] <- start_seq
    gc.data[["Stop"]] <- stop_seq
    return(gc.data)
  }
}

gc_filter <- function(gc_data,gc_min = 0.45,gc_high = 0.65,gc_target = 0.55,n_consecutive = 1){
  gc_filtered <- gc_data[which(gc_data[["Consecutive"]] < n_consecutive),]
  gc_filtered <- gc_filtered[which(gc_filtered[["GC"]] > gc_min & gc_filtered[["GC"]] < gc_high),]
  gc_filtered <- gc_filtered[order(abs(gc_filtered[["GC"]] - gc_target),decreasing = F),]
  return(gc_filtered)
}

blast_homolog <- function(gc_data,blast_db,BLAST_args = "-task blastn-short",parallel = FALSE){
  if(nrow(gc_data) == 0){
    warning("No data passed GC content.")
    return(data.frame(NULL))
  } else {
    if(parallel){
      bl.data <- future_lapply(gc_data[["Sequence"]],function(x){
        seq.temp <- DNAStringSet(x)
        cl <- rBLAST:::predict.BLAST(blast_db, seq.temp,BLAST_args = BLAST_args,check.names = FALSE)
        return(cl)
      })
    } else {
      bl.data <- pbmclapply(gc_data[["Sequence"]],function(x){
        seq.temp <- DNAStringSet(x)
        cl <- rBLAST:::predict.BLAST(blast_db, seq.temp,BLAST_args = BLAST_args,check.names = FALSE)
        return(cl)
      })
    }
    names(bl.data) <- gc_data[["Sequence"]]
    return(bl.data)
  }
}

blast_filter <- function(bl_data,transcriptID,genename,gc_data,alignment = 15,filter_trans = FALSE, probe_mode = FALSE, names.delim = "-", names.field = 1, parse = TRUE){
  if(nrow(gc_data) == 0){
    warning("No data passed GC content.")
    return(data.frame(NULL))
  } else {
    if(probe_mode){
      bl_filter <- lapply(bl_data,function(x){
        x <- x[-1,]
        temp <- x[which(x[["Alignment.Length"]] >= alignment),]
        temp[["SubjectID"]] <- as.character(temp[["SubjectID"]])
        temp[["ProbeSet"]] <- unlist(lapply(temp[["SubjectID"]], function(i){
          strsplit(i,names.delim)[[1]][names.field]
        }))
        return(temp)        
      })
      gc_data[["BLAST_Hits"]] <- unlist(lapply(bl_filter,nrow))
      gc_data[["BLAST_Hits_ProbeSet"]] <- unlist(lapply(bl_filter,function(x) paste(unique(x[["ProbeSet"]]),collapse = ",")))
    } else {
      if(parse == TRUE){
        bl_filter <- lapply(bl_data,function(x){
          temp <- x[which(x[["Alignment.Length"]] >= alignment),]
          temp[["SubjectID"]] <- as.character(temp[["SubjectID"]])
          temp[["TranscriptID"]] <- unlist(lapply(temp[["SubjectID"]], function(i){
            strsplit(strsplit(i,"[|]")[[1]][1],"[.]")[[1]][1]
          }))
          temp[["Genename"]] <- unlist(lapply(temp[["SubjectID"]], function(i){
            strsplit(i,"[|]")[[1]][6]
          }))
          if(filter_trans){
            temp <- temp[which(temp[["TranscriptID"]] != transcriptID),]
          } else {
            temp <- temp[which(temp[["Genename"]] != genename),]
          }
          return(temp)        
        })
      } else {
        bl_filter <- lapply(bl_data,function(x){
          temp <- x[which(x[["Alignment.Length"]] >= alignment),]
          temp[["TranscriptID"]] <- as.character(temp[["SubjectID"]])
          return(temp)        
        })
      }
      gc_data[["BLAST_Hits"]] <- unlist(lapply(bl_filter,nrow))
      gc_data[["BLAST_Hits_Transcripts"]] <- unlist(lapply(bl_filter,function(x) paste(unique(x[["TranscriptID"]]),collapse = ",")))
      gc_data[["BLAST_Hits_Gene"]] <- unlist(lapply(bl_filter,function(x) paste(unique(x[["Genename"]]),collapse = ",")))
    }
    return(gc_data)
  }
}

get_expression <- function(gc_data,gtex_expr,tissue_type = "All",trans_ex = NULL,genes_ex = NULL,gene_table, trans_table){
  if(nrow(gc_data) == 0){
    warning("No data passed GC content.")
    return(data.frame(NULL))
  } else {
    trans.hits <- lapply(gc_data[["BLAST_Hits_Transcripts"]],function(x) unlist(strsplit(x,",")))
    genes.hits <- lapply(gc_data[["BLAST_Hits_Gene"]],function(x) unlist(strsplit(x,",")))
    if(!is.null(genes_ex)){
      genes.hits <- lapply(genes.hits,function(x) x[which(!x %in% genes_ex)])
      genes.id <- gene_table[,1][which(gene_table[,2] %in% genes_ex)]
      trans_ex <- trans_table[,1][which(trans_table[,2] %in% genes.id)]
    }
    if(!is.null(trans_ex)) trans.hits <- lapply(trans.hits,function(x) x[which(!x %in% trans_ex)])
    gc_data[["BLAST_Hits_Transcripts"]] <- unlist(lapply(trans.hits,paste,collapse=","))
    gc_data[["BLAST_Hits_Gene"]] <- unlist(lapply(genes.hits,paste,collapse=","))
    if("All" %in% tissue_type){
      genes.expr <- lapply(trans.hits,function(x){
        genes.match <- which(rownames(gtex_expr) %in% x)
        base::colSums(log(gtex_expr[genes.match,,drop=FALSE]+1,2))
      })
    } else {
      genes.expr <- lapply(trans.hits,function(x){
        genes.match <- which(rownames(gtex_expr) %in% x)
        base::colSums(log(gtex_expr[genes.match,tissue_type,drop=FALSE]+1,2))
      })
    }
    genes.expr <- do.call(rbind,genes.expr)
    gc_data[["Off_Targets_Expr"]] <- rowMeans(genes.expr)
    return(gc_data)
  }
}

seq_select <- function(gc_data,TPM_cutoff = 15,n_space = 0,order_by_feature = NULL,feature_target = NULL,start_name = "Start", seqname = "Sequence"){
  if(nrow(gc_data) == 0){
    warning("No sequence passed QC.")
    return(data.frame(NULL))
  } else {
    gc_tpm <- gc_data[which(gc_data[["Off_Targets_Expr"]] <= TPM_cutoff),]
    if(nrow(gc_tpm) == 0){
      warning("No sequence passed QC.")
      return(data.frame(NULL))
    } else {
      if(!is.null(order_by_feature)){
        gc_ordered <- gc_tpm[order(abs(gc_tpm[[order_by_feature]]-feature_target),decreasing = F),]
      } else {
        gc_ordered <- gc_tpm[order(gc_tpm[[start_name]],decreasing = F),]
      }
      start_i <- gc_ordered[[start_name]][1]
      len_seq <- str_length(gc_data[[seqname]][1])
      seq_selected <- 1
      if(max(gc_ordered[[start_name]])-len_seq-n_space-1 > 0){
        while(start_i <= (max(gc_ordered[[start_name]])-len_seq-n_space-1)){
          selected_i <- which(gc_ordered[[start_name]] >= start_i+len_seq+n_space)[1]
          seq_selected <- c(seq_selected,selected_i)
          start_i <- gc_ordered[[start_name]][selected_i]
        }
        return(gc_ordered[unique(seq_selected),])
      } else {
        return(gc_ordered[1,])
      }
    }
  }
}

iter_select.old <- function(gc_data,TPM_cutoff = 15,n_space = 0,order_by_feature = NULL,feature_target = NULL, start_name = "Start", stop_name = "Stop", seqname = "Sequence"){
  gc_selected <- seq_select(gc_data,TPM_cutoff = TPM_cutoff,n_space = n_space,order_by_feature = order_by_feature,feature_target = feature_target, start_name = start_name, seqname = seqname)
  gc_selected_new <- gc_selected
  while(nrow(gc_selected_new) != 0){
    gc_intv <- Intervals(gc_data[,c(start_name,stop_name)])
    select_intv <- Intervals(gc_selected[,c(start_name,stop_name)])
    select_intv <- as.matrix(select_intv)
    select_intv[,1] <- select_intv[,1] - n_space 
    select_intv[,2] <- select_intv[,2] + n_space 
    select_intv <- Intervals(select_intv)
    gc_size <- unlist(pblapply(1:nrow(gc_intv@.Data), function(i){
      temp_size <- sum(size(interval_intersection(select_intv,gc_intv[i]))+1)
      if(length(temp_size) == 0) temp_size <- 0
      temp_size
    }))
    gc_data <- gc_data[which(gc_size == 0),]
    gc_selected_new <- seq_select(gc_data,TPM_cutoff = TPM_cutoff,n_space = n_space,order_by_feature = order_by_feature,feature_target = feature_target, start_name = start_name, seqname = seqname)
    gc_selected <- rbind(gc_selected,gc_selected_new)
  }
  return(gc_selected)
}

iter_select <- function(gc_data,TPM_cutoff = 15,n_space = 0,order_by_feature = NULL,feature_target = NULL, start_name = "Start", stop_name = "Stop", seqname = "Sequence"){
  gc_selected <- seq_select(gc_data,TPM_cutoff = TPM_cutoff,n_space = n_space,order_by_feature = order_by_feature,feature_target = feature_target, start_name = start_name, seqname = seqname)
  gc_selected_new <- gc_selected
  while(nrow(gc_selected_new) != 0){
    test1 <- gc_selected[[stop_name]][-length(gc_selected[[stop_name]])]
    test2 <- gc_selected[[start_name]][-1]
    test1 <- c(test1,gc_selected[[stop_name]][length(gc_selected[[stop_name]])]) + n_space
    test2 <- c(test2 - n_space,max(gc_data[[stop_name]]))
    if(gc_selected[[start_name]][1] != 1) {
      test1 <- c(0,test1)
      test2 <- c(gc_selected[[start_name]][1] - n_space, test2)
    }
    row_selected <- base::t(base::apply(gc_data,1,function(x){
      as.numeric(x[[start_name]]) > test1 & as.numeric(x[[stop_name]]) < test2
    }))
    row_selected <- apply(row_selected,1,any)
    gc_data <- gc_data[which(row_selected),]
    gc_selected_new <- seq_select(gc_data,TPM_cutoff = TPM_cutoff,n_space = n_space,order_by_feature = order_by_feature,feature_target = feature_target, start_name = start_name, seqname = seqname)
    gc_selected <- rbind(gc_selected,gc_selected_new)
    gc_selected <- gc_selected[order(gc_selected[[start_name]]),]
  }
  return(gc_selected)
}

find_rev_compl <- function(multi_data){
  conv_seq <- sapply(multi_data[["Sequence"]],function(s){
    paste(rev(break_str(chartr("ATGC","TACG",s))),collapse = "")
  })
  return(cbind(Probe=conv_seq,multi_data,stringsAsFactors = FALSE))
} 

plot_gtex <- function(genename,data_mean,data_pct,ensemble,gene_table,species = "human"){
  specie_name <- ifelse(species == "human","hgnc_symbol","mgi_symbol")
  ensembl_id <- ensemble[which(ensemble[[specie_name]] %in% genename),]
  gene.match <- which(gene_table[["gene_id"]] %in% unique(ensembl_id[["ensembl_gene_id"]]))
  data1 <- t(log(data_mean[gene.match,,drop=FALSE]+1,2))
  data2 <- t(data_pct[gene.match,,drop=FALSE])
  data1 <- reshape2::melt(data1)
  colnames(data1)[3] <- "Log_TPM"
  data2 <- reshape2::melt(data2)
  data.all <- cbind(data1,Percent_Expr = round(data2[,3] * 100))
  pl1 <- ggplot(data = as.data.frame(data.all), mapping = aes(x = Var2, y = Var1)) + geom_point(mapping = aes(size = Percent_Expr, color = Log_TPM)) + 
    scale_color_distiller(palette = "Spectral") + xlab("") + ggtitle("Tissue Specific Expression") + ylab("") + scale_size_continuous(range = c(1,5),breaks = seq(0,100,by = 25)) + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
  return(pl1)
}

get_gtex <- function(genename,data_mean,data_pct,ensemble,gene_table,species = "human",plot = TRUE){
  specie_name <- ifelse(species == "human","hgnc_symbol","mgi_symbol")
  ensembl_id <- ensemble[which(ensemble[[specie_name]] %in% genename),]
  gene.match <- which(gene_table[["gene_id"]] %in% unique(ensembl_id[["ensembl_gene_id"]]))
  data1 <- t(log(data_mean[gene.match,,drop=FALSE]+1,2))
  data2 <- t(data_pct[gene.match,,drop=FALSE])
  data1 <- reshape2::melt(data1)
  colnames(data1)[1:3] <- c("Tissue", "Transcript","Log_TPM")
  data2 <- reshape2::melt(data2)
  data.all <- cbind(data1,Percent_Expr = round(data2[,3] * 100))
  return(data.all)
}

getBed <- function(gene.name,trans_used,specie = "human",mart){
  message("Getting BED file informain...")
  getAtt <- c('chromosome_name', 'start_position', 'end_position',
              'strand','ensembl_transcript_id','exon_chrom_start','exon_chrom_end')
  if(specie == "human"){
    trans <- getBM(attributes=getAtt,filters="hgnc_symbol",values = gene.name,mart=mart)
  } else {
    trans <- getBM(attributes=getAtt,filters="mgi_symbol",values = gene.name,mart=mart)
  }
  message("Generating BED files...")
  trans <- trans[trans$ensembl_transcript_id %in% trans_used,]
  trans.split <- split(trans,trans$ensembl_transcript_id)
  trans.split <- lapply(trans.split,function(x) x[order(x[,6],decreasing = F),])
  blockCount <-  unlist(lapply(trans.split,nrow))
  blockStarts <- unlist(lapply(trans.split,function(x){
    if(x$strand[1] == -1){
      return(paste(x[,3]-x[,7],collapse = ','))
    } else {
      return(paste(x[,6]-x[,2],collapse = ','))
    }}))
  blockSize <- unlist(lapply(trans.split,function(x) paste(x[,7]-x[,6]+1,collapse = ',')))
  end <- unlist(lapply(trans.split,function(x) x[1,3]-x[1,2]+1))
  strand <- unlist(lapply(trans.split,function(x) x[1,4]))
  score <- names(trans.split)
  genesBED <- data.frame(chr=gene.name,start=0,end=end,name=gene.name,score=score,strand=strand,thickStart=end,thickEnd=end,itemRGB=NA,blockCount=blockCount,blockSize=blockSize,blockStarts=blockStarts)
  return(genesBED)
}

generateBrowser <- function(gene.bed,dir.name = "tempBrowser",specie = 'human',selected,mart){
  message("Generating genome browser...")
  if(specie == 'human'){
    seq.all <- biomaRt::getSequence(id=gene.bed[["name"]],type = "hgnc_symbol",seqType = "gene_exon_intron",mart = mart)
  } else {
    seq.all <- biomaRt::getSequence(id=gene.bed[["name"]],type = "mgi_symbol",seqType = "gene_exon_intron",mart = mart)
  }
  write.fasta(seq.all[1,1],names = gene.bed$chr[1],file.out = "temp.fa")
  selected <- which(rownames(gene.bed) == selected)
  gene.bed1 <- gene.bed[selected,]
  gene.bed2 <- gene.bed[-selected,]
  tracks <- list(Selected = gene.bed1, Isoforms = gene.bed2)
  types <- c("exons","exons")
  colors <- c("#B8860B","#000000")
  gb <- genomebrowser(getAssemblyFromFasta("temp.fa"), tracks, types, colors, dir = dir.name)
  genome_addSequence(gb, "temp.fa")
  invisible(file.remove("temp.fa"))
  return(gb)
}

generateTracks.blast <- function(probe_list,seq.all,gb,BLAST_args = "-task blastn-short",order = FALSE){
  if(length(probe_list) > 1000) stop("More than 1000 probes found, please first split them into multiple <=1000 probe set!")
  tempName <- paste("temp",seq.all[,2],".fa",sep = "")
  fileConn <- file(tempName,"w")
  writeLines(paste(">",seq.all[,2],sep = ""), fileConn)
  writeLines(seq.all[,1], fileConn)
  close(fileConn)
  makeblastdb(tempName)
  bl.used <- blast(db=tempName,type = "blastn")
  bl.used.homolog <- lapply(probe_list,function(x){
    seq.temp <- DNAStringSet(x)
    cl <- predict(bl.used, seq.temp,BLAST_args = BLAST_args,check.names = FALSE)
    return(cl)
  })
  pb_len <- str_length(as.character(probe_list[1]))
  names(bl.used.homolog) <- unlist(probe_list)
  bl.used.homolog <- lapply(bl.used.homolog,function(x){
    x[order(x[["E"]],x[["Q.start"]],decreasing = F),]
  })
  blockStarts <- unlist(lapply(bl.used.homolog,function(x){
    if(x[1,7] == 1 & x[1,8] == pb_len){
      return(x[1,9]-1)
    } else if (x[1,7] != 1 & x[1,8] == pb_len) {
      if(any(x[,7] == 1)){
        return(paste(x[1:which(x[,7] == 1)[1],9]- 1,collapse=","))
      } else {
        return(x[1,9]-1)
      }
    } else if (any(x[,8] == 25)) {
      return(x[1,9]-1)
    } else {
      return(x[1,9]-1)
    }
  }))
  blockSize <- unlist(lapply(bl.used.homolog,function(x){
    if(x[1,7] == 1 & x[1,8] == pb_len){
      return(x[1,8]-x[1,7]+1)
    } else if(x[1,7] != 1 & x[1,8] == pb_len){
      if(any(x[,7] == 1)){
        return(paste(x[1:which(x[,7] == 1)[1],8]-x[1:which(x[,7] == 1)[1],7] + 1,collapse=","))
      } else {
        return(x[1,8]-x[1,7]+1)
      }
    } else if  (any(x[,8] == 25)){
      return(paste(x[1:which(x[,8] == 25)[1],8]-x[1:which(x[,8] == 25)[1],7] + 1,collapse=","))
    } else {
      return(x[1,8]-x[1,7]+1)
    }
  }))
  blockCount <- unlist(lapply(blockSize, function(x) length(strsplit(as.character(x),"[,]")[[1]])))
  score <- paste(seq.all[,2],1:length(bl.used.homolog),sep = "_")
  probesBED <- data.frame(chr=seq.all[,2],start=0,end=str_length(seq.all[1,1]),name=score,score=score,strand=NA,thickStart=str_length(seq.all[1,1]),thickEnd=str_length(seq.all[1,1]),itemRGB=NA,blockCount=blockCount,blockSize=blockSize,blockStarts=blockStarts)
  if(order){
    probesBED <- probesBED[order(probesBED$blockStarts,decreasing = F),]
  }
  types <- c("exons")
  colors <- c("red")
  tracks <- list(Probes = probesBED)
  genome_addTrack(gb,track = probesBED,trackname = "Probes",type = types,color = colors)
  invisible(file.remove(list.files(pattern = tempName)))
}

generateTracks <- function(probes,transcript,mart,seq.all,gb,order = FALSE){
  if(nrow(probes) == 0) stop("No probe found.", call. = FALSE)
  if(nrow(probes) > 1000) stop("More than 1000 probes found, split into several <=1000 probe objects using SubsetProbes first!", call. = FALSE)
  seq.pos <- getBM(attributes = c("start_position","end_position","exon_chrom_start","exon_chrom_end"), filters = c("ensembl_transcript_id"), values = transcript, mart = mart)
  strand.use <- getBM(attributes = c("strand"), filters = c("ensembl_transcript_id"), values = transcript, mart = mart)
  if(strand.use == 1){
    seq.pos <- seq.pos[order(seq.pos[["exon_chrom_start"]]),]
    seq.adj <- seq.pos[["exon_chrom_start"]][1] - seq.pos[["start_position"]][1] 
  } else if(strand.use == -1){
    seq.pos <- seq.pos[order(seq.pos[["exon_chrom_end"]],decreasing = TRUE),]
    seq.adj <- seq.pos[["end_position"]][1] - seq.pos[["exon_chrom_end"]][1] 
  }
  seq.pos <- seq.pos[,c("exon_chrom_start","exon_chrom_end")]
  seq.pos <- seq.pos - seq.pos[["exon_chrom_start"]][1] + 1
  if(strand.use == -1) seq.pos[,2:1] <- max(seq.pos[,1:2]) - seq.pos[,1:2] + 1
  seq.pos <- seq.pos[order(seq.pos[["exon_chrom_start"]],decreasing = FALSE),]
  start.pos <- seq.pos[,2] - seq.pos[,1]
  start.pos <- cumsum(start.pos+1) + 1
  start.pos <- c(seq.pos[1,1],start.pos)
  seq.pos[["cdna_start"]] <- start.pos[-length(start.pos)]
  seq.pos[["cdna_end"]] <- start.pos[-1] - 1
  probe_len <- str_length(probes[["Sequence"]][1])
  blockStarts <- list()
  blockSizes <- list()
  for (i in 1:nrow(probes)){
    starting <- rev(which(probes[["Start"]][i] >= seq.pos[["cdna_start"]]))[1]
    ending <- which(probes[["Stop"]][i] <= seq.pos[["cdna_end"]])[1]
    blockstart <- c()
    blocksize <- c()
    j <- starting
    while(j <= ending){
      if(j == starting){
        blockstart <- c(blockstart, seq.pos[["exon_chrom_start"]][j] + probes[["Start"]][i] - seq.pos[["cdna_start"]][j] - 1)
        if(blockstart + probe_len <= seq.pos[["exon_chrom_end"]][j]){
          blocksize <- c(blocksize, probe_len)
        } else {
          blocksize <- c(blocksize, seq.pos[["exon_chrom_end"]][j] - blockstart)
        }
      } else {
        blockstart <- c(blockstart, seq.pos[["exon_chrom_start"]][j] - 1)
        if(seq.pos[["exon_chrom_start"]][j] + probe_len - cumsum(blocksize) <= seq.pos[["exon_chrom_end"]][j]){
          blocksize <- c(blocksize, probe_len - cumsum(blocksize))
        } else {
          blocksize <- c(blocksize, seq.pos[["exon_chrom_end"]][j] - seq.pos[["exon_chrom_start"]][j] + 1)
        }
      }
      j <- j + 1
    }
    blockStarts[[i]] <- paste(blockstart + seq.adj,collapse = ",")
    blockSizes[[i]] <- paste(blocksize,collapse = ",")
  }
  blockCount <- unlist(lapply(blockSizes, function(x) length(strsplit(as.character(x),"[,]")[[1]])))
  probe.names <- paste(seq.all[,2],1:nrow(probes),sep = "_")
  block.first <- lapply(blockStarts, function(x) as.numeric(strsplit(as.character(x),"[,]")[[1]][1]))
  blockStarts <- unlist(blockStarts)
  blockSizes <- unlist(blockSizes)
  probesBED <- data.frame(chr=seq.all[,2],start=0,end=str_length(seq.all[1,1]),name=probe.names,score=probe.names,strand=NA,thickStart=str_length(seq.all[1,1]),thickEnd=str_length(seq.all[1,1]),itemRGB=NA,blockCount=blockCount,blockSize=blockSizes,blockStarts=blockStarts)
  if(order){
    probesBED <- probesBED[order(probesBED$blockStarts,decreasing = F),]
  }
  types <- c("exons")
  colors <- c("red")
  tracks <- list(Probes = probesBED)
  genome_addTrack(gb,track = probesBED,trackname = "Probes",type = types,color = colors)
}

calc_2nd_struct <- function(sequence_data,cores = -1){
  temp.data <- t(run_RNAfold(sequence_data[["Sequence"]],parallel.cores = cores))
  sequence_data[["Structure"]] <- temp.data[,2]
  sequence_data[["Min_Free_Energy"]] <- as.numeric(temp.data[,3])
  return(sequence_data)
}

concat_probes <- function(sequences, indice, insetion = "T"){
  concated <- paste(sequences[indice],collapse = insetion)
  return(concated)
}

get_repmask <- function(trans_data, gc_data, genome, masks_used = NULL, exon_only = TRUE){
  chr <- genome[[trans_data[["Chromosome"]]]]
  if(!is.null(masks_used)) {
    if("All" %in% masks_used){
      IRanges::active(Biostrings::masks(chr)) <- TRUE
    } else {
      mask.names <- names(IRanges::active(Biostrings::masks(chr)))
      if(!masks_used %in% mask.names) stop("masks_used not found in all available built-in masks.", call. = FALSE)
      IRanges::active(Biostrings::masks(chr))[masks_used] <- TRUE
    }
  }
  chr <- as(chr, "XStringViews")
  chr <- intervals::Intervals(cbind(BiocGenerics::start(chr@ranges),BiocGenerics::end(chr@ranges)))
  if(exon_only){
    exons_starts <- as.numeric(strsplit(trans_data[["ExonStarts"]],",")[[1]])
    exons_ends <- as.numeric(strsplit(trans_data[["ExonEnds"]],",")[[1]])
    exon_data <- intervals::Intervals(cbind(exons_starts,exons_ends))
    convert_data <- exon_data - min(exon_data) + 1
    if(trans_data[["Strand"]] == -1){
      convert_data <- max(convert_data) - convert_data + 1
      if(nrow(convert_data) == 1)  convert_data <- as.matrix(convert_data)
      convert_data <- convert_data[,2:1,drop=FALSE]
    }
    seq_data <- convert_data
    if(nrow(exon_data) > 1){
      for(i in 1:(nrow(exon_data)-1)){
        seq_data[i+1,1] <- seq_data[i,2]+1
        seq_data[i+1,2] <- seq_data[i+1,1] + convert_data[i+1,2] - convert_data[i+1,1]
      } 
    }
    exon_intersect <- lapply(1:nrow(exon_data), function(x) intervals::interval_intersection(exon_data[x], chr))
    exon_position <- lapply(1:nrow(exon_data), function(x){
      if(nrow(exon_intersect[[x]] != 0)){
        temp <- exon_intersect[[x]] - exon_data[x][[1]] + convert_data[x][[1]]
      } else {
        temp <- exon_intersect[[x]]
      }
      return(temp)
    })
    exon_interval <- lapply(1:nrow(exon_data), function(x){
      return(exon_position[[x]] - convert_data[x][[1]] + seq_data[x][[1]])
    })
    exon_interval <- as(do.call(rbind, exon_interval),"Intervals")
    gc_data_interval <- as(as.matrix(gc_data[,c("Start","Stop")]),"Intervals")
    message("Calculating Repetitive Masks for ",trans_data$ensembl_transcript_id)
    gc_data_intersect <- pbapply::pblapply(1:nrow(gc_data_interval),function(x) intervals::interval_intersection(gc_data_interval[x],exon_interval))
    gc_len <- lapply(gc_data[["Sequence"]],stringr::str_length)
    gc_data_count <- lapply(1:length(gc_data_intersect), function(x){
      count <- gc_len[[x]]
      if(nrow(gc_data_intersect[[x]]) != 0) count <- gc_len[[x]] - (sum(gc_data_intersect[[x]][,2] - gc_data_intersect[[x]][,1]) + nrow(gc_data_intersect[[x]]))
      return(count)
    })
  } else {
    trans_strat <- as.numeric(trans_data[["TranscriptStart"]])
    trans_end <- as.numeric(trans_data[["TranscriptEnd"]])
    trans_intv <- intervals::Intervals(cbind(trans_strat,trans_end))
    gc_data_interval <- as(as.matrix(gc_data[,c("Start","Stop")] + trans_strat),"Intervals")
    gc_data_intersect <- pbapply::pblapply(1:nrow(gc_data_interval),function(x) intervals::interval_intersection(gc_data_interval[x],chr))
    gc_len <- lapply(gc_data[["Sequence"]],stringr::str_length)
    gc_data_count <- lapply(1:length(gc_data_intersect), function(x){
      count <- gc_len[[x]]
      if(nrow(gc_data_intersect[[x]]) != 0) count <- gc_len[[x]] - (sum(gc_data_intersect[[x]][,2] - gc_data_intersect[[x]][,1]) + nrow(gc_data_intersect[[x]]))
      return(count)
    })
  }
  gc_data[["RepMask_Length"]] <- unlist(gc_data_count)
  return(gc_data)
}

find_excutable <- function(exe){
  path <- Sys.which(exe)
  if(path == "") stop("Executable for ", exe, " not found!\nPlease make sure that the software is correctly installed and path variables are set.",call. = FALSE)
  return(path)
}

create_binding_data <- function(sequence, probe, prefix = NULL){
  if(is.null(prefix)) prefix = "temp"
  filepath <- paste(prefix, ".in", sep = "")
  fileConn <- file(filepath,"w")
  writeLines("2", fileConn)
  writeLines(as.character(sequence), fileConn)
  writeLines(as.character(probe), fileConn)
  writeLines("2", fileConn)
  close(fileConn)
}

create_binding_concentration <- function(sequence_conc = 1e-6, probe_conc = 1e-6, prefix = NULL){
  if(is.null(prefix)) prefix = "temp"
  filepath <- paste(prefix, ".con", sep = "")
  fileConn <- file(filepath,"w")
  writeLines(as.character(sequence_conc), fileConn)
  writeLines(as.character(probe_conc), fileConn)
  close(fileConn)
}

nupack_wrapper <- function(nupack_function, params = "", prefix = "temp"){
  find_excutable(nupack_function)
  cmd <- paste(find_excutable(nupack_function), params, prefix)
  system(paste(cmd))
}

calc_2nd_struct_nupack <- function(sequence, params = "", filepath = NULL, prefix = "temp"){
  if(!is.null(filepath)){
    filepath_prefix <- paste(filepath, prefix, sep = "/")
  } else {
    filepath <- "."
    filepath_prefix <- prefix
  }
  filename <- paste(filepath_prefix, ".in", sep = "")
  fileConn <- file(filename,"w")
  writeLines(as.character(sequence), fileConn)
  close(fileConn)
  nupack_wrapper("mfe", params, prefix = filepath_prefix)
  filename <- paste(filepath_prefix, ".mfe", sep = "")
  data <- read.csv(filename, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  data.temp <- data[which(!grepl("%", data[,1]))[2:3],]
  invisible(file.remove(list.files(path = filepath, pattern = prefix)))
  return(list(Structure = as.numeric(data.temp[1]), Min_Free_Energy = data.temp[2]))
}

duplex_binding <- function(sequence, probe, seq_M = 1e-6, prob_M = 1e-6, params = "", filepath = NULL, prefix = "temp"){
  if(!is.null(filepath)){
    filepath_prefix <- paste(filepath, prefix, sep = "/")
  } else {
    filepath <- "."
    filepath_prefix <- prefix
  }
  create_binding_data(sequence, probe, filepath_prefix)
  create_binding_concentration(seq_M, prob_M, filepath_prefix)
  nupack_wrapper("complexes", paste(params,"-quiet"), prefix = filepath_prefix)
  nupack_wrapper("concentrations", params = "-quiet", prefix = filepath_prefix)
  filename <- paste(filepath_prefix, ".eq", sep = "")
  out <- system(sprintf("wc -l %s",filename),intern=TRUE)
  n <- as.integer(sub(sprintf("[ ]*([0-9]+)[ ]%s",filename),"\\1",out))
  data <- read.csv(filename,skip = n-5,sep = "\t",header = FALSE)
  idx <- which(data[,3] == 1 & data[,4] == 1)
  pct <- data[idx,6] / sum(data[,6])
  invisible(file.remove(list.files(path = filepath, pattern = prefix)))
  return(list(Duplex_Energy = data[idx,5], Equilibrium_Percent = pct))
}

check_any_feature <- function(object, slot = "Target", feature){
  if(slot == "Readout"){
    feature_exist <- is.null(object@probes[[slot]][[feature]])
  } else if(slot == "Target"){
    feature_exist <- unlist(lapply(object@probes[[slot]],function(x){
      is.null(x[[feature]])
    }))
  }
  return(any(feature_exist))
}

grid_search <-  function(object, slot = "Target", transcript_id, features_list, target_number = 30, do_filter = FALSE, order_by_target = TRUE){
  filter.list <- lapply(features_list,function(x){
    if(length(x) == 0) stop("No min or max provided", call. = FALSE)
    if(!"min" %in% names(x)) x[["min"]] <- rep(-Inf, length(x[[1]]))
    if(!"max" %in% names(x)) x[["max"]] <- rep(Inf, length(x[[1]]))
    return(x)
  })
  filter.list.data <- data.frame(filter.list)
  data.new <- list()
  for(i in 1:(ncol(filter.list.data)/2)){
    data.new[[i]] <- apply(filter.list.data[,(2*i-1):(2*i)], 1, paste, collapse = " to ")
  }
  names(data.new) <- names(filter.list)
  data.new <- data.frame(data.new,stringsAsFactors = FALSE)
  param.grid <- expand.grid(data.new,stringsAsFactors = FALSE)
  filter.names <- names(filter.list)
  filtered.obj <- apply(param.grid, 1, function(x){
    grid.data <- data.frame(strsplit(x," to "), stringsAsFactors = FALSE)
    min.use <- as.numeric(grid.data[1,])
    max.use <- as.numeric(grid.data[2,])
    obj.temp <- FilterProbes(object, slot = slot, transcript_ids = transcript_id, filter_names = filter.names, min = min.use, max = max.use, do.filter = TRUE, warning = FALSE)
  })
  filtered.number <- lapply(filtered.obj, function(x) nrow(x@probes[[slot]][[transcript_id]]))
  data.res <- data.frame(param.grid)
  data.res <- cbind(data.res, filtered = unlist(filtered.number))
  if(order_by_target){
    data.order <- order(abs(data.res[["filtered"]] - target_number), decreasing = FALSE)
  } else {
    data.order <- order(data.res[["filtered"]], decreasing = TRUE)
  }
  filtered.obj <- filtered.obj[data.order]
  data.res <- data.res[data.order, ]
  if(do_filter){
    return(filtered.obj[[1]])
  } else {
    return(data.res)
  }
}

add_base <- function(..., base){
  paste(list(...), collapse = base)
}

codebook_index <- function(binary_codebook){
  mat_idx <- apply(binary_codebook, 1, function(x){
    which(x == 1)
  })
  return(t(mat_idx))
}

get_hamming_dist <- function(binary_codebook, min_hamming = NULL){
  binary_codebook <- as.matrix(binary_codebook)
  hamming_binary <- function(X){
    D <- X %*% t(1 - X)
    D + t(D)
  }
  ham_dist <- hamming_binary(binary_codebook)
  if(!is.null(min_hamming)){
    ham_dist[which(ham_dist < min_hamming)] = 0
    ham_dist[which(ham_dist >= min_hamming)] = 1
  }
  rownames(ham_dist) <- 1:nrow(ham_dist)
  colnames(ham_dist) <- 1:ncol(ham_dist)
  return(ham_dist)
}

codebook_convert <- function(idx_table, conversion_table, readout_probe){
  conversion_seq <- apply(conversion_table, 1, function(x){
    readout_probe[["Sequence"]][as.numeric(x[1])]
  })
  conversion_table_seq <- cbind(conversion_table, Sequence = conversion_seq)
  idx_conversion <- t(apply(idx_table, 1, function(x){
    conversion_table_seq[["Sequence"]][match(x, conversion_table_seq[["bit"]])]
  }))
  return(idx_conversion)
}

get_exon <- function(object, transcript_id){
  seq_data <- object@sequences[which(object@sequences[["ensembl_transcript_id"]] == transcript_id),]
  starts <- as.numeric(strsplit(seq_data[["ExonStarts"]],",")[[1]])
  ends <- as.numeric(strsplit(seq_data[["ExonEnds"]],",")[[1]])
  exons <- intervals::Intervals(cbind(starts,ends), type = "Z")
  return(exons)
}

get_junction <- function(object, transcript_id, do_exon = FALSE){
  probe_data <- intervals::Intervals(object[[transcript_id]][,c("Start","Stop")])
  trans_data <- object@sequences[which(object@sequences[["ensembl_transcript_id"]] %in% transcript_id),]
  exons_starts <- as.numeric(strsplit(trans_data[["ExonStarts"]],",")[[1]])
  exons_ends <- as.numeric(strsplit(trans_data[["ExonEnds"]],",")[[1]])
  exon_data <- intervals::Intervals(cbind(exons_starts,exons_ends))
  convert_data <- exon_data - min(exon_data) + 1
  if(trans_data[["Strand"]] == -1){
    exon_data <- exon_data[,2:1]
    convert_data <- max(convert_data) - convert_data + 1
    if(nrow(convert_data) == 1)  convert_data <- as.matrix(convert_data)
    convert_data <- convert_data[,2:1,drop=FALSE]
  }
  seq_data <- convert_data
  if(nrow(exon_data) > 1){
    for(i in 1:(nrow(exon_data)-1)){
      seq_data[i+1,1] <- seq_data[i,2]+1
      seq_data[i+1,2] <- seq_data[i+1,1] + convert_data[i+1,2] - convert_data[i+1,1]
    } 
  }
  if(do_exon){
    exon_data_new <- intervals::Intervals(seq_data)
    exon_intv <- intervals::interval_overlap(probe_data,exon_data_new)
    exon_info <- unlist(lapply(exon_intv, function(x){
      paste("Exon",x,sep = "_", collapse = "&")
    }))
    return(exon_info)
  } else {
    seq_data <- c(t(as.matrix(seq_data)))
    seq_data <- seq_data[-c(1,length(seq_data))]
    seq_data <- matrix(seq_data,ncol = 2,byrow = TRUE)
    seq_data <- intervals::Intervals(seq_data)
    probe_intv <- intervals::interval_included(probe_data, seq_data)
    idx <- unlist(lapply(probe_intv, length))
    probe_intv[which(idx == 0)] <- NA
    exon_data_intv <- unlist(lapply(probe_intv, function(x){
      if(is.na(x)){
        c("NA,NA")
      } else {
        data1 <- as.matrix(exon_data)[x,2]
        data2 <- as.matrix(exon_data)[x+1,1]
        data.temp <- cbind(data1, data2)
        data.temp <- apply(data.temp, 1, paste, collapse = ",")
        paste(data.temp, collapse = ";")
      }
    }))
    exon_data_intv[which(exon_data_intv == "NA,NA")] = "None"
    return(exon_data_intv)
  }
}

check_transcripts <- function(object, transcript_ids){
  if(!all(transcript_ids %in% Transcripts(object))) stop(paste("Cannot find transcript", paste(transcript_ids[which(!transcript_ids %in% Transcripts(object))])), call. = FALSE)
}

split_probes <- function(sequence, min, nick){
  len <- stringr::str_length(sequence)
  max <- len - min - nick
  probe_list_all <- c()
  for(i in min:max){
    probe_a <- substr(sequence,1,i)
    probe_b <- substr(sequence,i+1+nick,len)
    probe_list_a <- data.frame(Original = sequence, Part = "A", Sequence = probe_a, Start = 1, Stop = i, stringsAsFactors = FALSE)
    probe_list_b <- data.frame(Original = sequence, Part = "B", Sequence = probe_b, Start = i+1+nick, Stop = len, stringsAsFactors = FALSE)
    probe_list <- rbind.data.frame(probe_list_a, probe_list_b)
    probe_list_all <- rbind.data.frame(probe_list_all, probe_list, stringsAsFactors = FALSE)
  }
  return(probe_list_all)
}

find_rc <- function(sequence){
  seq_rc <- find_rev_compl(data.frame(Sequence = toupper(sequence)))
  return(seq_rc[1,1])
}

