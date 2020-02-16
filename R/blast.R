#' @include utils.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# BLAST Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Write to FASTA File
#'
#' Write sequences into FASTA format.
#' @param object FISHprobe object.
#' @param sequences Sequences to write into FASTA format.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param filepath  FASTA file path and name. Default is to write to temp.fa in current directory.
#' @param seq_name Sequence names or suffix. Default is adding Seq suffix.
#' @export
#' @examples \donttest{
#' probe.data = GetSlotData(object, slotname = "probes", subname = "Target")
#' CreateFASTAFile(probe.data[[1]][["Sequence"]], filepath = "probe.fa")
#' }
CreateFASTAFile <- function(object = NULL, sequences = NULL, transcript_ids = NULL, slot = "Target", filepath = "temp.fa", seq_name = "Seq"){
  if(is.null(object) & is.null(sequences)){
    stop("Either object or sequences must be provided", call. = FALSE)
  } else if(is.null(object)){
    fileConn = file(filepath,"w")
    if(length(sequences) == length(seq_name)){
      for(i in 1:length(sequences)){
        writeLines(paste(">",seq_name[i],sep = ""), fileConn)
        writeLines(as.character(sequences[i]), fileConn)
      }
    } else {
      for(i in 1:length(sequences)){
        writeLines(paste(">",as.character(seq_name),"-",i,sep = ""), fileConn)
        writeLines(as.character(sequences[i]), fileConn)
      }
    }
    close(fileConn)
  } else {
    data.use = GetSlotData(object, slotname = "probes", subname = slot)
    if(slot == "Readout"){
      data.use = list(data.use)
      seq.names = "readout"
    } else if(slot == "Target"){
      if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
      check_transcripts(object,transcript_ids)
      data.use = data.use[transcript_ids]
      seq.names = names(data.use)
    }
    fileConn = file(filepath,"w")
    for(i in 1:length(data.use)){
      sequences = data.use[[i]][["Sequence"]]
      for(j in 1:length(sequences)){
        writeLines(paste(">",as.character(sequences[j]),";",seq.names,"-",j,sep = ""), fileConn)
        writeLines(as.character(sequences[j]), fileConn)
      }
    }
    close(fileConn)
  }
}

#' Create BLAST Object
#'
#' Create BLAST object from a FASTA file or load pre-computed BLAST database. See \code{\link[rBLAST]{makeblastdb}} and \code{\link[rBLAST]{blast}} for detail.
#' 
#' Pre-computed human BLAST database object was built using \href{https://www.gencodegenes.org/human/}{GENCODE} v31 human transcript sequences.
#' Pre-computed mouse BLAST database object was built using \href{https://www.gencodegenes.org/mouse/}{GENCODE} M22 mouse transcript sequences.
#' Pre-computed repetitive masked sequence BLAST database object was built from \href{https://genome.ucsc.edu/cgi-bin/hgTables}{UCSC}. 
#' 
#' @param file FASTA file path and name.
#' @param specie Either \emph{human} or \emph{mouse}. Load pre-computed blast database object.
#' @param masked_repeats Whether to load pre-computed specie specific repetitive sequences. Default is FALSE.
#' @return Returns a BLAST database object.
#' @importFrom rBLAST makeblastdb blast
#' @export
#' @examples \donttest{
#' # generate customized BLAST database 
#' bl.used = CreateBLASTObject("customized.fa")
#' }
CreateBLASTObject <- function(file = NULL, specie = NULL, masked_repeats = FALSE){
  if(!is.null(specie)){
    if(specie == "human"){
      if(masked_repeats){
        bl.obj <- blast(db = system.file("extdata", "hg38_rpmk.fa", package = "FISHprobe"), type = "blastn")
      } else {
        bl.obj <- blast(db = system.file("extdata", "gencode.v32.transcripts.fa", package = "FISHprobe"), type = "blastn")
      }
    } else if(specie == "mouse"){
      if(masked_repeats){
        bl.obj <- blast(db = system.file("extdata", "mm10_rpmk.fa", package = "FISHprobe"), type = "blastn")
      } else {
        bl.obj <- blast(db = system.file("extdata", "gencode.vM23.transcripts.fa", package = "FISHprobe"), type = "blastn")
      }
    } else {
      stop("specie must be either human or mouse.", call. = F)
    }
  } else{
    makeblastdb(file)
    bl.obj <- blast(db=file, type = "blastn")
  }
  return(bl.obj)
}

#' Run BLAST on Probes
#'
#' Run BLAST on probe sets against a database. Please install NCBI \href{https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download}{BLAST} first.
#' Note the default blastn is performed on both strands, for short sequences (<50nt), with E value cut-off of 1. 
#' Parallelization is implemented by default using \code{\link[pbmcapply]{pbmclapply}} with process bars.
#' However when future_parallel is set to true, \code{\link[future]{future}} will be used to spead up, and the process bar will not be shown.
#' 
#' @param object FISHprobe object.
#' @param sequences Sequences to BLAST. Note that return_blast will be set to TRUE for this mode.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param gene_name Gene or transcript ID of the input sequences. Default is NULL.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param blast_db BLAST database object. See \code{\link{CreateBLASTObject}} for detail. Default is using precomputed human or mouse BLAST databases.
#' \itemize{
#'   \item pre-computed human BLAST database object using \href{https://www.gencodegenes.org/human/}{GENCODE} v32 human transcript sequences.
#'   \item pre-computed mouse BLAST database object using \href{https://www.gencodegenes.org/mouse/}{GENCODE} M23 mouse transcript sequences.
#' } 
#' @param repeat_mask Whether to use pre-computed specie specific mask. Default is FALSE. See \code{\link{CreateBLASTObject}} for detail.
#' @param BLAST_args Arguments passed to BLAST in command-line style. Default is \code{-task blastn-short -evalue 1},
#' BLAST for sequences shorter than 50 bases and E value cut-off of 10. See \code{\link[rBLAST]{blast}} for detail. 
#' @param future_parallel Parallelization with \code{\link[future]{future}} R package. Default is FALSE.
#' @param future_plan Parallelization \code{\link[future]{plan}} used. Default is multiprocess.
#' @param return_blast Whether to return the BLAST result instead. Default is FALSE.
#' @return Returns a FISHprobe object or BLAST results depends on return_blast choice.
#' @importFrom pbmcapply pbmclapply 
#' @importFrom Biostrings DNAStringSet
#' @import rBLAST
#' @import future
#' @import future.apply
#' @export
#' @examples \donttest{
#' object = BLASTProbes(object)
#' 
#' # use customized BLAST database 
#' bl.used = CreateBLASTObject("customized.fa")
#' object = BLASTProbes(object, blast_db = bl.used)
#' 
#' # return BLAST results instead of FISHprobe object
#' bl.results = BLASTProbes(object, blast_db = bl.used, return_blast = TRUE)
#' 
#' # BLAST on user defined sequences
#' bl.results = BLASTProbes(sequences = sequences.use, gene_name = "gene", blast_db = bl.used)
#' }
BLASTProbes <- function(object = NULL, sequences = NULL, gene_name = NULL, transcript_ids = NULL, slot = "Target", blast_db = NULL, repeat_mask = FALSE, BLAST_args = "-task blastn-short -evalue 10", future_parallel = FALSE, future_plan = "multiprocess", return_blast = FALSE){
  if(is.null(object) & is.null(sequences)){
    stop("Either object or sequences must be provided", call. = FALSE)
  } else if(is.null(object)){
    if(is.null(blast_db)) stop("blast_db cannot be empty when running customized sequence mode", call. = FALSE)
    bl.data = list(sequences = data.frame(Sequence = sequences,stringsAsFactors = FALSE))
    if(!is.null(gene_name)) names(bl.data) = gene_name
    return_blast = TRUE
  } else {
    if(slot == "Target"){
      if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
      check_transcripts(object,transcript_ids)
      bl.data = object@probes[[slot]][transcript_ids]
    } else if(slot == "Readout"){
      bl.data = list(object@probes[[slot]])
    }
  }
  if(is.null(blast_db)){
    bl.use = CreateBLASTObject(specie = object@specie, masked_repeats = repeat_mask)
  } else {
    bl.use = blast_db
  }
  time_start = Sys.time()
  if(future_parallel){
    if(length(bl.data) == 1){
      plan(future_plan)
      bl.res = list(blast_homolog(bl.data[[1]], blast_db = bl.use, BLAST_args = BLAST_args, parallel = future_parallel))
    } else {
      plan(list(sequential, future_plan))
      bl.res = future_lapply(bl.data, blast_homolog ,blast_db = bl.use, BLAST_args = BLAST_args, parallel = future_parallel)
    }
  } else {
    bl.res = lapply(bl.data, blast_homolog ,blast_db = bl.use, BLAST_args = BLAST_args)
  }
  time_end = Sys.time()
  time_diff = time_end - time_start
  message("Time Elapsed ", round(as.numeric(time_diff),2), " ", attributes(time_diff)[["units"]])
  if(slot == "Target") names(bl.res) = transcript_ids
  if(return_blast){
    return(bl.res)
  } else {
    if(slot == "Readout"){
      object@blast.data[[slot]] = bl.res[[1]]
    } else if(slot == "Target"){
      object@blast.data[[slot]] = bl.res
    }
    object@param.log[[paste("blastProbes",slot,sep = "-")]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                                                   , time = Sys.time()) 
    return(object)
  }
}

#' Filter BLAST Results
#'
#' Filter BLAST results to get off-target gene/transcript hits.
#' @param object FISHprobe object.
#' @param blast_result BLAST result. Note that return_blast will be set to TRUE for this mode.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param alignment_length BLAST alignment length cutoff. BLAST hits with alignment length longer than this value will be consider as off-target hits. Default is 11.
#' @param E_value_cutoff BLAST alignment E value cutoff. BLAST hits with E value bigger than this will not be removed. Default is None.
#' @param filter_by_transcript Whether to filter based on transcript level or gene level. Default is FALSE.
#' @param probe_mode Whether to activate probe mode to BLAST the probes against each other. Default is FALSE.
#' @param names_delim Choose which delimiter to use to filter the off-taget gene/transcript hits. Default is ; for example ACCTTGGCGGCCAAT;ENSMUST00000000001-1 where ; separates sequences and ID.
#' @param names_field Choose which field for gene/transcript ID when probe_mode is on. Eg. For ACCTTGGCGGCCAAT;ENSMUST00000000001-1. the names_field = 1 will give ACCTTGGCGGCCAAT as gene name to filter the off-target probe hits. Default is 1.
#' @param parse Whether to parse the BLAST results into genes/transcripts. Default is TRUE
#' @param return_blast Whether to return the BLAST result instead. Default is FALSE.
#' @return Returns a FISHprobe object or BLAST results depends on return_blast choice.
#' @importFrom pbmcapply pbmclapply 
#' @export
#' @examples \donttest{
#' object = FilterBLAST(object, alignment_length = 15, E_value_cutoff = 1)
#' 
#' # use customized BLAST results after running BLASTProbes
#' blast_res = BLASTProbes(sequences = sequences.use, gene_name = "gene", blast_db = bl.used)
#' blast_filtered = FilterBLAST(object, blast_results = blast_res, alignment_length = 15)
#' }
FilterBLAST <- function(object = NULL, blast_result = NULL, slot = "Target", transcript_ids = NULL, alignment_length = 11, E_value_cutoff = NULL, filter_by_transcript = FALSE, probe_mode = FALSE,
                        names_delim = "-", names_field = 1, parse = TRUE, return_blast = FALSE){
  if(is.null(object) & is.null(blast_result)){
    stop("Either object or blast_results must be provided", call. = FALSE)
  } else if(probe_mode){
    if(!is.null(object)){
      blast_results = object@blast.data[[slot]]
    } else {
      blast_results = blast_result
    }
    if(!is.null(E_value_cutoff)){
      names_blast = names(blast_results)
      blast_results = lapply(1:length(blast_results),function(i){
        lapply(blast_results[[i]],function(j){
          j[which(j[["E"]] <= E_value_cutoff),]
        })
      })
      names(blast_results) = names_blast
    }
    if(!is.null(transcript_ids)) blast_results = blast_results[transcript_ids]
    blast.filtered = pbmclapply(1:length(blast_results),function(i){
      blast_filter(blast_results[[i]],names(blast_results[[i]]),names(blast_results[[i]]),data.frame(Sequence=names(blast_results[[i]]),stringsAsFactors = FALSE),alignment_length,filter_by_transcript,probe_mode = TRUE,names_delim,names_field,parse)
    })
    names(blast.filtered) = names(blast_results)
  } else {
    if(slot == "Target"){
      if(is.null(transcript_ids)){
        transripts.use = names(object@blast.data[[slot]])
      } else {
        transripts.use = transcript_ids 
      }
      check_transcripts(object,transcript_ids)
      blast_results = object@blast.data[[slot]][transripts.use]
    } else if(slot == "Readout"){
      blast_results = object@blast.data[[slot]]
    }
    if(!is.null(E_value_cutoff)){
      if(slot == "Target"){
        blast_results = lapply(1:length(blast_results),function(i){
          lapply(blast_results[[i]],function(j){
            j[which(j[["E"]] <= E_value_cutoff),]
          })
        })
      } else if(slot == "Readout"){
        blast_results = lapply(blast_results,function(j){
          j[which(j[["E"]] <= E_value_cutoff),]
        })
      }
    }
    if(slot == "Target"){
      genes.use = object@sequences[["GeneName"]][which(object@sequences[["ensembl_transcript_id"]] %in% transripts.use)]
      blast.filtered = pbmclapply(1:length(transripts.use),function(i){
        blast_filter(blast_results[[i]],transripts.use[i],genes.use[i],object@probes[[slot]][[i]],alignment_length,filter_by_transcript,parse = parse)
      })
      names(blast.filtered) = transripts.use
    } else if(slot == "Readout"){
      blast.filtered = blast_filter(blast_results,"","",object@probes[[slot]],alignment_length,filter_by_transcript,parse = parse)
    }
  }
  if(return_blast){
    return(blast.filtered)
  } else {
    object@off.target.expr[[slot]] = blast.filtered
    object@param.log[[paste("FilterBLAST",slot,sep = "-")]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                                                   , time = Sys.time()) 
    return(object)
  }
}

#' Calculate Probe Off-Target Expressions
#'
#' Calculate off-target expression for each probe.
#' 
#' If both gene names and transcript IDs are provided to exclude, the gene names will override.
#' 
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param tissue Tissue of interest. Default is the mean of all available tissues. Check available tissues with \code{\link{ListTissues}}.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param transcripts_to_exclude Transcript IDs to exclude in the off-target BLAST filtering. Default is NULL.
#' @param genes_to_exclude Gene names to exclude in the off-target BLAST filtering. Default is NULL.
#' @return Returns a FISHprobe object.
#' @importFrom pbmcapply pbmclapply 
#' @export 
#' @examples \donttest{
#' object = OffTargetExpr(object)
#' }
OffTargetExpr <- function(object, slot = "Target", tissue = "All", transcript_ids = NULL, transcripts_to_exclude = NULL, genes_to_exclude = NULL){
  if(is.null(transcript_ids)) transcript_ids = Transcripts(object)
  check_transcripts(object,transcript_ids)
  if(object@specie == "human"){
    expr.use = FISHprobe::gtex_mean_human
    gene.table.use = FISHprobe::human_ensemble
    trans.table.use = FISHprobe::transcripts_table_human
  } else if(object@specie == "mouse"){
    expr.use = FISHprobe::gtex_mean_mouse
    gene.table.use = FISHprobe::mouse_ensemble
    trans.table.use = FISHprobe::transcripts_table_mouse
  }
  if(!"All" %in% tissue){
    if(any(! tissue %in% ListTissues(object))) stop(paste("Cannot find tissue", paste(tissue[which(!tissue %in% ListTissues(object))])),call. = FALSE)
  }
  if(slot == "Target"){
    data.use = object@off.target.expr[[slot]][transcript_ids]
    data.use = pbmclapply(1:length(data.use),function(i){
      get_expression(data.use[[i]],expr.use,tissue,transcripts_to_exclude,genes_to_exclude,gene.table.use,trans.table.use)
    })
    names(data.use) = transcript_ids
    object@probes[[slot]][transcript_ids] = lapply(1:length(object@probes[[slot]][transcript_ids]), function(i){
      object@probes[[slot]][transcript_ids][[i]] = object@probes[[slot]][transcript_ids][[i]][,!colnames(object@probes[[slot]][transcript_ids][[i]]) %in% c("BLAST_Hits","BLAST_Hits_Transcripts", "BLAST_Hits_Gene", "Off_Targets_Expr")]
      data.temp = cbind(object@probes[[slot]][transcript_ids][[i]],data.use[[i]][,c("BLAST_Hits","BLAST_Hits_Transcripts", "BLAST_Hits_Gene", "Off_Targets_Expr")])
    })
    names(object@probes[[slot]][transcript_ids]) = transcript_ids
  } else if(slot == "Readout"){
    data.use = object@off.target.expr[[slot]]
    data.use = get_expression(data.use,expr.use,tissue,transcripts_to_exclude,genes_to_exclude,gene.table.use,trans.table.use)
    data.use[["BLAST_Hits"]] = unlist(data.use[["BLAST_Hits"]])
    object@probes[[slot]] = object@probes[[slot]][,!colnames(object@probes[[slot]]) %in% c("BLAST_Hits","BLAST_Hits_Transcripts", "BLAST_Hits_Gene", "Off_Targets_Expr")]
    object@probes[[slot]] = cbind.data.frame(object@probes[[slot]],data.use[,c("BLAST_Hits","BLAST_Hits_Transcripts", "BLAST_Hits_Gene", "Off_Targets_Expr")],stringsAsFactors =FALSE)
  }
  object@param.log[[paste("OffTarget",slot,sep = "-")]] = list(call = list(transcript_ids = transcript_ids,  tissue = tissue, transcripts_to_exclude = transcripts_to_exclude, genes_to_exclude = genes_to_exclude)
                                                                  , time = Sys.time())
  return(object)
}

#' Select Probes
#'
#' Select probes based on off-target expression and probe GC contents.
#' @param object FISHprobe object.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param TPM_cutoff Off-target expression cutoff in [logTPM]. Default is 15.
#' @param n_space Spaces between probes. Negative value results in overlapping probes. Default is 0 (non-overlapping).
#' @param order_by_feature Whether to select probes based on a numerical feature in \code{\link{ListFeatures}}. Default is none.
#' Examples including GC content ordering (order probes by the distance to feature_target value). 
#' @param feature_target Targeting numerical values to order the probes before selection (order probes by the distance to this value). Default is none.
#' @param iterative Whether to use iterative mode. Default is TRUE.
#' @param recalculate Whether to reselect the probes from stash.probes. Default is FALSE.
#' @return Returns a FISHprobe object.
#' @importFrom pbapply pblapply
#' @export
#' @examples \donttest{
#' object = SelectProbes(object, TPM_cutoff = 10, n_space = -5, order_by_feature = "GC",
#'                       feature_target = 0.50, iterative = TRUE)
#' # order by GC high to low
#' object = SelectProbes(object, TPM_cutoff = 10, n_space = -5, order_by_feature = "GC",
#'                       feature_target = 0.50, iterative = TRUE)
#' # order by GC low to high
#' object = SelectProbes(object, TPM_cutoff = 10, n_space = -5, order_by_feature = "GC",
#'                       feature_target = 0.50, iterative = TRUE)
#' }
SelectProbes <- function(object, transcript_ids = NULL, TPM_cutoff = 15, n_space = 0, order_by_feature = NULL, feature_target = NULL, iterative = TRUE, recalculate = FALSE){
  slot = "Target"
  if(is.null(transcript_ids)) transcript_ids = names(object@probes[["Target"]])
  check_transcripts(object,transcript_ids)
  if(!is.null(order_by_feature)){
    if(is.null(feature_target)) stop("Please provide a numerical value for feature_target.", call. = FALSE)
    if(length(order_by_feature) > 1 | length(feature_target) > 1) stop("Only one input feature  is  allowed. (order_by_feature and feature_target)", call. = FALSE)
    if(check_any_feature(object, slot, order_by_feature)) stop("Feature not found.", call. = FALSE)
    if(!is.numeric(feature_target)) stop("feature_target must be numerical.", call. = FALSE)
  }
  if(recalculate) object@probes[[slot]][transcript_ids] = object@stash.probes[[slot]][transcript_ids]
  data.to.select = object@probes[[slot]][transcript_ids]
  if(!iterative){
    selected.probes = pblapply(1:length(data.to.select),function(i){
      seq_select(data.to.select[[i]],TPM_cutoff,n_space,order_by_feature,feature_target)
    })
  } else {
    selected.probes = pblapply(1:length(data.to.select),function(i){
      iter_select(data.to.select[[i]],TPM_cutoff,n_space,order_by_feature,feature_target)
    })
  }
  names(selected.probes) = transcript_ids
  object@probes[[slot]][transcript_ids] = selected.probes
  object@param.log[[paste("Select",slot,sep = "-")]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                                            , time = Sys.time()) 
  return(object)
}

