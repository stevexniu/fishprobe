#' @include utils.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class Definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The FISHprobe S4 Class
#'
#' The FISHprobe object stores information associated with the Target and Readout probes, including probe sequences, BLAST results and off-target expression etc. 
#' The slot information is listed as follow:
#'
#' @slot genes Gene names of interest.
#' @slot transcripts List of ENSEMBL transcript ids associated with genes slot.
#' @slot sequences Gene sequences associated with genes slot.
#' @slot sequence.type Type of gene sequences in sequences slot.
#' @slot probes List contains Target and Readout probes.
#' @slot specie Specie either \emph{mouse} or \emph{human}.
#' @slot expr.data List of expression data associated with genes slot.
#' @slot off.target.expr List contains off-target expression data.
#' @slot blast.data List contains BLAST results.
#' @slot filter.table List contains probe filtering conditions.
#' @slot stash.probes List contains stashed probe sets.
#' @slot codebook List of dataframes contains the binary codebook, indices, conversion table and readout probe sequences used.
#' \itemize{
#'    \item The $Binary codebook is a data frame with shape: N of total encoded genes by 1+M of total hybridization rounds.
#'    \item The $Indice table is a data frame with shape: N of total encoded genes by 1+K of total readout probe used per FISH probe construct.
#'    \item The $Conversion table is a data frame with shape: N of total encoded genes by 3 (Channel of used, Index of readout probes, Bit assigned to that channel).
#'    \item The $Sequence is a data frame with shape: N of total encoded genes by 1+K of total readout probe used per FISH probe construct.
#'    \item The $Assignment is a data frame with final pairs of readout flanks and their assignment to each of the transcript or gene.
#'    \item All the above codebook slots must have the first column reserved to fluorescent Channel. 
#'  }
#' @slot final.probes Data frame contains the final readout/target design of FISH probes.
#' @slot param.log List of command parameter used.
#' @name FISHprobe-class
#' @rdname FISHprobe-class
#' @exportClass FISHprobe
#'
setClass(Class = "FISHprobe", 
         slots = c(
           probes = "list",
           genes = "character",
           transcripts = "list",
           sequences = "data.frame",
           sequence.type = "character",
           specie = "character",
           expr.data = "list",
           off.target.expr = "list",
           blast.data = "list",
           filter.table = "list",
           stash.probes = "list",
           codebook = "list",
           final.probes = "list",
           param.log = "list"
         ),
         prototype = list(probes = list(Target = NULL, Readout = NULL),
                          specie = "human",
                          off.target.expr = list(Target = NULL, Readout = NULL),
                          blast.data = list(Target = NULL, Readout = NULL),
                          filter.table = list(Target = NULL, Readout = NULL),
                          stash.probes = list(Target = NULL, Readout = NULL),
                          codebook = list(Binary = data.frame(NULL), Indice = data.frame(NULL), Conversion = data.frame(NULL), Sequence = data.frame(NULL), Assignment = data.frame(NULL)),
                          final.probes = data.frame(NULL)
         ),
         validity = function(object){
           if(object@specie == "human" | object@specie == "mouse"){
             return(TRUE)
           } else{
             return("FISHprobe object specie slot must be either human or mouse")
           }
         }
)

#' The DiffProbe S4 Class
#'
#' The DiffProbe object is a subclass of FISHprobe to design probes for different isoforms.
#' It has additional slots as follow:
#'
#' @slot DiffProbe Data frame of differential or shared sequences between different transcript isoforms.
#' @slot DiffProbe.type Type of either sequences either \strong{differential} or \strong{shared}.
#' @name DiffProbe-class
#' @rdname DiffProbe-class
#' @exportClass DiffProbe
#'
setClass(Class = "DiffProbe",contains = "FISHprobe",
         slots = c(
           DiffProbe = "data.frame",
           DiffProbe.type = "character"
           ),
         prototype = list(DiffProbe.type = "None"),
         validity = function(object){
           if(all(object@DiffProbe.type %in% c("Differential","Shared","None"))){
             return(TRUE)
           } else{
             return("DiffProbe object DiffProbe.type slot must be either Differential, Shared or None.")
           }
         }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 Generic Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @importFrom Hmisc capitalize
setMethod(f = "show",signature = "FISHprobe",
          definition = function(object){
            if(validObject(object)){
              cat(Hmisc::capitalize(object@specie),"FISHProbes for:", object@genes, "\n")
              cat("Transcript IDs:", names(object@probes[["Target"]]), "\n")
              cat("Number of Targeting Probes:", unlist(lapply(object@probes[["Target"]],nrow)), "\n")
            }
          }
)

#' @importFrom Hmisc capitalize
setMethod(f = "show",signature = "DiffProbe",
          definition = function(object){
            if(validObject(object)){
              cat(Hmisc::capitalize(object@specie),"DiffProbe Probes for:", object@genes, "\n")
              cat("Transcript IDs:", names(object@probes[["Target"]]), "\n")
              cat("Number of Targeting Probes:", unlist(lapply(object@probes[["Target"]],nrow)), "\n")
              cat("Number of DiffProbe Regions:", nrow(object@DiffProbe),"\n")
              cat("DiffProbe Types:", object@DiffProbe.type, "\n")
            }
          }
)

#' @export
#' @method length FISHprobe
"length.FISHprobe" <- function(x){
  length(Transcripts(x))
}

#' @export
#' @method [[ FISHprobe
"[[.FISHprobe" <- function(x, i, ...){
  if(missing(i)){
    return(x)
  } else if(length(i) == 1){
    if(i == "Readout"){
      return(x@probes[[i]])
    } else {
      if(is.numeric(i)){
        return(x@probes[["Target"]][[i]])
      } else if(is.character(i)){
        if(!i %in% names(x@probes[["Target"]])) stop(paste("Cannot find", paste(i[which(!i %in% names(x@probes[["Target"]]))],collapse = " ")), call. = FALSE)
        return(x@probes[["Target"]][[i]])
      }
    }
  } else if(length(i) > 1){
    if(is.numeric(i)){
      return(x@probes[["Target"]][i])
    } else if(is.character(i)){
      if(!all(i %in% names(x@probes[["Target"]]))) stop(paste("Cannot find", paste(i[which(!i %in% names(x@probes[["Target"]]))],collapse = " ")), call. = FALSE)
      return(x@probes[["Target"]][i])
    }
  }
}


#' @export
#' @method [[<- FISHprobe
"[[<-.FISHprobe" <- function(x, i, value){
  if(length(i) > 1) stop("Only one name allowed.", call. = FALSE)
  if(i == "Readout"){
    if(!is.data.frame(value)) stop("Probe data must be a data.frame.", call. = FALSE)
    x@probes[["Readout"]] <- value
  } else {
    if(!i %in% names(x@probes[["Target"]])) stop(paste("Cannot find", i), call. = FALSE)
    if(!is.data.frame(value)) stop("Probe data must be a data.frame.", call. = FALSE)
    x@probes[["Target"]][[i]] <- value
  }
  return(x)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Object Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create FISH Probe Object
#'
#' Creat FISH probe object.
#' @param gene_names Genes of interest. Mouse gene symbols can be found in \code{\link{mousegenes}}. Human gene symbols can be found in \code{\link{humangenes}}.
#' @param specie_type Either \emph{human} or \emph{mouse}. Default is \emph{human}.
#' @param preload_readout Logical value whether to load pre-computed readout probes and codebook. Default is FALSE.
#' @param ensemble_transcript_id Ensemble transcript IDs used. Default is None.
#' @param ensemble_gene_id Ensemble gene IDs used. Default is None.
#' @param customized_name Customized gene name. Default is None.
#' @param customized_data Customized sequences used. Default is None.
#' @importFrom Hmisc capitalize 
#' @return Returns a \code{\link{FISHprobe-class}} object.
#' @export
#' @examples \donttest{
#' object = CreateProbeObject("Egfr", specie_type = "mouse")
#' object = CreateProbeObject(c("EGFR","ACTB"), specie_type = "human")
#' object = CreateProbeObject(c("EGFR","ACTB"), specie_type = "human")
#' }
CreateProbeObject <- function(gene_names = "", specie_type = "human", preload_readout = FALSE, ensemble_transcript_id = NULL, ensemble_gene_id = NULL, customized_name = NULL, customized_data = NULL){
  if(!is.null(customized_name) & !is.null(customized_data)){
    message("Customized Sequences Are Used.")
    trans_data = list(customized_name)
    names(trans_data) = customized_name
    data_use = data.frame(sequences = customized_data, GeneName = customized_name, ensembl_transcript_id = customized_name, stringsAsFactors = FALSE)
    object = new(Class = "FISHprobe", genes = customized_name, transcripts = trans_data, specie = specie_type, sequences = data_use) 
  } else {
    if (specie_type == "human"){
      if(!is.null(ensemble_transcript_id)){
        if(any(!ensemble_transcript_id %in% FISHprobe::transcripts_table_human[,1])) stop(paste("Cannot find human ensemble transcript id", paste(ensemble_transcript_id[which(!ensemble_transcript_id %in% FISHprobe::transcripts_table_human[,1])],collapse = " ")),call. = FALSE)
        gene_id_new = FISHprobe::transcripts_table_human[,2][match(ensemble_transcript_id, FISHprobe::transcripts_table_human[,1])] 
        gene_name_new = FISHprobe::human_ensemble[,2][match(gene_id_new, FISHprobe::human_ensemble[,1])] 
        trans_data = data.frame(gene = gene_name_new, trans = ensemble_transcript_id, stringsAsFactors = FALSE)
        trans_data = split(trans_data[["trans"]], f = trans_data[["gene"]])
        object = new(Class = "FISHprobe", genes = gene_name_new, specie = specie_type, transcripts = trans_data) 
      } else {
        if(!is.null(ensemble_gene_id)){
          if(any(!ensemble_gene_id %in% FISHprobe::transcripts_table_human[,2])) stop(paste("Cannot find human ensemble gene id", paste(ensemble_gene_id[which(!ensemble_gene_id %in% FISHprobe::transcripts_table_human[,2])],collapse = " ")),call. = FALSE)
          gene_name_new = FISHprobe::human_ensemble[,2][match(ensemble_gene_id, FISHprobe::human_ensemble[,1])] 
          print(gene_name_new)
          object = new(Class = "FISHprobe", genes = gene_name_new, specie = specie_type) 
        } else {
          if(any(!gene_names %in% FISHprobe::humangenes[,1])) stop(paste("Cannot find human gene", paste(gene_names[which(!gene_names %in% FISHprobe::humangenes[,1])],collapse = " ")),call. = FALSE)
          object = new(Class = "FISHprobe", genes = gene_names, specie = specie_type) 
        }
      }
    } else if (specie_type == "mouse"){
      if(!is.null(ensemble_transcript_id)){
        if(any(!ensemble_transcript_id %in% FISHprobe::transcripts_table_mouse[,1])) stop(paste("Cannot find mouse ensemble transcript id", paste(ensemble_transcript_id[which(!ensemble_transcript_id %in% FISHprobe::transcripts_table_mouse[,1])],collapse = " ")),call. = FALSE)
        gene_id_new = FISHprobe::transcripts_table_mouse[,2][match(ensemble_transcript_id, FISHprobe::transcripts_table_mouse[,1])] 
        gene_name_new = FISHprobe::mouse_ensemble[,2][match(gene_id_new, FISHprobe::mouse_ensemble[,1])] 
        trans_data = data.frame(gene = gene_name_new, trans = ensemble_transcript_id, stringsAsFactors = FALSE)
        trans_data = split(trans_data[["trans"]], f = trans_data[["gene"]])
        object = new(Class = "FISHprobe", genes = gene_name_new, specie = specie_type, transcripts = trans_data) 
      } else {
        if(!is.null(ensemble_gene_id)){
          if(any(!ensemble_gene_id %in% FISHprobe::transcripts_table_mouse[,2])) stop(paste("Cannot find mouse ensemble gene id", paste(ensemble_gene_id[which(!ensemble_gene_id %in% FISHprobe::transcripts_table_mouse[,2])],collapse = " ")),call. = FALSE)
          gene_name_new = FISHprobe::mouse_ensemble[,2][match(ensemble_gene_id, FISHprobe::mouse_ensemble[,1])] 
          object = new(Class = "FISHprobe", genes = gene_name_new, specie = specie_type)
        } else {
          if(any(!gene_names %in% FISHprobe::mousegenes[,1])) stop(paste("Cannot find mouse gene", paste(gene_names[which(!gene_names %in% FISHprobe::mousegenes[,1])],collapse = " ")),call. = FALSE)
          object = new(Class = "FISHprobe", genes = gene_names, specie = specie_type) 
        }
      }
    }
  }
  if(preload_readout){
    object = SetSlotData(object = object, data = FISHprobe::readout_v1, slotname = "probes", subname = "Readout") 
    object = SetSlotData(object = object, data = FISHprobe::codebook_v1, slotname = "codebook") 
  } 
  object@param.log[["Create"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe())))
                                      , time = Sys.time())
  return(object)
}

#' Get Slot Data
#'
#' Get data from a certain slot.
#' @param object FISHprobe object.
#' @param slotname Slot name.
#' @param subname Sub-slot name if there is any.
#' @return Returns the slot data.
#' @export
#' @examples \donttest{
#' # get the target probe data
#' GetSlotData(object, slotname = "probes", subname = "Target")
#' }
GetSlotData <- function(object, slotname, subname = NULL){
  data.temp = try(slot(object, name = slotname), silent = TRUE)
  if(class(data.temp) == "try-error") stop("slot_name not found", call. = FALSE) 
  if(!is.null(subname)){
    data.temp = try(data.temp[[subname]],  silent = TRUE)
    if(class(data.temp) == "try-error" | is.null(data.temp)) stop("subname not found", call. = FALSE) 
  }
  return(data.temp)
}

#' Get Transcript IDs
#'
#' Get transcript IDs.
#' @param object FISHprobe object.
#' @return Returns transcript IDs in the object.
#' @export
#' @examples \donttest{
#' Transcripts(object)
#' }
Transcripts <- function(object){
  return(names(object@probes[["Target"]]))
}

#' Set Slot Data
#'
#' Set data for a certain slot.
#' @param object FISHprobe object.
#' @param data Data to set with.
#' @param slotname Slot name.
#' @param subname Sub-slot name if there is any.
#' @return Returns a FISHprobe object.
#' @export
#' @examples \donttest{
#' object = SetSlotData(object, data = data.use, slotname = "probes", subname = "Target")
#' }
SetSlotData <- function(object, data, slotname, subname = NULL){
  data.temp = try(slot(object, name = slotname), silent = TRUE)
  if(class(data.temp) == "try-error") stop("slot_name not found", call. = FALSE) 
  if(!is.null(subname)){
    data.sub = try(data.temp[[subname]],  silent = TRUE)
    if(class(data.sub) == "try-error") stop("subname not found", call. = FALSE) 
    data.temp[[subname]] = data
    slot(object, name = slotname) <- data.temp
  } else {
    slot(object, name = slotname) <- data
  }
  return(object)
}

#' Select Transcripts
#'
#' Select transcripts of genes of interest by tissue specific expressions.
#' @param object FISHprobe object.
#' @param genes Genes of interest. Default is all the genes in the genes slot.
#' @param tissue Tissue of interest. Default is the mean of all available tissues. Check available tissues with \code{\link{ListTissues}}.
#' @param select_best Select the most expressed transcript id. Default is TRUE.
#' @param TPM_min Minimum logTPM expression level. 
#' @param TPM_max Maximum logTPM expression level.
#' @param percent_min Minimum percentage of expression among tissue replicates.
#' @param percent_max Maximum percentage of expression among tissue replicates.
#' @return Returns a FISHprobe object with the selected transcripts for each gene stored in object@transcripts$Selected.
#' @export
#' @examples \donttest{
#' # Automatic selection of the most expressed transcript
#' object = SelectTranscript(object)
#' # Select by at least 10 logTPM expression in lung tissue
#' object = SelectTranscript(object, genes = "ACTB", tissue = "lung", TPM_min = 10)
#' }
SelectTranscript <- function(object, genes = NULL, tissue = "All", select_best = TRUE, TPM_min = -Inf, TPM_max = Inf, percent_min = -Inf, percent_max = Inf){
  if(is.null(genes)){
    genes = object@genes
  } else {
    if(any(!genes %in% object@genes)) stop(cat("Cannot find gene", genes[which(!genes %in% object@genes)]),call. = FALSE)
    genes = intersect(genes,object@genes)
  }
  if(!"All" %in% tissue){
    if(any(! tissue %in% ListTissues(object))) stop(cat("Cannot find tissue", tissue[which(!tissue %in% ListTissues(object))]),call. = FALSE)
  }
  object@transcripts[["Selected"]] = list()
  for(g in genes){
    data.temp = object@expr.data[[g]]
    data.temp = split(data.temp, data.temp[["Transcript"]])
    if("All" %in% tissue){
      data.temp = lapply(data.temp, function(x){
        temp = x[1,]
        temp[1,1] = "All"
        temp[1,3:4] = colMeans(x[,3:4])
        return(temp)
      })
    } else {
      data.temp = lapply(data.temp, function(x){
        temp = x[1,]
        temp[1,1] = "Selected"
        temp[1,3:4] = colMeans(x[which(x[["Tissue"]] %in% unique(tissue)),3:4])
        return(temp)
      })
    }
    data.temp = do.call(rbind,data.temp)
    data.temp[["Transcript"]] = rownames(data.temp)
    if(select_best){
      transcript.selected = data.temp[["Transcript"]][which.max(data.temp[,3])]
    } else {
      transcript.selected = data.temp[["Transcript"]][which((data.temp[,3] < TPM_max & data.temp[,3] > TPM_min) & (data.temp[,4] < percent_max & data.temp[,4] > percent_min))]
    }
    object@transcripts[["Selected"]][[g]] = transcript.selected
  }
  object@param.log[["SelectTranscript"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                                , time = Sys.time())
  return(object)
}

#' Get Transcript Sequences
#'
#' Get transcript sequences.
#' I noted that the Ensembl server can be unstable for some reasons, which will result in returning a null gene sequence.
#' A warning will ring and manual curation may be needed in such cases.
#' @param object FISHprobe object.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param seqtype Type of sequences. See \code{\link[biomaRt]{getSequence}} for all possible sequence types. Default is cdna sequences.
#' @return Returns a FISHprobe object.
#' @importFrom pbapply pblapply
#' @export
#' @examples \donttest{
#' # Automatic selection of the most expressed transcript
#' object = GetSequence(object)
#' # Select by at least 10 logTPM expression in lung tissue
#' object = SelectTranscript(object, genes = "ACTB", tissue = "lung", TPM_min = 10)
#' }
GetSequence <- function(object, transcript_ids = NULL, seqtype ="cdna"){
  all.transcripts = object@transcripts
  all.transcripts[["Selected"]] = NULL
  all.transcripts = unlist(all.transcripts)
  if(is.null(transcript_ids)){
    transcript_ids = unlist(object@transcripts[["Selected"]])
    if(length(transcript_ids) == 0) transcript_ids = all.transcripts
  } else {
    if(any(!transcript_ids %in% all.transcripts)) stop(cat("Following Transcripts ID not found",transcript_ids[which(!transcript_ids %in% all.transcripts)]),call. = FALSE)
  }
  if(object@specie == "human") {
    seq.all = pblapply(transcript_ids,get_sequences,FISHprobe::human,FISHprobe::transcripts_table_human,FISHprobe::human_ensemble,"human",seqtype)
  } else if (object@specie == "mouse"){
    seq.all = pblapply(transcript_ids,get_sequences,FISHprobe::mouse,FISHprobe::transcripts_table_mouse,FISHprobe::mouse_ensemble,"mouse",seqtype)
  }
  seq.all = do.call(rbind.data.frame, c(seq.all, stringsAsFactors = FALSE))
  object@sequences = seq.all
  object@sequence.type = seqtype
  object@param.log[["GetSeq"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                   , time = Sys.time()) 
  return(object)
}

#' Subset Probe by Sequence
#'
#' Subset FISHprobe object with sequences
#' 
#' The order of subseting is defined as transcript_ids > gene_anmes > sequences, if all of them are provided.
#' 
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param gene_names Gene names to subset. Default is NULL.
#' @param sequences Probe targeting sequences (not the probe sequences) to subset.
#' @param indices Probe indices to subset.
#' @param reverse_select Whether to perform reverse/removal selection. Default is FALSE.
#' @return Returns a FISHprobe object.
#' @export
#' @examples \donttest{
#' object = SubsetProbes(object, transcript_ids = "ENST00000216037")
#' }
SubsetProbes <- function(object, slot = "Target", transcript_ids = NULL, gene_names = NULL, sequences = NULL, indices = NULL, reverse_select = FALSE){
  if(slot == "Readout"){
    if(is.null(sequences) & is.null(indices)) stop("Sequences or indices needed to subset Readout probes.", call. = FALSE)
    data.use = object@probes[[slot]]
    if(reverse_select){
      if(is.null(indices)){
        idx = which(!data.use[["Sequence"]] %in% sequences)
      } else {
        idx = setdiff(1:nrow(data.use),indices)
      }
    } else {
      if(is.null(indices)){
        idx = which(data.use[["Sequence"]] %in% sequences)
      } else {
        idx = indices
      }
    }
    object@probes[[slot]] = object@probes[[slot]][idx,,drop = FALSE]
    object@blast.data[[slot]] = object@blast.data[[slot]][idx,,drop = FALSE]
    object@off.target.expr[[slot]] = object@off.target.expr[[slot]][idx,,drop = FALSE]
    object@stash.probes[[slot]] = object@stash.probes[[slot]][idx,,drop = FALSE]
    object@filter.table[[slot]] = object@filter.table[[slot]][idx,,drop = FALSE]
  } else if(slot == "Target"){
    if(!is.null(transcript_ids)){
      if(!all(transcript_ids %in% names(object@probes[[slot]]))) stop(paste("Cannot find transcript", paste(transcript_ids[which(!transcript_ids %in% names(object@probes[[slot]]))],collapse = " ")),call. = FALSE)
      if(reverse_select) {
        trans_match = which(!Transcripts(object) %in% transcript_ids)
        transcript_ids = Transcripts(object)[trans_match]
      } else {
        trans_match = which(Transcripts(object) %in% transcript_ids)
        transcript_ids = Transcripts(object)[trans_match]
      }
      object@genes = unique(object@sequences[["GeneName"]][which(object@sequences[["ensembl_transcript_id"]] %in% transcript_ids)])
      object@probes[[slot]] = object@probes[[slot]][transcript_ids]
      object@blast.data[[slot]] = object@blast.data[[slot]][transcript_ids]
      object@off.target.expr[[slot]] = object@off.target.expr[[slot]][transcript_ids]
      object@stash.probes[[slot]] = object@stash.probes[[slot]][transcript_ids]
      object@filter.table[[slot]] = object@filter.table[[slot]][transcript_ids]
    } else if(!is.null(gene_names)){
      if(!all(gene_names %in% object@genes)) stop(paste("Cannot find gene", paste(gene_names[which(!gene_names %in% object@genes)],collapse = " ")),call. = FALSE)
      if(reverse_select) {
        gene_match = which(!object@sequences[["GeneName"]] %in% gene_names)
        object@genes = setdiff(object@genes, gene_names)
        transcripts_all = object@sequences[gene_match,"ensembl_transcript_id"]
        transcript_ids = intersect(Transcripts(object), transcripts_all)
      } else {
        gene_match = which(object@sequences[["GeneName"]] %in% gene_names)
        object@genes = intersect(object@genes, gene_names)
        transcripts_all = object@sequences[gene_match,"ensembl_transcript_id"]
        transcript_ids = intersect(Transcripts(object), transcripts_all)
      }
      object@probes[[slot]] = object@probes[[slot]][transcript_ids]
      object@blast.data[[slot]] = object@blast.data[[slot]][transcript_ids]
      object@off.target.expr[[slot]] = object@off.target.expr[[slot]][transcript_ids]
      object@stash.probes[[slot]] = object@stash.probes[[slot]][transcript_ids]
      object@filter.table[[slot]] = object@filter.table[[slot]][transcript_ids]
    } 
    if(!is.null(sequences)){
      for(i in 1:length(object@probes[[slot]])){
        data.use = object@probes[[slot]][[i]]
        if(reverse_select){
          idx = which(!data.use[["Sequence"]] %in% sequences)
        } else {
          idx = which(data.use[["Sequence"]] %in% sequences)
        }
        object@probes[[slot]][[i]] = object@probes[[slot]][[i]][idx,,drop = FALSE]
        object@blast.data[[slot]][[i]] = object@blast.data[[slot]][[i]][idx]
        object@off.target.expr[[slot]][[i]] = object@off.target.expr[[slot]][[i]][idx,,drop = FALSE]
        object@stash.probes[[slot]][[i]] = object@stash.probes[[slot]][[i]][idx,,drop = FALSE]
        object@filter.table[[slot]][[i]] = object@filter.table[[slot]][[i]][idx,,drop = FALSE]
      }
    } else if(!is.null(indices)){
      for(i in 1:length(object@probes[[slot]])){
        data.use = object@probes[[slot]][[i]]
        if(reverse_select){
          idx = setdiff(1:nrow(data.use), indices)
        } else {
          idx = indices
        }
        object@probes[[slot]][[i]] = object@probes[[slot]][[i]][idx,,drop = FALSE]
        object@blast.data[[slot]][[i]] = object@blast.data[[slot]][[i]][idx]
        object@off.target.expr[[slot]][[i]] = object@off.target.expr[[slot]][[i]][idx,,drop = FALSE]
        object@stash.probes[[slot]][[i]] = object@stash.probes[[slot]][[i]][idx,,drop = FALSE]
        object@filter.table[[slot]][[i]] = object@filter.table[[slot]][[i]][idx,,drop = FALSE]
      }
    }
  }
  object@param.log[[paste("Subset",slot,sep = "-")]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                   , time = Sys.time()) 
  return(object)
}

#' Merge Probe Objects
#'
#' Merge a list of FISHprobe objects. Please do not merge object with empty probe set.
#' @param object_list A list of FISHprobe objects.
#' @return Returns a merged FISHprobe object.
#' @export
#' @examples \donttest{
#' object_merged = MergeProbes(list(object1, object2))
#' }
MergeProbes <- function(object_list){
  merged.objects = object_list[[1]]
  diffprobe = FALSE
  if(any(unlist(lapply(object_list, class)) != class(merged.objects))){
    warning("Merging Objects from Different Classes", call. = FALSE)
    merged.objects = as(merged.objects, Class = "DiffProbe")
    object_list = lapply(object_list, as, Class = "DiffProbe")
    diffprobe = TRUE
  } 
  merged.objects@sequence.type = rep(merged.objects@sequence.type, length(merged.objects@probes[["Target"]]))
  merged.objects@specie = merged.objects@specie
  merged.objects@param.log = list(object_1 = merged.objects@param.log)
  trans.select = merged.objects@transcripts[["Selected"]] 
  merged.objects@transcripts = merged.objects@transcripts[which(names(merged.objects@transcripts) != "Selected")]
  for(i in 2:length(object_list)){
    if(object_list[[i]]@specie != merged.objects@specie) stop("Cannot merge probes across different species", call. = FALSE)
    if(diffprobe){
      merged.objects@DiffProbe = rbind(merged.objects@DiffProbe, object_list[[i]]@DiffProbe)
      merged.objects@DiffProbe.type = c(merged.objects@DiffProbe.type, object_list[[i]]@DiffProbe.type)
    }
    merged.objects@genes = c(merged.objects@genes, object_list[[i]]@genes)
    merged.objects@transcripts = c(merged.objects@transcripts, object_list[[i]]@transcripts)
    trans.select = c(trans.select, object_list[[i]]@transcripts[["Selected"]])
    merged.objects@sequences = rbind(merged.objects@sequences, object_list[[i]]@sequences[which(!object_list[[i]]@sequences[["ensembl_transcript_id"]] %in% merged.objects@sequences[["ensembl_transcript_id"]]),])
    merged.objects@sequence.type = c(merged.objects@sequence.type, rep(object_list[[i]]@sequence.type, length(object_list[[i]]@probes[["Target"]])))
    merged.objects@probes[["Target"]] = c(merged.objects@probes[["Target"]], object_list[[i]]@probes[["Target"]]) 
    merged.objects@probes[["Readout"]] = rbind(merged.objects@probes[["Readout"]], object_list[[i]]@probes[["Readout"]]) 
    merged.objects@expr.data = c(merged.objects@expr.data, object_list[[i]]@expr.data)
    merged.objects@off.target.expr[["Target"]] = c(merged.objects@off.target.expr[["Target"]], object_list[[i]]@off.target.expr[["Target"]]) 
    merged.objects@off.target.expr[["Readout"]] = rbind(merged.objects@off.target.expr[["Readout"]], object_list[[i]]@off.target.expr[["Readout"]]) 
    merged.objects@codebook[["Binary"]] = rbind(merged.objects@codebook[["Binary"]], object_list[[i]]@codebook[["Binary"]])
    merged.objects@codebook[["Indice"]] = rbind(merged.objects@codebook[["Indice"]], object_list[[i]]@codebook[["Indice"]])
    merged.objects@codebook[["Conversion"]] = rbind(merged.objects@codebook[["Conversion"]], object_list[[i]]@codebook[["Conversion"]])
    merged.objects@codebook[["Sequence"]] = rbind(merged.objects@codebook[["Sequence"]], object_list[[i]]@codebook[["Sequence"]])
    merged.objects@codebook[["Assignment"]] = rbind(merged.objects@codebook[["Assignment"]], object_list[[i]]@codebook[["Assignment"]])
    merged.objects@param.log[[paste("object",i,sep = "_")]] = object_list[[i]]@param.log
    merged.objects@blast.data[["Target"]] = c(merged.objects@blast.data[["Target"]], object_list[[i]]@blast.data[["Target"]]) 
    merged.objects@blast.data[["Readout"]] = rbind(merged.objects@blast.data[["Readout"]], object_list[[i]]@blast.data[["Readout"]])
    merged.objects@filter.table[["Target"]] = c(merged.objects@filter.table[["Target"]], object_list[[i]]@filter.table[["Target"]])
    merged.objects@filter.table[["Readout"]] = rbind(merged.objects@filter.table[["Readout"]], object_list[[i]]@filter.table[["Readout"]])
    merged.objects@stash.probes[["Target"]] = c(merged.objects@stash.probes[["Target"]], object_list[[i]]@stash.probes[["Target"]])
    merged.objects@stash.probes[["Readout"]] = rbind(merged.objects@stash.probes[["Readout"]], object_list[[i]]@stash.probes[["Readout"]])
    merged.objects@final.probes = c(merged.objects@final.probes, object_list[[i]]@final.probes) 
  }
  merged.objects@transcripts[["Selected"]] = trans.select
  return(merged.objects)
}

#' List All Possible Features to Filter Probes
#'
#' List all possible features to filter the probes for each transcript.
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param pattern Feature patterns to match, examples include GC, Tm etc. Default is none.
#' @param use.stash Whether to use stash.probes data. Default is FALSE.
#' @return Returns a list of possible features to filter on.
#' @export
#' @examples \donttest{
#' ListFeatures(object)
#' }
ListFeatures <- function(object, slot = "Target", transcript_ids = NULL, pattern = "", use.stash = FALSE){
  if(use.stash){
    data.use = object@stash.probes[[slot]]
  } else {
    data.use = object@probes[[slot]]
  }
  if(slot == "Readout"){
    features = list(Readout = grep(pattern, colnames(data.use), value = TRUE))
  } else if(slot == "Target"){
    if(is.null(transcript_ids)) transcript_ids = names(data.use)
    check_transcripts(object, transcript_ids)
    features = lapply(data.use[transcript_ids],function(x){
      grep(pattern, colnames(x), value = TRUE)
    })
    names(features) = transcript_ids
  }
  return(features)
}

#' Reset Probe Data 
#' 
#' Reset probe data to previously stashed version.
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @return Returns a FISHprobe object.
#' @export
#' @examples \donttest{
#' ResetProbes(object)
#' }
ResetProbes <- function(object, slot = "Target", transcript_ids = NULL){
  if(slot == "Target") 
    if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
    object@probes[[slot]][transcript_ids] = object@stash.probes[[slot]][transcript_ids]
  if(slot == "Readout")
    object@probes[[slot]] = object@stash.probes[[slot]]
  object@param.log[[paste("Reset",slot,sep = "-")]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                   , time = Sys.time())
  return(object)
}
