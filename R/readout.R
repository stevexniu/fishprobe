#' @include utils.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Readout Probes Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Add Readout Probes
#'
#' Add readout probes.
#' @param object FISHprobe object.
#' @param readout_data Data frame with a set of readout probes.
#' @param codebook_data Codebook object with correct format. See \code{\link{FISHprobe-class}} for details.
#' @return Returns a FISHprobe object.
#' @export
#' @examples \donttest{
#' object = AddReadout(object, readout.probes)
#' }
AddReadout <- function(object, readout_data, codebook_data){
  object = SetSlotData(object, readout_data, "probes", subname = "Readout")
  object = SetSlotData(object, codebook_data, "codebook")
  object@param.log[["AddReadout"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-c(1,length(mget(names(formals()),sys.frame(sys.nframe()))))])
                                                        , additional = match.call(),time = Sys.time()) 
  return(object)
}

#' Concatenate Readout Probe Sequences
#'
#' Concatenate readout probe sequences based on their indice.
#' @param object FISHprobe object.
#' @param concat_list A list of indice of concatenated probes.
#' @param insert Inserted based between each probes. Default is none.
#' @return Returns a FISHprobe object.
#' @export
#' @examples \donttest{
#' object = ConcatReadout(object, concat_list = list(c(1,2),c(1,3),c(2,3)))
#' }
ConcatReadout <- function(object, concat_list, insert = ""){
  concat_probes = lapply(concat_list, function(x){
    concat_probes(object@probes[["Readout"]][["Sequence"]], x, insertion = insert)
  })
  object = AddReadout(object, unlist(concat_probes))
  object@param.log[["ConcatReadout"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                             , time = Sys.time()) 
  return(object)
}

#' Generate Readout Probes
#'
#' Generate readout probes/flanks based on codebook and readout probes.
#' @param object FISHprobe object.
#' @param scheme Scheme for concatenate readout probes. A list contains 2 set of readout probes to be combined. Default 1-2 and 3-4 position combination.
#' @param flank_5end Sequence to be added to the 5' end flank regions. Default is none.
#' @param flank_3end Sequence to be added to the 3' end flank regions. Default is none.
#' @param concat Sequence to be added in between readout probes. Default is none.
#' @return Returns a FISHprobe object with generated readout probes/flanks in the codebook$Assignment slot.
#' @export
#' @examples \donttest{
#' object = GenerateReadout(object, flank = "AAAAAAAA", concat = "T")
#' }
GenerateReadout <- function(object, scheme = list(1:2,3:4), flank_5end = "", flank_3end = "", concat = ""){
  if(nrow(object@codebook[["Sequence"]]) == 0) stop("No readout sequence detected in codebool slot.", call. = FALSE)
  if(length(scheme) != 2) stop("The scheme must have length of 2.", call. = FALSE)
  data.use = object@codebook[["Sequence"]][,-1]
  channel.use = object@codebook[["Sequence"]][,1]
  concat.list = lapply(scheme, function(x){
    apply(data.use, 1, function(d){
      paste(d[x], collapse = concat)
    })
  })
  concat.list[[1]] = paste(flank_5end, concat.list[[1]], sep = "")
  concat.list[[2]] = paste(concat.list[[2]], flank_3end, sep = "")
  concat.data = do.call(cbind, concat.list)
  colnames(concat.data) = c("Flank1", "Flank2")
  concat.data = cbind(Channel = channel.use, concat.data)
  concat.data = data.frame(concat.data, check.names = FALSE, stringsAsFactors = FALSE)
  object@codebook[["Assignment"]] = concat.data
  object@param.log[["GenerateReadout"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                             , time = Sys.time()) 
  return(object)
}

#' Assign Transcripts/Genes to Readout Probes
#'
#' Assign transcripts/genes to the readout probes/flanks generated.
#' @param object FISHprobe object.
#' @param channel_transcript N channel list contains transcript ids in each channel. 
#' @param indices N channel list contains specific indices for each transcript in each channel. 
#' Default is NULL, the transcripts in each channel will be ordered the same as in channel_transcript list.
#' @return Returns a FISHprobe object with updated codebook$Assignment slot.
#' @export
#' @examples \donttest{
#' channel_list = list("488" = Transcript1, "561" = Transcript2, "647" = Transcript3)
#' object = AssignReadout(object, channel_transcript = channel_list)
#' }
AssignReadout <- function(object, channel_transcript, indices = NULL){
  if(nrow(object@codebook[["Assignment"]]) == 0) stop("No probes detected in codebook$Assignment slot.", call. = FALSE)
  channel.names = unique(object@codebook[["Assignment"]][["Channel"]])
  if(!all(names(channel_transcript) %in% channel.names)) stop(paste("Cannot find channel", paste(names(channel_transcript)[which(!names(channel_transcript) %in% channel.names)],collapse = " ")), call. = FALSE)
  object@codebook[["Assignment"]][["Transcript_ID"]] = NA
  object@codebook[["Binary"]][["Transcript_ID"]] = NA
  object@codebook[["Indice"]][["Transcript_ID"]] = NA
  object@codebook[["Sequence"]][["Transcript_ID"]] = NA
  for(i in names(channel_transcript)){
    trans.used = channel_transcript[[i]]
    idx.use = which(object@codebook[["Assignment"]][["Channel"]] == i)
    if(is.null(indices)){
      trans.len = length(trans.used)
      object@codebook[["Assignment"]][idx.use,"Transcript_ID"][1:trans.len] = trans.used
      object@codebook[["Binary"]][idx.use,"Transcript_ID"][1:trans.len] = trans.used
      object@codebook[["Indice"]][idx.use,"Transcript_ID"][1:trans.len] = trans.used
      object@codebook[["Sequence"]][idx.use,"Transcript_ID"][1:trans.len] = trans.used
      } else {
        indices.use = indices[[i]]
        object@codebook[["Assignment"]][idx.use,"Transcript_ID"][indices.use] = trans.used
        object@codebook[["Binary"]][idx.use,"Transcript_ID"][indices.use] = trans.used
        object@codebook[["Indice"]][idx.use,"Transcript_ID"][indices.use] = trans.used
        object@codebook[["Sequence"]][idx.use,"Transcript_ID"][indices.use] = trans.used
      }
  }
  object@codebook[c("Assignment","Binary","Indice","Sequence")] = lapply(object@codebook[c("Assignment","Binary","Indice","Sequence")], function(x){
    cbind.data.frame(x, Gene = object@sequences$GeneName[match(x[["Transcript_ID"]], object@sequences$ensembl_transcript_id)])
  })
  object@param.log[["AssignReadout"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                               , time = Sys.time())
  return(object)
}

#' Generate Final Probe Set
#'
#' Generate the final probe set.
#' @param object FISHprobe object.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param concat Sequence to be added in between readout and target probes. Default is none.
#' @return Returns a FISHprobe object with updated final.probes slot.
#' @export
#' @examples \donttest{
#' object = GenerateFinalProbes(object, concat = "T")
#' }
GenerateFinalProbes <- function(object, transcript_ids = NULL, concat = ""){
  if(is.null(object@codebook[["Assignment"]][["Transcript_ID"]])) stop("No transcript assignment found in codebook$Assignment.", call. = FALSE)
  trans.use = object@codebook[["Assignment"]][["Transcript_ID"]][which(!is.na(object@codebook[["Assignment"]][["Transcript_ID"]]))]
  if(is.null(transcript_ids)) transcript_ids = Transcripts(object)
  if(!all(transcript_ids %in% trans.use)) stop(paste("Cannot find assignment for transcript", paste(transcript_ids[which(!transcript_ids %in% trans.use)],collapse = " ")), call. = FALSE)
  data.list = lapply(transcript_ids, function(x){
    readout.data = object@codebook[["Assignment"]][which(object@codebook[["Assignment"]][["Transcript_ID"]] == x),]
    flk1 = readout.data[["Flank1"]]
    flk2 = readout.data[["Flank2"]]
    channel.use = readout.data[["Channel"]]
    probe.data = object@probes[["Target"]][[x]]
    probe.names =  object@probes[["Target"]][[x]][["Sequence"]]
    probes.final = apply(probe.data, 1, function(x){
      paste(flk1,x[["Probe"]],flk2, sep = concat)
    })
    codebook.data = object@codebook[["Indice"]][which(object@codebook[["Assignment"]][["Transcript_ID"]] == x),2:5]
    gene.name = object@sequences[["GeneName"]][which(object@sequences[["ensembl_transcript_id"]] == x)]
    probes.all = data.frame(Channel = channel.use, Transcript_ID = x, GeneName = gene.name, Final_Probe = probes.final, Hybrid_Round = codebook.data,
                            stringsAsFactors = FALSE, row.names = probe.names, check.names = FALSE, check.rows = FALSE)
  })
  names(data.list) = transcript_ids
  object@final.probes = data.list
  object@param.log[["GenerateFinalProbes"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                             , time = Sys.time())
  return(object)
}
