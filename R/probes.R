#' @include utils.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Target Probes Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Extract Possible Probes
#'
#' Extract possible probes of a specific length.
#' @param object FISHprobe object.
#' @param probe_length Probe length. Default is 25 nt.
#' @param consecutive_repeats Find consecutive repeats of N same nucleotides. Default is 5 nt.
#' @return Returns a FISHprobe object with all possible probes and information stored in object@probes$Target.
#' @importFrom stringr str_count
#' @export
#' @examples \donttest{
#' # Generate all possible probes of 30nt long 
#' # and find all consecutive repeats with 3 same nucleotides.
#' object = ExtractProbes(object, probe_length = 30, consecutive_repeats = 3)
#' }
ExtractProbes <- function(object, probe_length = 25, consecutive_repeats = 5){
  object@probes[["Target"]] = lapply(object@sequences[,1], prob_extract, probe_len = probe_length, consecutive_repeats = consecutive_repeats)
  names(object@probes[["Target"]]) = object@sequences[["ensembl_transcript_id"]]
  if(!is.null(object@sequences[["SeqStarts"]])){
    for(i in 1:length(object@probes[["Target"]])){
      object@probes[["Target"]][[i]][["Start"]] = object@probes[["Target"]][[i]][["Start"]] + object@sequences[["SeqStarts"]][i] - 1
      object@probes[["Target"]][[i]][["Stop"]] = object@probes[["Target"]][[i]][["Stop"]] + object@sequences[["SeqStarts"]][i] - 1
    }
  }
  object@param.log[["ExtractProbes"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                       , time = Sys.time()) 
  return(object)
}

#' Generate Probe Sequences
#'
#' Generate reverse complementary probe sequences against the targeted sequences.
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @return Returns a FISHprobe object.
#' @export
#' @examples \donttest{
#' object = GenerateProbes(object)
#' }
GenerateProbes <- function(object, slot = "Target", transcript_ids = NULL){
  if(slot == "Target"){
    if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
    check_transcripts(object,transcript_ids)
    object@probes[[slot]][transcript_ids] = lapply(object@probes[[slot]][transcript_ids], find_rev_compl)
    object@stash.probes[[slot]][transcript_ids] = lapply(object@stash.probes[[slot]][transcript_ids], find_rev_compl)
  } else if(slot == "Readout"){
    object@probes[[slot]] = lapply(list(object@probes[[slot]]), find_rev_compl)[[1]]
    object@stash.probes[[slot]] = lapply(object@stash.probes[[slot]], find_rev_compl)[[1]]
  }
  object@param.log[[paste("GenerateProbes",slot,sep = "-")]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                                              , time = Sys.time()) 
  return(object)
}

#' Calculate GC Content
#'
#' Calculate GC content of probes.
#' 
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param colname_use Column name used to store GC. Default is GC.
#' @return Returns a FISHprobe object.
#' @export
#' @examples \donttest{
#' object = CalcGC(object)
#' }
CalcGC <- function(object, slot = "Target", transcript_ids = NULL, colname_use = "GC"){
  if(slot == "Target"){
    if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
    check_transcripts(object,transcript_ids)
    gc.data = lapply(object@probes[[slot]][transcript_ids],function(data){
      unlist(lapply(data[["Sequence"]], function(x) gc_content(x)))
    })
    for(i in 1:length(gc.data)){
      object@probes[[slot]][transcript_ids][[i]][[colname_use]] = gc.data[[i]]
    }
  } else if(slot == "Readout"){
    gc.data = lapply(list(object@probes[[slot]]),function(data){
      unlist(lapply(data[["Sequence"]], function(x) gc_content(x)))
    })
    object@probes[[slot]][[colname_use]] = gc.data[[1]]
  }
  object@param.log[[paste("CalcGC",slot,sep = "-")]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-c(1,length(mget(names(formals()),sys.frame(sys.nframe()))))])
                                                        , additional = match.call(),time = Sys.time()) 
  return(object)
}

#' Calculate Consecutive Repeats
#'
#' Find consecutive repeats of N same nucleotides.
#' 
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param consecutive_repeats Find consecutive repeats of N same nucleotides. Default is 5 nt. See \code{\link{ExtractProbes}} for details.
#' @param colname_use Column name used to store consecutive repeat counts. Default is Consecutive.
#' @return Returns a FISHprobe object.
#' @export
#' @examples \donttest{
#' object = CalcConsecutive(object)
#' }
CalcConsecutive <- function(object, slot = "Target", transcript_ids = NULL, consecutive_repeats = 5, colname_use = "Consecutive"){
  if(slot == "Target"){
    if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
    check_transcripts(object,transcript_ids)
    gc.data = lapply(object@probes[[slot]][transcript_ids],function(data){
      unlist(lapply(data[["Sequence"]], function(x) find_consecut(x, times = consecutive_repeats)))
    })
    for(i in 1:length(gc.data)){
      object@probes[[slot]][transcript_ids][[i]][[colname_use]] = gc.data[[i]]
    }
  } else if(slot == "Readout"){
    gc.data = lapply(list(object@probes[[slot]]),function(data){
      unlist(lapply(data[["Sequence"]], function(x) find_consecut(x, times = consecutive_repeats)))
    })
    object@probes[[slot]][[colname_use]] = gc.data[[1]]
  }
  object@param.log[[paste("CalcConsecutive",slot,sep = "-")]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-c(1,length(mget(names(formals()),sys.frame(sys.nframe()))))])
                                                        , additional = match.call(),time = Sys.time()) 
  return(object)
}

#' Calculate Probe Secondary Structures
#'
#' Calculate secondary structures of the probes. This function utilizes run_RNAfold function from LncFinder package and complexes from NUPACK. 
#' See help page of \code{\link[LncFinder]{run_RNAfold}} and \href{http://www.nupack.org/downloads/documentation}{NUPACK} for details.
#' Please make sure you have installed \href{https://www.tbi.univie.ac.at/RNA/#download}{ViennaRNA} or \href{http://nupack.org/downloads/}{NUPACK}.
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param method Method used to calculate secondary structures, either viennaRNA or nupack. Default is viennaRNA.
#' @param pseudoknots Whether to consider pseudoknots. Default is FALSE.
#' @param Na_M Na concentration in unit molar [M]. Default is 0.33 (330mM, equal to 2XSSC).
#' @param temperature Temperature in degree celsius used. The default is 37 (with adjustment of formamide parameters below).
#' @param formamide_percent Formamide concentration in percentage. Default is 50.
#' @param formamide_factor Formamide factor for temperature calculation, it will adjust for an increase of temperature celsius in [formamide_factor] * [formamide_percent]. Default is 0.62.
#' @param nupack_params Additional parameters for NUPACK. Default is none.
#' @return Returns a FISHprobe object.
#' @importFrom LncFinder run_RNAfold
#' @export
#' @examples \donttest{
#' object = Calc2ndStruct(object)
#' }
Calc2ndStruct <- function(object, slot = "Target", transcript_ids = NULL, method = "viennaRNA", pseudoknots = FALSE, Na_M = 0.33, temperature = 37, formamide_percent = 50, formamide_factor = 0.62, nupack_params = ""){
  if(slot == "Readout"){
    if(check_any_feature(object, slot, "Probe")) stop("Probe not found. Run GenerateProbes to generate probes first.", call. = FALSE)
    if(method == "viennaRNA"){
      object@probes[[slot]] = lapply(list(object@probes[[slot]]), calc_2nd_struct)[[1]]
    } else if(method == "nupack"){
      temperature = temperature + formamide_percent * formamide_factor
      if(pseudoknots) nupack_params = paste("-pseudo",nupack_params)
      second.data = lapply(list(object@probes[[slot]]),function(data){
        data.temp = do.call(rbind,apply(data,1,function(x){
          calc_2nd_struct_nupack(sequence = x["Probe"], params = paste("-T",temperature,"-material dna","-sodium",Na_M,nupack_params))
          }))
        cbind(data, Min_Free_Energy = unlist(data.temp[,1]), Structure = unlist(data.temp[,2]), stringsAsFactors = FALSE)
        })
      object@probes[[slot]] = second.data[[1]]
      } 
    } else if (slot == "Target"){
    if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
    check_transcripts(object,transcript_ids)
    if(check_any_feature(object, slot, "Probe")) stop("Probe not found. Run GenerateProbes to generate probes first.", call. = FALSE)
    if(method == "viennaRNA"){
      object@probes[[slot]][transcript_ids] = lapply(object@probes[[slot]][transcript_ids], calc_2nd_struct)
    } else if(method == "nupack"){
      temperature = temperature + formamide_percent * formamide_factor
      if(pseudoknots) nupack_params = paste("-pseudo",nupack_params)
      second.data = lapply(object@probes[[slot]][transcript_ids] ,function(data){
        data.temp = do.call(rbind,apply(data,1,function(x){
          calc_2nd_struct_nupack(sequence = x["Probe"], params = paste("-T",temperature,"-material dna","-sodium",Na_M,nupack_params))
        }))
        cbind(data, Min_Free_Energy = unlist(data.temp[,1]), Structure = unlist(data.temp[,2]), stringsAsFactors = FALSE)
      })
      names(second.data) = transcript_ids
      object@probes[[slot]][transcript_ids] = second.data
    }
  }
  object@param.log[[paste("Calc2ndStruct",slot,sep = "-")]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                                                                                  , time = Sys.time()) 
  return(object)
}

#' Calculate Melting Temperatures 
#'
#' Calculate Tm melting temperatures of probes. This function relies on Tm_NN and chem_correction functions from TmCalculator package.
#' See \code{\link[TmCalculator]{Tm_NN}} and \code{\link[TmCalculator]{chem_correction}} help page for details.
#' Note that if Tm is negative after adjusted for formamide etc. it will be reset to zero.
#' 
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param Na Na concentration in [mM]. Default is 330mM (equal to 2XSSC).
#' @param probe_nM Concentration of the probes in [nM]. Default is 5nM.
#' @param target_nM Concentration of the targeting sequences in [nM].. Default is 5nM.
#' @param formamide_percent Formamide concentration in percentage. Default is 50.
#' @param formamide_factor Coeffecient of Tm decrease per percent formamide.. Default is 0.65.
#' @param colname_use Column name used to store Tm. Default is Tm.
#' @param ... Other arguments can be passed for Tm calculation. See \code{\link[TmCalculator]{Tm_NN}} for details.
#' @return Returns a FISHprobe object.
#' @importFrom TmCalculator Tm_NN chem_correction
#' @export
#' @examples \donttest{
#' object = CalcTm(object, Na = 350, probe_nM = 10, target_nM = 10, formamide_percent = 10)
#' }
CalcTm <- function(object, slot = "Target", transcript_ids = NULL, Na = 330, probe_nM = 5, target_nM = 5, formamide_percent = 50, formamide_factor = 0.65, colname_use = "Tm", ...){
  if(slot == "Target"){
    if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
    check_transcripts(object,transcript_ids)
    if(check_any_feature(object, slot, "Probe")) stop("Probe not found. Run GenerateProbes to generate probes first.", call. = FALSE)
    tm.data = lapply(object@probes[[slot]][transcript_ids],function(data){
      unlist(lapply(data[["Probe"]], function(x) Tm_calc(x, Na = Na, dnac1 = probe_nM, dnac2 = target_nM, nn_table = "R_DNA_NN1", formamide_percent, fmdfactor = formamide_factor, ...)))
    })
    for(i in 1:length(tm.data)){
      object@probes[[slot]][transcript_ids][[i]][[colname_use]] = tm.data[[i]]
    }
  } else if(slot == "Readout"){
    if(check_any_feature(object, slot, "Probe")) stop("Probe not found. Run GenerateProbes to generate probes first.", call. = FALSE)
    tm.data = lapply(list(object@probes[[slot]]),function(data){
      unlist(lapply(data[["Probe"]], function(x) Tm_calc(x, Na = Na, dnac1 = probe_nM, dnac2 = target_nM, formamide_percent, fmdfactor = formamide_factor, ...)))
    })
    object@probes[[slot]][[colname_use]] = tm.data[[1]]
  }
  object@param.log[[paste("CalcTm",slot,sep = "-")]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-c(1,length(mget(names(formals()),sys.frame(sys.nframe()))))])
                                                        , additional = match.call(),time = Sys.time()) 
  return(object)
}

#' Calculate Repetitive Masks
#'
#' Find the length of overlapping betweem probes and repetitive masked region in the genome.
#' This length varies between 0 (non-overlapping) to the length of the probe (full-overlapping).
#' 
#' @param object FISHprobe object.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param masks_names Names of built-in masks. Default is NULL (the default active build-in masks, see below).
#' \itemize{
#'   \item AGAPS, assembly gaps (big N-blocks), active by default.
#'   \item AMB, intra-contig ambiguities (IUPAC ambiguity letter), active by default.
#'   \item RM, RepeatMasker, not active by default. 
#'   \item TRF, Tandem Repeats Finder, not active by default. 
#'   \item All, will activate all above masks.
#'   }
#' @return Returns a FISHprobe object.
#' @importFrom intervals Intervals interval_intersection
#' @importFrom Biostrings masks
#' @importFrom stringr str_length
#' @importFrom IRanges active
#' @importFrom BiocGenerics start end
#' @importFrom pbapply pblapply
#' @export 
#' @examples \donttest{
#' object = CalcRepMask(object)
#' }
CalcRepMask <- function(object, transcript_ids = NULL, masks_names = NULL){
  if(is.null(transcript_ids)) transcript_ids = names(object@probes[["Target"]])
  check_transcripts(object,transcript_ids)
  if(object@specie == "human"){
    genome.use = BSgenome.Hsapiens.UCSC.hg38.masked::BSgenome.Hsapiens.UCSC.hg38.masked
  } else if(object@specie == "mouse"){
    genome.use = BSgenome.Mmusculus.UCSC.mm10.masked::BSgenome.Mmusculus.UCSC.mm10.masked
  }
  if(length(object@sequence.type) > 1){
    data.use = lapply(1:length(transcript_ids), function(i){
      trans.data = object@sequences[which(object@sequences[["ensembl_transcript_id"]] == transcript_ids[i]),]
      probe.data = object@probes[["Target"]][[transcript_ids[i]]]
      get_repmask(trans.data,probe.data,genome.use, masks_used = masks_names, exon_only = !grepl("intron",object@sequence.type[i]))
    })
  } else {
    data.use = lapply(transcript_ids, function(x){
      trans.data = object@sequences[which(object@sequences[["ensembl_transcript_id"]] == x),]
      probe.data = object@probes[["Target"]][[x]]
      get_repmask(trans.data,probe.data,genome.use, masks_used = masks_names, exon_only = !grepl("intron",object@sequence.type))
    })
  }
  

  names(data.use) = transcript_ids
  object@probes[["Target"]][transcript_ids] = data.use
  object@param.log[[paste("CalcRepMask",paste(transcript_ids, collapse = "_"),sep = "-")]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                                                     , time = Sys.time()) 
  return(object)
}

#' Calculate Duplexing Free Energy and Percentage (Probability)
#'
#' Calculate duplex binding free energy and duplex precetage (can be treated as duplex probability) between probes and targets.
#' 
#' See \href{http://www.nupack.org/downloads/serve_public_file/nupack_user_guide_3.2.2.pdf?type=pdf}{NUPACK 3.2.2 User Guide} for details.
#' Also see \href{https://www.pnas.org/content/115/10/E2183}{OligoMiner} for details.
#' 
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param Na_M Na concentration in unit molar [M]. Default is 0.33 (330mM, equal to 2XSSC).
#' @param probe_M Concentration of the probes in unit molar [M]. Default is 1e-6 (1uM).
#' @param target_M Concentration of the targeting sequences in unit molar [M]. Default is 1e-6 (1uM).
#' @param temperature Temperature in degree celsius used. The default is 37 (with adjustment of formamide parameters below).
#' @param formamide_percent Formamide concentration in percentage. Default is 50.
#' @param formamide_factor Formamide factor for temperature calculation, it will adjust for an increase of temperature celsius in [formamide_factor] * [formamide_percent]. Default is 0.62.
#' @param nupack_params Additional parameters for NUPACK. Default is none.
#' @return Returns a FISHprobe object.
#' @export
#' @examples \donttest{
#' object = CalcDuplex(object, formamide_percent = 20)
#' }
CalcDuplex <- function(object, slot = "Target", transcript_ids = NULL, Na_M = 0.33, probe_M = 1e-6, target_M = 1e-6, temperature = 37, formamide_percent = 50, formamide_factor = 0.62, nupack_params = ""){
  temperature = temperature + formamide_percent * formamide_factor
  if(slot == "Target"){
    if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
    check_transcripts(object,transcript_ids)
    if(check_any_feature(object, slot, "Probe")) stop("Probe not found. Run GenerateProbes to generate probes first.", call. = FALSE)
    duplex.data = lapply(object@probes[[slot]][transcript_ids],function(data){
      if(length(data[["Probe"]]) == 0) stop("Reverse Complementary Probes Not Found.", call. = FALSE)
      data.temp = do.call(rbind,apply(data,1,function(x){
        duplex_binding(sequence = x["Sequence"], probe = x["Probe"], seq_M = target_M, prob_M = probe_M, params = paste("-T",temperature,"-sodium",Na_M,"-material dna",nupack_params))
      }))
      cbind(data, Duplex_Energy = unlist(data.temp[,1]), Equilibrium_Percent = unlist(data.temp[,2]))
    })
    names(duplex.data) = transcript_ids
    object@probes[[slot]][transcript_ids] = duplex.data
  } else if(slot == "Readout"){
    if(check_any_feature(object, slot, "Probe")) stop("Probe not found. Run GenerateProbes to generate probes first.", call. = FALSE)
    duplex.data = lapply(list(object@probes[[slot]]),function(data){
      if(length(data[["Probe"]]) == 0) stop("Reverse Complementary Probes Not Found.", call. = FALSE)
      data.temp = do.call(rbind,apply(data,1,function(x){
        duplex_binding(sequence = x["Sequence"], probe = x["Probe"], seq_M = target_M, prob_M = probe_M, params = paste("-T",temperature,"-sodium",Na_M,"-material dna",nupack_params))
      }))
      cbind(data, Duplex_Energy = unlist(data.temp[,1]), Equilibrium_Percent = unlist(data.temp[,2]))
    })
    object@probes[[slot]] = duplex.data[[1]]
  }
  object@param.log[[paste("CalcDuplex",slot,sep = "-")]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-c(1,length(mget(names(formals()),sys.frame(sys.nframe()))))])
                                                        , additional = match.call(),time = Sys.time()) 
  return(object)
}

#' Calculate Exon Junction Probes
#'
#' Calculate exon junction reads for each probe.
#' Please note the numbering of the exons may be different from transcript to transcript of the same gene. 
#' For example Exon-5 can have different sequences between two isoforms.
#' @param object FISHprobe object.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @return Returns a FISHprobe object.
#' @importFrom pbmcapply pbmclapply
#' @import intervals
#' @export
#' @examples \donttest{
#' object = CalcJunction(object)
#' }
CalcJunction <- function(object, transcript_ids = NULL){
  if(is.null(transcript_ids)) transcript_ids = names(object@probes[["Target"]])
  check_transcripts(object,transcript_ids)
  junction_data = pbmclapply(transcript_ids, function(x) get_junction(object, x))
  for(i in 1:length(transcript_ids)){
    object[[transcript_ids[i]]] = cbind.data.frame(object[[transcript_ids[i]]], Junction = junction_data[[i]], stringsAsFactors = FALSE)
  }
  object@param.log[["CalcJunction"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-c(1,length(mget(names(formals()),sys.frame(sys.nframe()))))])
                                                               , additional = match.call(),time = Sys.time()) 
  return(object)
}

#' Calculate Exon Information 
#'
#' Calculate exon information of the probes.
#' Please note the numbering of the exons may be different from transcript to transcript of the same gene. 
#' For example Exon-5 can have different sequences between two isoforms.
#' @param object FISHprobe object.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @return Returns a FISHprobe object.
#' @importFrom pbmcapply pbmclapply
#' @import intervals
#' @export
#' @examples \donttest{
#' object = CalcJunction(object)
#' }
CalcExon <- function(object, transcript_ids = NULL){
  if(is.null(transcript_ids)) transcript_ids = names(object@probes[["Target"]])
  check_transcripts(object,transcript_ids)
  junction_data = pbmclapply(transcript_ids, function(x) get_junction(object, x, do_exon = TRUE))
  for(i in 1:length(transcript_ids)){
    object[[transcript_ids[i]]] = cbind.data.frame(object[[transcript_ids[i]]], Exon = junction_data[[i]], stringsAsFactors = FALSE)
  }
  object@param.log[["CalcExon"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-c(1,length(mget(names(formals()),sys.frame(sys.nframe()))))])
                                        , additional = match.call(),time = Sys.time()) 
  return(object)
}

#' Filter Probes 
#'
#' Filter probes based on given pramameters. 
#' Be cautious for this step, use it as the last step.
#' If accidental filtering is done please do object@probe=object@stash.probes if probes are stashed first. Then redo filtering again.
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param filter_names Names to filtering the probes on. See \code{\link{ListFeatures}} for features names for filtering.
#' @param min Minimum values used for filtering. Default is -Inf (one feature).
#' @param max Maximum values used for filtering. Default is Inf (one feature).
#' @param equal_to Filter probes with features equal to this value. Default is None (if provided will override min and max filters).
#' @param reverse_select Whether to perform reverse selection on equal_to, if set TRUE will do "not equal to" operation. Default is FALSE.
#' @param do.filter Wether to perform filtering on the probes. Default is FALSE.
#' @param warning Wether to show warning message. Default is FALSE.
#' @return Returns a FISHprobe object.
#' @export
#' @examples \donttest{
#' object = FilterProbes(object, filter_names = "GC", min = 50, max = 60)
#' }
FilterProbes <- function(object, slot = "Target", transcript_ids = NULL, filter_names, min = -Inf, max = Inf, equal_to = NULL, reverse_select = FALSE, do.filter = FALSE, warning = FALSE){
  if(slot == "Readout"){
    feature.all = ListFeatures(object, slot = slot)
    if(length(feature.all) == 0) stop("No feature to filter", call. = FALSE)
    for(i in 1:length(feature.all)){
      if(!all(filter_names %in% feature.all[[i]]))
        stop(paste("Cannot find filter_names", filter_names[which(!filter_names %in% feature.all[[i]])],"in readout probes"),call. = FALSE) 
    }
    if(length(unique(length(filter_names), length(min), length(max))) != 1) stop("filter_names, min and max must be in equal length!", call. = FALSE)
    filter_names.used = paste(paste(filter_names, min, max, sep = "_"),collapse = " & ")
    filter.all = lapply(list(object@filter.table[[slot]]),colnames)
    if(length(filter.all) != 0){
      for(i in 1:length(filter.all)){
        if(filter_names.used %in% filter.all[[i]]){
          if(warning) warning("Duplicated filtering is found in readout probe",call. = FALSE) 
          object@filter.table[[slot]] = object@filter.table[[slot]][,which(colnames(object@filter.table[[slot]]) != filter_names.used),drop=FALSE]
        }
      }
    }
    if(is.null(equal_to)){
      probe.temp = lapply(list(object@probes[[slot]]), function(x){
        for(i in 1:length(filter_names)){
          x = x[which(x[[filter_names[i]]] >= min[i] & x[[filter_names[i]]] <= max[i]),]
        }
        return(x)
      })
      filter_names.used = paste(paste(filter_names, min, max, sep = "_"),collapse = " & ")
    } else {
      if(reverse_select){
        probe.temp = lapply(list(object@probes[[slot]]), function(x){
          for(i in 1:length(filter_names)){
            x = x[which(x[[filter_names[i]]] != equal_to[i]),]
          }
          return(x)
        })
      } else {
        probe.temp = lapply(list(object@probes[[slot]]), function(x){
          for(i in 1:length(filter_names)){
            x = x[which(x[[filter_names[i]]] == equal_to[i]),]
          }
          return(x)
        })
      }
      filter_names.used = paste(paste(filter_names, equal_to, sep = "_"),collapse = " & ")
    }
    filter.temp = matrix(0, nrow = nrow(object@probes[[slot]]), ncol = 1)
    colnames(filter.temp) = filter_names.used
    rownames(filter.temp) = rownames(object@filter.table[[slot]])
    filter.temp[which(rownames(filter.temp) %in% probe.temp[[1]][["Sequence"]]),] = 1
    object@filter.table[[slot]] = cbind(object@filter.table[[slot]], filter.temp)
    if(do.filter){
      if(warning) warning("Redout probes are filtered!", call. = FALSE)
      object@probes[[slot]] = probe.temp[[1]]
    }
  } else if(slot == "Target"){
    if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
    feature.all = ListFeatures(object, slot = slot, transcript_ids = transcript_ids)
    if(length(feature.all) == 0) stop("No feature to filter", call. = FALSE)
    for(i in 1:length(feature.all)){
      if(!all(filter_names %in% feature.all[[i]]))
        stop(paste("Cannot find filter_names", filter_names[which(!filter_names %in% feature.all[[i]])],"in transcript", names(feature.all)[i]),call. = FALSE) 
    }
    if(length(unique(length(filter_names), length(min), length(max))) != 1) stop("filter_names, min and max must be in equal length!", call. = FALSE)
    filter_names.used = paste(paste(filter_names, min, max, sep = "_"),collapse = " & ")
    filter.all = lapply(object@filter.table[[slot]][transcript_ids],colnames)
    if(length(filter.all) != 0){
      for(i in 1:length(filter.all)){
        if(filter_names.used %in% filter.all[[i]]){
          if(warning) warning(paste("Duplicated filtering is found in","in transcript", transcript_ids[i]),call. = FALSE) 
          object@filter.table[[slot]][transcript_ids][[i]] = object@filter.table[[slot]][transcript_ids][[i]][,which(colnames(object@filter.table[[slot]][transcript_ids][[i]]) != filter_names.used),drop=FALSE]
        }
      }
    }
    if(is.null(equal_to)){
      probe.temp = lapply(object@probes[[slot]][transcript_ids], function(x){
        for(i in 1:length(filter_names)){
          x = x[which(x[[filter_names[i]]] >= min[i] & x[[filter_names[i]]] <= max[i]),]
        }
        return(x)
      })
      filter_names.used = paste(paste(filter_names, min, max, sep = "_"),collapse = " & ")
    } else {
      if(reverse_select){
        probe.temp = lapply(object@probes[[slot]][transcript_ids], function(x){
          for(i in 1:length(filter_names)){
            x = x[which(x[[filter_names[i]]] != equal_to[i]),]
          }
          return(x)
        })
      } else {
        probe.temp = lapply(object@probes[[slot]][transcript_ids], function(x){
          for(i in 1:length(filter_names)){
            x = x[which(x[[filter_names[i]]] == equal_to[i]),]
          }
          return(x)
        })
      }
      filter_names.used = paste(paste(filter_names, equal_to, sep = "_"),collapse = " & ")
    }
    names(probe.temp) = transcript_ids
    if(length(object@filter.table[[slot]][transcript_ids]) == 0){
      object@filter.table[[slot]][transcript_ids] = lapply(1:length(probe.temp),function(i){
        filter.temp = matrix(0, nrow = nrow(object@probes[[slot]][transcript_ids][[i]]), ncol = 1)
        colnames(filter.temp) = filter_names.used
        rownames(filter.temp) = object@probes[[slot]][transcript_ids][[i]][["Sequence"]]
        filter.temp[which(rownames(filter.temp) %in% probe.temp[[i]][["Sequence"]]),] = 1
        object@filter.table[[slot]][[transcript_ids[i]]] = cbind(object@filter.table[[slot]][[transcript_ids[i]]], filter.temp)
      })
    } else {
      object@filter.table[[slot]][transcript_ids] = lapply(1:length(probe.temp),function(i){
        filter.temp = matrix(0, nrow = nrow(object@filter.table[[slot]][transcript_ids][[i]]), ncol = 1)
        colnames(filter.temp) = filter_names.used
        rownames(filter.temp) = rownames(object@filter.table[[slot]][transcript_ids][[i]])
        filter.temp[which(rownames(filter.temp) %in% probe.temp[[i]][["Sequence"]]),] = 1
        object@filter.table[[slot]][[transcript_ids[i]]] = cbind(object@filter.table[[slot]][[transcript_ids[i]]], filter.temp)
      })
    }
    if(do.filter){
      if(warning) warning(paste(paste(transcript_ids, collapse = " "),"probes are filtered!"), call. = FALSE)
      object@probes[[slot]][transcript_ids] = probe.temp
    }
  }
  object@param.log[[paste("FilterProbes",slot,sep = "-")]] = list(call = data.frame(filter_names, max, min), do.filter = do.filter
                                                                     , time = Sys.time())
  return(object)
}

#' Order Probes
#'
#' Order probes by a given parameter (starting/end positions, GC content etc).
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param order_by Parameter to order the probes. Any columns names in object@probes and object@off.target.expr can be used such as GC, Start etc. Default is -5.
#' @param decrease Whether order by decreasing order. Default is FALSE.
#' @return Returns a FISHprobe object.
#' @export
#' @examples \donttest{
#' # order probes by their starting positions in increasing order
#' object = OrderProbe(object, order_by = "Start")
#' # order probes by their GC contents in decreasing order
#' object = OrderProbe(object, order_by = "GC", decrease = TRUE)
#' }
OrderProbe <- function(object, slot = "Target", transcript_ids = NULL, order_by = "Start", decrease = FALSE){
  if(slot == "Readout"){
    name.1 = colnames(object@probes[[slot]])
    name.2 = colnames(object@off.target.expr[[slot]])
    if(order_by %in% name.1){
      data.use = object@probes[[slot]]
    } else if(order_by %in% name.2){
      data.use = object@off.target.expr[[slot]]
    } else {
      stop("order_by option not found in @probes or @off.target.expr in readout probes", call. = FALSE)
    }
    probe_order = order(unlist(data.use[[order_by]]), decreasing = decrease)
    object@probes[[slot]] = object@probes[[slot]][probe_order,]
  } else if(slot == "Target"){
    if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
    check_transcripts(object,transcript_ids)
    for(i in 1:length(transcript_ids)){
      name.1 = colnames(object@probes[[slot]][[transcript_ids[i]]])
      name.2 = colnames(object@off.target.expr[[slot]][[transcript_ids[i]]])
      if(order_by %in% name.1){
        data.use = object@probes[[slot]][[transcript_ids[i]]]
      } else if(order_by %in% name.2){
        data.use = object@off.target.expr[[slot]][[transcript_ids[i]]]
      } else {
        stop("order_by option not found in @probes or @off.target.expr in ", transcript_ids[i], call. = FALSE)
      }
      probe_order = order(unlist(data.use[[order_by]]), decreasing = decrease)
      object@probes[[slot]][[transcript_ids[i]]] = object@probes[[slot]][[transcript_ids[i]]][probe_order,]
    }
  }
  object@param.log[[paste("OrderProbes",slot,sep = "-")]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-c(1,length(mget(names(formals()),sys.frame(sys.nframe()))))])
                                                               , additional = match.call(),time = Sys.time())                                                               
  return(object)
}

#' Generate Genome Browser with Probe Tracks
#'
#' Generate genome browsers with probe tracks. For probe set over 1000 probes please split it with \code{\link{SubsetProbes}} first.
#' 
#' If a FISHprobe object is given, it will generate genome browser files for all the genes and transcripts in the object. 
#' A folder contains the genome browser files for each gene will be created in the current directory with the name of \emph{GeneName-TranscriptID}.
#' 
#' If no FISHprobe object is given, instead gene, transcript_id and seqeuences are given, it will generate genome browser for the input gene/transcript and probe sequences.
#' A genome browser folder named as \emph{GeneName} will be created.
#' Note that under this mode the probe tracks are generated using BLAST and the probes contain exon-exon junction may show multiple bands (>2) instead of 2.
#' 
#' The genome browser is created using D3GB R package, in order to visualize it with your local browser please refer to this \href{ http://d3gb.usal.es/help/index.html}{page}.
#' 
#' @param object FISHprobe object.
#' @param gene Gene symbols. Either Target or Readout. Default is NULL.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param seqeuences Probe sequences. Default is NULL.
#' @param specie Either human or mouse. Default is human.
#' @param BLAST_args Arguments passed to BLAST in command-line style. Default is -task blastn-short. See  \code{\link{BLASTProbes}} for detail.
#' @param file_name Additional file name for genome browser folder. Default is none.
#' @param sort Sort probes by targeting positions in increasing order.
#' @return Returns a FISHprobe object.
#' @importFrom biomaRt getBM getSequence
#' @importFrom D3GB genome_addTrack genomebrowser getAssemblyFromFasta genome_addSequence
#' @importFrom stringr str_length 
#' @importFrom seqinr write.fasta
#' @importFrom rBLAST makeblastdb blast
#' @importFrom Biostrings DNAStringSet
#' @export
#' @examples \donttest{
#' GenerateProbeTracks(object)
#' # generate customized genome browser 
#' GenerateProbeTracks(gene = "Actb",transcript_id = "ENSMUSG00000029580",
#'                     specie = "mouse",seqeuences = ProbeSequences)
#'}
GenerateProbeTracks <- function(object = NULL, gene = NULL, transcript_ids = NULL, seqeuences = NULL, specie = "human", BLAST_args = "-task blastn-short", file_name = NULL, sort = FALSE){
  if(is.null(object)){
    if(!is.null(gene) & !is.null(seqeuences) & !is.null(transcript_ids)){
      gene.symbol = ifelse(specie == "human", "hgnc_symbol","mgi_symbol")
      if(specie == "human" ){
        mart.use = FISHprobe::human
      } else if(specie == "mouse"){
        mart.use = FISHprobe::mouse
      }
      gene.bed = getBed(gene.name = gene, trans_used = transcript_ids, specie = specie, mart = mart.use)
      gb = generateBrowser(gene.bed = gene.bed, dir.name = paste(c(gene,file_name), collapse = "-"), selected = transcript_ids, specie = specie, mart = mart.use)
      seq.use = biomaRt::getSequence(id = gene.bed[["name"]][1], type = gene.symbol, seqType = "gene_exon_intron", mart = mart.use)
      if(nrow(seq.use) > 1){
        cat("More than one transcript sequences found, only the LONGEST one were used. Please check further for gene:", genes.use)
        gene.len = str_length(seq.use[["gene_exon_intron"]])
        seq.use = seq.use[which.max(gene.len),]
      }
      generateTracks.blast(seqeuences, seq.use, gb, BLAST_args = BLAST_args, order = sort)
      invisible(file.remove(list.files(pattern = "probe_blast.fa")))
    }
  } else {
    if(object@specie == "human" ){
      mart.use = FISHprobe::human
    } else if(object@specie == "mouse"){
      mart.use = FISHprobe::mouse
    }
    if(is.null(transcript_ids)){
      transripts.use = names(object@probes[["Target"]])
    } else {
      check_transcripts(object,transcript_ids)
      transripts.use = transcript_ids
    }
    for(i in 1:length(transripts.use)){
      gene.symbol = ifelse(object@specie == "human", "hgnc_symbol","mgi_symbol")
      genes.use = object@sequences[which(object@sequences[["ensembl_transcript_id"]] == transripts.use[i]), "GeneName"]
      all.trans = object@sequences[which(object@sequences[["GeneName"]] == genes.use), "ensembl_transcript_id"]
      gene.bed = getBed(gene.name = genes.use, trans_used = all.trans, specie = object@specie,mart = mart.use)
      gb = generateBrowser(gene.bed = gene.bed, dir.name = paste(c(genes.use,transripts.use[i],file_name),collapse = "-"), selected = transripts.use[i], specie = object@specie, mart = mart.use)
      seq.use = biomaRt::getSequence(id = gene.bed[["name"]][1], type = gene.symbol, seqType = "gene_exon_intron", mart = mart.use)
      if(nrow(seq.use) > 1){
        cat("More than one transcript sequences found, only the LONGEST one were used. Please check further for gene:", genes.use)
        gene.len = str_length(seq.use[["gene_exon_intron"]])
        seq.use = seq.use[which.max(gene.len),]
      }
      generateTracks(object@probes[["Target"]][[transripts.use[i]]], transripts.use[i], mart.use, seq.use, gb, sort)
    }
  }
}

#' Save Probe Sets
#' 
#' Save probe sets as .csv files with file names as <TranscriptID>_probes.csv.
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param suffix Additional suffix to be added to the probe file name. Default is none. 
#' @export
#' @examples \donttest{
#' SaveProbes(objcet)
#' }
SaveProbes <- function(object, slot = "Target", transcript_ids = NULL, suffix = ""){
  if(slot == "Readout"){
    write.csv(object@probes[[slot]], file = paste("Final.Probes",suffix,"csv",sep = "."))
  } else if(slot == "Target"){
    if(is.null(transcript_ids)){
      transcripts.use = names(object@probes[[slot]])
    } else {
      check_transcripts(object,transcript_ids)
      transcripts.use = transcript_ids
    }
    for(i in 1:length(transcripts.use)){
      data.save = object@probes[[slot]][[transcripts.use[i]]]
      data.save = apply(data.save, 2, as.character)
      gene.name = object@sequences[["GeneName"]][which(object@sequences[["ensembl_transcript_id"]] == transcripts.use[i])]
      write.csv(data.save, file = paste(gene.name,transcripts.use[i],suffix,"probes.csv",sep = "_"))
    }
    if(length(object@final.probes) > 0){
      final.probes.data = do.call(rbind,object@final.probes)
      write.csv(final.probes.data, file = paste("Final.Probes",suffix,"csv",sep = "."))
    }
  }
}

#' Stash Probe Sets
#' 
#' Stash probe sets into stash.probes slot.
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @export
#' @examples \donttest{
#' StashProbes(object)
#' }
StashProbes <- function(object, slot = "Target", transcript_ids = NULL){
  if(is.null(transcript_ids)){
    transcripts.use = names(object@probes[[slot]])
  } else {
    check_transcripts(object,transcript_ids)
    transcripts.use = transcript_ids
  }
  object@stash.probes[[slot]][transcripts.use] = object@probes[[slot]][transcripts.use]
  object@param.log[[paste("StashProbes",slot,sep = "-")]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-c(1,length(mget(names(formals()),sys.frame(sys.nframe()))))])
                                                                 , additional = match.call(),time = Sys.time())  
  return(object)
}

#' Grid Search for Optimal Probe Filtering
#'
#' Grid search for best parameters to filter the probes.
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param features_list List of feature names with min and max to filter, similar to \code{\link{FilterProbes}}
#' @param target_number Targeted number of probes. Default is 30 probes.
#' @param do.filter Wether to perform filtering on the probes. Default is FALSE.
#' @param order_by_target Order according to the closeness to target_number. If set FALSE, the filtering will be order by number from high to low. Default is TRUE.
#' @return Returns a data frame of grid search results or a FISHprobe object (depends on do.filter).
#' @export
#' @examples \donttest{
#' # Create a feature list with feature names and min/max values
#' feature.list = list(BLAST_Hits = list(min = c(5,10)),GC = list(min = c(.4,.5), max = c(.5,.7)))
#' search.res = GridSelect(object, features_list = feature.list, target_number = 25)
#' print(search.res)
#' 
#' # Perform filtering on the object
#' object.new = GridSelect(object, features_list = feature.list, target_number = 25, do.filter = TRUE)
#' }
GridSelect <- function(object, slot = "Target", transcript_ids = NULL, features_list, target_number = 30, do.filter = FALSE, order_by_target = TRUE){
  if(slot == "Readout"){
    return(grid_search(object = object, slot = slot, transcript_id = transcript_ids, features_list = features_list, target_number = target_number, do_filter = do.filter, order_by_target = order_by_target))
  } else if(slot == "Target"){
    if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
    check_transcripts(object,transcript_ids)
    grid.list = lapply(transcript_ids, function(trans){
      grid_search(object = object, slot = slot, transcript_id = trans, features_list = features_list, target_number = target_number, do_filter = do.filter, order_by_target = order_by_target)
    })
    if(do.filter){
      grid.list = lapply(grid.list, function(obj) obj@probes[[slot]][[1]])
      names(grid.list) = transcript_ids
      object@probes[[slot]] = grid.list
      return(object)
    } else {
      names(grid.list) = transcript_ids
      return(grid.list)
    }
  }
}

#' Differenctial Probes Between Transcript Isoforms
#' 
#' Get differential or shared Probes between different transcript isoforms.
#' Note this won't work if the gene spans multiple locations on different chromosomes.
#' @param object FISHprobe object.
#' @param transcripts_selected Transcript IDs selected to get the differential sequences. 
#' @param transcripts_to_compare Transcript IDs to compare to. Default is all other isoforms than transcript_selected from the same gene.
#' @param find_shared Whether to find common sequences between isoforms.
#' @param probe_length Probe length. Default is 25 nt. See \code{\link{ExtractProbes}} for details.
#' @param consecutive_repeats Find consecutive repeats of N same nucleotides. Default is 5 nt. See \code{\link{ExtractProbes}} for details.
#' @return Returns a \code{\link{DiffProbe-class}} object (subclass of \code{\link{FISHprobe-class}}) with differential or shared probes between isoforms.
#' @import stringr 
#' @import intervals 
#' @import biomaRt
#' @export
#' @examples \donttest{
#' # ENST00000216037 and ENST00000344347 are isoforms of human XBP1 genes. 
#' xbp1 = CreateProbeObject("XBP1")
#' xbp1 = TissueExpr(xbp1)
#' xbp1 = GetSequence(xbp1)
#' xbp1_diff = DiffProbes(xbp1, transcripts_selected = "ENST00000216037", 
#'                        transcripts_to_compare = "ENST00000344347")
#' }
DiffProbes <- function(object, transcripts_selected, transcripts_to_compare = NULL, find_shared = FALSE, probe_length = 25, consecutive_repeats = 5){
  if(object@sequence.type != "cdna") stop("Please use cdna as input sequences in @sequences slot.", call. = FALSE)
  seqs_data = object@sequences
  if(!all(transcripts_selected %in% seqs_data[["ensembl_transcript_id"]])) stop(paste("Cannot find selected ensemble transcript id", paste(transcripts_selected[which(!transcripts_selected %in% seqs_data[["ensembl_transcript_id"]])],collapse = " ")),call. = FALSE)
  gene_name = unique(seqs_data[["GeneName"]][which(seqs_data[["ensembl_transcript_id"]] %in% transcripts_selected)])
  if(length(gene_name) > 1) stop("Isoforms provided are selected from different genes: ", paste(gene_name, collapse = ","), call. = FALSE)
  seqs_data = seqs_data[which(seqs_data[["GeneName"]] == gene_name),]
  if(length(unique(seqs_data[["Chromosome"]])) > 1) stop("Gene spans multiple chromosomes.")
  if(is.null(transcripts_to_compare)){
    transcripts_to_compare = seqs_data[["ensembl_transcript_id"]][which(!seqs_data[["ensembl_transcript_id"]] %in% transcripts_selected)]
  } else {
    if(length(intersect(transcripts_selected,transcripts_to_compare)) > 0) stop(paste("Overlap between selected and compared transcript IDs", paste(intersect(transcripts_selected,transcripts_to_compare), collapse = " ")), call. = FALSE)
    if(!all(transcripts_to_compare %in% seqs_data[["ensembl_transcript_id"]])) stop(paste("Cannot find ensemble transcript id to compare", paste(transcripts_to_compare[which(!transcripts_to_compare %in% seqs_data[["ensembl_transcript_id"]])],collapse = " ")),call. = FALSE)
  }
  exon_selected = do.call(c, lapply(transcripts_selected, function(x) get_exon(object, x)))
  exon_to_compare = do.call(c, lapply(transcripts_to_compare, function(x) get_exon(object, x)))
  if(find_shared){
    exon_intv = interval_intersection(exon_selected, exon_to_compare)
    if(length(exon_intv) == 0) stop("No Common Sequences found between given isoforms.", call. = FALSE)
  } else {
    exon_intv = close_intervals(interval_difference(exon_selected, exon_to_compare))
    if(length(exon_intv) == 0) stop("No Differential Sequences found between given isoforms.", call. = FALSE)
  }
  exon_intv = as.matrix(exon_intv)
  strand_use = unique(seqs_data[["Strand"]])
  if(object@specie == "human") {
    message("Getting Sequences...")
    seq_all = biomaRt::getSequence(id = unique(seqs_data[["ensembl_gene_id"]]), type = "ensembl_gene_id", seqType = "gene_exon_intron", mart = FISHprobe::human)
    seq_start = biomaRt::getBM(attributes=c('start_position'),
                               filters = 'ensembl_gene_id',
                               values = unique(seqs_data[["ensembl_gene_id"]]),
                               mart = FISHprobe::human)[[1]]
    trans_start_end = biomaRt::getBM(attributes=c('transcript_start','transcript_end'),
                                     filters = 'ensembl_transcript_id',
                                     values = transcripts_selected,
                                     mart = FISHprobe::human)[1,]
    message("Calculating...")
    exon_intv_new = exon_intv - seq_start + 1
    if(strand_use == -1) exon_intv_new = str_length(seq_all[,1]) - exon_intv_new[,2:1,drop=FALSE] 
    intv_seq = apply(exon_intv_new, 1, function(x) substr(seq_all[,1], start = x[1], stop = x[2]))
  } else if (object@specie == "mouse"){
    message("Getting Sequences...")
    seq_all = biomaRt::getSequence(id = unique(seqs_data[["ensembl_gene_id"]]), type = "ensembl_gene_id", seqType = "gene_exon_intron", mart = FISHprobe::mouse)
    seq_start = biomaRt::getBM(attributes=c('start_position'),
                               filters = 'ensembl_gene_id',
                               values = unique(seqs_data[["ensembl_gene_id"]]),
                               mart = FISHprobe::mouse)[[1]]
    trans_start_end = biomaRt::getBM(attributes=c('transcript_start','transcript_end'),
                                     filters = 'ensembl_transcript_id',
                                     values = unique(seqs_data[["ensembl_transcript_id"]]),
                                     mart = FISHprobe::mouse)[1,]
    message("Calculating...")
    exon_intv_new = exon_intv - seq_start + 1
    if(strand_use == -1) exon_intv_new = str_length(seq_all[,1]) - exon_intv_new[,2:1,drop=FALSE] 
    intv_seq = apply(exon_intv_new, 1, function(x) substr(seq_all[,1], start = x[1], stop = x[2]))
  }
  intv_seq = data.frame(DiffProbe = intv_seq, stringsAsFactors = FALSE)
  intv_seq[["GeneName"]] = gene_name
  intv_seq[["ExonStarts"]] = as.character(exon_intv[,1])
  intv_seq[["ExonEnds"]] = as.character(exon_intv[,2])
  if(strand_use == -1){
    intv_seq[["FullSeqStarts"]] = as.numeric(trans_start_end[1,2]) - as.numeric(intv_seq[["ExonEnds"]]) + 1
    intv_seq[["FullSeqEnds"]] = as.numeric(trans_start_end[1,2]) - as.numeric(intv_seq[["ExonStarts"]]) + 1
  } else {
    intv_seq[["FullSeqStarts"]] = as.numeric(intv_seq[["ExonStarts"]]) - trans_start_end[1,1] + 1
    intv_seq[["FullSeqEnds"]] = as.numeric(intv_seq[["ExonEnds"]]) - trans_start_end[1,1] + 1
  }
  intv_seq[["Transcript_Selected"]] = paste(transcripts_selected, collapse = ",")
  intv_seq[["Transcript_Compared"]] = paste(transcripts_to_compare, collapse = ",")
  intv_seq[["Chromosome"]] = unique(seqs_data[["Chromosome"]])
  intv_seq[["Strand"]] = unique(seqs_data[["Strand"]])
  intv_seq[["ensembl_transcript_id"]] = transcripts_selected
  intv_seq = intv_seq[order(intv_seq[["FullSeqStarts"]], decreasing = FALSE),]
  rownames(intv_seq) = paste(paste(transcripts_selected, collapse = "+"), 1:nrow(intv_seq), sep = "_DiffProbe_")
  object_new = CreateProbeObject(customized_name = gene_name, specie_type = object@specie, customized_data = " ")
  object_new@sequences = intv_seq
  
  trans_data = object@sequences[which(object@sequences[["ensembl_transcript_id"]] %in% transcripts_selected),]
  exons_starts = as.numeric(strsplit(trans_data[["ExonStarts"]],",")[[1]])
  exons_ends = as.numeric(strsplit(trans_data[["ExonEnds"]],",")[[1]])
  exon_data = intervals::Intervals(cbind(exons_starts,exons_ends))
  intv_data = intervals::Intervals(as.numeric(as.matrix(intv_seq[,c("ExonStarts","ExonEnds")])))
  intv_include = apply(intv_data,1,function(x) which(interval_included(exon_data,intervals::Intervals(x)) != 0))
  if(trans_data[["Strand"]] == -1){
    intv_diff = exon_data[intv_include] - intv_data 
    intv_diff = intv_diff[,2:1,drop=FALSE]
  } else {
    intv_diff = intv_data - exon_data[intv_include]
  }
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
  intv_new = seq_data[intv_include,] + intv_diff
  object_new@sequences[["SeqStarts"]] = intv_new[,1]
  object_new@sequences[["SeqEnds"]] = intv_new[,2]
  intv_seq[["SeqStarts"]] = intv_new[,1]
  intv_seq[["SeqEnds"]] = intv_new[,2]
  object_new@sequences[["Diff_Exon"]] = paste("Exon",intv_include)
  
  object_new = ExtractProbes(object_new, probe_length = probe_length, consecutive_repeats = consecutive_repeats)
  probe.data = do.call(rbind,object_new@probes[["Target"]])
  probe.data[["Diff_Region"]] = rep(object_new@sequences[["Diff_Exon"]], unlist(lapply(object_new@probes[["Target"]],nrow)))
  rownames(probe.data) = NULL
  object_new@probes[["Target"]] = list(probe.data)
  names(object_new@probes[["Target"]]) = paste(transcripts_selected,collapse = "+") 
  object_new@sequences = object@sequences
  object_new@transcripts[[1]] = c(transcripts_selected,transcripts_to_compare)
  object_new@sequence.type = object@sequence.type
  object_new = as(object_new, "DiffProbe")
  object_new@DiffProbe = intv_seq
  object_new@DiffProbe.type = ifelse(find_shared, "Shared",  "Differential")
  object_new@param.log[["DiffProbes"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-c(1,length(mget(names(formals()),sys.frame(sys.nframe()))))])
                                                                 , additional = match.call(),time = Sys.time())  
  return(object_new)
}

#' Differenctial Exon Junction Probes Between Transcript Isoforms
#' 
#' Get differential or shared Probes that targeting exon junctions between different transcript isoforms.
#' Note this won't work if the gene spans multiple locations on different chromosomes.
#' @param object FISHprobe object.
#' @param transcripts_selected Transcript IDs selected to get the differential sequences. 
#' @param transcripts_to_compare Transcript IDs to compare to. Default is all other isoforms than transcript_selected from the same gene.
#' @return Returns a FISHprobe object, with "DiffJunction" added to @probes slot with differential, shared or none labels indicating exon junction types.
#' @importFrom pbmcapply pbmclapply
#' @export
#' @examples \donttest{
#' # ENST00000216037 and ENST00000344347 are isoforms of human XBP1 genes. 
#' xbp1 = CreateProbeObject("XBP1")
#' xbp1 = TissueExpr(xbp1)
#' xbp1 = GetSequence(xbp1)
#' xbp1_diff = ExtractProbes(xbp1)
#' object = CalcJunction(object)
#' object = DiffJunction(object, transcripts_selected = "ENST00000216037", 
#'                       transcripts_to_compare = "ENST00000344347")
#' }
DiffJunction <- function(object, transcripts_selected, transcripts_to_compare = NULL){
  check_transcripts(object,transcripts_selected)
  if(is.null(transcripts_to_compare)){
    transcripts_to_compare = setdiff(Transcripts(object),transcripts_selected)
  } else {
    check_transcripts(object,transcripts_to_compare)
    if(length(intersect(transcripts_selected,transcripts_to_compare)) > 0) stop(paste("Overlap between selected and compared transcript IDs", paste(intersect(transcripts_selected,transcripts_to_compare), collapse = " ")), call. = FALSE)
  }
  data_select = do.call(rbind.data.frame,object@probes[["Target"]][transcripts_selected])[["Junction"]]
  data_compare = do.call(rbind.data.frame,object@probes[["Target"]][transcripts_to_compare])[["Junction"]]
  data_select = unique(data_select[which(data_select != "None")])
  data_compare = unique(data_compare[which(data_compare != "None")])
  junc_share = intersect(data_select,data_compare)
  junc_diff1 = setdiff(data_select,data_compare)
  junc_diff2 = setdiff(data_compare,data_select)
  for(i in 1:length(transcripts_selected)){
    object[[transcripts_selected[i]]][["DiffJunction"]] = "None"
    object[[transcripts_selected[i]]][["DiffJunction"]][which(object[[transcripts_selected[i]]][["Junction"]] %in% junc_diff1)] = "Differential"
    object[[transcripts_selected[i]]][["DiffJunction"]][which(object[[transcripts_selected[i]]][["Junction"]] %in% junc_share)] = "Shared"
  }
  for(i in 1:length(transcripts_to_compare)){
    object[[transcripts_to_compare[i]]][["DiffJunction"]] = "None"
    object[[transcripts_to_compare[i]]][["DiffJunction"]][which(object[[transcripts_to_compare[i]]][["Junction"]] %in% junc_diff2)] = "Differential"
    object[[transcripts_to_compare[i]]][["DiffJunction"]][which(object[[transcripts_to_compare[i]]][["Junction"]] %in% junc_share)] = "Shared"
  }
  object@param.log[["DiffJunction"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-c(1,length(mget(names(formals()),sys.frame(sys.nframe()))))])
                                                                 , additional = match.call(),time = Sys.time())  
  return(object)
}
 
#' Split Probes into Two Sub Sequences
#' 
#' Split the current transcript specific targeting probes into two small sub-sequences with anticipated lengths.
#' New probes will have a new feature named "Part" either A or B indicating the first and second sub-sequences. 
#' Note that this will add to the previously calculated variable names with an "Original_" prefix such as "Original_GC".
#' Therefore you have to recalculate GC content, Tm etc. for each of the newly splitted probes.
#' @param object FISHprobe object.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param min_length Minimal sub-probe length in nucleotides. Default is 10 nt. 
#' @param nick_length Nick region length in nucleotides. Default is 0 nt.
#' @return Returns a new FISHprobe object with splitted probes.
#' @importFrom stringr str_length 
#' @importFrom pbapply pblapply 
#' @export
#' @examples \donttest{
#' object = SplitProbes(object, min_length = 15, nick_length = 1)
#' }
SplitProbes <- function(object, transcript_ids = NULL, min_length = 10, nick_length = 0){
  if(is.null(transcript_ids)) transcript_ids = names(object@probes[["Target"]])
  check_transcripts(object,transcript_ids)
  object@probes[["Target"]][transcript_ids] = pblapply(object@probes[["Target"]][transcript_ids], function(trans){
    colnames(trans) = paste("Original", colnames(trans), sep = "_") 
    data.temp = lapply(trans[["Original_Sequence"]], split_probes, min = min_length, nick = nick_length)
    for(i in 1:length(data.temp)){
      data.temp[[i]][,c("Start","Stop")] = data.temp[[i]][,c("Start","Stop")] + trans[["Original_Start"]][i] - 1
      data.temp[[i]] = cbind(data.temp[[i]], trans[i,-1], row.names = NULL)
    }
    data.temp = do.call(rbind, data.temp)
  })
  object@filter.table[["Target"]] = NULL
  object@param.log[["SplitProbes"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                      , time = Sys.time()) 
  return(object)
}

#' Feature Differences Between Two Probe Groups 
#' 
#' Calculate the absolute numerical value differences of selected feature between two different probe groups.
#' @param object FISHprobe object.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param feature Feature used to compare. 
#' @param group_by Group of probes. Default is splitted probe groups ("Part").
#' @return Returns a FISHprobe object.
#' @importFrom pbapply pblapply 
#' @export
#' @examples \donttest{
#' object = DiffFeature(object, feature = "Tm")
#' }
DiffFeature <- function(object, transcript_ids = NULL, feature, group_by = "Part"){
  if(is.null(transcript_ids)) transcript_ids = names(object@probes[["Target"]])
  check_transcripts(object,transcript_ids)
  for(i in 1:length(transcript_ids)){
    if(!group_by %in% ListFeatures(object, transcript_ids = transcript_ids[i])[[1]]) stop(group_by, " groups cannot be found in ", transcript_ids[i], call. = FALSE)
    if(!feature %in% ListFeatures(object, transcript_ids = transcript_ids[i])[[1]]) stop(feature, " feature cannot be found in ", transcript_ids[i], call. = FALSE)
    groups = as.factor(object@probes[["Target"]][[transcript_ids[i]]][[group_by]])
    feature.use = object@probes[["Target"]][[transcript_ids[i]]][[feature]]
    feature.values = split(feature.use, groups)
    feature.diff = rep(abs(as.numeric(feature.values[[1]]) - as.numeric(feature.values[[2]])), each = 2)
    object@probes[["Target"]][[transcript_ids[i]]][[paste(group_by,feature,"Diff",sep = "_")]] = feature.diff
  }
  object@param.log[["DiffFeature"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                           , time = Sys.time()) 
  return(object)
}

#' Filter Probe Pairs
#' 
#' Filter pairs of splitted probes to remove singlet probes after filtering.
#' To make use of this feature \code{\link{FilterProbes}} has to be run first with \code{do.filter = FALSE}.
#' @param object FISHprobe object.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param filter_name Filter table feature in @filter.taable slot used for filtering.  
#' @return Returns a FISHprobe object.
#' @export
#' @examples \donttest{
#' object = FilterProbes(object,filter_names = "GC",max = 65,do.filter = FALSE)
#' object = FilterPairs(object, feature = "Tm")
#' }
#' 
FilterPairs <- function(object, transcript_ids = NULL, filter_name){
  if(is.null(transcript_ids)) transcript_ids = names(object@probes[["Target"]])
  check_transcripts(object, transcript_ids)
  for(i in 1:length(transcript_ids)){
    data.use = object@filter.table[["Target"]][[transcript_ids[i]]][, filter_name]
    if((length(data.use) %% 2) != 0) stop("Singlet probe found in ", transcript_ids[i], call. = FALSE)
    data.use = split(data.use, rep(1:(length(data.use)/2), each = 2))
    data.use = which(unlist(lapply(data.use, function(x) all(x == 1))))
    data.use = sort(c(data.use*2,data.use*2-1))
    object@probes[["Target"]][[transcript_ids[i]]] = object@probes[["Target"]][[transcript_ids[i]]][data.use,,drop = FALSE]
    object@blast.data[["Target"]][[i]] = object@blast.data[["Target"]][[i]][data.use]
    object@off.target.expr[["Target"]][[i]] = object@off.target.expr[["Target"]][[i]][data.use,,drop = FALSE]
    object@stash.probes[["Target"]][[i]] = object@stash.probes[["Target"]][[i]][data.use,,drop = FALSE]
    object@filter.table[["Target"]][[i]] = object@filter.table[["Target"]][[i]][data.use,,drop = FALSE]
  }
  object@param.log[["FilterPairs"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                           , time = Sys.time()) 
  return(object)
}

#' Order Probe Pairs
#' 
#' Order pairs of splitted probes based on feature values.
#' @param object FISHprobe object.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param feature Feature name used for ordering also see \code{\link{ListFeatures}}.
#' @param feature_target Targeting numerical values to order the probes before selection (order probes by the distance to this value). Default is none.
#' @param decrease Whether to order probes in decreasing order. Default is FALSE.
#' @return Returns a FISHprobe object.
#' @export
#' @examples \donttest{
#' object = OrderPairs(object, feature = "Tm")
#' }
#' 
OrderPairs <- function(object, transcript_ids = NULL, feature, feature_target = NULL,decrease = FALSE){
  if(is.null(transcript_ids)) transcript_ids = names(object@probes[["Target"]])
  check_transcripts(object, transcript_ids)
  if(length(feature) > 1) stop("Only one input feature is allowed.", call. = FALSE)
  if(!is.null(feature_target)){
    if(!is.numeric(feature_target)) stop("feature_target must be numerical.", call. = FALSE)
  }
  for(i in 1:length(transcript_ids)){
    data.use = object@probes[["Target"]][[transcript_ids[i]]][,feature]
    if(!is.null(feature_target)) data.use = abs(data.use - feature_target)
    if((length(data.use) %% 2) != 0) stop("Singlet probe found in ", transcript_ids[i], call. = FALSE)
    data.use = split(data.use, rep(1:(length(data.use)/2), each = 2))
    idx = order(unlist(lapply(data.use, mean)),decreasing = decrease)
    idx = unlist(lapply(idx, function(x) c(2*x-1,2*x)))
    object@probes[["Target"]][[transcript_ids[i]]] = object@probes[["Target"]][[transcript_ids[i]]][idx,,drop = FALSE]
    object@blast.data[["Target"]][[i]] = object@blast.data[["Target"]][[i]][idx]
    object@off.target.expr[["Target"]][[i]] = object@off.target.expr[["Target"]][[i]][idx,,drop = FALSE]
    object@stash.probes[["Target"]][[i]] = object@stash.probes[["Target"]][[i]][idx,,drop = FALSE]
    object@filter.table[["Target"]][[i]] = object@filter.table[["Target"]][[i]][idx,,drop = FALSE]
  }
  object@param.log[["OrderPairs"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                              , time = Sys.time()) 
  return(object)
}

#' Select Pair Probes
#' 
#' Select probe prairs after \code{\link{SplitProbes}}, similar as the implementation of \code{\link{SelectProbes}}.
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
#' object = PairSelect(object, n_space = 10)
#' }
#' 
PairSelect <- function(object, transcript_ids = NULL, TPM_cutoff = 15, n_space = 0, order_by_feature = NULL, feature_target = NULL, iterative = TRUE, recalculate = FALSE){
  slot = "Target"
  if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
  check_transcripts(object,transcript_ids)
  if(!is.null(order_by_feature)){
    if(is.null(feature_target)) stop("Please provide a numerical value for feature_target.", call. = FALSE)
    if(length(order_by_feature) > 1 | length(feature_target) > 1) stop("Only one input feature is allowed. (order_by_feature and feature_target)", call. = FALSE)
    if(check_any_feature(object, slot, order_by_feature)) stop("Feature not found.", call. = FALSE)
    if(!is.numeric(feature_target)) stop("feature_target must be numerical.", call. = FALSE)
    object = OrderPairs(object,transcript_ids,order_by_feature,feature_target)
  }
  if(recalculate) object@probes[[slot]][transcript_ids] = object@stash.probes[[slot]][transcript_ids]
  data.to.select = object@probes[[slot]][transcript_ids]
  if(!iterative){
    selected.probes = pblapply(1:length(data.to.select),function(i){
      seqs.selected = seq_select(data.to.select[[i]],TPM_cutoff,n_space,start_name = "Original_Start", seqname = 'Original')
      seqs.match = sort(match(rownames(seqs.selected), rownames(data.to.select[[i]])))
      seq.match.all = c()
      for(j in seqs.match){
        if(j %% 2 == 1){
          seq.match.all = c(seq.match.all, j, j+1)
        } else {
          seq.match.all = c(seq.match.all, j, j-1)
        }
      }
      seq.match.all = sort(unique(seq.match.all))
      data.to.select[[i]][seq.match.all,]
    })
  } else {
    selected.probes = lapply(1:length(data.to.select),function(i){
      seqs.selected = iter_select(data.to.select[[i]],TPM_cutoff,n_space,start_name = "Original_Start", stop_name = "Original_Stop", seqname = 'Original')
      seqs.match = sort(match(rownames(seqs.selected), rownames(data.to.select[[i]])))
      seq.match.all = c()
      for(j in seqs.match){
        if(j %% 2 == 1){
          seq.match.all = c(seq.match.all, j, j+1)
        } else {
          seq.match.all = c(seq.match.all, j, j-1)
        }
      }
      seq.match.all = sort(unique(seq.match.all))
      data.to.select[[i]][seq.match.all,]
    })
  }
  names(selected.probes) = transcript_ids
  object@probes[[slot]][transcript_ids] = selected.probes
  object@param.log[[paste("PairSelect",slot,sep = "-")]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                                            , time = Sys.time()) 
  return(object)
}
