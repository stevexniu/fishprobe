#' @include object.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Plot Filtering Binary Heatmap
#'
#' Plot a binary heatmap of filtering conditions for probe seb lection.
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param color Color used for heatmap. Default is red.
#' @param show_probe_names Whether to show probe names in the row of heatmap. Default is FALSE. 
#' Note if set to TRUE when the object contains many probes, the heatmap can be well-covered by the fonnts of the probes. 
#' @importFrom Ringo plotBM 
#' @export
#' @examples \donttest{
#' FilterHeatmap(object, cols_use = "blue")
#' }
FilterHeatmap <- function(object, slot = "Target", transcript_ids = NULL, color = "red", show_probe_names = FALSE){
  if(slot == "Readout"){
    plotBM(as.matrix(object@filter.table[[slot]]), boxCol = cols.use, main = "Readout Probe") 
  } else {
    if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
    if(length(transcript_ids) < 3){
      par(mfrow = c(1,length(transcript_ids)))
    } else {
      par(mfrow = c(ceiling(length(transcript_ids)/3),3))
    }
    for(trans in transcript_ids){
      data.to.plot = as.matrix(object@filter.table[[slot]][[trans]])
      if(ncol(data.to.plot) == 1) data.to.plot = cbind(data.to.plot, data.to.plot)
      colnames.use = colnames(data.to.plot)
      rownames.use = rownames(data.to.plot)
      dimnames(data.to.plot) = NULL
      plotBM(data.to.plot, boxCol = color, main = trans)
      if(show_probe_names) axis(side = 2, at = seq(nrow(data.to.plot)) - 0.5, labels = NULL, tick = FALSE, line = 0, las = 2)
      axis(side = 1, at = seq(ncol(data.to.plot)) - 0.5, labels = colnames.use,
           tick = FALSE, line = 0, las = 2, hadj = 0.6)
    }
    par(mfrow=c(1,1))
  }
}

#' Plot Filtering Barplot
#'
#' Plot a binary heatmap of filtering conditions for probe seb lection.
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @export
#' @examples \donttest{
#' FilterBarplot(object)
#' }
FilterBarplot <- function(object, slot = "Target", transcript_ids = NULL){
  if(slot == "Readout"){
    barplot(as.matrix(object@filter.table[[slot]]), ylab = "Counts", main = "Readout Probe") 
  } else {
    if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
    if(length(transcript_ids) < 3){
      par(mfrow = c(1,length(transcript_ids)))
    } else {
      par(mfrow = c(ceiling(length(transcript_ids)/3),3))
    }
    for(trans in transcript_ids){
      barplot(as.matrix(object@filter.table[[slot]][[trans]]), ylab = "Counts", main = trans)
    }
    par(mfrow=c(1,1))
  }
}

#' Print Number of Probes After Filtering 
#'
#' Print number of survival probes after different filtering conditions.
#' @param object FISHprobe object.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @return Returns numbers of probes survived after each filtering condition.
#' @export
#' @examples \donttest{
#' PrintFilter(object)
#' }
PrintFilter <- function(object, slot = "Target", transcript_ids = NULL){
  if(slot == "Readout"){
    numbers = list(Readout = base::colSums(object@filter.table[[slot]]))
  } else {
    if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
    numbers = lapply(transcript_ids, function(x){
      base::colSums(object@filter.table[[slot]][[x]])
    })
    names(numbers) = transcript_ids
  }
  return(numbers)
}

#' Barplot for Selected Features
#'
#' Generate barplots for selected features
#' @param object FISHprobe object.
#' @param features Features to plot on. See \code{\link{ListFeatures}} for features names for filtering.
#' @param slot Probe slot. Either Target or Readout. Default is Target probe set.
#' @param transcript_ids Transcript IDs. Default is NULL (using all the available transcripts in the object).
#' @param plot_type Plot type: "dot" for dot plot, "bar" for bar plot, "hist" for histogram. Default is dot plot.
#' @export
#' @examples \donttest{
#' PlotFeatures(object,"GC", plot_type = "bar")
#' }
PlotFeatures <- function(object, features, slot = "Target", transcript_ids = NULL, plot_type = "dot"){
  feature.all = ListFeatures(object, slot = slot, transcript_ids = transcript_ids)
  if(length(feature.all) == 0) stop("No features found", call. = FALSE)
  if(slot == "Readout"){
    if(!all(features %in% feature.all)) stop(paste("Cannot find features", features[which(!features %in% feature.all)],"in slot", slot),call. = FALSE) 
    if(length(features) < 3){
      par(mfrow = c(1,length(features)))
    } else {
      par(mfrow = c(ceiling(length(features)/3),3))
    }
    for(i in features){
      if(plot_type == "dot"){
        plot(object@probes[[slot]][[i]], ylab = "Values", main = "Readout Probes")
      } else if(plot_type == "bar"){
        barplot(object@probes[[slot]][[i]], ylab = "Values", main = "Readout Probes")
      } else if(plot_type == "hist"){
        hist(object@probes[[slot]][[i]], ylab = "Values", main = "Readout Probes")
      }
    }
  } else if(slot == "Target") {
    for(i in 1:length(feature.all)){
      if(!all(features %in% feature.all[[i]]))
        stop(paste("Cannot find features", features[which(!features %in% feature.all[[i]])],"in transcript", names(feature.all)[i]),call. = FALSE) 
    }
    if(is.null(transcript_ids)) transcript_ids = names(object@probes[[slot]])
    par(mfrow = c(length(transcript_ids),length(features)))
    for(i in transcript_ids){
      for(j in features){
        if(plot_type == "dot"){
          plot(object@probes[[slot]][[i]][[j]], ylab = "Values", main = paste(i,j))
        } else if(plot_type == "bar"){
          barplot(object@probes[[slot]][[i]][[j]], ylab = "Values", main = paste(i,j))
        } else if(plot_type == "hist"){
          hist(object@probes[[slot]][[i]][[j]], ylab = "Values", main = paste(i,j))
        }
      }
    }
  }
  par(mfrow = c(1, 1))
}
 