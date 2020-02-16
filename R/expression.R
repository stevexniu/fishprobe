#' @include utils.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get Tissue Expression of Targeting Genes
#'
#' This function allows you to get tissue specific expression of your genes of interests.
#' 
#' For Mouse genes, the expression profile is ENCODE RNA-seq \code{\link{gtex_mean_mouse}}.
#' For Human genes, the expression profile is GTEX RNA-seq \code{\link{gtex_mean_human}}.
#' @param object FISHprobe object.
#' @param genes Which gene to use. Default is all the genes in the object.
#' @return Returns a FISHprobe object.
#' @importFrom reshape2 melt
#' @export
#' @examples \donttest{
#' TissueExpr(object)
#' }
TissueExpr <- function(object, genes = NULL){
  if(is.null(genes)) genes = object@genes
  if (object@specie == "mouse"){
    for(gene in genes){
      data.expr = get_gtex(genename = gene,data_mean = FISHprobe::gtex_mean_mouse,data_pct = FISHprobe::gtex_pct_mouse,ensemble = FISHprobe::mouse_ensemble,gene_table = FISHprobe::transcripts_table_mouse,species = "mouse",plot = plot)
      object@expr.data[[gene]] = data.expr
      object@expr.data[[gene]][,1] = as.character(object@expr.data[[gene]][,1])
      object@expr.data[[gene]][,2] = as.character(object@expr.data[[gene]][,2])
      object@transcripts[[gene]] = as.character(unique(data.expr[["Transcript"]]))
    }
  } else if (object@specie == "human"){
    for(gene in object@genes){
      data.expr = get_gtex(genename = gene,data_mean = FISHprobe::gtex_mean_human,data_pct = FISHprobe::gtex_pct_human,ensemble = FISHprobe::human_ensemble,gene_table = FISHprobe::transcripts_table_human,species = "human",plot = plot)
      object@expr.data[[gene]] = data.expr
      object@expr.data[[gene]][,1] = as.character(object@expr.data[[gene]][,1])
      object@expr.data[[gene]][,2] = as.character(object@expr.data[[gene]][,2])
      object@transcripts[[gene]] = as.character(unique(data.expr[["Transcript"]]))
    }
  }
  object@param.log[["TissueExpr"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                      , time = Sys.time()) 
  return(object)
}

#' Plot Tissue Expression of Targeting Genes
#'
#' Plot tissue expressions data after calling \emph{TissueExpr}.
#' @param object FISHprobe object.
#' @param gene_plot Genes of interest to plot. Default is all the genes in the object.
#' @param return_plot Return ggplot object. Default is FALSE.
#' @import ggplot2
#' @import cowplot
#' @export
#' @examples \donttest{
#' PlotTissueExpr(object)
#' PlotTissueExpr(object, gene_plot = "ACTB")
#' }
PlotTissueExpr <- function(object, gene_plot = NULL, return_plot = FALSE){
  if(is.null(gene_plot)) gene_plot = object@genes
  plot.list = list()
  for(gene in gene_plot){
    plot.list[[gene]] = ggplot(data = object@expr.data[[gene]], mapping = aes(x = Transcript, y = Tissue)) + geom_point(mapping = aes(size = Percent_Expr, color = Log_TPM)) +
      scale_color_distiller(palette = "Spectral") + xlab("") + ggtitle(gene) + ylab("") + scale_size_continuous(range = c(1,5),breaks = seq(0,100,by = 25)) +
      theme_cowplot() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) 
  }
  if(return_plot){
    return(plot_grid(plotlist = plot.list))
  } else {
    plot_grid(plotlist = plot.list)
  }
}

#' List All Possible Tissues
#'
#' List All possible tissues in the object.
#' @param object FISHprobe object.
#' @param pattern Feature patterns to match, examples include Brain, Colon for human etc. Default is none.

#' @return Returns a list of possible tissues.
#' @export
#' @examples \donttest{
#' ListTissues(object)
#' ListTissues(object, pattern = "Brain")
#' }
ListTissues <- function(object, pattern = ""){
  if(object@specie == "human"){
    return(grep(pattern, colnames(FISHprobe::gtex_mean_human), value = TRUE))
  } else {
    return(grep(pattern, colnames(FISHprobe::gtex_mean_mouse), value = TRUE))
  }
}

#' Add Additional Gene Expression Data
#' 
#' Add gene expression data for the genes/transcripts used in the object.
#' @param object FISHprobe object.
#' @param expr_data Expression matrix with rownames as gene/transcripts names.
#' @param name_match Names used to match the expr_data with object@sequences slot data. Default is GeneName.
#' @return Returns a FISHprobe object.
#' @export
#' @examples \donttest{
#' AddGeneExpressions(object, expression_data, name_match = "ensembl_transcript_id")
#' }
AddGeneExpressions <- function(object, expr_data, name_match = "GeneName"){
  if(!all(name_match %in% object@sequences))  stop(paste("Cannot find names to match", paste(name_match[which(!name_match %in% colnames(object@sequences))],collapse = " ")),call. = FALSE) 
  object@sequences[["expression"]] = NULL
  object@sequences[,"expression"] = expr_data[match(rownames(expr_data),object@sequences[[name_match]]),]
  return(object)
  object@param.log[["AddExpr"]] = list(call = do.call(cbind, mget(names(formals()),sys.frame(sys.nframe()))[-1])
                                                           , time = Sys.time())
}
