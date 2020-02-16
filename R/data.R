#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Datasets
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Human Gene Names
#'
#' A data frame of human gene names.
#'
"humangenes"

#' Mouse Gene Names
#'
#' A data frame of mouse gene names.
#'
"mousegenes"

#' Human Gene Names and ENSEMBLE Gene ID Conversion GRCh38(hg38)
#'
#' A data frame of human gene names and their ENSEMBL gene IDs (GR).
#'
"human_ensemble"

#' Mouse Gene Names and ENSEMBLE Gene ID Conversion GRCh38(mm10)
#'
#' A data frame of mouse gene names and their ENSEMBL gene IDs.
#'
"mouse_ensemble"

#' Human ENSEMBLE Gene ID and Transcript ID GRCh38(hg38)
#'
#' A data frame of human ENSEMBL gene ID and transcript IDs.
#'
"transcripts_table_human"

#' Mouse ENSEMBLE Gene ID and Transcript ID GRCh38(mm10)
#'
#' A data frame of mouse ENSEMBL gene ID and transcript IDs.
#'
"transcripts_table_mouse"

#' Human BioMart Object GRCh38(hg38)
#'
#' A human biomart object. See \code{\link[biomaRt]{useMart}} for details
#'
"human"

#' Mouse BioMart Object GRCh38(mm10)
#'
#' A mouse biomart object. See \code{\link[biomaRt]{useMart}} for details
#'
"mouse"

#' Human Transcript Mean Tissue Specific Expression (GTEx Analysis V8)
#'
#' An expression matrix of human ENSEMBL transcript IDs and their mean GTEX tissue specific expressions in [log2TPM].
#' Rows are transcript IDs and Columns are GTEX tissues.
#' Check \href{https://gtexportal.org/home/datasets}{GTEX} datasets for more details.
#'
"gtex_mean_human"

#' Mouse Transcript Mean Tissue Specific Expression
#'
#' An expression matrix of mouse ENSEMBL transcript IDs and their mean ENCODE tissue specific RNA-seq expressions in [log2TPM].
#' Rows are transcript IDs and Columns are ENCODE tissues.
#' Check \href{https://www.encodeproject.org/matrix/?type=Experiment&status=released&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&assay_slims=Transcription}{ENCODE} datasets for more details.
#' Check meta data using \code{\link{gtex_mouse_meta}}
#' 
"gtex_mean_mouse"

#' Mouse ENCODE RNA-seq Metadata
#'
#' A data frame of mouse ENCODE RNA-seq data used for \code{\link{gtex_mean_mouse}}.
#' Check \href{https://www.encodeproject.org/matrix/?type=Experiment&status=released&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&assay_slims=Transcription}{ENCODE} datasets for more details.
#' 
"gtex_mouse_meta" 

#' Human Transcript Tissue Specific Expression Percentage (GTEx Analysis V8)
#'
#' A matrix of human transcript tissue expression percentage using cutoff of 0 logTPM as expressed.
#' Rows are transcript IDs and Columns are GTEX tissues.
#' Check \href{https://gtexportal.org/home/datasets}{GTEX} datasets for more details.
#'
"gtex_pct_human"

#' Mouse Transcript Tissue Specific Expression Percentage
#'
#' A matrix of mouse transcript tissue expression percentage using cutoff of 0 logTPM as expressed.
#' Rows are transcript IDs and Columns are ENCODE tissues.
#' Check \href{https://www.encodeproject.org/matrix/?type=Experiment&status=released&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&assay_slims=Transcription}{ENCODE} datasets for more details.
#' Check meta data using \code{\link{gtex_mouse_meta}}
#' 
"gtex_pct_mouse"

#' Pre-Computed Codebook
#'
#' A codebook object with 480 pre-computed 4-bit code in 32 hybridization rounds in 3 channels(488nm, 561nm, 647nm).
#'  
"codebook_v1"

#' Pre-Computed Readout Probe
#'
#' A set of 226 pre-computed 20mer readout probes
#'  
"readout_v1"
