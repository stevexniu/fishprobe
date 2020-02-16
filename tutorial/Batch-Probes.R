library(FISHprobe)
genenames = c("H2-Eb1","Arg1","Adgre1","Pdgfra","Aif1","Cdh1","Mki67","Il1b") # Genes of interest
tissue = "mammary gland" # Tissue of interest
for(i in 1:length(genenames)){
  genename = genenames[i]
  print(genename)
  test = CreateProbeObject(genename,"mouse")
  test = TissueExpr(test)
  test = SelectTranscript(test,tissue = tissue)
  test = GetSequence(test)
  test = ExtractProbes(test,probe_length = 25,consecutive_repeats = 5)
  test = StashProbes(test, slot = "Target")  
  test = BLASTProbes(test)
  test = FilterBLAST(test, alignment_length = 11, E_value_cutoff = 1)
  test = OffTargetExpr(test,tissue = tissue)
  test = SelectProbes(test, n_space = -5, TPM_cutoff = 15, order_by_feature = "GC", feature_target = 0.55, iterative = T)
  test = GenerateProbes(test)
  test = CalcTm(test)
  test = CalcExon(test)
  test = Calc2ndStruct(test, method = "nupack", pseudoknots = TRUE, Na_M = 0.39, formamide_percent = 50)
  test = CalcDuplex(test, Na_M = 0.39, formamide_percent = 50)
  test = FilterProbes(test,filter_names = "Equilibrium_Percent",min = 0.2,max = 1,do.filter = TRUE)
  test = FilterProbes(test,filter_names = "Min_Free_Energy",min = -0.5,max = 0,do.filter = TRUE)
  test = FilterProbes(test,filter_names = "Consecutive",min = -Inf,max = 0.99,do.filter = TRUE)
  test = FilterProbes(test,filter_names = "GC",min = 0.45,max = 0.75,do.filter = TRUE)
  test = OrderProbe(test,order_by = "Off_Targets_Expr",decrease = FALSE)
  if(nrow(test@probes$Target[[1]]) <= 30){
    message("Not sufficient number of probes: ",nrow(test@probes$Target[[1]]))
    GenerateProbeTracks(test)
    saveRDS(test,file = paste(paste(genename,tissue,sep = "_"),".rds",sep = ""))
  } else {
    test = SubsetProbes(test,sequences = test@probes$Target[[1]]$Sequence[1:30])
    GenerateProbeTracks(test)
    saveRDS(test,file = paste(paste(genename,tissue,sep = "_"),".rds",sep = ""))
  }
}

############################################################
#  Merge Probes
############################################################
object.list = list()
object.names = list.files(pattern = ".rds")
for(i in 1:length(object.names)){
  object.list[[i]] = readRDS(object.names[i])
}
object.merge = MergeProbes(object.list)
object.merge = SetSlotData(object.merge, data = FISHprobe::codebook_v1, slotname = "codebook")
object.merge = GenerateReadout(object.merge, flank_5end ="GTGATTACAATTCGGCAGGCTT ", flank_3end = " CCCTATAGTGAGTCGTATTAGATCGAGT", concat = " ") # Add forward and reverse flank and concatenated nucleotides
object.merge = AssignReadout(object.merge, channel_transcript = list("561" = Transcripts(object.merge)))
object.merge = GenerateFinalProbes(object.merge, concat = " ")
SaveProbes(object.merge)
saveRDS(object.merge,file = "MergedProbes.rds")


