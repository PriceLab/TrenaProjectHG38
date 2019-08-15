# test_TrenaProjectHG38
#------------------------------------------------------------------------------------------------------------------------
library(TrenaProjectHG38)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
# TODO: test getEnhancers and getEncodeDHS with genes which have no GeneHancer regions (pshannon, 10may 2019)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("trenaProj")){

   genes <- c("TREM2", "INPP5D")
   genomeName <- "hg38"
   geneInfoTable.path <- system.file(package="TrenaProjectHG38", "extdata", "geneInfoTable.RData")
   stopifnot(file.exists(geneInfoTable.path))

   footprintDatabaseHost <- "khaleesi.systemsbiology.net"
   footprintDatabaseNames <- c("brain_hint_20, brain_wellington_16")

   packageDataDirectory <- system.file(package="TrenaProjectHG38", "extdata")

   trenaProj <- TrenaProjectHG38(projectName="HG38 test",
                                 supportedGenes=genes,
                                 geneInfoTable.path=geneInfoTable.path,
                                 footprintDatabaseHost=footprintDatabaseHost,
                                 footprintDatabaseNames=footprintDatabaseNames,
                                 packageDataDirectory=packageDataDirectory,
                                 quiet=TRUE)
   } # creating trenaProj for use in multiple functions below

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
   test_ctor_withFootprintDatabasePortSpecified()
   # test_getEnhancers()
   test_getPrimaryTranscriptInfo()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))

   tbl.geneInfo <- getGeneInfoTable(trenaProj)
   checkTrue(nrow(tbl.geneInfo) > 50000)
   geneInfo.columns <- c( "ensg", "chrom", "start", "end", "tss", "strand", "geneSymbol",
                         "entrez", "appris", "tsl", "transcript", "type")
   checkEquals(colnames(tbl.geneInfo), geneInfo.columns)

   checkEquals(getSupportedGenes(trenaProj), character(0))
   checkEquals(getFootprintDatabaseHost(trenaProj), footprintDatabaseHost)
   checkEquals(getFootprintDatabaseNames(trenaProj), footprintDatabaseNames)

   message(sprintf("--- testing get/setTargetGene"))
   setTargetGene(trenaProj, genes[1], curatedGenesOnly=FALSE)
   checkEquals(getTargetGene(trenaProj), genes[1])

   message(sprintf("--- getting transcript info for %s", genes[1]))

   tbl.transcripts <- getTranscriptsTable(trenaProj)
   checkTrue(nrow(tbl.transcripts) == 1)
   checkEquals(colnames(tbl.transcripts), geneInfo.columns)

   checkEquals(tbl.transcripts$geneSymbol, genes[1])

   message(sprintf("--- testing get/setTargetGene"))
   setTargetGene(trenaProj, "PIGF")           # a placental gene, not in the IGAP project
   checkEquals(getTargetGene(trenaProj), "PIGF")

   message(sprintf("--- getting transcript info for %s", "PIGF"))
   tbl.transcripts <- getTranscriptsTable(trenaProj)
   checkTrue(nrow(tbl.transcripts) == 1)
   checkEquals(tbl.transcripts$geneSymbol, "PIGF")

     # return to TREM2, whose coordinates we check below
   setTargetGene(trenaProj, genes[1])

   checkEquals(getExpressionMatrixNames(trenaProj), c("dummyExpressionSet_1", "dummyExpressionSet_2"))

   expected <- c("someGene.region.vcf", "tbl.snp.gwas.minimal")
   file.list <- getVariantDatasetNames(trenaProj)
   checkTrue(all(expected %in% file.list))
   #checkTrue(file.exists("someGene.region.vcf"))

     # most variant files - other than vcfs - are serialized into .RData files,
     # with that suffix stripped off for human readers (in a presumed Shiny UI)

   #checkTrue(file.exists(sprintf("%s.RData", file.list[["tbl.snp.gwas.minimal"]])))

   checkEquals(getCovariateDatasetNames(trenaProj), list())

   # tbl.enhancers <- getEnhancers(trenaProj)
   # checkEquals(colnames(tbl.enhancers), c("chrom", "start", "end", "type", "combinedScore", "geneSymbol"))

   # checkTrue(nrow(tbl.enhancers) >= 5)
   # checkEquals(unique(tbl.enhancers$geneSymbol), getTargetGene(trenaProj))

   #tbl.dhs <- getEncodeDHS(trenaProj)
   #checkEquals(colnames(tbl.dhs), c("chrom", "chromStart", "chromEnd", "count", "score"))
      # these open chromatin regions bear no necessary relation to the targetGene
      # instead, they span the region of the targetGene-associated enhancers.  check that
   #loc.min <- min(tbl.enhancers$start)
   #loc.max <- max(tbl.enhancers$end)
   #chromosome <- unique(tbl.enhancers$chrom)
   #checkEquals(length(chromosome), 1)
   #checkTrue(all(tbl.dhs$chromStart >= loc.min))
   #checkTrue(all(tbl.dhs$chromStart <= loc.max))
   #checkTrue(all(tbl.dhs$chromEnd >= loc.min))
   #checkTrue(all(tbl.dhs$chromEnd <= loc.max))
   #checkTrue(all(tbl.dhs$chrom == chromosome))

   #tbl.chipseq <- getChipSeq(trenaProj, chrom=chromosome, start=loc.min, end=loc.max, tfs=NA)
   #checkTrue(nrow(tbl.chipseq) > 2000)

   #checkEquals(colnames(tbl.chipseq), c("chrom", "start", "endpos", "tf", "name", "strand", "peakStart", "peakEnd"))

   checkEquals(getGeneRegion(trenaProj)$chromLocString,      "chr6:41158506-41163186")
   checkEquals(getGeneRegion(trenaProj, flankingPercent=20)$chromLocString, "chr6:41157570-41164122")

   #checkEquals(getGeneEnhancersRegion(trenaProj)$chromLocString,                     "chr6:41154324-41210533")
   #checkEquals(getGeneEnhancersRegion(trenaProj, flankingPercent=10)$chromLocString, "chr6:41148703-41216154")

   vf <- getVariantDatasetNames(trenaProj)

   checkTrue(nrow(getGeneInfoTable(trenaProj)) > 15000)
   checkTrue(ncol(getGeneInfoTable(trenaProj)) >= 10)

   checkEquals(getFootprintDatabasePort(trenaProj), 5432)


   checkTrue(!recognizedGene(trenaProj, "bogusGene"))      # only genes in the tbl.geneInfo are recognized

} # test_ctor
#------------------------------------------------------------------------------------------------------------------------
test_ctor_withFootprintDatabasePortSpecified <- function()
{
   message(sprintf("--- test_ctor_withFootprintDatabasePortSpecified"))

   trenaProj <- TrenaProjectHG38(projectName="HG38 test",
                                 supportedGenes=genes,
                                 geneInfoTable.path=geneInfoTable.path,
                                 footprintDatabaseHost=footprintDatabaseHost,
                                 footprintDatabasePort=5433,
                                 footprintDatabaseNames=footprintDatabaseNames,
                                 packageDataDirectory=packageDataDirectory,
                                 quiet=TRUE)

   checkEquals(getFootprintDatabasePort(trenaProj), 5433)

} # test_ctor_withFootprintDatabsePortSpecified
#------------------------------------------------------------------------------------------------------------------------
# test_getEnhancers <- function()
# {
#    message(sprintf("--- test_getEnhancers"))
#
#    setTargetGene(trenaProj, "TREM2")
#    tbl.trem2 <- getEnhancers(trenaProj)
#    checkTrue(all(tbl.trem2$geneSymbol == "TREM2"))
#
#    tbl.mef2c <- getEnhancers(trenaProj, "MEF2C")
#    checkTrue(all(tbl.mef2c$geneSymbol == "MEF2C"))
#
#    tbl.trem2.again <- getEnhancers(trenaProj, "TREM2")
#    checkEquals(tbl.trem2, tbl.trem2.again)
#
#    tbl.bogus <- getEnhancers(trenaProj, "bogus99")
#    checkEquals(nrow(tbl.bogus), 0)
#
# } # test_getEnhancers
#------------------------------------------------------------------------------------------------------------------------
test_getPrimaryTranscriptInfo <- function()
{
   message(sprintf("--- test_getPrimaryTranscriptInfo"))

   checkEquals(getTranscriptsTable(trenaProj, "CRH")$tss, 66178725)
   # checkEquals(getPrimaryTranscriptInfo(trenaProj, "CRH")$tss, 66178725)

   setTargetGene(trenaProj, "TREM2")
   checkEquals(getTranscriptsTable(trenaProj)$tss, 41163176)
   # checkEquals(getPrimaryTranscriptInfo(trenaProj)$tss, 41163176)

} # test_getPrimaryTranscriptInfo
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
