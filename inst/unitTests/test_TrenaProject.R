# test_TrenaProject
#------------------------------------------------------------------------------------------------------------------------
library(TrenaProject)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   printf("--- test_ctor")

   genes <- c("TREM2", "INPP5D")
   genomeName <- "hg38"

   footprintDatabaseHost <- "khaleesi.systemsbiology.net"
   footprintDatabaseNames <- c("brain_hint_20, brain_wellington_16")

   expressionDirectory <- system.file(package="TrenaProject", "extdata", "expression")
   variantsDirectory <- system.file(package="TrenaProject", "extdata", "variants")
   covariatesFile <- system.file(package="TrenaProject", "extdata", "covariates", "dummyCovariates.RData")
   checkTrue(file.exists(expressionDirectory))
   checkTrue(file.exists(variantsDirectory))
   checkTrue(file.exists(covariatesFile))

   proj <- TrenaProject(supportedGenes=genes,
                        genomeName=genomeName,
                        footprintDatabaseHost=footprintDatabaseHost,
                        footprintDatabaseNames=footprintDatabaseNames,
                        expressionDirectory=expressionDirectory,
                        variantsDirectory=variantsDirectory,
                        covariatesFile=covariatesFile,
                        quiet=TRUE)

   checkEquals(getSupportedGenes(proj), genes)
   checkEquals(getFootprintDatabaseHost(proj), footprintDatabaseHost)
   checkEquals(getFootprintDatabaseNames(proj), footprintDatabaseNames)

   checkTrue(is.null(getTargetGene(proj)))
   setTargetGene(proj, genes[1])
   checkEquals(getTargetGene(proj), genes[1])

   tbl.transcripts <- getTranscriptsTable(proj)
   checkTrue(nrow(tbl.transcripts) >= 3)

   checkEquals(getExpressionMatrixNames(proj), c("dummyExpressionSet_1", "dummyExpressionSet_2"))

   checkTrue(is.matrix(getExpressionMatrix(proj, "dummyExpressionSet_1")))
   checkTrue(is.matrix(getExpressionMatrix(proj, "dummyExpressionSet_2")))

   expected <- c("someGene.region.vcf", "tbl.snp.gwas.minimal")
   file.list <- getVariantDatasetNames(proj)
   checkTrue(all(expected %in% names(file.list)))
   checkTrue(file.exists(file.list[["someGene.region.vcf"]]))

     # most variant files - other than vcfs - are serialized into .RData files,
     # with that suffix stripped off for human readers (in a presumed Shiny UI)

   checkTrue(file.exists(sprintf("%s.RData", file.list[["tbl.snp.gwas.minimal"]])))

   tbl.covariates <- getCovariatesTable(proj)

   tbl.enhancers <- getEnhancers(proj)
   checkEquals(colnames(tbl.enhancers), c("chrom", "start", "end", "type", "combinedScore", "geneSymbol"))

   checkTrue(nrow(tbl.enhancers) >= 5)
   checkEquals(unique(tbl.enhancers$geneSymbol), getTargetGene(proj))

   tbl.dhs <- getEncodeDHS(proj)
   checkEquals(colnames(tbl.dhs), c("chrom", "chromStart", "chromEnd", "count", "score"))
      # these open chromatin regions bear no necessary relation to the targetGene
      # instead, they span the region of the targetGene-associated enhancers.  check that
   loc.min <- min(tbl.enhancers$start)
   loc.max <- max(tbl.enhancers$end)
   chromosome <- unique(tbl.enhancers$chrom)
   checkEquals(length(chromosome), 1)
   checkTrue(all(tbl.dhs$chromStart >= loc.min))
   checkTrue(all(tbl.dhs$chromStart <= loc.max))
   checkTrue(all(tbl.dhs$chromEnd >= loc.min))
   checkTrue(all(tbl.dhs$chromEnd <= loc.max))
   checkTrue(all(tbl.dhs$chrom == chromosome))

   tbl.chipseq <- getChipSeq(proj, chrom=chromosome, start=loc.min, end=loc.max, tfs=NA)
   checkTrue(nrow(tbl.chipseq) > 2000)
   checkEquals(colnames(tbl.chipseq), c("chr", "start", "end", "tf", "name", "peakStart", "peakEnd"))

   checkEquals(getGeneRegion(proj),                     "chr6:41158506-41163186")
   checkEquals(getGeneRegion(proj, flankingPercent=20), "chr6:41157570-41164122")

   checkEquals(getGeneEnhancersRegion(proj),                     "chr6:41154324-41210533")
   checkEquals(getGeneEnhancersRegion(proj, flankingPercent=10), "chr6:41148703-41216154")

   browser()
   vf <- getVariantDatasetNames(proj)

} # test_ctor
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
