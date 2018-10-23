# test_psg1.R
#------------------------------------------------------------------------------------------------------------------------
library(TrenaGene)
library(RUnit)
library(TrenaGenePlacentaData)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_basic()
   test_withExpressionData()
   test_getEnhancers()
   test_getEncodeDHS()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_basic <- function()
{
   printf("--- test_basic")

   psg1 <- TrenaGene("PSG1", "hg38")
   checkTrue("TrenaGene" %in% is(psg1))

   tbl <- getTranscriptsTable(psg1)
   checkEquals(dim(tbl), c(12, 30))
   gene.min <- subset(tbl, moleculetype=="gene")$start
   gene.max <- subset(tbl, moleculetype=="gene")$endpos
   checkTrue(min(tbl$start) == gene.min)
   checkTrue(max(tbl$endpos) == gene.max)
   transcript.lengths <- sort(with(tbl, 1 + endpos-start))
   checkEquals(transcript.lengths, c(3084, 4180, 5579, 6624, 10003, 12475, 12514, 13177, 13189, 13235, 13245, 13359))

   checkEquals(getExpressionMatrixNames(psg1), list())

} # test_basic
#------------------------------------------------------------------------------------------------------------------------
test_withExpressionData <- function()
{
   printf("--- test_withExpressionData")

   dataDirectory <- system.file(package="TrenaGenePlacentaData", "extdata")
   checkTrue(file.exists(dataDirectory))

   psg1 <- TrenaGene("PSG1", "hg38", expressionDataDirectory=dataDirectory)
   checkTrue("TrenaGene" %in% is(psg1))

   matrix.names <- sort(getExpressionMatrixNames(psg1))
   checkTrue(length(matrix.names) >= 2)
   checkTrue("SRP094910" %in% matrix.names)

   checkTrue("FilteredCountData8282018" %in% matrix.names)
   checkTrue("FilteredLengthScaledTPM8282018" %in% matrix.names)
   checkTrue("SRP094910" %in% matrix.names)
   checkTrue("UnfilteredCountData8282018" %in% matrix.names)
   checkTrue("UnfilteredLengthScaledTPM8282018" %in% matrix.names)
   checkTrue("FilteredCountData8282018-vsn"   %in% matrix.names)
   checkTrue("FilteredLengthScaledTPM8282018-vsn" %in% matrix.names)

   mtx <- loadExpressionData(psg1, "FilteredLengthScaledTPM8282018")
   checkEquals(dim(mtx), c(16664, 112))

} # test_withExpressionData
#------------------------------------------------------------------------------------------------------------------------
test_getEnhancers <- function()
{
   printf("--- test_getEnhancers")

   dataDirectory <- system.file(package="TrenaGenePlacentaData", "extdata")

   psg1 <- TrenaGene("PSG1", "hg38", expressionDataDirectory=dataDirectory)
   checkTrue("TrenaGene" %in% is(psg1))
   tbl.enhancers <- getEnhancers(psg1)
   checkEquals(dim(tbl.enhancers), c(8, 6))
   checkEquals(colnames(tbl.enhancers), c("chrom", "start", "end", "type", "combinedScore", "geneSymbol"))

} # test_withExpressionData
#------------------------------------------------------------------------------------------------------------------------
test_getEncodeDHS <- function()
{
   printf("--- test_getEncodeDHS")

   dataDirectory <- system.file(package="TrenaGenePlacentaData", "extdata")
   psg1 <- TrenaGene("PSG1", "hg38", expressionDataDirectory=dataDirectory)
   checkTrue("TrenaGene" %in% is(psg1))

   tbl.dhs <- getEncodeDHS(psg1)
   checkTrue(nrow(tbl.dhs) > 30)   # 36 rows (10 oct 2018)
   checkEquals(colnames(tbl.dhs), c("chrom", "chromStart", "chromEnd", "count", "score"))

} # test_withExpressionData
#------------------------------------------------------------------------------------------------------------------------
if(!interactive()) runTests()

