# test_col1a1.R
#------------------------------------------------------------------------------------------------------------------------
library(TrenaGene)
library(RUnit)
library(TrenaGeneSkinData)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_basic()
   test_withExpressionData()
   test_getEnhancers()
   test_getEncodeDHS()
   test_getChipSeq()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_basic <- function()
{
   printf("--- test_basic")

   col1a1 <- TrenaGene("COL1A1", "hg38")
   checkTrue("TrenaGene" %in% is(col1a1))

   tbl <- getTranscriptsTable(col1a1)
   checkEquals(dim(tbl), c(14, 30))
   gene.min <- subset(tbl, moleculetype=="gene")$start
   gene.max <- subset(tbl, moleculetype=="gene")$endpos
   checkTrue(min(tbl$start) == gene.min)
   checkTrue(max(tbl$endpos) == gene.max)
   transcript.lengths <- sort(with(tbl, 1 + endpos-start))
   checkEquals(transcript.lengths, c(517, 574, 616, 867, 893, 903, 948, 1279, 2116, 2183, 2355, 3150, 18344, 18344))

   checkEquals(getExpressionMatrixNames(col1a1), list())

} # test_basic
#------------------------------------------------------------------------------------------------------------------------
test_withExpressionData <- function()
{
   printf("--- test_withExpressionData")

   dataDirectory <- system.file(package="TrenaGeneSkinData", "extdata")
   checkTrue(file.exists(dataDirectory))

   col1a1 <- TrenaGene("COL1A1", "hg38", expressionDataDirectory=dataDirectory)
   checkTrue("TrenaGene" %in% is(col1a1))

   matrix.names <- getExpressionMatrixNames(col1a1)
   checkTrue(length(matrix.names) >= 4)
   checkTrue("mtx.protectedAndExposed" %in% matrix.names)

} # test_withExpressionData
#------------------------------------------------------------------------------------------------------------------------
test_getEnhancers <- function()
{
   printf("--- test_getEnhancers")

   dataDirectory <- system.file(package="TrenaGeneSkinData", "extdata")

   col1a1 <- TrenaGene("COL1A1", "hg38", expressionDataDirectory=dataDirectory)
   checkTrue("TrenaGene" %in% is(col1a1))
   tbl.enhancers <- getEnhancers(col1a1)
   checkEquals(dim(tbl.enhancers), c(32, 6))
   checkEquals(colnames(tbl.enhancers), c("chrom", "start", "end", "type", "combinedScore", "geneSymbol"))

} # test_getEnhancers
#------------------------------------------------------------------------------------------------------------------------
test_getEncodeDHS <- function()
{
   printf("--- test_getEncodeDHS")

   dataDirectory <- system.file(package="TrenaGeneSkinData", "extdata")
   col1a1 <- TrenaGene("COL1A1", "hg38", expressionDataDirectory=dataDirectory)
   checkTrue("TrenaGene" %in% is(col1a1))

   tbl.dhs <- getEncodeDHS(col1a1)
   checkTrue(nrow(tbl.dhs) > 2000)
   checkEquals(colnames(tbl.dhs), c("chrom", "chromStart", "chromEnd", "count", "score"))

} # test_getEncodeDHS
#------------------------------------------------------------------------------------------------------------------------
test_getChipSeq <- function()
{
   printf("--- test_getChipSeq")

   dataDirectory <- system.file(package="TrenaGeneSkinData", "extdata")
   col1a1 <- TrenaGene("COL1A1", "hg38", expressionDataDirectory=dataDirectory)
   checkTrue("TrenaGene" %in% is(col1a1))

   tbl.enhancers <- getEnhancers(col1a1)
   tbl.strong <- subset(tbl.enhancers, combinedScore > 500)
   chrom <- tbl.strong$chrom[1]
   loc.min <- tbl.strong$start[1]
   loc.max <- tbl.strong$end[1]


   tbl.chipSeq <- getChipSeq(col1a1, "chr17", 50200750, 50201000)
   checkTrue(nrow(tbl.chipSeq) > 30)
   checkEquals(colnames(tbl.chipSeq),
               c("chr", "start", "end", "tf", "name", "peakStart", "peakEnd"))

} # test_getChipSeq
#------------------------------------------------------------------------------------------------------------------------
if(!interactive()) runTests()

