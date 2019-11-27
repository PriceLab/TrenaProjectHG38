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

   footprintDatabaseHost <- "khaleesi.systemsbiology.net"
   footprintDatabaseNames <- c("brain_hint_20, brain_wellington_16")

   packageDataDirectory <- system.file(package="TrenaProjectHG38", "extdata")

   trenaProj <- TrenaProjectHG38(projectName="HG38 test",
                                 supportedGenes=genes,
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
   test_getEnhancers()
   test_getPrimaryTranscriptInfo()
   test_getGeneRegulatoryRegions_varyingTissue()
   test_getGeneRegulatoryRegions_supplementcombineGeneHancerWithGenericPromoterOrProximalPromoter()

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

   checkEquals(getSupportedGenes(trenaProj), c("TREM2", "INPP5D"))
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

   checkEquals(getCovariateDatasetNames(trenaProj), list())


   checkEquals(getGeneRegion(trenaProj)$chromLocString,      "chr6:41158506-41163186")
   checkEquals(getGeneRegion(trenaProj, flankingPercent=20)$chromLocString, "chr6:41157570-41164122")
   checkEquals(getProximalPromoter(trenaProj, 0, 0)$chromLocString, "chr6:41163176-41163176")
   checkEquals(getProximalPromoter(trenaProj, 2500, 500)$chromLocString, "chr6:41162676-41165676")

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
                                 footprintDatabaseHost=footprintDatabaseHost,
                                 footprintDatabasePort=5433,
                                 footprintDatabaseNames=footprintDatabaseNames,
                                 packageDataDirectory=packageDataDirectory,
                                 quiet=TRUE)

   checkEquals(getFootprintDatabasePort(trenaProj), 5433)

} # test_ctor_withFootprintDatabsePortSpecified
#------------------------------------------------------------------------------------------------------------------------
test_getEnhancers <- function()
{
   message(sprintf("--- test_getEnhancers"))

   all.tissues <- getEnhancerTissues(trenaProj)
   setTargetGene(trenaProj, "TREM2")
   tbl.trem2.all <- getEnhancers(trenaProj)
   checkTrue(all(tbl.trem2.all$geneSymbol == "TREM2"))

   tbl.mef2c <- getEnhancers(trenaProj, "MEF2C")
   checkTrue(all(tbl.mef2c$geneSymbol == "MEF2C"))

   tbl.trem2.again <- getEnhancers(trenaProj, "TREM2")
   checkEquals(tbl.trem2.all, tbl.trem2.again)

   suppressWarnings(tbl.bogus <- getEnhancers(trenaProj, "bogus99"))
   checkEquals(nrow(tbl.bogus), 0)

} # test_getEnhancers
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
test_getGeneRegulatoryRegions_varyingTissue <- function()
{
   message(sprintf("--- test_getGeneRegulatoryRegions_varyingTissue"))
   all.tissues <- getEnhancerTissues(trenaProj);

   tbl <- getGeneRegulatoryRegions(trenaProj, "TREM2", "all") # tissues)
   checkTrue(nrow(tbl) >= 9)

   tbl.transcript <- getTranscriptsTable(trenaProj, targetGene="TREM2")[1,]
   checkEquals(tbl.transcript$strand, -1)
   checkEquals(tbl.transcript$tss, 41163176)

   brain.tissues <-  c("brain", "Brain")
   checkTrue(!any(brain.tissues %in% tbl$tissue))

   suppressWarnings({
      tbl.genericPromoter.10kb <- getGeneRegulatoryRegions(trenaProj, "TREM2", brain.tissues,
                                                           geneHancerMissing.promoter.upstream=5000,
                                                           geneHancerMissing.promoter.downstream=5000)
      with(tbl.genericPromoter.10kb, {
              checkEquals(start, 41163176 - 5000);
              checkEquals(end,   41163176 + 5000);
              })
      tbl.genericPromoter.small <- getGeneRegulatoryRegions(trenaProj, "TREM2", brain.tissues,
                                                            geneHancerMissing.promoter.upstream=10,
                                                            geneHancerMissing.promoter.downstream=3)
      with(tbl.genericPromoter.small, {
              checkEquals(start, 41163176 - 3);
              checkEquals(end,   41163176 + 10);
              })
       }) # suppressWarnings

      # INPP5D is on the positive strand.  check the fallback calculation for that case

   tbl.transcript <- getTranscriptsTable(trenaProj, targetGene="INPP5D")[1,]
   checkEquals(tbl.transcript$strand, 1)
   checkEquals(tbl.transcript$tss, 233060398)

   suppressWarnings({
      tbl <- getGeneRegulatoryRegions(trenaProj, "INPP5D", "bogus")
      checkEquals(nrow(tbl), 0)
      tbl <- getGeneRegulatoryRegions(trenaProj, "INPP5D", "bogus",
                                      geneHancerMissing.promoter.upstream=3,
                                      geneHancerMissing.promoter.downstream=10)
      with(tbl, {
          checkEquals(start, 233060398-3);
          checkEquals(end,   233060398+10)
          })
       }) # suppressWarnings

} # test_getGeneRegulatoryRegions_varyingTissue
#------------------------------------------------------------------------------------------------------------------------
# request the genehancer regions for a gene, in a tissue, where we know from prior examination
# that there ARE enhancers, but that they do not include anything like a traditional promoter around the TSS
# then make the request again, this time providing the +/- coordinates of that putative traditional promoter
test_getGeneRegulatoryRegions_supplementcombineGeneHancerWithGenericPromoterOrProximalPromoter <- function()
{
   message(sprintf("--- test_getGeneRegulatoryRegions_combineGeneHancerWithGenericPromoterOrProximalPromoter"))

     # first, a gene for which we have transcript info, but with no entry in GeneHancer
   tbl.transcript <- getTranscriptsTable(trenaProj, targetGene="TREM2")[1,]
   tss <- tbl.transcript$tss
   tbl <- getGeneRegulatoryRegions(trenaProj, "TREM2", "placenta")
   checkEquals(nrow(tbl), 2)

   tbl.ghSupplemented <- getGeneRegulatoryRegions(trenaProj, "TREM2", "placenta",
                                      geneHancerSupplemental.promoter.upstream=2000,
                                      geneHancerSupplemental.promoter.downstream=100)
   checkEquals(nrow(tbl.ghSupplemented), 3)

     # does the new row have the right coordinates?  tss +/- 2000, 100?
     # the new row should be first, since the tbl is sorted just before return
   checkEquals(tbl.ghSupplemented$start[1], (tss-100))
   checkEquals(tbl.ghSupplemented$end[1], (tss+2000))

     # the getTeneRegulatoryRegions method supports a second strategy to
     # supplement GeneHancer: if the target gene is unknown to GeneHancer
     # then the "geneHancerMissing.promoter.[up|down]Stream coordinates
     # are used to create a region around the TSS - which we expect the
     # user will want to be something like +/- 5kb
   suppressWarnings(
      checkEquals(nrow(getGeneRegulatoryRegions(trenaProj, "TREM2", "ear")), 0)
      )

   suppressWarnings(
       tbl.ghMissing <- getGeneRegulatoryRegions(trenaProj, "TREM2", "ear",
                                             geneHancerMissing.promoter.upstream=6000,
                                             geneHancerMissing.promoter.downstream=4000)
       )
   checkEquals(nrow(tbl.ghMissing), 1)
   checkEquals(tbl.ghMissing$start, tss-4000)
   checkEquals(tbl.ghMissing$end, tss+6000)
   checkEquals(tbl.ghMissing$source, "TrenaProjectHG38.gh.missing")

      # the geneHancerMising and geneHancerSupplemental arguments
      # should never be used at the same time.   ensure that this is so

   suppressWarnings(
       tbl.ghMissing.2 <- getGeneRegulatoryRegions(trenaProj, "TREM2", "ear",
                                             geneHancerMissing.promoter.upstream=6000,
                                             geneHancerMissing.promoter.downstream=4000,
                                             geneHancerSupplemental.promoter.upstream=2000,
                                             geneHancerSupplemental.promoter.downstream=200)
       )
   checkEquals(tbl.ghMissing, tbl.ghMissing.2)

      # now do the complementary check: gh regions present, sup & missing both
      # requested, only supplemental should be used

   tbl.ghSupplemented.2 <- getGeneRegulatoryRegions(trenaProj, "TREM2", "placenta",
                                      geneHancerSupplemental.promoter.upstream=2000,
                                      geneHancerSupplemental.promoter.downstream=100,
                                      geneHancerMissing.promoter.upstream=6000,
                                      geneHancerMissing.promoter.downstream=4000)

   checkEquals(tbl.ghSupplemented, tbl.ghSupplemented.2)


} # test_getGeneRegulatoryRegions_combineGeneHancerWithGenericPromoterOrProximalPromoter
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
