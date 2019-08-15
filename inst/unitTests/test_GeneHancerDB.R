# test_GeneHancerDB.R
#------------------------------------------------------------------------------------------------------------------------
library(TrenaProjectHG38)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("ghdb")){
   ghdb <- GeneHancerDB()
   } # creating for several tests below

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
   test_listTissues()
   test_foxo6()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))
   checkTrue("GeneHancerDB" %in% is(ghdb))

} # test_ctor
#------------------------------------------------------------------------------------------------------------------------
test_listTissues <- function()
{
   message(sprintf("--- test_listTissues"))
   tissues <- listTissues(ghdb)
   checkTrue(length(tissues) > 300)
   checkTrue(all(c("placenta", "Placenta") %in% tissues))

} # test_listTissues
#------------------------------------------------------------------------------------------------------------------------
# a simple test case, chosen for historical (and no longer important) reasons
test_foxo6 <- function()
{
   message(sprintf("--- test_foxo6"))

      # ENCODE and Ensembl use different capitalization
   tissues <-  c("placenta", "Placenta")
   tbl <- retrieveEnhancersFromDatabase(ghdb, "FOXO6", tissues)
   checkTrue(is.data.frame(tbl))
   checkEquals(nrow(tbl), 13)
   checkTrue(all(tissues %in% tbl$tissue))
   checkEquals(length(which(duplicated(tbl$sig))), 0)  # no duplicated regions

   tbl.all <- retrieveEnhancersFromDatabase(ghdb, "FOXO6", tissues="all")
   checkTrue(is.data.frame(tbl.all))
   checkEquals(nrow(tbl.all), 43)
   checkEquals(length(which(duplicated(tbl.all$sig))), 0)  # no duplicated regions

} # test_gata6
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
