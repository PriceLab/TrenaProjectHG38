#' @import trena
#' @importFrom DBI dbConnect dbListTables dbGetQuery dbListConnections dbDisconnect
#' @importFrom RPostgreSQL dbConnect dbListTables dbGetQuery dbListConnections dbDisconnect
#' @import RPostgreSQL
#' @import org.Hs.eg.db
#' @importFrom methods new
#'
#' @title TrenaProjectHG38
#------------------------------------------------------------------------------------------------------------------------
#' @name TrenaProjectHG38-class
#' @rdname TrenaProjectHG38-class
#' @aliases TrenaProjectHG38
#'
## @import methods

.TrenaProjectHG38 <- setClass ("TrenaProjectHG38",
                        representation = representation(
                           supportedGenes="character",
                           geneInfoTable="data.frame",
                           genomeName="character",
                           footprintDatabaseHost="character",
                           footprintDatabaseNames="character",
                           footprintDatabasePort="numeric",
                           expressionDirectory="character",
                           variantsDirectory="character",
                           covariatesFile="character",
                           state="environment",
                           quiet="logical"
                           )
                         )


#------------------------------------------------------------------------------------------------------------------------
setGeneric('getSupportedGenes',         signature='obj', function(obj) standardGeneric('getSupportedGenes'))
setGeneric('setTargetGene',             signature='obj', function(obj, targetGene, curatedGenesOnly=FALSE)
              standardGeneric('setTargetGene'))
setGeneric('getGenome',                 signature='obj', function(obj) standardGeneric('getGenome'))
setGeneric('getTargetGene',             signature='obj', function(obj) standardGeneric('getTargetGene'))
setGeneric('getGeneInfoTable',          signature='obj', function(obj) standardGeneric('getGeneInfoTable'))
setGeneric('getFootprintDatabaseHost',  signature='obj', function(obj) standardGeneric ('getFootprintDatabaseHost'))
setGeneric('getFootprintDatabasePort',  signature='obj', function(obj) standardGeneric ('getFootprintDatabasePort'))
setGeneric('getFootprintDatabaseNames', signature='obj', function(obj) standardGeneric ('getFootprintDatabaseNames'))
setGeneric('getTranscriptsTable',       signature='obj', function(obj) standardGeneric ('getTranscriptsTable'))
setGeneric('getPrimaryTranscriptInfo',  signature='obj', function(obj, targetGene=NA) standardGeneric ('getPrimaryTranscriptInfo'))
setGeneric('getExpressionMatrixNames',  signature='obj', function(obj) standardGeneric ('getExpressionMatrixNames'))
setGeneric('getExpressionMatrix',       signature='obj', function(obj, matrixName) standardGeneric ('getExpressionMatrix'))
setGeneric('getVariantDatasetNames',    signature='obj', function(obj) standardGeneric ('getVariantDatasetNames'))
setGeneric('getVariantDataset',         signature='obj', function(obj, datasetName) standardGeneric ('getVariantDataset'))
setGeneric('getEnhancers',              signature='obj', function(obj, targetGene=NA) standardGeneric ('getEnhancers'))
setGeneric('getEncodeDHS',              signature='obj', function(obj, targetGene) standardGeneric ('getEncodeDHS'))
setGeneric('getChipSeq',                signature='obj', function(obj, chrom, start, end, tfs=NA) standardGeneric ('getChipSeq'))
setGeneric('getCovariatesTable',        signature='obj', function(obj) standardGeneric ('getCovariatesTable'))
setGeneric('getGeneRegion',             signature='obj', function(obj, flankingPercent=0) standardGeneric ('getGeneRegion'))
setGeneric('getGeneEnhancersRegion',    signature='obj', function(obj, flankingPercent=0) standardGeneric ('getGeneEnhancersRegion'))
setGeneric('recognizedGene',            signature='obj', function(obj, geneName) standardGeneric ('recognizedGene'))
#------------------------------------------------------------------------------------------------------------------------
#' Define an object of class Trena
#'
#' @description
#' TrenaProjectHG38 and its (projected) subclasses provide convenient containers in which to collect
#'  trena-related aggregation of a gene's (a hierarchy of classes) including expression data,
#' transcript and variant info, genomic and epigenomic context, trena models and/or the means to create them
#'
#' @rdname TrenaProjectHG38-class
#'
#' @param supportedGenes a vector of character strings
#' @param footprintDatabaseHost Character string (e.g., "khaleesi.systemsbiology.net")
#' @param footprintDatabaseNames Character string (e.g., "hint_brain_20")
#' @param expressionDirectory A string pointing to a collection of RData expression matrices
#' @param variantsDirectory A string pointing to a collection of RData variant files
#' @param covariatesFile  the (optional) name of a covariates files
#' @param quiet A logical indicating whether or not the Trena object should print output
#'
#' @return An object of the TrenaProjectHG38 class
#'
#' @export
#'
#'
TrenaProjectHG38 <- function(supportedGenes,
                             geneInfoTable.path,
                             footprintDatabaseHost,
                             footprintDatabaseNames,
                             footprintDatabasePort=5432,
                             expressionDirectory,
                             variantsDirectory,
                             covariatesFile,
                             quiet)
{

   state <- new.env(parent=emptyenv())
      # gene-specific information, freshly assigned with every call to setTargetGene
   state$targetGene <- NULL
   state$tbl.transcripts <- NULL

   geneInfoTable.filepath <- system.file(package="TrenaProjectHG38", "extdata", "geneInfoTable_hg38.RData")


   tbl.name <-load(geneInfoTable.filepath)
   stopifnot(tbl.name == "tbl.geneInfo")

   .TrenaProjectHG38(supportedGenes=supportedGenes,
                 geneInfoTable=tbl.geneInfo,
                 genomeName="hg38",
                 footprintDatabaseHost=footprintDatabaseHost,
                 footprintDatabaseNames=footprintDatabaseNames,
                 footprintDatabasePort=footprintDatabasePort,
                 expressionDirectory=expressionDirectory,
                 variantsDirectory=variantsDirectory,
                 covariatesFile=covariatesFile,
                 state=state,
                 quiet=quiet)


} # ctor
#------------------------------------------------------------------------------------------------------------------------
#' get a summary of the objects
#'
#' @rdname show
#' @aliases show
#'
#' @param obj An object of class TrenaProjectHG38
#'
#' @export

setMethod('show', 'TrenaProjectHG38',

    function(object) {
       msg <- sprintf ('--- TrenaProjectHG38')
       cat(msg, '\n', sep='')

       msg <- sprintf("expression matrix directory: %s", object@expressionDirectory)
       cat(msg, "\n", sep='')
       })
#------------------------------------------------------------------------------------------------------------------------
#' get the list of genes supported in this project
#'
#' @rdname getSupportedGenes
#' @aliases getSupportedGenes
#'
#' @param obj An object of class TrenaProjectHG38
#'
#' @export

setMethod('getSupportedGenes', 'TrenaProjectHG38',

   function(obj) {
      obj@supportedGenes
      })

#------------------------------------------------------------------------------------------------------------------------
#' Set a single gene for analysis
#'
#' @rdname setTargetGene
#' @aliases setTargetGene
#'
#' @param obj An object of class TrenaProjectHG38
#' @param targetGene a characteor string, the name of the gene
#'
#' @export

setMethod('setTargetGene', 'TrenaProjectHG38',

   function(obj, targetGene, curatedGenesOnly=FALSE) {
      if(curatedGenesOnly){
         if(!all(is.na(getSupportedGenes(obj))))
            stopifnot(targetGene %in% getSupportedGenes(obj))
         }
      obj@state$targetGene <- targetGene
      xyz <- "about to set tbl.transcripts"
      targetGene.regex <- sprintf("^%s$", targetGene)
      index <- grep(toupper(targetGene.regex), toupper(obj@geneInfoTable$geneSymbol))
      if(length(index) == 0)
         return(FALSE)
      tbl.tmp <- obj@geneInfoTable[index,]
      obj@state$tbl.transcripts  <- tbl.tmp
      return(TRUE)
      })

#------------------------------------------------------------------------------------------------------------------------
#' get the single gene currently designated for analysis
#'
#' @rdname getTargetGene
#' @aliases getTargetGene
#'
#' @param obj An object of class TrenaProjectHG38
#'
#' @export

setMethod('getTargetGene', 'TrenaProjectHG38',

   function(obj) {
      obj@state$targetGene
      })

#------------------------------------------------------------------------------------------------------------------------
#' get the standard short name (e.g., "hg38", "mm10") of the project under study
#'
#' @rdname getGenome
#' @aliases getGenome
#'
#' @param obj An object of class TrenaProjectHG38
#'
#' @export

setMethod('getGenome', 'TrenaProjectHG38',

   function(obj) {
      obj@genomeName
      })

#------------------------------------------------------------------------------------------------------------------------
#' get the host on which the footprints Postgres database is running
#'
#' @rdname getFootprintDatabaseHost
#' @aliases getFootprintDatabaseHost
#'
#' @param obj An object of class TrenaProjectHG38
#'
#' @export

setMethod('getFootprintDatabaseHost', 'TrenaProjectHG38',

   function(obj) {
      obj@footprintDatabaseHost
      })

#------------------------------------------------------------------------------------------------------------------------
#' get the port  on which the footprints Postgres database is running
#'
#' @rdname getFootprintDatabasePort
#' @aliases getFootprintDatabasePort
#'
#' @param obj An object of class TrenaProjectHG38
#'
#' @export

setMethod('getFootprintDatabasePort', 'TrenaProjectHG38',

   function(obj) {
      obj@footprintDatabasePort
      })

#------------------------------------------------------------------------------------------------------------------------
#' get the names of the database tables holding footprints
#'
#' @rdname getFootprintDatabaseNames
#' @aliases getFootprintDatabaseNames
#'
#' @param obj An object of class TrenaProjectHG38
#'
#' @export

setMethod('getFootprintDatabaseNames', 'TrenaProjectHG38',

   function(obj) {
      obj@footprintDatabaseNames
      })

#------------------------------------------------------------------------------------------------------------------------
#' Get ENSEMBL-derived genomic information on what they judget to be the primary transcript of this gene
#'
#' @rdname getPrimaryTranscriptInfo
#' @aliases getPrimarytTranscriptInfo
#'
#' @param obj An object of class TrenaProjectHG38
#' @param targetGene default NA, in which case the previously assigned targetGene is assumed
#'
#' @seealso setTargetGene

#' @export
setMethod('getPrimaryTranscriptInfo',  'TrenaProjectHG38',

   function(obj, targetGene=NA){

      if(!is.na(targetGene)){
         return(subset(getGeneInfoTable(obj), geneSymbol==targetGene))
         }

      targetGene <- getTargetGene(obj)
      if(is.null(targetGene)){
         message("no targetGene set for this project, none supplied as argument to this function")
         return(data.frame())
         }

      return(subset(getGeneInfoTable(obj), geneSymbol==targetGene))
      })

#------------------------------------------------------------------------------------------------------------------------
#' Get the transcripts for the gene
#'
#' @rdname getTranscriptsTable
#' @aliases getTranscriptsTable
#'
#' @param obj An object of class TrenaProjectHG38
#'
#' @export

setMethod('getTranscriptsTable',  'TrenaProjectHG38',

   function(obj) {
      return(obj@state$tbl.transcripts)
      })

#------------------------------------------------------------------------------------------------------------------------
#' Get the names of the expression matrices - their names with directory and .RData suffix stripped out
#'
#' @rdname getExpressionMatrixNames
#' @aliases getExpressionMatrixNames
#'
#' @param obj An object of class TrenaProjectHG38
#'
#' @export

setMethod('getExpressionMatrixNames',  'TrenaProjectHG38',

   function(obj) {
      if(is.na(obj@expressionDirectory))
         return(list())
      all.files <- list.files(obj@expressionDirectory)
      rdata.filenames <- grep(".RData$", all.files, value=TRUE)
      sub(".RData", "", rdata.filenames, fixed=TRUE)
      })

#------------------------------------------------------------------------------------------------------------------------
#' Get the a specifc expression matrix
#'
#' @rdname getExpressionMatrix
#' @aliases getExpressionMatrix
#'
#' @param obj An object of class TrenaProjectHG38
#' @param matrixName A numeric matrix
#'
#' @export

setMethod('getExpressionMatrix',  'TrenaProjectHG38',

    function(obj, matrixName){
       if(is.na(obj@expressionDirectory)){
          return(NA)
          }
       if(!matrixName %in% getExpressionMatrixNames(obj)){
          return(NA)
          }
       filename <- sprintf("%s.RData", matrixName)
       full.path <- file.path(obj@expressionDirectory, filename)
       stopifnot(file.exists(full.path))
       mtx <- NULL
       eval(parse(text=paste("mtx <- ", load(full.path))))
       invisible(mtx)
        })

#------------------------------------------------------------------------------------------------------------------------
#' List the RData files in the variants directory
#'
#' @rdname getVariantDatasetNames
#' @aliases getVariantDatasetNames
#'
#' @param obj An object of class TrenaProjectHG38
#'
#' @export

setMethod('getVariantDatasetNames', 'TrenaProjectHG38',

      function(obj){
          if(obj@variantsDirectory == "/dev/null")
             return(list())
          filenames <- sub(".RData", "", list.files(obj@variantsDirectory), fixed=TRUE)
          #full.paths <- file.path(obj@variantsDirectory, filenames)
          #names(full.paths) <- filenames
          #return(as.list(full.paths))
          return(filenames)
          })

#------------------------------------------------------------------------------------------------------------------------
#' Get the specified variants table
#'
#' @rdname getVariantDataset
#' @aliases getVariantDataset
#'
#' @param obj An object of class TrenaProjectHG38
#' @param datasetName character string, the variant dataset of interest
#'
#' @export

setMethod('getVariantDataset', 'TrenaProjectHG38',

    function(obj, datasetName){
        stopifnot(!obj@variantsDirectory == "/dev/null")
        file.name <- sprintf("%s.RData", file.path(obj@variantsDirectory, datasetName))
        tbl <- NULL
        eval(parse(text=sprintf("tbl <- %s", load(file.name))))
        tbl
        })

#------------------------------------------------------------------------------------------------------------------------
#' Get all the enhancer regions for the gene
#'
#' @rdname getEnhancers
#' @aliases getEnhancers
#'
#' @param obj An object of class TrenaProjectHG38
#' @param targetGene default NA, in which case the current object's targetGene is used.
#'
#' @seealso setTargetGene
#'
#' @export

setMethod('getEnhancers',  'TrenaProjectHG38',

     function(obj, targetGene=NA_character_){
        if(is.na(targetGene))
           targetGene <- getTargetGene(obj)
        stopifnot(!is.null(targetGene))
        tbl.enhancers <- data.frame() # suppress R CMD CHECK NOTE
        full.path <- system.file(package="TrenaProjectHG38", "extdata", "genomeAnnotation", "geneHancer.v4.7.allGenes.RData")
        stopifnot(file.exists(full.path))
        load(full.path)
        #geneSymbol <- NULL
        tbl.out <- subset(tbl.enhancers, toupper(geneSymbol) == toupper(targetGene))
        printf("found %d enhancers for %s", nrow(tbl.out), targetGene)
        tbl.out
        })

#------------------------------------------------------------------------------------------------------------------------
#' Get all the dnase hypersensitivity regions in the expansive region covered by the enhancer
#'
#' @rdname getEncodeDHS
#' @aliases getEncodeDHS
#'
#' @param obj An object of class TrenaProjectHG38
#'
#' @export

setMethod('getEncodeDHS',   'TrenaProjectHG38',

    function(obj){
       hdf <- HumanDHSFilter("hg38", "wgEncodeRegDnaseClustered", pwmMatchPercentageThreshold=0,
                             geneInfoDatabase.uri="bogus", regions=data.frame(), pfms=list())
       tbl.enhancers <- getEnhancers(obj)
       #browser()
       chrom <- tbl.enhancers$chrom[1]
       loc.min <- min(tbl.enhancers$start)
       loc.max <- max(tbl.enhancers$end)
       tbl.dhs <- getRegulatoryRegions(hdf, "wgEncodeRegDnaseClustered", chrom, loc.min, loc.max)
           # todo: close the MySQL connection used above.
       return(tbl.dhs)
       })

#------------------------------------------------------------------------------------------------------------------------
#' Get all the ReMap 2018 ChIP-seq binding sites in the extended (i.e., enhancers-extended) region of the gene
#'
#' @rdname getChipSeq
#' @aliases getChipSeq
#'
#' @param obj An object of class TrenaProjectHG38
#' @param chrom string
#' @param start numeric
#' @param end numeric
#' @param tfs one of more tfs - limit the hits to those for this transcription factor/s
#'
#' @export

setMethod('getChipSeq',  'TrenaProjectHG38',

    function(obj, chrom, start, end, tfs=NA){
       db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host="khaleesi")
       query <- sprintf("select * from chipseq where chrom='%s' and start >= %d and endpos <= %d", chrom, start, end)
       tbl.chipSeq <- dbGetQuery(db, query)
       tf <- NULL  # quiet R CMD CHECK NOTE
       if(!obj@quiet) message(sprintf("tfs before filtering: %d", length (tbl.chipSeq$tf)))
       if(!(all(is.na(tfs))))
         tbl.chipSeq <- subset(tbl.chipSeq, tf %in% tfs)
       if(!obj@quiet){
          message(sprintf("incoming tf filter count: %d", length(tfs)))
          message(sprintf("tfs after filtering: %d", length(unique((tbl.chipSeq$tf)))))
          }
       dbDisconnect(db)
       return(tbl.chipSeq)
       })

#------------------------------------------------------------------------------------------------------------------------
#' Get the covariates table in which each sample, for which expression data is available, is described
#'
#' @rdname getCovariatesTable
#' @aliases getCovariatesTable
#'
#' @param obj An object of class TrenaProjectHG38
#'
#' @return A data.frame
#' @export

setMethod('getCovariatesTable', 'TrenaProjectHG38',

      function(obj){
         tbl <- data.frame()
         if(!is.na(obj@covariatesFile))
            eval(parse(text=sprintf("tbl <- %s", load(obj@covariatesFile))))
         tbl
         })

#------------------------------------------------------------------------------------------------------------------------
#' Get the chromosomal region surrounding the current targetGene, with a flanking percentage added up and downstream
#'
#' @rdname getGeneRegion
#' @aliases getGeneRegion
#'
#' @param obj An object of class TrenaProjectHG38
#' @param flankingPercent a numeric percentage of the gene's total span
#'
#' @return a chrom.loc (chrom:start-end) string
#' @export

setMethod('getGeneRegion',  'TrenaProjectHG38',
          function(obj, flankingPercent=0){
             tbl.transcripts <- getTranscriptsTable(obj)[1,]  # currently always nrow of 1
             chrom <- tbl.transcripts$chr
             start <- tbl.transcripts$start
             end   <- tbl.transcripts$end
             span <- 1 + end - start
             flank <- round(span * (flankingPercent/100))
             chromLocString <- sprintf("%s:%d-%d", chrom, start - flank, end + flank)
             list(chrom=chrom, start=start, end=end, chromLocString=chromLocString)
             })

#------------------------------------------------------------------------------------------------------------------------
#' Get the chromosomal region enclosing the enhancers of the current targetGene, with a flanking percentage added
#'
#' @rdname getGeneEnhancersRegion
#' @aliases getGeneEnhancersRegion
#'
#' @param obj An object of class TrenaProjectHG38
#' @param flankingPercent a numeric percentage of the gene's total span
#'
#' @return a chrom.loc (chrom:start-end) string
#'
#' @export

setMethod('getGeneEnhancersRegion',  'TrenaProjectHG38',
          function(obj, flankingPercent=0){
             tbl.enhancers <- getEnhancers(obj)
             if(nrow(tbl.enhancers) == 0){  # fake it
                message(sprintf("no enhancers for this TrenaProjectHG38"))
                return(getGeneRegion(obj, flankingPercent=100))
                }

             chrom <- tbl.enhancers$chrom[1]
             start <- min(tbl.enhancers$start)
             end <- max(tbl.enhancers$end)
             span <- 1 + end - start
             flank <- round(span * (flankingPercent/100))
             chromLocString <- sprintf("%s:%d-%d", chrom, start - flank, end + flank)
             list(chrom=chrom, start=start, end=end, chromLocString=chromLocString)
             })

#------------------------------------------------------------------------------------------------------------------------
#' return the data.frame with gene ids, chromosome, tss of the primary transcript, and strand for all genes in the project
#'
#' @rdname getGeneInfoTable
#' @aliases getGeneInfoTable
#'
#' @param obj An object of class TrenaProjectHG38
#'
#' @return a data.frame with gene names as rownames
#'
#' @export

setMethod('getGeneInfoTable',  'TrenaProjectHG38',

   function(obj){
      return(obj@geneInfoTable)
      })

#------------------------------------------------------------------------------------------------------------------------
#' Do we have expression data for the suggested gene? genomic and epigenetic information?
#'
#' @rdname recognizedGene
#' @aliases recognizedGene
#'
#' @param obj An object of class TrenaProjectHG38
#' @param geneName A character string in the same protocol as the project's expression matrices
#'
#' @return a chrom.loc (chrom:start-end) string
#'
#' @export

setMethod('recognizedGene',  'TrenaProjectHG38',

   function(obj, geneName){
      geneName.regex <- sprintf("^%s$", geneName)
      index <- grep(toupper(geneName.regex), toupper(obj@geneInfoTable$geneSymbol))
      return(length(index) > 0)
      #tbl.geneInfo <-getGeneInfoTable(obj)
      #return(geneName %in% tbl.geneInfo$geneSymbol)
      })

#------------------------------------------------------------------------------------------------------------------------

