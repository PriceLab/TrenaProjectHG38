#' @import trena
#' @import TrenaProject
#' @importMethodsFrom TrenaProject getEncodeDHS
#' @importMethodsFrom TrenaProject getChipSeq
#' @importMethodsFrom TrenaProject getGeneRegion
#' @importMethodsFrom TrenaProject getGeneEnhancersRegion
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
#' @import methods

.TrenaProjectHG38 <- setClass("TrenaProjectHG38", contains="TrenaProject")

#------------------------------------------------------------------------------------------------------------------------
#' Define an object of class TrenaProjectHG38
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
#' @param packageDataDirectory A string pointing to the parent of a more-or-less standard set of data subdirectories
#' @param quiet A logical indicating whether or not the Trena object should print output
#'
#' @return An object of the TrenaProjectHG38 class
#'
#' @export
#'
#'
TrenaProjectHG38 <- function(projectName,
                             supportedGenes,
                             geneInfoTable.path,
                             footprintDatabaseHost,
                             footprintDatabaseNames,
                             footprintDatabasePort=5432,
                             packageDataDirectory,
                             quiet)
{

      # gene-specific information, freshly assigned with every call to setTargetGene

   geneInfoTable.path <- system.file(package="TrenaProjectHG38", "extdata", "geneInfoTable_hg38.RData")
   stopifnot(file.exists(geneInfoTable.path))

   .TrenaProjectHG38(TrenaProject(projectName,
                                  supportedGenes=character(0),
                                  genomeName="hg38",
                                  geneInfoTable.path=geneInfoTable.path,
                                  footprintDatabaseHost=footprintDatabaseHost,
                                  footprintDatabasePort=footprintDatabasePort,
                                  footprintDatabaseNames=footprintDatabaseNames,
                                  packageDataDirectory=packageDataDirectory,
                                  quiet=quiet
                                  ))


} # ctor
#------------------------------------------------------------------------------------------------------------------------
## #' Get all the enhancer regions for the gene
## #'
## #' @rdname getEnhancers
## #' @aliases getEnhancers
## #'
## #' @param obj An object of class TrenaProjectHG38
## #' @param targetGene default NA, in which case the current object's targetGene is used.
## #'
## #' @seealso setTargetGene
## #'
## #' @export
##
## setMethod('getEnhancers',  'TrenaProjectHG38',
##
##      function(obj, targetGene=NA_character_){
##         if(is.na(targetGene))
##            targetGene <- getTargetGene(obj)
##         if(is.null(targetGene)) return(data.frame())
##         tbl.enhancers <- data.frame() # suppress R CMD CHECK NOTE
##         full.path <- system.file(package="TrenaProjectHG38", "extdata", "genomeAnnotation", "geneHancer.v4.7.allGenes.RData")
##         stopifnot(file.exists(full.path))
##         load(full.path)
##         subset(tbl.enhancers, toupper(geneSymbol) == toupper(targetGene))
##         })
##
#------------------------------------------------------------------------------------------------------------------------
#' Get all the dnase hypersensitivity regions in the expansive region covered by the enhancer
#'
#' @rdname getEncodeDHS
#' @aliases getEncodeDHS
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getEncodeDHS',   'TrenaProject',

    function(obj){

       hdf <- HumanDHSFilter("hg38", "wgEncodeRegDnaseClustered", pwmMatchPercentageThreshold=0,
                             geneInfoDatabase.uri="bogus", regions=data.frame(), pfms=list())
       tbl.enhancers <- getEnhancers(obj)
       if(nrow(tbl.enhancers) == 0) return(data.frame())
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
#' @param obj An object of class TrenaProject
#' @param chrom string
#' @param start numeric
#' @param end numeric
#' @param tfs one of more tfs - limit the hits to those for this transcription factor/s
#'
#' @export

setMethod('getChipSeq',  'TrenaProject',

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
#' Get the chromosomal region enclosing the enhancers of the current targetGene, with a flanking percentage added
#'
#' @rdname getGeneEnhancersRegion
#' @aliases getGeneEnhancersRegion
#'
#' @param obj An object of class TrenaProject
#' @param flankingPercent a numeric percentage of the gene's total span
#'
#' @return a chrom.loc (chrom:start-end) string
#'
#' @export

setMethod('getGeneEnhancersRegion',  'TrenaProject',
          function(obj, flankingPercent=0){
             tbl.enhancers <- getEnhancers(obj)
             if(nrow(tbl.enhancers) == 0){  # fake it
                message(sprintf("no enhancers for this TrenaProject"))
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
