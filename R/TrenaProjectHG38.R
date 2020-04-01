#' @import trena
#' @import TrenaProject
#' @importMethodsFrom TrenaProject getEncodeDHS
#' @importMethodsFrom TrenaProject getChipSeq
#' @importMethodsFrom TrenaProject getGeneRegion
#' @importMethodsFrom TrenaProject getGeneEnhancersRegion
## ' @importMethodsFrom TrenaProject getGeneRegulatoryTissues
#' @importMethodsFrom TrenaProject getGeneRegulatoryRegions
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

.TrenaProjectHG38 <- setClass("TrenaProjectHG38",
                              contains="TrenaProject",
                              representation=representation(
                                 genehancer="GeneHancerDB"
                                 ))

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
                             footprintDatabaseHost,
                             footprintDatabaseNames,
                             footprintDatabasePort=5432,
                             packageDataDirectory,
                             quiet)
{

      # gene-specific information, freshly assigned with every call to setTargetGene

   geneInfoTable.path <- system.file(package="TrenaProjectHG38", "extdata", "geneInfoTable_hg38.RData")
   stopifnot(file.exists(geneInfoTable.path))
   ghdb <- GeneHancerDB()

   .TrenaProjectHG38(TrenaProject(projectName,
                                  supportedGenes=supportedGenes,
                                  genomeName="hg38",
                                  geneInfoTable.path=geneInfoTable.path,
                                  footprintDatabaseHost=footprintDatabaseHost,
                                  footprintDatabasePort=footprintDatabasePort,
                                  footprintDatabaseNames=footprintDatabaseNames,
                                  packageDataDirectory=packageDataDirectory,
                                  quiet=quiet
                                  ),
                     genehancer=ghdb)


} # ctor
#------------------------------------------------------------------------------------------------------------------------
#' Get gene regulatory regions for the current target gene, or one explicitly named, using the specified strategy
#'
#' @rdname getGeneRegulatoryRegions
#' @aliases getGeneRegulatoryRegions
#'
#' @param obj An object of class TrenaProjectHG38
#' @param targetGene default NA, in which case the current object's targetGene is used.
#'
#' @seealso setTargetGene
#'
#' @export

setMethod('getGeneRegulatoryRegions',  'TrenaProjectHG38',

    function(obj, targetGene=NA, tissues="all",
             geneHancerMissing.promoter.upstream=0,
             geneHancerMissing.promoter.downstream=0,
             geneHancerSupplemental.promoter.upstream=0,
             geneHancerSupplemental.promoter.downstream=0){

       if(is.na(targetGene))
          targetGene <- getTargetGene(obj)

       tbl <- retrieveEnhancersFromDatabase(obj@genehancer, targetGene, tissues)
       tbl.ghMissing <- data.frame()
       tbl.ghSupplemental <- data.frame()
          # if no regions from GeneHancer, and geneHancerMissing locations are provided:
       gh.missing <-!(geneHancerMissing.promoter.upstream  == 0 && geneHancerMissing.promoter.downstream == 0)
       gh.missing <- gh.missing & nrow(tbl) == 0
          # if > 0 regions obtained from GeneHancer, and supplemental promoter locations provided:
       gh.supplementing <- !(geneHancerSupplemental.promoter.upstream == 0 && geneHancerSupplemental.promoter.downstream == 0)
       gh.supplementing <- gh.supplementing & nrow(tbl) > 0

       #if(all(c(gh.missing, gh.supplementing))){
       #   msg <- sprintf("error! you requested both a 'geneHancerMissing' and a 'geneHancerSupplemental' promoter")
       #   stop(msg)
       #   }

          # get vital genomic info on the targetGene
       tbl.transcripts <- getTranscriptsTable(obj, targetGene)  # just one row, by convention
       chrom <- tbl.transcripts$chrom[1]
       tss <-  tbl.transcripts$tss[1]
       strand <- tbl.transcripts$strand[1]

       if(gh.missing){ # add a generic promoter, typically +/- 5kb
          start.loc <- tss - geneHancerMissing.promoter.upstream
          end.loc   <- tss + geneHancerMissing.promoter.downstream
          if(strand == -1){
             start.loc <- tss - geneHancerMissing.promoter.downstream
             end.loc   <- tss + geneHancerMissing.promoter.upstream
             }
          tbl.ghMissing <- data.frame(chrom=chrom, start=start.loc, end=end.loc, gene=targetGene,
                                      eqtl=NA, hic=NA, erna=NA, coexpression=NA, distancescore=NA, tssproximity=NA,
                                      combinedscore=NA, elite=NA, source="TrenaProjectHG38.gh.missing",
                                      type=NA, ghid=NA, tissue=NA, sig=NA, stringsAsFactors=FALSE)
          tbl <- tbl.ghMissing
          } # if gh.missing

       if(gh.supplementing){ # add a very conservative classical promoter, typically +/- 2kb
          start.loc <- tss - geneHancerSupplemental.promoter.upstream
          end.loc   <- tss + geneHancerSupplemental.promoter.downstream
          if(strand == -1){
             start.loc <- tss - geneHancerSupplemental.promoter.downstream
             end.loc   <- tss + geneHancerSupplemental.promoter.upstream
             }
          tbl.ghSupplemental <- data.frame(chrom=chrom, start=start.loc, end=end.loc, gene=targetGene,
                                      eqtl=NA, hic=NA, erna=NA, coexpression=NA, distancescore=NA, tssproximity=NA,
                                      combinedscore=NA, elite=NA, source="TrenaProjectHG38.gh.supplemental",
                                      type=NA, ghid=NA, tissue=NA,  stringsAsFactors=FALSE)
          tbl <- rbind(tbl, tbl.ghSupplemental)
          } # gh.supplementing

          # make sure the tbl is in sort order, ascending by start location
       if(nrow(tbl) > 0){
          new.order <- order(tbl$start, decreasing=FALSE)
          tbl <- tbl[new.order,]
          }
       tbl
       })

#------------------------------------------------------------------------------------------------------------------------
#' Get all the enhancer tissues included in the current genehancer
#'
#' @rdname getEnhancerTissues
#' @aliases getEnhancerTissues
#'
#' @param obj An object of class TrenaProjectHG38
#'
#' @seealso getEnhancers
#'
#' @export

setMethod('getEnhancerTissues',  'TrenaProjectHG38',

     function(obj){
        listTissues(obj@genehancer)
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

     function(obj, targetGene=NA_character_, tissues="all", maxSize=10000){
        if(is.na(targetGene))
           targetGene <- getTargetGene(obj)
        if(is.null(targetGene)) return(data.frame())
        tbl <- retrieveEnhancersFromDatabase(obj@genehancer, targetGene, tissues)
        if(nrow(tbl) == 0)
           return(data.frame())
        size <- with(tbl, 1 + end - start)
        deleters <- which(size > maxSize)
        if(length(deleters) > 0)
           tbl <- tbl[-deleters,]
        tbl
        })

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
#' @param chrom string vector
#' @param start numeric vector
#' @param end numeric vector
#' @param tfs one of more tfs - limit the hits to those for this transcription factor/s
#'
#' @return data.frame
#'
#' @export

setMethod('getChipSeq',  'TrenaProject',

    function(obj, chrom, start, end, tfs=NA){

       db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host="khaleesi")
       # query <- sprintf("select * from chipseq where chrom='%s' and start >= %d and endpos <= %d", chrom, start, end)
       # tbl.chipSeq <- dbGetQuery(db, query)
       runQuery <- function(i){
          query <- sprintf("select * from chipseq where chrom='%s' and start >= %d and endpos <= %d", chrom[i], start[i], end[i])
          printf(query)
          tbl.hits <- dbGetQuery(db, query)
          printf("ChIP-seq hits: %d", nrow(tbl.hits))
          tbl.hits
          }
       x <- lapply(seq_len(length(chrom)), runQuery)
       dbDisconnect(db)
       tbl.chipSeq <- do.call(rbind, x)
       tf <- NULL  # quiet R CMD CHECK NOTE
       if(!obj@quiet) message(sprintf("tfs before filtering: %d", length (tbl.chipSeq$tf)))
       if(nrow(tbl.chipSeq) == 0){
          message(sprintf("TrenaProjectHG38:getChipSeq, no ChIP for %s in %s:%d-%d\n",
                          paste(tfs, collapse=","), chrom[1], min(start), max(end)))
          return(data.frame())
          }
       if(!(all(is.na(tfs))))
         tbl.chipSeq <- subset(tbl.chipSeq, tf %in% tfs)
       if(!obj@quiet){
          message(sprintf("incoming tf filter count: %d", length(tfs)))
          message(sprintf("tfs after filtering: %d", length(unique((tbl.chipSeq$tf)))))
          }
       colnames(tbl.chipSeq)[grep("endpos", colnames(tbl.chipSeq))] <- "end"
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
