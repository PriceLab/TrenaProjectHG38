#' @importFrom RPostgreSQL dbConnect dbListTables dbGetQuery dbListConnections dbDisconnect
#' @import RPostgreSQL
#' @import org.Hs.eg.db
#' @importFrom methods new
#'
#' @title GeneHancerDB
#------------------------------------------------------------------------------------------------------------------------
#' @name GeneHancerDB-class
#' @rdname GeneHancerDB-class
#' @aliases GeneHancerDB
#'
#' @import methods

.GeneHancerDB <- setClass("GeneHancerDB",
                          representation = representation(
                             db="DBIConnection",
                             state="environment"
                             )
                          )
#------------------------------------------------------------------------------------------------------------------------
setGeneric('retrieveEnhancersFromDatabase', signature='obj', function(obj, targetGene, tissues)
              standardGeneric('retrieveEnhancersFromDatabase'))
setGeneric('listTissues', signature='obj', function(obj) standardGeneric('listTissues'))
#------------------------------------------------------------------------------------------------------------------------
#' Create a GeneHancerDB connection
#'
#' @rdname GeneHancerDB-class
#'
#' @return An object of the GeneHancerDB class
#'
#' @export
#'
#'
GeneHancerDB <- function()
{

   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="gh411", host="khaleesi")
   state <- new.env(parent=emptyenv())

   .GeneHancerDB(db=db, state=state)

} # ctor
#------------------------------------------------------------------------------------------------------------------------
#' Get all the enhancer regions for the gene
#'
#' @rdname retrieveEnhancersFromDatabase
#' @aliases retrieveEnhancersFromDatabase
#'
#' @param obj An object of class TrenaProjectHG38
#' @param targetGene a HUGO gene symbol
#' @param tissues "all" or a vector of case-agnostic tissue names
#'
#' @seealso listTissues
#'
#' @export

setMethod('retrieveEnhancersFromDatabase',  'GeneHancerDB',

     function(obj, targetGene, tissues){

        if(length(tissues) == 1 & tissues[1] == "all")
           tissueClause <- ""
        else {
           tissueSet <- paste(tissues, collapse="','")
           tissueClause <- sprintf("AND t.tissue in ('%s') ", tissueSet)
           }

        query <- paste0("select e.chr as chrom, ",
                        "e.element_start as start, ",
                        "e.element_end as end, ",
                        "a.symbol as gene, ",
                        "a.eqtl_score as eqtl, ",
                        "a.chic_score as HiC, ",
                        "a.erna_score as erna, ",
                        "a.expression_score as coexpression, ",
                        "a.distance_score as distanceScore, ",
                        "a.tss_proximity as tssProximity, ",
                        "a.combined_score as combinedScore, ",
                        "a.is_elite as elite, ",
                        "t.source as source, ",
                        "t.tissue as tissue, ",
                        "e.type as type, ",
                        "a.ghid as ghid ",
                        "from associations AS a, ",
                        "tissues AS t, elements as e ",
                        "where a.symbol='%s' ",
                        "%s",
                        "AND a.ghid=t.ghid ",
                        "AND e.ghid=a.ghid")
        query <- sprintf(query, targetGene, tissueClause)

        tbl <- dbGetQuery(obj@db, query)

           #------------------------------------------------------------
           # todo: when eliminating duplicates, collecte the sometimes
           #       multiple tissues in which the same enhancer is found
           #       pshannon (14 aug 2019)
           #------------------------------------------------------------

        tbl$sig <- with(tbl, sprintf("%s:%d-%d", chrom, start, end))
        dups <- which(duplicated(tbl$sig))
        tbl.1 <- tbl
        if(length(dups) > 0)
           tbl.1 <- tbl[-dups, ]

          # our current best guess is that eQTL, Hi-C, and enhancer RNA are credible indicators
          # of enhancer/gene association.  so keep only the rows with a value in one or more
          # of these columns, or with a combined score > 5.
          # combinedscore is some unstated function of all the scores.  we include as a fallback
          # an alternative threshold, just in case.

        tbl.2 <- subset(tbl.1, !(is.nan(eqtl) & is.nan(hic) & is.nan(erna)) | combinedscore >= 5)

        return(tbl.2)
        })

#------------------------------------------------------------------------------------------------------------------------
#' Return a character vector containing all of the tissues known to GeneHancer
#'
#' @rdname listTissues
#' @aliases listTissues
#'
#' @param obj An object of class GeneHancerDB
#'
#' @export

setMethod('listTissues', 'GeneHancerDB',

    function(obj){
      dbGetQuery(obj@db, "select distinct tissue from tissues")$tissue
      })

#------------------------------------------------------------------------------------------------------------------------
