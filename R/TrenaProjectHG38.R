#' @import trena
#' @import TrenaProject
#' @importMethodsFrom TrenaProject getEnhancers
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
TrenaProjectHG38 <- function(projectName,
                             supportedGenes,
                             geneInfoTable.path,
                             footprintDatabaseHost,
                             footprintDatabaseNames,
                             footprintDatabasePort=5432,
                             expressionDirectory,
                             variantsDirectory,
                             covariatesFile,
                             quiet)
{

      # gene-specific information, freshly assigned with every call to setTargetGene

   geneInfoTable.path <- system.file(package="TrenaProjectHG38", "extdata", "geneInfoTable_hg38.RData")
   stopifnot(file.exists(geneInfoTable.path))

   .TrenaProjectHG38(TrenaProject(projectName,
                                  supportedGenes=supportedGenes,
                                  genomeName="hg38",
                                  geneInfoTable.path=geneInfoTable.path,
                                  footprintDatabaseHost=footprintDatabaseHost,
                                  footprintDatabaseNames=footprintDatabaseNames,
                                  footprintDatabasePort=footprintDatabasePort,
                                  expressionDirectory=expressionDirectory,
                                  variantsDirectory=variantsDirectory,
                                  covariatesFile=covariatesFile,
                                  quiet=quiet))


} # ctor
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
        printf("---- TrenaProjectHG38.R::getEnhancers")

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

