#----------------------------------------------------------------------------------------------------
#' @import methods
#' @import TrenaProject
#' @importFrom AnnotationDbi select
#' @import org.Mm.eg.db
#'
#' @title TrenaProjectMM10-class
#'
#' @name TrenaProjectMM10-class
#' @rdname TrenaProjectMM10-class
#' @aliases TrenaProjectMM10
#' @exportClass TrenaProjectMM10
#'

.TrenaProjectMM10 <- setClass("TrenaProjectMM10", contains="TrenaProject")

#----------------------------------------------------------------------------------------------------
#' Define an object of class TrenaProjectMM10
#'
#' @description
#' Expression, variant and covariate data for the genes of interest (perhaps unbounded) for pre-term birth studies
#'
#' @rdname TrenaProjectMM10-class
#'
#' @export
#'
#' @return An object of the TrenaProjectMM10 class
#'

TrenaProjectMM10 <- function(quiet=TRUE)

{
   genomeName <- "mm10"

   directory <- system.file(package="TrenaProjectMM10", "extdata", "geneSets")
   geneSet.files <- list.files(directory)
   geneSets <- list()
   for(file in geneSet.files){
      full.path <- file.path(directory, file)
      genes <- scan(full.path, sep="\t", what=character(0), quiet=TRUE)
      geneSet.name <- sub(".txt", "", file)
      geneSets[[geneSet.name]] <- genes
      }

   footprintDatabaseNames <- NA_character_;
   expressionDirectory <- system.file(package="TrenaProjectMM10", "extdata", "expression")
   variantsDirectory <- system.file(package="TrenaProjectMM10", "extdata", "variants")
   footprintDatabaseHost <- NA_character_;
   geneInfoTable.path <- NA_character_;
   genomicRegionsDirectory <- NA_character_;
   covariatesFile <- NA_character_;

   stopifnot(file.exists(expressionDirectory))

   .TrenaProjectMM10(TrenaProject("TrenaProjectMM10",
                                  supportedGenes=geneSets[[1]],
                                  genomeName=genomeName,
                                  geneInfoTable.path = geneInfoTable.path,
                                  footprintDatabaseHost=footprintDatabaseHost,
                                  footprintDatabaseNames=footprintDatabaseNames,
                                  footprintDatabasePort=5432,
                                  expressionDirectory=expressionDirectory,
                                  genomicRegionsDirectory = genomicRegionsDirectory,
                                  variantsDirectory=variantsDirectory,
                                  covariatesFile=covariatesFile,
                                  quiet=quiet
                                  ))

} # TrenaProjectMM10, the constructor
#----------------------------------------------------------------------------------------------------
