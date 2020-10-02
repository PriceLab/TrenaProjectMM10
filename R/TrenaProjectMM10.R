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
   dataDirectory <- system.file(package="TrenaProjectMM10", "extdata")
   footprintDatabaseHost <- NA_character_;
   geneInfoTable.path <- system.file(package="TrenaProjectMM10", "extdata", "geneInfoTable.RData")

   .TrenaProjectMM10(TrenaProject("TrenaProjectMM10",
                                  supportedGenes=geneSets[[1]],
                                  genomeName=genomeName,
                                  geneInfoTable.path = geneInfoTable.path,
                                  footprintDatabaseHost=footprintDatabaseHost,
                                  footprintDatabaseNames=footprintDatabaseNames,
                                  footprintDatabasePort=5432,
                                  packageDataDirectory=dataDirectory,
                                  quiet=quiet
                                  ))

} # TrenaProjectMM10, the constructor
#----------------------------------------------------------------------------------------------------
#' Get all the enhancer regions for the gene
#'
#' @rdname getEnhancers
#' @aliases getEnhancers
#'
#' @param obj An object of class TrenaProjectMM10
#' @param targetGene default NA, in which case the current object's targetGene is used.
#'
#' @seealso setTargetGene
#'
#' @export

setMethod('getEnhancers',  'TrenaProjectMM10',

     function(obj, targetGene=NA_character_, tissues="all", maxSize=10000){

        if(is.na(targetGene))
           targetGene <- getTargetGene(obj)

        if(!exists("tp.hg38")){
           require("TrenaProjectHG38.generic")
           tp.hg38 <- TrenaProjectHG38.generic()
           }
        if(!exists("tbl.gha")){
           f <- "~/github/TrenaProjectMM10/inst/extdata/hg38.genehancer.411.geneAssociations.RData"
           #f <- system.file(package="TrenaProjectMM10", "extdata", "hg38.genehancer.411.geneAssociations.RData")
           message(sprintf("loading hg38 gh gene associations: %s", f))
           tbl.gha <- get(load(f))
           }
        if(!exists("tbl.psl")){
           f <- "~/github/TrenaProjectMM10/inst/extdata/hg38.mm10.genehancer.411.map.RData"
           #f <- system.file(package="TrenaProjectMM10", "extdata", "hg38.mm10.genehancer.411.map.RData")
           message(sprintf("loading human-to-mouse mapping file: %s", f))
           tbl.psl <- get(load(f))
           }

        targetGene.hg38 <- toupper(targetGene)
        targetGene.mm10 <- paste(substring(targetGene.hg38, 1, 1),
                                 tolower(substring(targetGene.hg38, 2, nchar(targetGene.hg38))), sep="")
        tbl.hg38.gene <- getTranscriptsTable(tp.hg38, targetGene.hg38)
        chrom.hg38  <- tbl.hg38.gene$chrom
        tss.hg38    <- tbl.hg38.gene$tss
        strand.hg38 <- tbl.hg38.gene$strand
        printf("hg38 %s %s %d %s", targetGene.hg38, chrom.hg38, tss.hg38, strand.hg38)

        tbl.mm10.gene <- getTranscriptsTable(obj, targetGene.mm10)
        tbl.mm10.gene <- subset(tbl.mm10.gene, nchar(chrom) < 8)
        chrom.mm10  <- tbl.mm10.gene$chrom
        tss.mm10    <- tbl.mm10.gene$tss
        strand.mm10 <- tbl.mm10.gene$strand
        printf("mm10 %s %s %d %s", targetGene.mm10, chrom.mm10, tss.mm10, strand.mm10)

             # find all enhancers.  calling client will do any thresholding needed
        ghids <- subset(tbl.gha, symbol==targetGene.hg38)$GHid
        length(ghids)

        tbl.mouse <- subset(tbl.psl, gh %in% ghids)
        dim(tbl.mouse)
        tbl.mouse <- tbl.mouse[, c("tName", "tStart", "tEnd", "gh")]
        colnames(tbl.mouse) <- c("chrom", "start", "end", "GHid")
        tbl.mouse.merged <- merge(tbl.mouse, tbl.gha, by="GHid")
        dim(tbl.mouse.merged)
        tbl.mouse <- subset(merge(tbl.mouse, tbl.gha, by="GHid"), toupper(symbol)==toupper(targetGene))
        tbl.mouse$targetGene <- targetGene

        minScore <- 2.0
        maxWidth <- 10000

        tbl.mouse <- subset(tbl.mouse, combined_score >= minScore)
        widths <- with(tbl.mouse, 1 + end - start)
        keepers <- which(widths <= 10000)
        length(widths); length(keepers)
        tbl.mouse <- tbl.mouse[keepers,]

        dim(tbl.mouse)
        chrom.matches <- length(grep(chrom.mm10, tbl.mouse$chrom))
        printf("credible enhancers: %d/%d %5.2f%%", chrom.matches, nrow(tbl.mouse), 100 * chrom.matches/nrow(tbl.mouse))
        tbl.mouse
        })

#------------------------------------------------------------------------------------------------------------------------
