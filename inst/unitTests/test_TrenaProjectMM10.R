library(TrenaProjectMM10)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tProj")) {
   message(sprintf("--- creating instance of TrenaProjectMM10"))
   tProj <- TrenaProjectMM10();
   }
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_getTranscriptTable()
   test_supportedGenes()
   test_variants()
   test_expressionMatrices()
   test_setTargetGene()

   test_getEnhancers()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue(all(c("TrenaProjectMM10", "TrenaProject") %in% is(tProj)))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_getTranscriptTable <- function()
{
   message(sprintf("--- test_getTranscriptable"))

   expected.colnames <- c("ensg", "chrom", "start", "end", "tss", "strand", "geneSymbol", "entrez", "appris",
                          "tsl", "transcript", "type")
   tbl.apoh <- getTranscriptsTable(tProj, "Apoh")
   checkEquals(dim(tbl.apoh), c(1, 12))
   checkEquals(colnames(tbl.apoh), expected.colnames)

   tbl.all <-  getTranscriptsTable(tProj, all=TRUE)
   checkTrue(nrow(tbl.all) > 25000)
   checkEquals(colnames(tbl.all), expected.colnames)

} # test_getTranscriptTable
#------------------------------------------------------------------------------------------------------------------------
test_supportedGenes <- function()
{
   message(sprintf("--- test_supportedGenes"))

   subset.expected <- c("Abca1")
   checkTrue(all(subset.expected %in% getSupportedGenes(tProj)))

} # test_supportedGenes
#------------------------------------------------------------------------------------------------------------------------
test_variants <- function()
{
   message(sprintf("--- test_variants"))

   checkEquals(length(getVariantDatasetNames(tProj)), 0)

} # test_variants
#------------------------------------------------------------------------------------------------------------------------
test_expressionMatrices <- function()
{
   expected <- c("thioglycollate-elicited-peritoneal-macrophages")
   checkTrue(all(expected %in% getExpressionMatrixNames(tProj)))

   # mtx <- getExpressionMatrix(tProj, expected[1])
   # checkEquals(dim(mtx), c(6890, 14))

} # test_expressionMatrices
#------------------------------------------------------------------------------------------------------------------------
# setting the target gene implies a few other assignements, all tested here:
#   geneInfo (temporarily also masquerading at tbl.transcripts
#   geneRegion
#   geneEnhancersRegion (when avaialable, defaults to geneRegion)
#
test_setTargetGene <- function()
{
   message(sprintf("--- test_setTargetGene"))

   setTargetGene(tProj, "Abca1")
   checkEquals(getTargetGene(tProj), "Abca1")


} # test_setTargetGene
#------------------------------------------------------------------------------------------------------------------------
test_getEnhancers <- function()
{
   message(sprintf("--- test_getEnhancers"))

   tbl.gh.trem2 <- getEnhancers(tProj, targetGene="Trem2")
   checkEquals(dim(tbl.gh.trem2), c(5, 14))
   checkTrue(with(tbl.gh.trem2, all(end-start) < 10000))
   checkTrue(all(tbl.gh.trem2$combined_score >= 2.0))

} # test_getEnhancers
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
