\name{CAT}
\alias{CAT}
\title{
Conditional association test
}
\description{
Conditional association test for microbiome data using a permutation-based approach
}
\usage{CAT(testList,otutable,taxonomy,
metric="Weighted UniFrac",metaData,
outcomeVar,tree=NULL,method="PERMANOVA",
numBS=1000,parallel=TRUE,nCore=2)
}
\arguments{
  \item{testList}{A list of taxa to be tested.}
  \item{otutable}{An OTU table, rows correspond to OTUs, while columns correspond to observations. Row names must match the row names in taxonomy, while the column names must match the row names in metaData.}
  \item{taxonomy}{A taxonomy table, rows correspond to OTUs, while columns correspond to taxonomy assignment at different levels. Row names must match the row names in the OTU table.}
  \item{metric}{Can be one single metric or a vector of metrics(for MiRKAT only). It accepts "Weighted UniFrac", "Unweighted UniFrac" and all metrics that can be calculated from the "vegdist"" function in the R package "vegan". The default is "Weighted UniFrac".}
  \item{metaData}{ The dataset containing the information of the study objects. Row names must match the column names in the OTU table.}
  \item{outcomeVar}{The outcome variable in metaData. If it is a vector of two elements, the function will consider the outcome is time-to-event (for now only compatible with MiRKAT) and use the first element in the list as the time variable and the second element in the list as the event indicator. For now it only accepts right censored data.}
  \item{tree}{An object of the class "tree", which is needed for calculating weighted and unweighted UniFrac distances.}
  \item{method}{Can be "PERMANOVA" or "MiRKAT". The default choice is "PERMANOVA".}
  \item{numBS}{Number of bootstrap samples used in the test. The default is 1000}
  \item{parallel}{Whether or not use parallel computing.}
  \item{nCore}{How many cores to use for parallel computing. This is only relevant when parallel=TRUE.}
}

\value{
A vector of p-values for the taxa list given.
}
\details{
The CAT function performs conditional association test for the list of the taxa provided.
}
\source{
Shi Y, Zhang L, Do KA, Jenq RR, Peterson CB (2022) \emph{CAT: a conditional association test for microbiome data using a permutation-based approach}}

\examples{
# Example 1 with PERMANOVA
otuPath<-system.file("extdata","GopalakrishnanOTUtable.csv", 
package = "CATMicrobiome")
otutable<-read.csv(otuPath,header=TRUE,row.names = 1)
taxonomyPath<-system.file("extdata","GopalakrishnanTaxonomy.csv", 
package = "CATMicrobiome")
taxonomy<-read.csv(taxonomyPath,header=TRUE,row.names = 1)
metaPath<-system.file("extdata","GopalakrishnanMeta.csv", 
package = "CATMicrobiome")
metaData<-read.csv(metaPath,header=TRUE,row.names = 1)
treePath<-system.file("extdata","GopalakrishnanTree.tree", 
package = "CATMicrobiome")
tree<-read.tree(treePath)

testList<-c("p__Firmicutes",
            "c__Clostridia",
            "o__Clostridiales",
            "f__Ruminococcaceae",
            "s__prausnitzii",
            "g__Faecalibacterium",
            "s__bromii",
            "g__Ruminococcus")

testResult<-CAT(testList,otutable,taxonomy,
metric="Weighted UniFrac",metaData,outcomeVar="ContOutcomes",
tree,numBS=1000,parallel=TRUE,nCore=2)
testResult
testResult<-CAT(testList,otutable,taxonomy,
metric="Weighted UniFrac",metaData,outcomeVar="ContOutcomes",
tree,numBS=1000,parallel=FALSE)
testResult

testResult<-CAT(testList,otutable,taxonomy,
metric="Weighted UniFrac",metaData,outcomeVar="BinOutcomes",
tree,numBS=1000,parallel=TRUE,nCore=2)
testResult
testResult<-CAT(testList,otutable,taxonomy,
metric="Weighted UniFrac",metaData,outcomeVar="BinOutcomes",
tree,numBS=1000,parallel=FALSE)
testResult

#Example 2 with survival outcomes
otuPath<-system.file("extdata","RiquelmeOTUtable.csv", 
package = "CATMicrobiome")
otutable<-read.csv(otuPath,header=TRUE,row.names = 1)
taxonomyPath<-system.file("extdata","RiquelmeTaxonomy.csv", 
package = "CATMicrobiome")
taxonomy<-read.csv(taxonomyPath,header=TRUE,row.names = 1)
metaPath<-system.file("extdata","RiquelmeMeta.csv", 
package = "CATMicrobiome")
metaData<-read.csv(metaPath,header=TRUE,row.names = 1)
treePath<-system.file("extdata","RiquelmeTree.tree", 
package = "CATMicrobiome")
tree<-read.tree(treePath)

testList<-c("Clostridia",
            "Clostridiales",
            "Lachnospiraceae",
            "Corynebacteriales",
            "Bacteroidia",
            "Bacteroidales",
            "Corynebacteriaceae",
            "Corynebacterium",
            "Pseudomonadales",
            "Pseudomonas")


testResult<-CAT(testList,otutable,taxonomy,
metric=c("bray","Unweighted UniFrac"),metaData,
outcomeVar=c("OS.Months","event"),tree=tree,
method="MiRKAT",parallel=TRUE,nCore=2)
testResult

testResult<-CAT(testList,otutable,taxonomy,
metric=c("bray","Weighted UniFrac","Unweighted UniFrac"),
metaData,outcomeVar=c("OS.Months","event"),
tree=tree,method="MiRKAT",parallel=FALSE)
testResult
}

