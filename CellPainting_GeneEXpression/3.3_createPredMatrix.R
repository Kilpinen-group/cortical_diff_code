setwd(".")


library(oob)

inputDir <- "zenodoArchive/"
dir.create("out3.3", showWarnings = FALSE)

gex2cpModel<-fastRead(paste0(inputDir,"/CPGEX_GeXmodellLasso.tsv"))

computeTheoMat <- function(mat, model) {
    nfeatures <- nrow(model) - 1
    res <- mat %*% model[2:(nfeatures + 1), ]
    sweep(res, 1, model[1, ], "+")
}
#create cp theo mat from gex model
scData<-readRDS("countTableCorrected.rds") |> t() #available upon request
cpTheoMat<-computeTheoMat(scData, gex2cpModel)

fastWrite(cpTheoMat, "out3.3/CPGEX_predictedCPfeatureFromGeX.tsv")


