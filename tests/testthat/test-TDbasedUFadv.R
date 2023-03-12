test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
require(RTCGA.rnaseq)
LIST <- list(ACC.rnaseq,
             BLCA.rnaseq,
             BRCA.rnaseq)
dummy <- prepareexpDrugandDisease(LIST)
expDrug <- dummy[[1]]
expDisease <- dummy[[2]]
require(Biobase)
Z <- prepareTensorfromMatrix(exprs(expDrug[seq_len(100),seq_len(10)]),exprs(expDisease[seq_len(100),seq_len(10)]))
sample<- outer(colnames(expDrug)[seq_len(10)],colnames(expDisease)[seq_len(10)],function(x,y){paste(x,y)})
require(TDbasedUFE)
Z <- PrepareSummarizedExperimentTensor(sample=sample,feature=rownames(expDrug)[seq_len(100)],value=Z)
HOSVD <- computeHosvd(Z)
cond <- list(NULL,rep(1:2,each=50),rep(1:2,each=50))
