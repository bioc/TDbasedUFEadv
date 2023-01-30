test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
require(RTCGA.rnaseq)
LIST <- list(ACC.rnaseq,
             BLCA.rnaseq,
             BRCA.rnaseq,
             CESC.rnaseq,
             COAD.rnaseq,
             ESCA.rnaseq,
             GBM.rnaseq,
             HNSC.rnaseq,
             KICH.rnaseq,
             KIRC.rnaseq,
             KIRP.rnaseq,
             LGG.rnaseq,
             LIHC.rnaseq,
             LUAD.rnaseq,
             LUSC.rnaseq,
             OV.rnaseq,
             PAAD.rnaseq,
             PCPG.rnaseq,
             PRAD.rnaseq,
             READ.rnaseq,
             SARC.rnaseq,
             SKCM.rnaseq,
             STAD.rnaseq,
             TGCT.rnaseq,
             THCA.rnaseq,
             UCEC.rnaseq,
             UCS.rnaseq)
dummy <- prepareexpDrugandDisease(LIST)
expDrug <- dummy[[1]]
expDisease <- dummy[[2]]
require(Biobase)
Z <- prepareTensorfromMatrix(exprs(expDrug[seq_len(1000),seq_len(100)]),exprs(expDisease[seq_len(1000),seq_len(100)]))
sample<- outer(colnames(expDrug)[seq_len(100)],colnames(expDisease)[seq_len(100)],function(x,y){paste(x,y)})
require(TDbasedUFE)
Z <- PrepareSummarizedExperimentTensor(sample=sample,feature=rownames(expDrug)[seq_len(1000)],value=Z)
HOSVD <- computeHosvd(Z)
cond <- list(NULL,rep(1:2,each=50),rep(1:2,each=50))
