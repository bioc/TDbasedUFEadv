require(RTCGA.rnaseq)
Cancer_cell_lines <- list(ACC.rnaseq,BLCA.rnaseq,BRCA.rnaseq)
Drug_and_Disease <- prepareexpDrugandDisease(Cancer_cell_lines)
expDrug <- Drug_and_Disease$expDrug
expDisease <- Drug_and_Disease$expDisease
require(Biobase)
Z <- prepareTensorfromMatrix(exprs(expDrug[seq_len(20),seq_len(10)]),
                             exprs(expDisease[seq_len(20),seq_len(10)]))
sample<- outer(colnames(expDrug)[seq_len(10)],
               colnames(expDisease)[seq_len(10)],function(x,y){paste(x,y)})
require(TDbasedUFE)
Z <- PrepareSummarizedExperimentTensor(
    sample=sample,feature=rownames(expDrug)[seq_len(20)],value=Z)
HOSVD <- computeHosvd(Z)
Cond <- prepareCondDrugandDisease(expDrug)

require(MOFAdata)
data("CLL_data")
data("CLL_covariates")
Z <- prepareTensorfromMatrix(t(CLL_data$Drugs[seq_len(200),seq_len(50)]),
                             t(CLL_data$Methylation[seq_len(200),seq_len(50)]))
Z <- prepareSummarizedExperimentTensorRect(
    sample=colnames(CLL_data$Drugs)[seq_len(50)],
    feature=list(Drugs=rownames(CLL_data$Drugs)[seq_len(200)],
                 Methylatiion=rownames(CLL_data$Methylation)[seq_len(200)]),
    sampleData=list(CLL_covariates[,1][seq_len(50)]),
    value=Z)
HOSVD <- computeHosvd(Z)
cond <- list(attr(Z,"sampleData")[[1]],NULL,NULL)
index_all <- selectFeatureTransRect(HOSVD,cond,de=c(0.01,0.01),
                                    input_all=8) #batch mode


SVD <- computeSVD(t(CLL_data$Drugs),t(CLL_data$Methylation))
Z <-CLL_data$Drugs %*%t(CLL_data$Methylation)
sample<- colnames(CLL_data$Methylation)
Z <- prepareSummarizedExperimentTensorRect(sample=sample,
    feature=list(rownames(CLL_data$Drugs),rownames(CLL_data$Methylation)),
    value=array(NA,dim(Z)),sampleData=list(CLL_covariates[,1]))
cond <- list(NULL,attr(Z,"sampleData")[[1]],attr(Z,"sampleData")[[1]])
SVD <- transSVD(SVD)
index_all <- selectFeatureRect(SVD,cond,de=c(0.5,0.5),input_all=6) #batch mode

Z <- prepareTensorfromList(CLL_data,as.integer(10))

Multi <- list(BLCA.rnaseq[seq_len(100),1+seq_len(1000)],
              BRCA.rnaseq[seq_len(100),1+seq_len(1000)],
              CESC.rnaseq[seq_len(100),1+seq_len(1000)],
              COAD.rnaseq[seq_len(100),1+seq_len(1000)])

Z <- prepareTensorfromList(Multi,as.integer(10))
Z <- aperm(Z,c(2,1,3))

require(RTCGA.clinical)
Clinical <- list(BLCA.clinical,BRCA.clinical,CESC.clinical,COAD.clinical)
Multi_sample <- list(BLCA.rnaseq[seq_len(100),1,drop=FALSE],
                     BRCA.rnaseq[seq_len(100),1,drop=FALSE],
                     CESC.rnaseq[seq_len(100),1,drop=FALSE],
                     COAD.rnaseq[seq_len(100),1,drop=FALSE])
#patient.stage_event.tnm_categories.pathologic_categories.pathologic_m
ID_column_of_Multi_sample  <- c(770,1482,773,791)
#patient.bcr_patient_barcode
ID_column_of_Clinical  <- c(20,20,12,14)
Z <- PrepareSummarizedExperimentTensor(
    feature =colnames(ACC.rnaseq)[1+seq_len(1000)],
    sample=array("",1),value=Z,
    sampleData=prepareCondTCGA(Multi_sample,
                               Clinical,ID_column_of_Multi_sample,
                               ID_column_of_Clinical))
HOSVD <- computeHosvd(Z)
cond<- attr(Z,"sampleData")
index <- selectFeatureProj(HOSVD,Multi,cond,de=1e-3,input_all=3) #Batch mode
