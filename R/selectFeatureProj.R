#' Title Select feature when projection strategy is employed for the
#'  case where features are shared with multiple omics profiles
#'
#' @param HOSVD HOSVD 
#' @param Multi list of omics profiles, row: sample, column: feature
#' @param cond list of conditions for individual omics profiles
#' @param de initial value for optimization of statdard deviation
#' @param p0 Threshold P-value
#' @param breaks The number of bins of histogram of P-values
#' @param input_all The number of selected feature. if null, intearactive mode
#' is activated 
#'
#' @return list composed of logical vector that represent which features are selected and p-values
#' @export
#'
#' @examples
#' require(TDbasedUFE)
#' Multi <- list(matrix(runif(1000),10),matrix(runif(1000),10),
#' matrix(runif(1000),10),matrix(runif(1000),10))
#' Z <- prepareTensorfromList(Multi,10)
#' Z <- aperm(Z,c(2,1,3))
#' Z <- PrepareSummarizedExperimentTensor(feature =as.character(1:10),
#'                                       sample=array("",1),value=Z)
#' HOSVD <- computeHosvd(Z)
#' cond <- rep(list(rep(1:2,each=5)),4)
#' index <- selectFeatureProj(HOSVD,Multi,cond,de=0.1,input_all=2)
selectFeatureProj <-
    function(HOSVD,Multi,cond,de=1e-4,p0=0.01,breaks=100,input_all=NULL){
    if (is.null(input_all))
    {
    LIST1 <- lapply(Multi,function(x){data.matrix(x)%*%data.matrix(HOSVD$U[[1]])})
    par(mfrow=c(length(cond),1))
    j<-1
    while(j %in% seq_len(dim(HOSVD$U[[1]])[2]))
    {
    for (i in seq_len(length(cond)))
    {
        boxplot(LIST1[[i]][,j]~cond[[i]],main=paste(j,i,sep="-"))
        abline(0,0,col=2,lty=2)
    }
        input <- menu(c("NEXT","PREV","SELCT"))
        if (input==2){if (j!=1){j<-j-1}} 
        else if (input==3){break} 
        else {if (j<dim(HOSVD$U[[1]])[2])j<-j+1}
    }
    input_all <- j
    }
    th <- function(sd,breaks,p0){
        P2 <- pchisq((u/sd)^2,1,lower.tail=FALSE)
        hc<- hist(1-P2,breaks=breaks,plot=FALSE)
        return(sd(hc$count[seq_len(sum(hc$breaks
                                       <1-min(P2[p.adjust(P2,"BH")>p0])))]))
    }
        u<- HOSVD$U[[1]][,input_all]
        sd <- optim(de,function(x){th(x,breaks,p0)},
                    control=list(warn.1d.NelderMead=FALSE))$par
        sd1 <- seq(0.1*sd,2*sd,by=0.1*sd)
        th0 <- apply(matrix(sd1,ncol=1),1,function(x){th(x,breaks,p0)})
        P2 <- pchisq((u/sd)^2,1,lower.tail=FALSE)
        plot.new()
        par(mfrow=c(1,2))
        plot(sd1,th0,type="o")
        arrows(sd,max(th0),sd,min(th0),col=2)
        hist(1-P2,breaks=breaks)
        par(mfrow=c(1,1))
        index <- p.adjust(P2,"BH")<p0
        index_all <- list(index=index,p.value=P2)
        return(index_all)
}


