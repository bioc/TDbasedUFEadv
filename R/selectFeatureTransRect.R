#' Title Select features for a tensor generated from two matricies that
#'  share samples.
#'
#' @param HOSVD  HOSVD 
#' @param cond list of conditions
#' @param de initial values for optimization of standard deviation
#' @param p0 threshold value for the significance
#' @param breaks number of bins of the histgram of P-values
#' @param input_all The selected singular value vectors attributed to samples.
#' if NULL, interactive mode
#'
#' @return list of logical vector that represent if the individual features
#' are selected and P-values. 
#' @export 
#'
#' @examples
#' require(TDbasedUFE)
#' matrix1 <- matrix(runif(10000),200) #row features, column samples
#' matrix2 <- matrix(runif(20000),400) #row features, column samples
#' Z <- prepareTensorfromMatrix(t(matrix1),t(matrix2))
#' Z <- PrepareSummarizedExperimentTensorRect(sample=as.character(seq_len(50)),
#' feature=list(as.character(seq_len(200)),as.character(seq_len(400))),
#' sampleData=list(rep(seq_len(2),each=25)),value=Z)
#' HOSVD <- computeHosvd(Z)
#' cond <- list(attr(Z,"sampleData")[[1]],NULL,NULL)
#' index_all <- selectFeatureTransRect(HOSVD,cond,de=c(0.01,0.01),input_all=2)
selectFeatureTransRect <- function(HOSVD,cond,de=rep(1e-4,2),p0=0.01,
                                   breaks=100,input_all=NULL){
    {
        interact<-FALSE
        if(is.null(input_all))
        {
            interact<-TRUE
            j<-1
            while(j %in% seq_len(dim(HOSVD$U[[1]])[2]))
            {
                boxplot(HOSVD$U[[1]][,j]~cond[[1]],main=j)
                abline(0,0,col=2,lty=2)
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
        index_all <- rep(list(NA))
        for (i in seq_len(2))
        {
            u<- HOSVD$U[[i+1]][,input_all]
            sd <- optim(de[i],function(x){th(x,breaks,p0)},
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
            if (interact) readline("Press Enter to proceed:")
            index <- p.adjust(P2,"BH")<p0
            index_all[[i]] <- list(index=index,p.value=P2)
        }
        return(index_all)
    }
}