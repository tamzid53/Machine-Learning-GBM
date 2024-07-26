library(UCSCXenaTools)
library(data.table)
library(R.utils)
library(tidyverse)
library(readr)
library(stringr)
library(skater)
library(caret)

ExprSubsetBySamp <- read_csv("ExprSubsetBySamp.csv")

## case-control classification analysis-----------


sample_vec_cc <- sapply(colnames(ExprSubsetBySamp), function(x) unlist(strsplit(x, "-"))[[1]][1]) 
table(sample_vec_cc)

caseControlResponse <- rep(0, (dim(ExprSubsetBySamp)[2]-1))
caseControlResponse[which(sample_vec_cc == "TCGA")] <- 1


#ExprDD <- as.data.frame(ExprSubsetBySamp[,1:(dim(ExprSubsetBySamp)[2]-1)])
#rownames(ExprDD) <- ExprSubsetBySamp$sample
#dim(ExprSubsetBySamp)

## map transcripts to gene symbols-----------


NCBI_sample_cc = ExprSubsetBySamp$sample  


NCBI_sample_cc <- str_split_fixed(as.character(NCBI_sample_cc), "\\.", 2)[,1]


# You can see that the subset data is listed as a integer, we now need to convert
# this to a vector to pass it into the annotation mapping

NCBI_sample_cc = as.vector(NCBI_sample_cc)
#names(data) <- c("ENSEMBL")

ExprSubsetBySamp$ENSEMBL <- NCBI_sample_cc

# Using the org.Hs.eg.db we set up mapping info - if you look at the documentation you
# can also obtain other keytypes

library(org.Hs.eg.db)

annots_cc <- AnnotationDbi::select(org.Hs.eg.db, keys=NCBI_sample_cc, columns="SYMBOL", keytype="ENSEMBL")

#head(annots)

result_cc <- merge(ExprSubsetBySamp, annots_cc, by.x="ENSEMBL", by.y="ENSEMBL")

result1_cc <- result_cc[which(!is.na(result_cc$SYMBOL)) , ]
ExpDD_geneSym_cc <- as.matrix(result1_cc[, c(2:1302)])
rownames(ExpDD_geneSym_cc) <- as.vector(result1_cc$SYMBOL)

dim(ExpDD_geneSym_cc)


X_cc <- as.matrix(t(ExpDD_geneSym_cc))
View(X_cc)


#Modeling Phase started-----


## Outcome-1: case-control classification analysis ----

#### single gene based DE analysis results as strong gene set
#### 

##EdgeR---------
library(edgeR)


group_cc <- caseControlResponse
y_cc <- DGEList(counts=ExpDD_geneSym_cc,group=group_cc)

keep_cc <- filterByExpr(y_cc)  #### https://rdrr.io/bioc/edgeR/man/filterByExpr.html
y_cc <- y_cc[keep_cc,,keep.lib.sizes=FALSE]
y_cc <- calcNormFactors(y_cc)
design_cc <- model.matrix(~group_cc)
y_cc <- estimateDisp(y_cc,design_cc)


fit_cc <- glmQLFit(y_cc,design_cc)
qlf_cc <- glmQLFTest(fit_cc,coef=2)
topTags(qlf_cc)


qlfT_cc <- qlf_cc$table
qlfT_sort_cc <- qlfT_cc[order(qlfT_cc[,4]),]

#qlfT_sort$fdr <- p.adjust(qlfT_sort[,4], method ="fdr")


#write.csv(qlfT_sort_cc,"3.Outputs/selected_genes_cc.csv")

edgeR.genes_cc <- rownames(qlfT_sort_cc)
length(edgeR.genes_cc)


edgeR.genes_cc <- sapply(rownames(qlfT_sort_cc), 
                         function(x) unlist(strsplit(x, "\\."))[[1]][1]) 

edgeR.set_cc <- match(edgeR.genes_cc, colnames(X_cc))

which(is.na(edgeR.set_cc))

#X_edgeR_keep <- X[,edgeR.set[-which(is.na(edgeR.set))]]
X_edgeR_keep_cc <- X_cc[,edgeR.set_cc]      ###???
dim(X_edgeR_keep_cc)


Xsd_cc <- apply(X_edgeR_keep_cc, 2, sd)
length(which(Xsd_cc==0))


##Mlda model---------

Xnew_cc <- X_edgeR_keep_cc ##????
Xnew_cc <- apply(Xnew_cc, 2, function(x) (x-mean(x))/sd(x))
#View(Xnew_cc)

dim(Xnew_cc)

Y_cc <- as.vector(caseControlResponse)

source("library_CIS_imagingData_01042023.R")

###Cross validation-----------

#k=1
#for (tau in c(20,30,50)) {

# for (nu in c(150,200,250)) {

tau=50
nu=300
system.time(try.cv_cc <- mLDA.cv(Xnew_cc, Y_cc, Z=NULL, fold=5, seed=1, 
                                 strong.X.list=list(c(1:100)), 
                                 tau=tau, alpha=0.9, nu=nu, d=3, nb=10))
try.cv_cc$err.tot
Y.cv.orig_cc <- as.vector(unlist(try.cv_cc$Y.test.arr))
Y.cv.pred_cc <- as.vector(unlist(try.cv_cc$Y.pred.arr))
(cm_cc<-confusionMatrix(as.factor(Y.cv.pred_cc), as.factor(Y.cv.orig_cc)))

#  accu[k]<-cm[[3]][1]
#  k=k+1
# }
#}  


cm_cc[[3]][1]

###final model------

system.time(try2_new_cc  <- mLDA.pair(Xnew_cc, Y_cc, Xnew_cc, 
                                      Z=NULL, Z.new=NULL, 
                                      pair=c(1,0), 
                                      strong.X.set=c(1:100), 
                                      tau=50, alpha=0.9, nu=300, d=3, nb=10))

pred_cc <- as.numeric(try2_new_cc$PredClass)
length(which(pred_cc!=Y_cc))

confusionMatrix(as.factor(pred_cc), as.factor(Y_cc))

`saveRDS(try2_new_cc,file=paste("casecontrol_mLDA", ".rds", sep=""))
#try2 <- readRDS(file=paste("casecontrol_mLDA", ".rds", sep=""))

#write.csv(Xnew_cc,"./3.Outputs/Tables/Xnew_cc.csv")

#ROC curve---------
fisher_cc <- try2_new_cc$FisherDR
pred_prob_cc <- (fisher_cc-min(fisher_cc))/(max(fisher_cc)-min(fisher_cc))


library(ROCR)

pred_cc <- ROCR::prediction(pred_prob_cc, Y_cc)
perf1 <- ROCR::performance( pred_cc, "tpr", "fpr" )
auc_ROCR1 <- ROCR::performance(pred_cc, measure="auc")
auc_cc <- auc_ROCR1@y.values[[1]]
auc_cc



par(mar=c(5,5,2,1) + .1)
plot(1, type="n", col="red", ylim=c(0.6, 1), xlim=c(0, 0.8), xlab="False positive rate", ylab="True positive rate", cex.lab=1.8, cex.axis=1.5)
plot(perf1, add = TRUE, col=ROC_col[1], lwd=3, lty=1, alpha=0.6)
     
     
library(pROC)     
roc_obj <- roc(Y_cc, pred_prob_cc)
plot(roc_obj, main="ROC Curve", col="#1c61b6")
auc_value <- auc(roc_obj)
print(paste("AUC:", auc_value))
