library(UCSCXenaTools)
library(data.table)
library(R.utils)
library(readr)
library(stringr)
library(skater)
library(caret)
library(tidyverse)

# Survival loading Data codes and prep------------
TCGA_survival_data<-fread("TCGA_survival_data") #Surv data
ExprSubsetBySamp <- read_csv("ExprSubsetBySamp.csv")

Surv_TCGA_use = subset(TCGA_survival_data, sample %in% colnames(ExprSubsetBySamp))


##One yaer---------
Surv_TCGA_use1<- Surv_TCGA_use %>% 
  mutate(one_year=ifelse(OS.time>=365,1,
                         ifelse(OS.time<365 & OS==1,0,99))) %>% 
  filter(one_year!=99) %>% dplyr::select(sample,one_year)

##Two yaer---------
#Surv_TCGA_use2<- Surv_TCGA_use %>% 
#  mutate(two_year=ifelse(OS.time>=730,1,
 #                        ifelse(OS.time<730 & OS==1,0,99))) %>% 
#  filter(two_year!=99)


prop.table(table(Surv_TCGA_use1$one_year))
#prop.table(table(Surv_TCGA_use2$two_year))


##Preparing the expr data-----
ExprSubsetBySamp_surv<-as.data.frame(subset(ExprSubsetBySamp,select =Surv_TCGA_use1$sample ))

rownames(ExprSubsetBySamp_surv) <- ExprSubsetBySamp$sample


ExprSubsetBySamp_surv$sample=ExprSubsetBySamp$sample

## map transcripts to gene symbols-----------------


NCBI_sample_surv = ExprSubsetBySamp_surv$sample  


NCBI_sample_surv <- str_split_fixed(as.character(NCBI_sample_surv), "\\.", 2)[,1]


# You can see that the subset data is listed as a integer, we now need to convert
# this to a vector to pass it into the annotation mapping

NCBI_sample_surv = as.vector(NCBI_sample_surv)
#names(data) <- c("ENSEMBL")

ExprSubsetBySamp_surv$ENSEMBL <- NCBI_sample_surv

# Using the org.Hs.eg.db we set up mapping info - if you look at the documentation you
# can also obtain other keytypes

library(org.Hs.eg.db)

annots_surv <- AnnotationDbi::select(org.Hs.eg.db, keys=NCBI_sample_surv, columns="SYMBOL", keytype="ENSEMBL")

head(annots_surv)

result_surv <- merge(ExprSubsetBySamp_surv, annots_surv, by.x="ENSEMBL", by.y="ENSEMBL")

result1_surv <- result_surv[which(!is.na(result_surv$SYMBOL)) , ]
ExpDD_geneSym_surv <- as.matrix(result1_surv[, c(2:132)])
rownames(ExpDD_geneSym_surv) <- as.vector(result1_surv$SYMBOL)

dim(ExpDD_geneSym_surv)


X_surv <- as.matrix(t(ExpDD_geneSym_surv))

#Inside mlda model-------


##Edge R library-------
library(edgeR)

group_surv <- Surv_TCGA_use1$one_year
y_surv <- DGEList(counts=ExpDD_geneSym_surv,group=group_surv)

keep_surv <- filterByExpr(y_surv)
y_surv <- y_surv[keep_surv,,keep.lib.sizes=FALSE]
y_surv <- calcNormFactors(y_surv)
design_surv <- model.matrix(~group_surv)
y_surv <- estimateDisp(y_surv,design_surv)


fit_surv <- glmQLFit(y_surv,design_surv)
qlf_surv <- glmQLFTest(fit_surv,coef=2)
topTags(qlf_surv)


qlfT_surv <- qlf_surv$table
qlfT_sort_surv <- qlfT_surv[order(qlfT_surv[,4]),]

#qlfT_sort[1:10,]

#write.csv(qlfT_sort_surv,"3.Outputs/selected_genes_surv.csv")

edgeR.genes_surv <- rownames(qlfT_sort_surv)
length(edgeR.genes_surv)


edgeR.genes_surv <- sapply(rownames(qlfT_sort_surv), function(x) unlist(strsplit(x, "\\."))[[1]][1]) 

edgeR.set_surv <- match(edgeR.genes_surv, colnames(X_surv))

which(is.na(edgeR.set_surv))

#X_edgeR_keep <- X[,edgeR.set[-which(is.na(edgeR.set))]]
X_edgeR_keep_surv <- X_surv[,edgeR.set_surv]      ###???
dim(X_edgeR_keep_surv)
#write.csv(X_edgeR_keep_surv,"X_edgeR_keep_surv.csv",row.names = F)

Xsd_surv <- apply(X_edgeR_keep_surv, 2, sd)
length(which(Xsd_surv==0))

##Mlda model---------


Xnew_surv <- X_edgeR_keep_surv  ##????
Xnew_surv <- apply(Xnew_surv, 2, function(x) (x-mean(x))/sd(x))
#View(Xnew)

dim(Xnew_surv)

### Comment: Need to add in covariates such as age, gender,etc.

Y_surv <- as.vector(Surv_TCGA_use1$one_year)

source("library_CIS_imagingData_01042023.R")


###Cross validation-----------

#k=1
#for (tau in c(20,30,50)) {

# for (nu in c(150,200,250)) {

tau=50
nu=300
system.time(try.cv_surv <- mLDA.cv(Xnew_surv, Y_surv, Z=NULL, fold=5,seed=1 ,
                                     
                                   tau=tau, alpha=0.9, nu=nu, d=3, nb=10))
try.cv_surv$err.tot
Y.cv.orig_surv <- as.vector(unlist(try.cv_surv$Y.test.arr))
Y.cv.pred_surv <- as.vector(unlist(try.cv_surv$Y.pred.arr))
(cm_surv<-confusionMatrix(as.factor(Y.cv.pred_surv), as.factor(Y.cv.orig_surv)))

#  asurvu[k]<-cm[[3]][1]
#  k=k+1
# }
#}  


cm_surv[[3]][1]

###final model------

system.time(try2_new_surv  <- mLDA.pair(Xnew_surv, Y_surv, 
                                        Xnew_surv, Z=NULL, Z.new=NULL, 
                                        pair=c(1,0), strong.X.set=c(1:100), tau=tau, 
                                        alpha=0.9, nu=nu, d=3, nb=10))

pred_surv <- as.numeric(try2_new_surv$PredClass)
length(which(pred_surv!=Y_surv))

table(pred_surv,Y_surv)
confusionMatrix(as.factor(pred_surv), as.factor(Y_surv))


saveRDS(try2_new_surv,file=paste("survival_mLDA", ".rds", sep=""))
#try2 <- readRDS(file=paste("casecontrol_mLDA", ".rds", sep=""))

pred_data<-data.frame(pred_surv,sample=rownames(Xnew_surv))
#write.csv(pred_data,"./3.Outputs/Tables/pred_table_surv.csv")
#write.csv(Xnew_surv,"./3.Outputs/Tables/X_surv.csv")


#ROC----------
library(pROC)
fisher_surv <- try2_new_surv$FisherDR
pred_prob_surv <- (fisher_surv-min(fisher_surv))/(max(fisher_surv)-min(fisher_surv))

roc_obj <- roc(Y_surv, pred_prob_surv)
plot(roc_obj, main="ROC Curve", col="#1c61b6")
auc_value <- auc(roc_obj)
print(paste("AUC:", auc_value))


