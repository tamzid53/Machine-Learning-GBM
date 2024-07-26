library(pathview)
library(readr)

##case control------
Xnew_cc <- read_csv("3.Outputs/Tables/Xnew_cc.csv")
total_genes_cc <- read_csv("3.Outputs/total_genes_cc.csv")

Xnew_cc_2<-Xnew_cc[,c("sample",total_genes_cc$genes)]

library(data.table)

X_path<-transpose(Xnew_cc_2)
View(X_path)
#redefine row and column names
rownames(X_path) <- colnames(Xnew_cc_2)
colnames(X_path) <- Xnew_cc_2$sample
X_path<-X_path[-1,]
X_path<-total_genes_cc[,c("genes","logFC")]
X_path<-total_genes_cc$logFC

names(X_path)<-total_genes_cc$genes

View(X_path)

pv.out <- pathview(gene.data = X_path, pathway.id="hsaT30436" , gene.idtype = "SYMBOL",
                   species = "hsa", out.suffix = "cc", kegg.native = T)


##survival------

X_surv <- read_csv("3.Outputs/Tables/X_surv.csv")
total_surv <- read_csv("3.Outputs/total_surv.csv")
X_surv_2<-X_surv[,c("sample",total_surv$genes)]

library(data.table)

X_path<-transpose(X_surv_2)
View(X_path)
#redefine row and column names
rownames(X_path) <- colnames(X_surv_2)
colnames(X_path) <- X_surv_2$sample
X_path<-X_path[-1,]
X_path<-total_surv[,c("genes","logFC")]
X_path<-total_surv$logFC

names(X_path)<-total_surv$genes

View(X_path)

pv.out <- pathview(gene.data = X_path, pathway.id="hsa03040" , gene.idtype = "SYMBOL",
                   species = "hsa", out.suffix = "surv", kegg.native = T)



