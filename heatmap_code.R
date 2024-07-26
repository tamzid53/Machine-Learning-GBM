
library(readr)
library(tidyverse)
library(pheatmap)
##surv_heatmap--------------
X_edgeR_keep_surv<-read.csv("./3.Outputs/X_edgeR_keep_surv.csv")
Surv_TCGA_use_kp<-read.csv("./3.Outputs/Tables/Surv_TCGA_use_kp.csv")

X_edgeR_keep_surv<-X_edgeR_keep_surv %>% 
  inner_join(Surv_TCGA_use_kp[,2:3],by=c("X"="sample"))

X_edgeR_keep_surv_new<-X_edgeR_keep_surv[,-1]

X_edgeR_keep_surv_new<-X_edgeR_keep_surv_new %>% arrange(one_year)

surv_mLDA_network0<-read.csv("./3.Outputs/surv_mLDA_network.csv")

surv_mLDA_network0$genes <- gsub("-", ".", surv_mLDA_network0$genes)

surv_mLDA_network<-surv_mLDA_network %>% filter(net %in% c("Network41","Network3"))

X_edgeR_keep_surv_new0<-X_edgeR_keep_surv_new %>% dplyr::select(c(surv_mLDA_network$genes))

X_edgeR_keep_surv_new2 <- t(X_edgeR_keep_surv_new0)

# Convert the transposed matrix back to a dataframe
X_edgeR_keep_surv_new2 <- as.data.frame(X_edgeR_keep_surv_new2, stringsAsFactors = FALSE)

X_edgeR_keep_surv_new2$genes<-rownames(X_edgeR_keep_surv_new2)
rownames(X_edgeR_keep_surv_new2)<-c()

X_edgeR_keep_surv_new2<-X_edgeR_keep_surv_new2 %>% inner_join(surv_mLDA_network,by=c("genes"))

X_edgeR_keep_surv_new2$Label <-
  paste(X_edgeR_keep_surv_new2$genes, X_edgeR_keep_surv_new2$net, sep = "-")

ann_color<-list(
  "type"=c('strong'='turquoise','weak'="#8B008B"),
  "status"=c('>=1- year'='green','< 1-year'='red')
  )



ann_row2<-data.frame(X_edgeR_keep_surv_new2[,c("type")])
ann_col2<-data.frame(X_edgeR_keep_surv_new$one_year)


ann_row2[,1]<-as.factor(ann_row2[,1])
ann_col2[,1]<-as.factor(ann_col2[,1])

colnames(ann_row2)<-c("type")
rownames(ann_row2)<-X_edgeR_keep_surv_new2$Label

colnames(ann_col2)<-c("status")
ann_col2$status<-ifelse(ann_col2$status=='0',"< 1-year",">=1- year")
rownames(ann_col2)<-colnames(X_edgeR_keep_surv_new2[,1:131])

#mat_work<-as.matrix(data_Salman_Coded[,c(11:38)])
mat_work<-as.matrix(X_edgeR_keep_surv_new2[,1:131])
rownames(mat_work)<-X_edgeR_keep_surv_new2$Label


(p_surv<-pheatmap(mat=mat_work,
         scale="none",
         cluster_cols=F,
         cluster_rows =F,
         clustering_method = "ward.D2",
         col = colorRampPalette(c("steelblue", "orange"))(1024),
         annotation_row   = ann_row2,
         annotation_col = ann_col2,
         show_colnames = F,
         annotation_colors = ann_color,
         legend = T,cellwidth = 6,cellheight =9,margins=c(40,20),fontsize=10))

ggsave("./1.Manuscript/Images/heatmap_surv2.tiff", plot = p_surv, width =15 , height =10 , 
       dpi = 300)

##casecontrol_heatmap--------------

X_edgeR_keep_new<-read.csv("./3.Outputs/X_edgeR_keep.csv")
casecontrol_mLDA_network0<-read.csv("./3.Outputs/casecontrol_mLDA_network.csv")

casecontrol_mLDA_network0$genes <- gsub("-", ".", casecontrol_mLDA_network0$genes)

casecontrol_mLDA_network<-casecontrol_mLDA_network0 %>% 
  filter(net %in% c("Network2","Network26","Network15"))

X_edgeR_keep_new2<-X_edgeR_keep_new %>% dplyr:: select(casecontrol_mLDA_network$genes)

X_edgeR_keep_new2 <- t(X_edgeR_keep_new2)

# Convert the transposed matrix back to a dataframe
X_edgeR_keep_new2 <- as.data.frame(X_edgeR_keep_new2, stringsAsFactors = FALSE)

X_edgeR_keep_new2$genes<-rownames(X_edgeR_keep_new2)
rownames(X_edgeR_keep_new2)<-c()

X_edgeR_keep_new2<-X_edgeR_keep_new2 %>% inner_join(casecontrol_mLDA_network,by=c("genes"))

X_edgeR_keep_new2$Label <- paste(X_edgeR_keep_new2$genes, X_edgeR_keep_new2$net, sep = "-")

ann_color<-list(
  "type"=c('strong'='turquoise','weak'="#8B008B"),
  "status"=c('control'='green','case'='red')
  
)


ann_row2<-data.frame(X_edgeR_keep_new2[,c("type")])
ann_col2<-data.frame(status=c(rep("control",times=1148),rep("case",times=153)))


ann_row2[,1]<-as.factor(ann_row2[,1])
ann_col2[,1]<-as.factor(ann_col2[,1])
#ann_row2[,3]<-as.factor(ann_row2[,3])

rownames(ann_row2)<-X_edgeR_keep_new2$Label
colnames(ann_row2)<-c("type")

colnames(ann_col2)<-c("status")
rownames(ann_col2)<-colnames(X_edgeR_keep_new2[,1:1301])




mat_work2<-as.matrix(X_edgeR_keep_new2[,1:1301])
rownames(mat_work2)<-X_edgeR_keep_new2$Label

(p_cc<-pheatmap(mat=mat_work2,
         scale="none",
         cluster_cols=F,
         cluster_rows =F,
         show_colnames = F,
         clustering_method = "ward.D2",
         col = colorRampPalette(c("steelblue", "orange"))(1024),
         annotation_row = ann_row2,
         annotation_col  = ann_col2 ,
         annotation_colors = ann_color,
         legend = T,cellwidth = 1,cellheight =7,margins=c(60,60),fontsize = 8))

ggsave("./1.Manuscript/Images/heatmap_cc2.tiff", plot = p_cc, width =22 , height =12 , 
       dpi = 300)




##casecontrol_heatmap--------------

X_edgeR_keep_new<-read.csv("./3.Outputs/X_edgeR_keep.csv")
casecontrol_mLDA_network0<-read.csv("./3.Outputs/casecontrol_mLDA_network.csv")

casecontrol_mLDA_network0$genes <- gsub("-", ".", casecontrol_mLDA_network0$genes)

casecontrol_mLDA_network<-casecontrol_mLDA_network0 %>% 
  filter(net %in% c("Network21","Network12","Networ28",
                    "Network20","Network17","Network18"))

X_edgeR_keep_new2<-X_edgeR_keep_new %>% dplyr:: select(casecontrol_mLDA_network$genes)

X_edgeR_keep_new2 <- t(X_edgeR_keep_new2)

# Convert the transposed matrix back to a dataframe
X_edgeR_keep_new2 <- as.data.frame(X_edgeR_keep_new2, stringsAsFactors = FALSE)

X_edgeR_keep_new2$genes<-rownames(X_edgeR_keep_new2)
rownames(X_edgeR_keep_new2)<-c()

X_edgeR_keep_new2<-X_edgeR_keep_new2 %>% inner_join(casecontrol_mLDA_network,by=c("genes"))

X_edgeR_keep_new2$Label <- paste(X_edgeR_keep_new2$genes, X_edgeR_keep_new2$net, sep = "-")

ann_color<-list(
  "type"=c('strong'='turquoise','weak'="#8B008B"),
  "status"=c('control'='green','case'='red')
  
)


ann_row2<-data.frame(X_edgeR_keep_new2[,c("type")])
ann_col2<-data.frame(status=c(rep("control",times=1148),rep("case",times=153)))


ann_row2[,1]<-as.factor(ann_row2[,1])
ann_col2[,1]<-as.factor(ann_col2[,1])
#ann_row2[,3]<-as.factor(ann_row2[,3])

rownames(ann_row2)<-X_edgeR_keep_new2$Label
colnames(ann_row2)<-c("type")

colnames(ann_col2)<-c("status")
rownames(ann_col2)<-colnames(X_edgeR_keep_new2[,1:1301])




mat_work2<-as.matrix(X_edgeR_keep_new2[,1:1301])
rownames(mat_work2)<-X_edgeR_keep_new2$Label

(p_cc<-pheatmap(mat=mat_work2,
                scale="none",
                cluster_cols=F,
                cluster_rows =F,
                show_colnames = F,
                clustering_method = "ward.D2",
                col = colorRampPalette(c("steelblue", "orange"))(1024),
                annotation_row = ann_row2,
                annotation_col  = ann_col2 ,
                annotation_colors = ann_color,
                legend = T,cellwidth = 1,cellheight =7,margins=c(60,60),fontsize = 8))

ggsave("./1.Manuscript/Images/heatmap_cc3.tiff", plot = p_cc, width =22 , height =12 , 
       dpi = 300)
