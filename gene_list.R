library(tidyverse)
library(readr)
library(readxl)
#casecontrol----------
try_casecontrol <- readRDS(file=paste("casecontrol_mLDA", ".rds", sep=""))
View(try_casecontrol)


##Differential Effect

##gene analysis---------
network_nodes_cc<-as.data.frame(try_casecontrol$network.nodes)
names(network_nodes_cc)[1]<-"nodes"
network_nodes_cc$genes<-rownames(network_nodes_cc)

totGene_cc <- data.frame(nodes=try_casecontrol$screenset)
totGene_cc<-totGene_cc %>% inner_join(network_nodes_cc,by="nodes")

strongGene_cc <- data.frame(nodes=try_casecontrol$MIset)
strongGene_cc<-strongGene_cc %>% inner_join(network_nodes_cc,by="nodes") 
strongGene_cc$type<-rep("strong",times=nrow(strongGene_cc))

weakGene_cc <- data.frame(nodes=setdiff(totGene_cc$nodes, strongGene_cc$nodes),
  genes=setdiff(totGene_cc$genes, strongGene_cc$genes))
weakGene_cc$type<-rep("weak",times=nrow(weakGene_cc))

##log fold change and p-value accumulation (from case control script)----------
log_fold_genes<-data.frame(genes=rownames(qlfT_sort_cc), qlfT_sort_cc)
rownames(log_fold_genes)<-NULL

strong_final<-log_fold_genes %>% inner_join(strongGene_cc[,-3],by="genes") 

weak_final<-log_fold_genes %>% inner_join(weakGene_cc[,-3],by="genes")

totGene_cc_00<- rbind(strongGene_cc,weakGene_cc)
totGene_cc_final<- totGene_cc %>% inner_join(totGene_cc_00,by=c("nodes","genes")) %>% 
                            inner_join(log_fold_genes,by="genes")

#write.csv(strong_final,"./3.Outputs/strong_cc.csv",row.names=F)
#write.csv(weak_final,"./3.Outputs/weak_cc.csv",row.names=F)
#write.csv(totGene_cc_final,"./3.Outputs/total_genes_cc.csv",row.names=F)

##Extracting network nodes-----------

casecontrol_mLDA_network <- NULL

for (i in 1:29) {
  temp <- as.data.frame(try_casecontrol[["local.network_list"]][[i]])
  temp$genes <- rownames(temp)
  temp$net <- paste("Network",i, sep = "")
  temp <- temp[,-1]
  casecontrol_mLDA_network <- rbind(casecontrol_mLDA_network,temp)
}
casecontrol_mLDA_network<-casecontrol_mLDA_network %>% 
  inner_join(totGene_cc_final[,c(2,3)],by="genes")

View(casecontrol_mLDA_network %>% group_by(net) %>% count() %>% arrange(desc(n)))

c(casecontrol_mLDA_network %>% filter(net=="Network2") %>% dplyr::select('genes'))

#write.csv(casecontrol_mLDA_network,"./3.Outputs/casecontrol_mLDA_network.csv",row.names=F)

##extracting DE values------
DE_CC0<-data.frame(DE=try_casecontrol$iffcond, 
                  nodes=as.numeric(rownames(data.frame(try_casecontrol$iffcond))))
DE_CC<-DE_CC0 %>% inner_join(network_nodes_cc,by="nodes") %>% 
  inner_join(totGene_cc_final[,c(2,3)],by="genes")

#write.csv(DE_CC,"./3.Outputs/DE_CC.csv",row.names=F)

##top 10 DE gene------

(de_1<-DE_CC[1:20,] %>%  ggplot() +
  geom_bar(aes(x=reorder(genes, abs(DE)),y=DE,fill=type),stat = 'identity')+
  coord_flip()+
  theme_minimal() +
  labs(x = "Genes", y = "DE",fill="Gene type") +
  theme(text = element_text(size = 12, face = "bold"),
    legend.position = 'bottom')+ggtitle("(A)"))

#surv----------
try_surv <- readRDS(file=paste("survival_mLDA", ".rds", sep=""))


##gene analysis---------
network_nodes_surv<-as.data.frame(try_surv$network.nodes)
names(network_nodes_surv)[1]<-"nodes"
network_nodes_surv$genes<-rownames(network_nodes_surv)

totGene_surv <- data.frame(nodes=try_surv$screenset)
totGene_surv<-totGene_surv %>% full_join(network_nodes_surv,by="nodes")


strongGene_surv <- data.frame(nodes=try_surv$MIset)
strongGene_surv<-strongGene_surv %>% inner_join(network_nodes_surv,by="nodes")
strongGene_surv$type<-rep("strong",times=nrow(strongGene_surv))


weakGene_surv <- data.frame(nodes=setdiff(totGene_surv$nodes, strongGene_surv$nodes),
                          genes=setdiff(totGene_surv$genes, strongGene_surv$genes))
weakGene_surv$type<-rep("weak",times=nrow(weakGene_surv))



##log fold change and p-value accumulation (from case control script)----------
log_fold_genes<-data.frame(genes=rownames(qlfT_sort_surv), qlfT_sort_surv)
rownames(log_fold_genes)<-NULL

strong_final<-log_fold_genes %>% inner_join(strongGene_surv,by="genes") 
weak_final<-log_fold_genes %>% inner_join(weakGene_surv,by="genes")


totGene_surv_00<- rbind(strongGene_surv,weakGene_surv)
totGene_surv_final<- totGene_surv_00 %>% 
  inner_join(log_fold_genes,by="genes")



#write.csv(strong_final,"./3.Outputs/strong_surv.csv",row.names=F)
#write.csv(weak_final,"./3.Outputs/weak_surv.csv",row.names=F)
#write.csv(totGene_surv_final,"./3.Outputs/total_surv.csv",row.names=F)

##Extracting network nodes-----------

surv_mLDA_network <- NULL

for (i in 1:48) {
  temp <- as.data.frame(try_surv[["local.network_list"]][[i]])
  temp$genes <- rownames(temp)
  temp$net <- paste("Network",i, sep = "")
  temp <- temp[,-1]
  surv_mLDA_network <- rbind(surv_mLDA_network,temp)
}
surv_mLDA_network<-surv_mLDA_network %>% 
  inner_join(totGene_surv_final[,c(2,3)],by="genes")

View(surv_mLDA_network %>% group_by(net) %>% count() %>% arrange(desc(n)))

write.csv(surv_mLDA_network,"./3.Outputs/surv_mLDA_network.csv",row.names=F)

##extracting DE values------
DE_surv0<-data.frame(DE=try_surv$iffcond, 
                   nodes=as.numeric(rownames(data.frame(try_surv$iffcond))))
DE_surv<-DE_surv0 %>% full_join(network_nodes_surv,by="nodes") %>% 
  full_join(totGene_surv_final[,c(2,3)],by="genes")

#write.csv(DE_surv,"./3.Outputs/DE_surv.csv",row.names=F)
mergede<-totGene_surv_final %>% full_join(DE_surv[c(1,3)],by="genes")

##top 10 DE gene------

(de_2<-DE_surv[1:20,] %>%  ggplot() +
  geom_bar(aes(x=reorder(genes, abs(DE)),y=DE,fill=type),stat = 'identity')+
  coord_flip()+
  theme_minimal() +
  labs(x = "Genes", y = "DE",fill="Gene type") +
  theme(text = element_text(size = 12, face = "bold"),
        legend.position = 'bottom')+ggtitle("(B)"))


ggarrange(de_1,de_2,common.legend = T)
#subtype----------
try_subtype <- readRDS(file=paste("subtype_mLDA", ".rds", sep=""))


##gene analysis---------
network_nodes_subtype_new<-data.frame()
for (i in 1:6) {
network_nodes_subtype<-as.data.frame(try_subtype$network.nodes_list[i])
names(network_nodes_subtype)[1]<-"nodes"
network_nodes_subtype$genes<-rownames(network_nodes_subtype)
rownames(network_nodes_subtype)<-c()
network_nodes_subtype_new<- rbind(network_nodes_subtype_new,network_nodes_subtype)
}
network_nodes_subtype_new2<-unique(network_nodes_subtype_new)

totGene_subtype <- data.frame(nodes=try_subtype$screenset)
totGene_subtype<-totGene_subtype %>% inner_join(network_nodes_subtype_new2,by="nodes")


strongGene_subtype <- data.frame(nodes=c(try_subtype$MImatrix.cord))
strongGene_subtype2 <- unique(strongGene_subtype)

strongGene_subtype2<-strongGene_subtype2 %>% 
  inner_join(network_nodes_subtype_new2,by="nodes")


weakGene_subtype <- data.frame(nodes=setdiff(totGene_subtype$nodes, strongGene_subtype2$nodes),
                            genes=setdiff(totGene_subtype$genes, strongGene_subtype2$genes))


for (j in 1:19) {
  temp <- as.data.frame(try_subtype[["local.network_list"]][[j]])
  temp$gene <- rownames(temp)
  temp$net <- paste("Network",j, sep = "")
  temp <- temp[,-1]
  subtypeival_mLDA_network <- rbind(subtypeival_mLDA_network,temp)
}

##ven diagram-----------

# Install and load the VennDiagram package
install.packages("VennDiagram")
library(ggVennDiagram)

# Create list


x <- list(
  "Strong gene disease/non-disease" = strongGene_cc$genes, 
  "Strong gene survival" = strongGene_surv$genes)
   
intersect(strongGene_cc$genes,strongGene_surv$genes)
# Create a Venn diagram_strong gene
(vennc_V<-ggVennDiagram(x, label_alpha = 0,label_color  = "white",
              label_size=10,set_size =0)+
  theme(legend.position = ""))

ggsave(vennc_V,"cc_vs_surv.png")

x <- list(
  strong_casecontrol = strongGene_cc$genes, 
  strong_surv = strongGene_surv$genes,
  strong_subtype = strongGene_subtype2$genes) 

# Create a Venn diagram_strong gene
ggVennDiagram(x, label_alpha = 0,label_color  = "white",set_size = 2.5)+
  theme(legend.position = "")


strog_gene<-data.frame(case_control=strongGene_cc$genes,survival= strongGene_surv$genes)
write.csv(strog_gene,"gene_strong22.csv")



# Create list weak


y <- list(
  "Weak gene disease/non-disease" = weakGene_cc$genes, 
  "Weak gene survival" = weakGene_surv$genes)

intersect(weakGene_cc$genes,weakGene_surv$genes)
# Create a Venn diagram_strong gene
(vennc_weak<-ggVennDiagram(y, label_alpha = 0,label_color  = "white",
                        label_size=10,set_size =0)+
    theme(legend.position = ""))




#pathview
#make subtype binary from literature
#miset 1st column category 1 and 2, second column 1 and 3....



##sig_genes_both_cc_surv

sig_surv<-read.csv("./3.Outputs/p_values_df_sig.csv")




intersect(strongGene_cc$genes,sig_surv$genes)
intersect(weakGene_cc$genes,sig_surv$genes)

intersect(strongGene_surv$genes,sig_surv$genes)

intersect(weakGene_surv$genes,sig_surv$genes)

