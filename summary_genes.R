library(readr)
library(tidyverse)
# Create a data frame with edgeR results (logFC and FDR columns)

selected_genes_cc0 <- read_csv("3.Outputs/selected_genes_cc.csv")

total_genes_cc<- read_csv("3.Outputs/total_genes_cc.csv")

selected_genes_cc<-selected_genes_cc0 %>% left_join(total_genes_cc[,1:2],by="genes")
# Create a volcano plot using ggplot2
##cc----------
selected_genes_cc<-selected_genes_cc %>% 
  mutate(cc=ifelse(logFC>1 & PValue<0.05 ,"Upregulated",
                     ifelse(logFC< -1 & PValue<0.05,"Downregulated",
                     
                             "Not Significant")))

selected_genes_cc<-selected_genes_cc %>% 
  mutate(select=ifelse(is.na(nodes),"not selected","selected by model"))
cc_reg<-ggplot(selected_genes_cc, aes(x = logFC, y = -log10(PValue),color=select)) +
  geom_point(alpha=0.5) +
  #geom_text(data = subset(selected_genes_cc,cc %in% c("Upregulated","Downregulated")),
  #          aes(label=genes),vjust=1.5,size=3)+
  scale_color_manual(values = c("steelblue", "red"), name = "selection status") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange",linewidth=1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "orange",linewidth=1) +
  theme_bw() +
  labs(x = "Log2 fold change", y = "-log10 (Pvalue)", title = "(A)")

ggsave("./3.Outputs/Images/volcano_cc.tiff",width = 7,height = 4)

selected_genes_cc %>% group_by(cc) %>% count()


##surv----------

selected_genes_surv0 <- read_csv("3.Outputs/selected_genes_surv.csv")
total_genes_surv<- read_csv("3.Outputs/total_surv.csv")

selected_genes_surv<-selected_genes_surv0 %>% left_join(total_genes_surv[,1:2],by="genes")

selected_genes_surv<-selected_genes_surv %>% 
  mutate(surv=ifelse(logFC>1 & PValue<0.05,"Upregulated",
                   ifelse(logFC< -1 & PValue<0.05,"Downregulated","Not Significant")))

selected_genes_surv<-selected_genes_surv %>% 
  mutate(select=ifelse(is.na(nodes),"not selected","selected by model"))

surv_reg<-ggplot(selected_genes_surv, aes(x = logFC, y = -log10(PValue),color=select)) +
  geom_point(alpha=0.5) +
 # geom_text(data = subset(selected_genes_surv,surv %in% c("Upregulated","Downregulated")),
  #          aes(label=genes),vjust=1.5,size=3)+
  scale_color_manual(values = c("steelblue", "red"), name = "selection status") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange",linewidth=1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "orange",linewidth=1) +
  theme_bw() +
  labs(x = "Log2 fold change", y = "-log10 (Pvalue)", title = "(B)")

ggsave("./3.Outputs/Images/volcano_surv_1.tiff",width = 7,height = 4)

library(ggpubr)

ggarrange(cc_reg,surv_reg,common.legend = T)
####more trimmed-----

selected_genes_surv<-selected_genes_surv %>% 
  mutate(surv2=ifelse(logFC>0.1 & PValue<0.05,"Upregulated",
                     ifelse(logFC< -0.1 & PValue<0.05,"Downregulated",
                            
                            "Not Significant")))

selected_genes_surv %>% group_by(surv2) %>% count()

surv2_reg<-ggplot(selected_genes_surv, aes(x = logFC, y = -log10(PValue),color=surv2)) +
  geom_point(alpha=0.5) +
  geom_text(data = subset(selected_genes_surv,surv2 %in% c("Upregulated","Downregulated")),
            aes(label=genes),vjust=1.5,size=3)+
  scale_color_manual(values = c("steelblue", "black","red"), name = "Significant") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray40") +
  theme_bw() +
  labs(x = "Log fold change", y = "-log10 (Pvalue)", title = "Survival")

ggsave("./3.Outputs/Images/volcano_surv_2.tiff",width = 7,height = 4)

ggarrange(cc_reg,surv2_reg,common.legend = T)

library(EnhancedVolcano)




EnhancedVolcano(selected_genes_cc,
                lab = selected_genes_cc$genes,
                x = 'logFC',
                y = 'PValue',
                title = 'case control',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0)

##survival------------

selected_genes_surv <- read_csv("3.Outputs/selected_genes_surv.csv")


EnhancedVolcano(selected_genes_surv,
                lab = selected_genes_surv$genes,
                x = 'logFC',
                y = 'PValue',
                title = 'survival',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0)
