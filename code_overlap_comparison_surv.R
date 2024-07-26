library(readxl)
library(tidyverse)
library(ggpubr)
library(ggVennDiagram)
library(clipr)

total_surv<-read.csv("./3.Outputs/total_surv.csv")
pahtway_comparison_all_surv <-  read_excel("3.Outputs/Tables/Consensus_all_genes_surv.xlsx")
pahtway_comparison_total_surv <- read_excel("3.Outputs/Tables/Consensus_total_genes_surv.xlsx")

##pathways list----------
pahtway_comparison_all_surv %>% group_by(source) %>% count() %>% 
  ggplot(aes(x=n,y=reorder(`source`,n)))+
  geom_bar(stat="identity",fill="steelblue")+
  geom_text(aes(label=n),hjust=1)+labs(y="source")

pahtway_comparison_total_surv$target_n<-  as.numeric(substr(pahtway_comparison_total_surv$overlap, 1, 2))

pahtway_comparison_all_surv_with_total<-pahtway_comparison_all_surv %>% left_join(pahtway_comparison_total_surv[,c(3,10,11)],by="pathway")


#creating matching column
pahtway_comparison_all_surv_total_matches <- pahtway_comparison_all_surv_with_total %>%
  mutate(unlist_words = strsplit(members_input_overlap, "; ")) %>%
  mutate(MatchedWords = Map(intersect, unlist_words, list(total_surv$genes))) %>%
  mutate(MatchedCount = lengths(Map(intersect, unlist_words, list(total_surv$genes))))

pahtway_comparison_all_surv_total_matches$percent_total<- 
  round(pahtway_comparison_all_surv_total_matches$MatchedCount/pahtway_comparison_all_surv_total_matches$effective_size*100,2)

pahtway_comparison_all_surv_total_matches<-pahtway_comparison_all_surv_total_matches %>% 
  mutate(overlap_new=paste0(MatchedCount,"/",effective_size,"(",percent_total,"%)"))

pahtway_comparison_all_surv_total_matches$MatchedWords <- 
  sapply(pahtway_comparison_all_surv_total_matches$MatchedWords, function(x) paste(x, collapse = ";"))

library(writexl)
write_xlsx(pahtway_comparison_all_surv_total_matches,
          "3.Outputs/Tables/pahtway_comparison_all_surv_total_matches.xlsx")

##top 10 plot-----------

(top10_all_cc_kegg<-pahtway_comparison_all_cc_total_matches %>% filter(source=="KEGG") %>% slice(1:10) %>% 
   ggplot(aes(x= -log10(`p-value`),y=reorder(`pathway`,-log10(`p-value`))))+
   geom_bar(stat = "identity",fill="steelblue")+
   geom_text(aes(label=`overlap_new`),hjust=1.25)+
   labs(x="-log10(PValue)",y="Pathway")+
   ggtitle("all genes cc (KEGG)")+
   theme_minimal())

(top10_all_cc_reactom<-pahtway_comparison_all_cc_total_matches %>% filter(source=="Reactome") %>% slice(1:10) %>% 
    ggplot(aes(x= -log10(`p-value`),y=reorder(`pathway`,-log10(`p-value`))))+
    geom_bar(stat = "identity",fill="steelblue")+
    geom_text(aes(label=`overlap_new`),hjust=1.25)+
    labs(x="-log10(PValue)",y="Pathway")+
    ggtitle("all genes cc (Reactom)")+
    theme_minimal())


(top10_all_cc_Wikipathways<-pahtway_comparison_all_cc_total_matches %>% filter(source=="Wikipathways") %>% slice(1:10) %>% 
    ggplot(aes(x= -log10(`p-value`),y=reorder(`pathway`,-log10(`p-value`))))+
    geom_bar(stat = "identity",fill="steelblue")+
    geom_text(aes(label=`overlap_new`),hjust=1.25)+
    labs(x="-log10(PValue)",y="Pathway")+
    ggtitle("all genes cc (Wikipathways)")+
    theme_minimal())

ggarrange(top10_all_cc_kegg,top10_all_cc_reactom,top10_all_cc_Wikipathways,nrow = 3)
ggsave("./3.Outputs/Images/all_pahtway_top_10.tiff")



##top 15 plot-----------

(top15_all_surv_kegg<-pahtway_comparison_all_surv_total_matches %>% filter(source=="KEGG") %>% 
  dplyr:: slice(1:15) %>% 
   ggplot(aes(x= -log10(`p-value`),y=reorder(`pathway`,-log10(`p-value`))))+
   geom_bar(stat = "identity",fill="steelblue")+
   geom_text(aes(label=`overlap_new`),hjust=1.25,size=4)+
   labs(x="-log10(PValue)",y="Pathway")+
   ggtitle("all genes surv (KEGG)")+
   theme_minimal())

ggsave("./3.Outputs/Images/all_pahtway_top_15_KEGG_surv.tiff",height = 7,width = 8)

(top15_all_surv_reactom<-pahtway_comparison_all_surv_total_matches %>% filter(source=="Reactome") %>% 
    dplyr::slice(1:15) %>% 
    ggplot(aes(x= -log10(`p-value`),y=reorder(`pathway`,-log10(`p-value`))))+
    geom_bar(stat = "identity",fill="steelblue")+
    geom_text(aes(label=`overlap_new`),hjust=1.25,size=4)+
    labs(x="-log10(PValue)",y="Pathway")+
    ggtitle("all genes surv (Reactom)")+
    theme_minimal())

ggsave("./3.Outputs/Images/all_pahtway_top_15_reactom_surv.tiff",height = 7,width = 8)


(top15_all_surv_Wikipathways<-pahtway_comparison_all_surv_total_matches %>% filter(source=="Wikipathways") %>% slice(1:15) %>% 
    ggplot(aes(x= -log10(`p-value`),y=reorder(`pathway`,-log10(`p-value`))))+
    geom_bar(stat = "identity",fill="steelblue")+
    geom_text(aes(label=`overlap_new`),hjust=1.25,size=4)+
    labs(x="-log10(PValue)",y="Pathway")+
    ggtitle("all genes surv (Wikipathways)")+
    theme_minimal())

ggarrange(top15_all_surv_kegg,top15_all_surv_reactom,nrow = 2)
ggsave("./3.Outputs/Images/all_pahtway_top_15_surv.tiff",height = 9,width = 8)

###code_total_enriched------------
(top10_total<-pahtway_comparison_total %>% 
   ggplot(aes(x= log_pvalue,y=reorder(`pathway/source`,log_pvalue)))+
   geom_bar(stat = "identity",fill="steelblue")+
   geom_text(aes(label=`total(covered)`),hjust=1)+
   labs(x="-log10(PValue)",y="Pathway/Source")+
   theme_minimal())