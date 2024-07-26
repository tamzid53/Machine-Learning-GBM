##GABRA6
library(readr)
library(survival)
library(ggsurvfit)
library("survminer")
library(tidyverse)
library(data.table)
library(WriteXLS)
# Survival loading Data codes and prep------------
TCGA_survival_data<-fread("TCGA_survival_data") #Surv data
ExprSubsetBySamp <- read_csv("ExprSubsetBySamp.csv")

pred_data <- read_csv("3.Outputs/Tables/pred_table_surv.csv")

Surv_TCGA_use = subset(TCGA_survival_data, sample %in% colnames(ExprSubsetBySamp))


##One yaer---------
Surv_TCGA_use1<- Surv_TCGA_use %>% 
  mutate(one_year=ifelse(OS.time>=365,1,
                         ifelse(OS.time<365 & OS==1,0,99))) %>% 
  filter(one_year!=99) %>% dplyr::select(sample,one_year)


Surv_TCGA_use_kp<- Surv_TCGA_use1 %>% inner_join(Surv_TCGA_use[,1:3],by="sample") %>% 
                      inner_join(pred_data[,-1],by="sample")

#write.csv(Surv_TCGA_use_kp,"3.Outputs/Tables/Surv_TCGA_use_kp.csv")

##test model 1----
(surv1<-Surv(Surv_TCGA_use$OS.time,Surv_TCGA_use$OS))
(km_fit <- survfit(surv1 ~ 1))

survfit2(surv1 ~ 1, data = Surv_TCGA_use) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  )



##test model 12----
(surv2<-Surv(Surv_TCGA_use_kp$OS.time,Surv_TCGA_use_kp$one_year))
(km_fit2 <- survfit(surv1 ~ 1))

survfit2(surv2 ~ 1, data = Surv_TCGA_use_kp) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  )



##test model 12----

km_fit1<-survfit(Surv(OS.time,one_year) ~ 1,data=Surv_TCGA_use_kp)
km_fit2<-survfit(Surv(OS.time,pred_surv) ~ 1,data=Surv_TCGA_use_kp)
fit <- list(Original = km_fit1, Predicted = km_fit2)

ggsurvplot(fit, data = Surv_TCGA_use_kp, combine = TRUE, # Combine curves
          risk.table = F,                  # Add risk table
          conf.int = F,                    # Add confidence interval
          conf.int.style = "step",            # CI style, use "step" or "ribbon"
          censor = FALSE,                     # Remove censor points
          tables.theme = theme_cleantable(),  # Clean risk table
          palette = "jco")


##blend with genes------
all_genes_cc <- read_csv("3.Outputs/Tables/Xnew_cc.csv")


Surv_TCGA_use_kp_all_gene<-Surv_TCGA_use_kp %>% left_join(all_genes_cc,by=("sample"))



# Initialize a list to store p-values for each column
p_values_list <- list()

# Loop through each column and calculate p-value
for (col_name in colnames(Surv_TCGA_use_kp_all_gene)[6:4667]) {
  # Create two groups based on the median of the column
  median_value <- round(median(Surv_TCGA_use_kp_all_gene[[col_name]]),5)
  expression <- ifelse(Surv_TCGA_use_kp_all_gene[[col_name]] < median_value, "Low", "High")
  
  # Create a survival object
  surv_obj <- with(Surv_TCGA_use_kp_all_gene, Surv(OS.time, OS))
  
  # Perform log-rank test
  log_rank_test <- survdiff(Surv(OS.time,OS) ~ expression,data=Surv_TCGA_use_kp_all_gene)
  
  # Get the p-value from the test result
  p_value <- log_rank_test$pvalue
  
  # Store the p-value in the list
  p_values_list[[col_name]] <- p_value
}



# Convert the list of p-values to a dataframe
p_values_df <- data.frame(genes = colnames(Surv_TCGA_use_kp_all_gene)[6:4667], 
                          P_Value = unlist(p_values_list))
p_values_df$sig<-ifelse(p_values_df$P_Value<0.05,"significant","not significant")
# Print the resulting dataframe
print(p_values_df)
#write.csv(p_values_df,"p_values_log_rank.csv",row.names = F)


p_values_df_sig<- p_values_df %>% filter(sig=="significant") %>% dplyr:: select(genes)
#write.csv(p_values_df_sig,"3.Outputs/p_values_df_sig.csv",row.names = F)

#significant gene in case control model



#creating matching column with signficant genes
pahtway_comp_surv_mathced <- pahtway_comparison_all_cc_total_matches %>%
  #mutate(unlist_words_surv = strsplit(MatchedWords, ", ")) %>%
  mutate(MatchedWords_surv = lapply(MatchedWords, function(lst) intersect(lst, p_values_df_sig$genes))) %>% 
  mutate(MatchedCount_surv = lengths(lapply(MatchedWords, function(lst) intersect(lst, p_values_df_sig$genes))))

  pahtway_comp_surv_mathced_save<- pahtway_comp_surv_mathced 
  
  pahtway_comp_surv_mathced_save$MatchedWords <- 
    sapply(pahtway_comp_surv_mathced_save$MatchedWords, function(x) paste(x, collapse = ";"))
  
  pahtway_comp_surv_mathced_save$MatchedWords_surv <- 
    sapply(pahtway_comp_surv_mathced_save$MatchedWords_surv, function(x) paste(x, collapse = ";"))

 library(writexl) 
  write_xlsx (pahtway_comp_surv_mathced_save,"./3.Outputs/Tables/pahtway_comp_surv_mathced_save.xlsx")  

   write_clip(pahtway_comp_surv_mathced)
# fitting genes--------------------------------


summary(km_fit_GABRA6)$surv

survdiff(Surv(OS.time,OS) ~ expression,data=Surv_TCGA_use_kp_gene)$pvalue


###GABRA6-----

Surv_TCGA_use_kp_GABRA6<-Surv_TCGA_use_kp_all_gene %>% 
  select(sample,one_year, OS.time,OS,GABRA6)
                                                              )

med<-median(Surv_TCGA_use_kp_GABRA6$GABRA6,na.rm = T)

Surv_TCGA_use_kp_GABRA6$expression<-ifelse(Surv_TCGA_use_kp_GABRA6$GABRA6 < med,
                                        "low","high")

km_fit_GABRA6<-survfit(Surv(OS.time,OS) ~ expression,data=Surv_TCGA_use_kp_gene)

# Calculate the log-rank p-value
logrank_p_value <- surv_pvalue(km_fit_GABRA6, data = Surv_TCGA_use_kp_gene, method = "log-rank")$pval




(GABRA6<-ggsurvplot(km_fit_GABRA6, data = Surv_TCGA_use_kp_gene,# Combine curves
              pval=TRUE, pval.method = T,
           risk.table = F,                  # Add risk table
           conf.int = T,                    # Add confidence interval
           #conf.int.style = "step",            # CI style, use "step" or "ribbon"
           censor = TRUE,                     # Remove censor points
           tables.theme = theme_cleantable(),  # Clean risk table
           palette = "jco")+ggtitle("GABRA6"))


(GABRA6<-ggsave("./3.Outputs/Images/gene_GABRA6.tiff"))

###CAMK2A----

CAMK2A<- as.tibble(all_genes_cc) %>% dplyr:: select("sample","CAMK2A")

Surv_TCGA_use_kp_gene<-Surv_TCGA_use_kp %>% inner_join(CAMK2A,by="sample")

Surv_TCGA_use_kp_gene$expression<-ifelse(Surv_TCGA_use_kp_gene$CAMK2A< -1.162102,
                                        "low","high")

km_fit_CAMK2A<-survfit(Surv(OS.time,OS) ~ expression,data=Surv_TCGA_use_kp_gene)

# Calculate the log-rank p-value
logrank_p_value <- surv_pvalue(km_fit_CAMK2A, data = Surv_TCGA_use_kp_gene, method = "log-rank")$pval



(CAMK2A<-ggsurvplot(km_fit_CAMK2A, data = Surv_TCGA_use_kp_gene,# Combine curves
           pval=TRUE, pval.method = T,
           risk.table = F,                  # Add risk table
           conf.int = T,                    # Add confidence interval
           #conf.int.style = "step",            # CI style, use "step" or "ribbon"
           censor = TRUE,                     # Remove censor points
           tables.theme = theme_cleantable(),  # Clean risk table
           palette = "jco")+ggtitle("CAMK2A"))


ggsave("./3.Outputs/Images/gene_CAMK2A.tiff")

####POGK------------

POGK<- as.tibble(all_genes_cc) %>% dplyr:: select("sample","POGK")

Surv_TCGA_use_kp_gene<-Surv_TCGA_use_kp %>% inner_join(POGK,by="sample")
med<-median(Surv_TCGA_use_kp_gene$POGK)
Surv_TCGA_use_kp_gene$expression<-ifelse(Surv_TCGA_use_kp_gene$POGK< med,
                                        "low","high")

km_fit_POGK<-survfit(Surv(OS.time,OS) ~ expression,data=Surv_TCGA_use_kp_gene)

(POGK<-ggsurvplot(km_fit_POGK, data = Surv_TCGA_use_kp_gene,# Combine curves
           risk.table = F,                  # Add risk table
           conf.int = T,                    # Add confidence interval
           #conf.int.style = "step",            # CI style, use "step" or "ribbon"
           censor = FALSE,                     # Remove censor points
           tables.theme = theme_cleantable(),  # Clean risk table
           palette = "jco") +ggtitle('POGK'))

ggsave("./3.Outputs/Images/gene_POGK.tiff")

#######GLS2-----------

GLS2_use<- as.tibble(all_genes_cc) %>% dplyr:: select("sample","GLS2")

Surv_TCGA_use_kp_gene<-Surv_TCGA_use_kp %>% inner_join(GLS2_use,by="sample")
med<-median(Surv_TCGA_use_kp_gene$GLS2)
Surv_TCGA_use_kp_gene$expression<-ifelse(Surv_TCGA_use_kp_gene$GLS2< med,
                                        "low","high")

km_fit_GLS2<-survfit(Surv(OS.time,OS) ~ expression,data=Surv_TCGA_use_kp_gene)

(GLS2<-ggsurvplot(km_fit_GLS2, data = Surv_TCGA_use_kp_gene,# Combine curves
           pval=TRUE, pval.method = T,
           risk.table = F,                  # Add risk table
           conf.int = T,                    # Add confidence interval
           #conf.int.style = "step",            # CI style, use "step" or "ribbon"
           censor = FALSE,                     # Remove censor points
           tables.theme = theme_cleantable(),  # Clean risk table
           palette = c("red","green4")) + ggtitle("GLS2"))

ggsave("./3.Outputs/Images/gene_GLS2.tiff")


###Using prediction model--------
splots <- list()
km_fit_pred_mod<-survfit(Surv(OS.time,OS) ~ pred_surv,data=Surv_TCGA_use_kp)

(splots[[1]]<-ggsurvplot(km_fit_pred_mod, data = Surv_TCGA_use_kp,# Combine curves
           pval=TRUE, pval.method = T,
           risk.table = F,                  # Add risk table
           conf.int = T,                    # Add confidence interval
           #conf.int.style = "step",            # CI style, use "step" or "ribbon"
           censor = FALSE,                     # Remove censor points
           tables.theme = theme_cleantable(),  # Clean risk table
           palette = c("red","green4")) + ggtitle("MLDA Model"))

ggsave("./3.Outputs/Images/pred_model.tiff")


library(ggpubr)

ggarrange(GABRA6,CAMK2A,GLS2,pred, nrow = 2,ncol=2)

par(mfrow=c(2,2))


##LAIR1-------

Surv_TCGA_use_kp_sg<-Surv_TCGA_use_kp_all_gene %>% 
  select(sample,one_year, OS.time,OS,LAIR1)


med<-median(Surv_TCGA_use_kp_sg$LAIR1,na.rm = T)

med2<-quantile(Surv_TCGA_use_kp_sg$LAIR1,na.rm = T,probs = 0.45)


Surv_TCGA_use_kp_sg$expression<-ifelse(Surv_TCGA_use_kp_sg$LAIR1 < med2,
                                          "low","high")
km_fit<-survfit(Surv(OS.time,OS) ~ expression,data=Surv_TCGA_use_kp_sg)

(splots[[2]]<-ggsurvplot(km_fit, data = Surv_TCGA_use_kp_sg,# Combine curves
                  pval=TRUE, pval.method = T,
                  risk.table = F,                  # Add risk table
                  conf.int = T,                    # Add confidence interval
                  #conf.int.style = "step",            # CI style, use "step" or "ribbon"
                  censor = FALSE,                     # Remove censor points
                  tables.theme = theme_cleantable(),  # Clean risk table
                  palette = c("red","green4")) + ggtitle("LAIR1"))

##CD68------

Surv_TCGA_use_kp_sg<-Surv_TCGA_use_kp_all_gene %>% 
  select(sample,one_year, OS.time,OS,CD68)


med<-median(Surv_TCGA_use_kp_sg$CD68,na.rm = T)
med2<-quantile(Surv_TCGA_use_kp_sg$CD68,na.rm = T,probs = 0.45)


Surv_TCGA_use_kp_sg$expression<-ifelse(Surv_TCGA_use_kp_sg$CD68 < med2,
                                      "low","high")
km_fit<-survfit(Surv(OS.time,OS) ~ expression,data=Surv_TCGA_use_kp_sg)

(splots[[3]]<-ggsurvplot(km_fit, data = Surv_TCGA_use_kp_sg,# Combine curves
                   pval=TRUE, pval.method = T,
                   risk.table = F,                  # Add risk table
                   conf.int = T,                    # Add confidence interval
                   #conf.int.style = "step",            # CI style, use "step" or "ribbon"
                   censor = FALSE,                     # Remove censor points
                   tables.theme = theme_cleantable(),  # Clean risk table
                   palette = c("red","green4")) + ggtitle("CD68"))

##MOBP------

Surv_TCGA_use_kp_sg<-Surv_TCGA_use_kp_all_gene %>% 
  select(sample,one_year, OS.time,OS,MOBP)


med<-median(Surv_TCGA_use_kp_sg$MOBP,na.rm = T)

Surv_TCGA_use_kp_sg$expression<-ifelse(Surv_TCGA_use_kp_sg$MOBP < med,
                                      "low","high")
km_fit<-survfit(Surv(OS.time,OS) ~ expression,data=Surv_TCGA_use_kp_sg)

(splots[[4]]<-ggsurvplot(km_fit, data = Surv_TCGA_use_kp_sg,# Combine curves
                  pval=TRUE, pval.method = T,
                  risk.table = F,                  # Add risk table
                  conf.int = T,                    # Add confidence interval
                  #conf.int.style = "step",            # CI style, use "step" or "ribbon"
                  censor = FALSE,                     # Remove censor points
                  tables.theme = theme_cleantable(),  # Clean risk table
                  palette = c("red","green4")   ) + 
                    ggtitle("MOBP"))


###Using prediction model--------

splots <- list()
km_fit_pred_mod<-survfit(Surv(OS.time,OS) ~ pred_surv,data=Surv_TCGA_use_kp)

(splots[[4]]<-ggsurvplot(km_fit_pred_mod, data = Surv_TCGA_use_kp,# Combine curves
                         pval=TRUE, pval.method = T,
                         risk.table = F,                  # Add risk table
                         conf.int = T,                    # Add confidence interval
                         #conf.int.style = "step",            # CI style, use "step" or "ribbon"
                         censor = FALSE,                     # Remove censor points
                         tables.theme = theme_cleantable(),  # Clean risk table
                         palette = c("blue4","orange2"),
                         legend.title = "Survival Status",
                         legend.labs =
                           c("<1-year survival", "> 1-year survival"))+
                          ggtitle("netLDA Model"))

##MAG------

Surv_TCGA_use_kp_sg<-Surv_TCGA_use_kp_all_gene %>% 
  select(sample,one_year, OS.time,OS,MAG)


med<-median(Surv_TCGA_use_kp_sg$MAG,na.rm = T)

Surv_TCGA_use_kp_sg$expression<-ifelse(Surv_TCGA_use_kp_sg$MAG < med,
                                       "low","high")
km_fit<-survfit(Surv(OS.time,OS) ~ expression,data=Surv_TCGA_use_kp_sg)

(splots[[1]]<-ggsurvplot(km_fit, data = Surv_TCGA_use_kp_sg,# Combine curves
                         pval=TRUE, pval.method = T,
                         risk.table = F,                  # Add risk table
                         conf.int = T,                    # Add confidence interval
                         #conf.int.style = "step",            # CI style, use "step" or "ribbon"
                         censor = FALSE,                     # Remove censor points
                         tables.theme = theme_cleantable(),  # Clean risk table
                         palette = c("blue4","orange2")   ) + 
    ggtitle("MAG"))

##CNDP1-----------
Surv_TCGA_use_kp_sg<-Surv_TCGA_use_kp_all_gene %>% 
  select(sample,one_year, OS.time,OS,CNDP1)

med<-median(Surv_TCGA_use_kp_sg$CNDP1,na.rm = T)

Surv_TCGA_use_kp_sg$expression<-ifelse(Surv_TCGA_use_kp_sg$CNDP1 < med,
                                       "low","high")
km_fit<-survfit(Surv(OS.time,OS) ~ expression,data=Surv_TCGA_use_kp_sg)

(splots[[2]]<-ggsurvplot(km_fit, data = Surv_TCGA_use_kp_sg,# Combine curves
                         pval=TRUE, pval.method = T,
                         risk.table = F,                  # Add risk table
                         conf.int = T,                    # Add confidence interval
                         #conf.int.style = "step",            # CI style, use "step" or "ribbon"
                         censor = FALSE,                     # Remove censor points
                         tables.theme = theme_cleantable(),  # Clean risk table
                         palette = c("blue4","orange2")   ) + 
    ggtitle("CNDP1"))



##PTPRN-----------
Surv_TCGA_use_kp_sg<-Surv_TCGA_use_kp_all_gene %>% 
  select(sample,one_year, OS.time,OS,ANGPTL4)

med<-median(Surv_TCGA_use_kp_sg$ANGPTL4,na.rm = T)

Surv_TCGA_use_kp_sg$expression<-ifelse(Surv_TCGA_use_kp_sg$ANGPTL4 < med,
                                       "low","high")
km_fit<-survfit(Surv(OS.time,OS) ~ expression,data=Surv_TCGA_use_kp_sg)

(splots[[3]]<-ggsurvplot(km_fit, data = Surv_TCGA_use_kp_sg,# Combine curves
                         pval=TRUE, pval.method = T,
                         risk.table = F,                  # Add risk table
                         conf.int = T,                    # Add confidence interval
                         #conf.int.style = "step",            # CI style, use "step" or "ribbon"
                         censor = FALSE,                     # Remove censor points
                         tables.theme = theme_cleantable(),  # Clean risk table
                         palette = c("blue4","orange2")   ) + 
    ggtitle("ANGPTL4"))

library(cowplot)

arrange_ggsurvplots(splots, print = TRUE, nrow = 2,ncol=2)

