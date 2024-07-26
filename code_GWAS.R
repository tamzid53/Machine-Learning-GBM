library(readxl)
library(tidyverse)
gwas_image <- read_excel("gwas_image.xlsx")

gwas_image %>% 
  ggplot(aes(x=`Reported Trait`,y=-log10(`p-value`),color=variant))+
  geom_point(size=3)+
  ylim(0,60)+
  coord_flip()+
  geom_text(aes(label=variant),hjust=-0.05)+
  facet_wrap(~gene,scales = "free_y",ncol = 2)+
  theme_light()+
  theme(legend.position = "",
        axis.title = element_text(color='black',size=12),
        axis.text = element_text(color='black',size=12),
        strip.text = element_text(color = "black",face ='bold',size=13), # Title color
        strip.background = element_rect(fill = "white", color = "black") )

        