#Load required libraries
library(tidyverse)
library(pals)

########################################################################
## Figure 2: Number of isolates per ST identified for each species. 
########################################################################
mlst_meta <- read_csv("data/metadata_mlst.csv")

#ST counts per species
species_cts<-mlst_meta %>% group_by(species,ST) %>% count(ST)
str(species_cts)
spec2 <-  species_cts %>%
  mutate_at("species", str_replace, "1", "Human") %>%
  mutate_at("species", str_replace, "2", "Bovine")
#keep order for species
spec2$species <- factor(spec2$species, levels = c("Bovine", "Human"))

p <-ggplot(spec2, 
       aes(fill=fct_reorder(as.factor(ST),n), y=n, x=species)) + 
  geom_bar(position="stack", stat="identity")+
  #geom_text(size = 3, position = position_stack(vjust = 0.5))+
  coord_flip()+
  labs (x="Species",
        y="N isolates", 
        fill="STs")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=as.vector(polychrome(26)))
ggsave(filename="figures/Figure2.pdf", plot=p, width=12,height=6)

########################################################################
## Figure 3: Number of isolates per ST observed for each sampling location and sampling round. 
########################################################################
#ST counts per farm
farm_counts<-mlst_meta %>% group_by(farm) %>% count(ST)
ST_counts<-mlst_meta %>% group_by(ST,farm,round) %>% count(ST)

# Stacked
farm_counts$farm<-addNA(farm_counts$farm)
farm_counts$farm<-as.factor(farm_counts$farm)
farm_counts$ST<-as.factor(farm_counts$ST)
farm_counts %>% mutate(farm_counts = fct_reorder(ST,n)) %>% pull(ST)

#facet grid for every farm with rounds
levels(ST_counts$farm)<-c("Farm 1", "Farm 2","Farm 3", "Other")
ST_counts$farm<-addNA(ST_counts$farm)
ST_counts$farm<-as.factor(ST_counts$farm)
ST_counts$ST<-as.factor(ST_counts$ST)
ST_counts$round<-as.factor(ST_counts$round)

cols <- c("#009E73","#D55E00", "#56B4E9")

ST_countsW <- ST_counts %>%
  mutate_at("farm", str_replace, "1", "Farm 1") %>%
  mutate_at("farm", str_replace, "2", "Farm 2") %>%
  mutate_at("farm", str_replace, "3", "Farm 3")
ST_countsW$farm[is.na(ST_countsW$farm)] <- "School"

p2<-ggplot(ST_countsW,
          aes(fill=fct_reorder(round,n), y=n, x=ST)) +
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  labs (x="Sequence Type",
        y="Number of Isolates",
        fill="Round")+
  scale_fill_manual(values=cols)+
  facet_wrap(vars(farm), ncol = 4)
ggsave(filename="figures/Figure3.pdf", plot=p2, width=12,height=6)