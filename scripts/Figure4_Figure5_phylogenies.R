
library("ggtree")
library("treeio")
library("tidyverse")
library("ggnewscale")
library("pals")

######################################################
########### Load trees ##############
######################################################
treefile<-"input/gubbins_Staph_all_clean_variant.fasta.contree"   
treefile2<-"input/gubbins_ST1_clean_variant.fasta.contree"
all_tree<-read.iqtree(treefile)
ST1_tree<-read.iqtree(treefile2)

######################################################
########### add traits ##############
######################################################

#convert tree to tibble
tree_all_tib <- as_tibble(all_tree)
tree_ST1_tib <- as_tibble(ST1_tree)
#load metadata
traits <- read_csv("results/metadata_mlst.csv")
#rename isolate to label for join
names(traits) <- c("label", "age", "farm", "round", "school", "species", "milk","ST")
#merge tree with traits
cols <- colnames(traits)
traits %<>%  mutate_at(cols, factor)
str(traits)
levels(traits$farm)<- c("Farm 1", "Farm 2", "Farm 3")
levels(traits$round)<-c("Round 1","Round 2", "Round 3")


y <- full_join(tree_all_tib, traits, by = 'label')
y2 <- left_join(tree_ST1_tib, traits, by = 'label')
#convert back to tree object
all_tree2<-as.treedata(y)
ST1_tree2<-as.treedata(y2)

############################################################
################### Figure 4 - Phylogeny ###################
############################################################

#Species as colored tip nodes, farm as colored tip labels, sampling round as text, 
cols<-c("1"='black', "2"='black')
p<-ggtree(all_tree2)+    
  geom_tiplab(aes(fill=farm),
              color="black",
              geom="label",
              label.padding = unit(0.15, "lines"), 
              label.size = 0.1,
              align=T,
              linesize=0.1,
              size=2,
              show.legend=FALSE) + 
  scale_fill_brewer("Blues")+

  new_scale_fill()+
  geom_tippoint(aes(color=species, shape=species), size=2)+
  scale_color_manual(values=cols)+
  scale_shape_manual(values=c(16, 5))+
  geom_tiplab(aes(label=round), align=T, linetype=NA, size=2, offset=0.02)+
  theme(legend.position="right")+
  geom_treescale(x=0, y=80)

p

### add heatmap with STs 
ST<-data.frame("ST" = traits[,c("ST")])
rownames(ST)<-traits$label
#use same colors as for ST graph - struggling to get this working
ST_list<-levels(ST$ST)
str(ST_list)

st_cols<-polychrome(26)

h1<-gheatmap(p, ST, 
             offset = 0.03, 
             color=NULL, 
             colnames=FALSE,
             width=0.015)+
  scale_fill_manual(breaks=c("1",    "5",    "6" ,   "7"   , "30"  , "34" ,  "39",   "45"  , "72",   "151",  "188"  ,"398" , "425" , "582"  ,"705" , "1247" ,"1598", "1649", "2000" ,"2001" ,"2002", "2003", "2004" ,"2005" ,"2006" ,"2007"), values=c( "#782AB6","#1CBE4F","#5A5156","#AA0DFE","#F8A19F","#FA0087",  "#E4E1E3","#1C8356", "#F6222E", "#F7E1A0" ,"#FC1CBF" , "#FE00FA","#16FF32","#3283FE","#C075A6"  ,"#FEAF16","#85660D"   ,"#B00068","#B10DA1"   ,  "#1CFFCE",  "#90AD1C" ,"#325A9B", "#C4451C"  , "#2ED9FF", "#FBE426" , "#DEA0FD" ), name="ST")
h1

ggsave(filename="figures/Staph_all_n.pdf", plot=h1, width=8,height=13)

#26 colours needed

############################################################
################## Figure 5 - Phylogeny ST1   ##############
############################################################
cols<-c("1"='black', "2"='black')

#drop 2 isolates
to_drop<-c("20MR0155", "20MR0134")
ST1_tree2_red <- drop.tip(ST1_tree2, to_drop)

p2<-ggtree(ST1_tree2_red)+    
  geom_tiplab(aes(fill=farm),
               color="black",
               geom="label",
               label.padding = unit(0.15, "lines"), # amount of padding around the labels
               label.size = 0.1,
               align=TRUE,
               linesize=0.1,
               size=2,
               show.legend=FALSE) + # size of label border
  scale_fill_brewer("Blues")+
  new_scale_fill()+
  geom_tippoint(aes(color=species, shape=species), size=2)+
  scale_color_manual(values=cols)+
  scale_shape_manual(values=c(16, 5))+
  geom_tiplab(aes(label=round), align=T, linetype=NA, size=2, offset=0.001)+
  theme(legend.position="right")+
  geom_treescale()
p2

#ggsave(filename="figures/Staph_ST1.pdf", plot=p2, width=8,height=6)
