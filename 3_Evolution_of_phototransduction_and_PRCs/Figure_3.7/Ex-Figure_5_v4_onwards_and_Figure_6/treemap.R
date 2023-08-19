# library
library(treemap)

# Based on data updated up to 6 July 2023.
# Build Dataset
groups <- c(rep("Helix-turn-helix",8),rep("Basic Domains group",4),rep("Unclassified Structure",3),
           rep("Zinc-Coordinating Group",3),rep("Beta-Scaffold Factors",2),rep("Other Alpha-Helix Group",1))
subgroup <- paste("subgroup" , c(1,2,3,4,1,2,1,2,3), sep="-")

subgroups <- c("Homeobox","ARID","Fork_head","ETS","MYB","SRF","Pou","COE",
               "TF_bZIP","HLH","bHLH","AP-2",
               "Others","T-box","MH1",
               "zf-C2H2","RXR-like","zf-GATA",
               "CSD","STAT",
               "HMG")
values <- c(5,3,1,1,1,1,1,1,
           7,3,3,1,
           12,1,1,
           5,1,1,
           1,1,
           4)
data <- data.frame(groups,subgroups,values)

# treemap
treemap(data,
        index=c("groups","subgroups"),
        vSize="values",
        type="index",
        fontsize.labels=c(15,12),
        fontcolor.labels=c("white","black"),
        fontface.labels=c(2,2),
        bg.labels=c("transparent"),
        align.labels=list(
          c("center", "center"), 
          c("left", "top")
        ),
        overlap.labels=1,
        inflate.labels=F,
        border.col=c("black","white"),
        border.lwds=c(4,1),
        title="Transcription Factor Families",
        fontsize.title=18
) 
