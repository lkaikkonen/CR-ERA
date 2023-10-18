### Bayesian network model for estimating the impacs of seabed mining on the Chatham Rise, SW Pacific Ocean
### Model script (lines 5-646) and visualisations (lines 645-2276)
### Kaikkonen L. 2023

#--------------------------------------------------------------------------#

# PACKAGES
# load packages from bioconductor

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("gRain", "gRbase","Rgraphviz"))

# Install and load R packages 
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages<-c("bnlearn", "bnviewer", "Rgraphviz", "DiagrammeR", "ggplot2", "gRain", "gRbase",
            "rhandsontable",
            "plotrix",
            "grid",
            "igraph",
            "shape",
            "RColorBrewer",
            "lattice",
            "latticeExtra",
            "bnmonitor", "reshape2", "ggpubr")
ipak(packages)

#--------------------------------------------------------------------------#

# DEFINING THE BN STRUCTURE 

# Read network nodes and connection
nodes_t<-read.table("Model structure/Node_CR-ERA.csv", header=FALSE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
nodes2<-as.matrix(nodes_t[,1])
states<-nodes_t[,2:6]
colnames(states)<-NULL

# Import network connections as matrix

data2<- read.table("Model structure/Edge_CR-ERA.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE) #read edges
arc.set =as.matrix(data2)
dagCR<-empty.graph(nodes2) # create empty DAG with nodes in the data
arcs(dagCR) = arc.set # create DAG

plot(dagCR) #check DAG structure, simple

graphviz.plot(dagCR, layout = "dot")

g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(dagCR))
graph::nodeRenderInfo(g) <- list(fontsize=140)
Rgraphviz::renderGraph(g)

# Set levels/possible outcome states for discrete variables

Mined <-states[1,1:2]
NRem<-states[2,1:2] #nodule removal
ST <-states[3,1:3] #sediment type
SC<-states[4,1:3] #sediment contaminants
DE<-states[5,1:3] #depth of extraction
PR<-states[6,1:2] #Plume release
DistMS<-states[7,1:3] # distance from mining site
MI<-states[8,1:3] # Mining intensity

VE<-states[9,1:3] #volume extraction
SS<-states[10,1:3] #suspended sediment
CRel<-states[11,1:2] #contaminant release
SDep<-states[12,1:3] # sediment deposition
SChan<-states[13,1:2] #sediment changes

SSI_dir<-states[14,1:5]
SSI_indir<-states[14,1:5]
SSI_t2<-states[14,1:5]
SSI_t3<-states[14,1:5]

SMI_dir<-states[14,1:5]
SMI_indir<-states[14,1:5]
SMI_t2<-states[14,1:5]
SMI_t3<-states[14,1:5]

DM_dir<-states[14,1:5]
DM_indir<-states[14,1:5]
DM_t2<-states[14,1:5]
DM_t3<-states[14,1:5]

SM_dir<-states[14,1:5]
SM_indir<-states[14,1:5]
SM_t2<-states[14,1:5]
SM_t3<-states[14,1:5]

SSBM_dir<-states[14,1:5]
SSBM_indir<-states[14,1:5]
SSBM_t2<-states[14,1:5]
SSBM_t3<-states[14,1:5]

LMM_dir<-states[14,1:5]
LMM_indir<-states[14,1:5]
LMM_t2<-states[14,1:5]
LMM_t3<-states[14,1:5]

LSI_dir<-states[14,1:5]
LSI_indir<-states[14,1:5]
LSI_t2<-states[14,1:5]
LSI_t3<-states[14,1:5]

MPE_dir<-states[14,1:5]
MPE_indir<-states[14,1:5]
MPE_t2<-states[14,1:5]
MPE_t3<-states[14,1:5]

GH_dir<-states[14,1:5]
GH_indir<-states[14,1:5]
GH_t2<-states[14,1:5]
GH_t3<-states[14,1:5]

MGE_dir<-states[14,1:5]
MGE_indir<-states[14,1:5]
MGE_t2<-states[14,1:5]
MGE_t3<-states[14,1:5]

PH_dir<-states[14,1:5]
PH_indir<-states[14,1:5]
PH_t2<-states[14,1:5]
PH_t3<-states[14,1:5]

SESF_dir<-states[14,1:5]
SESF_indir<-states[14,1:5]
SESF_t2<-states[14,1:5]
SESF_t3<-states[14,1:5]

SEFF_dir<-states[14,1:5]
SEFF_indir<-states[14,1:5]
SEFF_t2<-states[14,1:5]
SEFF_t3<-states[14,1:5]

SCFF_dir<-states[14,1:5]
SCFF_indir<-states[14,1:5]
SCFF_t2<-states[14,1:5]
SCFF_t3<-states[14,1:5]

SCSF_dir<-states[14,1:5]
SCSF_indir<-states[14,1:5]
SCSF_t2<-states[14,1:5]
SCSF_t3<-states[14,1:5]


#--------------------------------------------------------------------------#

# PARAMETERISATION OF BN NODES

# CPTs of independent variables
MIp<-array(c(0.4,0.4,0.2), dim=c(3), dimnames=list(Mining_intensity = MI)) #Mining intensity 
DEp<- array(c(0.10,0.45, 0.45), dim = 3,dimnames=list(Depth_extraction = DE)) # Depth of extraction
STp <- array(c(0.03, 0.97,0), dim = 3,dimnames=list(Sediment_type = ST))  #Sediment type
PRp <- array(c(0.5,0.5), dim = 2,dimnames=list(Plume_release = PR)) # Plume release
DistMS.pr<- array(c(0.5,0.25, 0.25), dim = 3,dimnames=list(Distance_mining_site = DistMS)) # Distance mining site

# Dependent physicochemical variables

#tip: argument dim corresponds to the max extent of each of the variables
# add parent nodes in reverse order than in the table

# Mined
Mp<-array(c(0,1, 1,0,1,0), dim=c(2,3), dimnames=list(Mined = Mined, Distance_mining_site=DistMS))
#Nodule removal

NRp<-array(c(1,0, 0,1), dim=c(2,2), dimnames=list(Nodule_removal= NRem, Mined=Mined)) #Nodule removal y/n

# Volume extraction
VEp <-read.csv("Final CPTS/Volume_extraction.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
VEp<-as.matrix(VEp[,3:5])
colnames(VEp) <- NULL
VE.pr<-array(t(VEp),dim=c(3,3,3),dimnames = list(Volume_extraction=VE,Depth_extraction=DE,Mining_intensity=MI)) 

# Sediment contaminants
SCp <-read.csv("Final CPTS/Sediment_contaminants.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SCp<-as.matrix(SCp[,2:4])
colnames(SCp) <- NULL
SC.pr<-array(t(SCp),dim=c(3,3),dimnames = list(Sediment_contaminants=SC,Sediment_Type=ST))   #Sediment contaminants

# Suspended sediment 
SSed<-read.csv("Final CPTS/Suspended_sediment-v2.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SSed<-as.matrix(SSed[,5:7])
colnames(SSed) <- NULL
SS.pr<-array(t(SSed),dim=c(3,3,2,3,3),dimnames = list(Suspended_sediment=SS, Volume_extraction=VE, Plume_release=PR, Sediment_type=ST, Distance_mining_site=DistMS)) 

# Sediment deposition
SDepp<-read.csv("Final CPTS/Sediment_deposition-v2.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SDepp<-as.matrix(SDepp[,5:7])
colnames(SDepp) <- NULL
SD.pr<-array(t(SDepp),dim=c(3,3,2,3,3),dimnames = list(Sediment_deposition=SDep,Suspended_sediment = SS, Plume_release=PR, Sediment_type=ST, Distance_mining_site=DistMS)) 

# Contaminant release

CRelp<-read.csv("Final CPTS/Contaminant_release.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
CRelp<-as.matrix(CRelp[,5:6])
colnames(CRelp) <- NULL
CRel.pr<-array(t(CRelp),dim=c(2,3,3,3,3),dimnames = list(Contaminant_release = CRel, Sediment_contaminants=SC,Volume_extraction = VE, Sediment_type=ST,Distance_mining_site=DistMS )) 

# Sediment changes

SChp<-read.csv("Final CPTS/Sediment_changes.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SChp<-as.matrix(SChp[,5:6])
colnames(SChp) <- NULL
SCh.pr<-array(t(SChp),dim=c(2,3,2,3,3),dimnames = list(Sediment_changes=SChan,Depth_extraction=DE, Nodule_removal=NRem, Mining_intensity=MI, Sediment_deposition=SDep)) 

#Read CPTs for faunal groups

SSI_dir.p<-read.csv("Final CPTS/SSI_dir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SSI_dir.p<-as.matrix(SSI_dir.p[,3:7])
colnames(SSI_dir.p) <- NULL
SSI_dir.pr<-array(t(SSI_dir.p),dim=c(5,3,2),dimnames = list(SSI_dir=SSI_dir, Mining_intensity=MI, Mined=Mined)) 

SSI_indir.p<-read.csv("Final CPTS/SSI_indir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SSI_indir.p<-as.matrix(SSI_indir.p[,4:8])
colnames(SSI_indir.p) <- NULL
SSI_indir.pr<-array(t(SSI_indir.p),dim=c(5,2,3,3),dimnames = list(SSI_indir=SSI_indir, Contaminant_release=CRel,Sediment_deposition=SDep, Suspended_sediment=SS)) 

SSI_t2.p<-read.csv("Final CPTS/TotalMortality.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SSI_t2.p<-as.matrix(SSI_t2.p[,3:7])
colnames(SSI_t2.p) <- NULL
SSI_t2.pr<-array(t(SSI_t2.p),dim=c(5,5,5),dimnames = list(SSI_t2=SSI_indir, SSI_dir=SSI_dir, SSI_indir=SSI_indir))

SSI_t3.p<-read.csv("Final CPTS/SSI_t3.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SSI_t3.p<-as.matrix(SSI_t3.p[,3:7])
colnames(SSI_t3.p) <- NULL
SSI_t3.pr<-array(t(SSI_t3.p),dim=c(5,2,5),dimnames = list(SSI_t3=SSI_t3, Sediment_changes=SChan, SSI_t2=SSI_t2))

# Small mobile infauna SMI

SMI_dir.p<-read.csv("Final CPTS/SSI_dir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SMI_dir.p<-as.matrix(SMI_dir.p[,3:7])
colnames(SMI_dir.p) <- NULL
SMI_dir.pr<-array(t(SMI_dir.p),dim=c(5,3,2),dimnames = list(SMI_dir=SMI_dir, Mining_intensity=MI,Mined=Mined)) 

SMI_indir.p<-read.csv("Final CPTS/SMI_indir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SMI_indir.p<-as.matrix(SMI_indir.p[,4:8])
colnames(SMI_indir.p) <- NULL
SMI_indir.pr<-array(t(SMI_indir.p),dim=c(5,2,3,3),dimnames = list(SMI_indir=SMI_indir, Contaminant_release=CRel,Sediment_deposition=SDep, Suspended_sediment=SS)) 

SMI_t2.p<-read.csv("Final CPTS/TotalMortality.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SMI_t2.p<-as.matrix(SMI_t2.p[,3:7])
colnames(SMI_t2.p) <- NULL
SMI_t2.pr<-array(t(SMI_t2.p),dim=c(5,5,5),dimnames = list(SMI_t2=SMI_indir, SMI_dir=SMI_dir, SMI_indir=SMI_indir))

SMI_t3.p<-read.csv("Final CPTS/SMI_t3.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SMI_t3.p<-as.matrix(SMI_t3.p[,3:7])
colnames(SMI_t3.p) <- NULL
SMI_t3.pr<-array(t(SMI_t3.p),dim=c(5,2,5),dimnames = list(SMI_t3=SMI_t3, Sediment_changes=SChan, SMI_t2=SMI_t2))

# Surface meiofauna

SM_dir.p<-read.csv("Final CPTS/SSI_dir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SM_dir.p<-as.matrix(SM_dir.p[,3:7])
colnames(SM_dir.p) <- NULL
SM_dir.pr<-array(t(SM_dir.p),dim=c(5,3,2),dimnames = list(SM_dir=SM_dir, Mining_intensity=MI,Mined=Mined)) 

SM_indir.p<-read.csv("Final CPTS/SM_indir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SM_indir.p<-as.matrix(SM_indir.p[,4:8])
colnames(SM_indir.p) <- NULL
SM_indir.pr<-array(t(SM_indir.p),dim=c(5,2,3,3),dimnames = list(SM_indir=SM_indir, Contaminant_release=CRel,Sediment_deposition=SDep, Suspended_sediment=SS)) 

SM_t2.p<-read.csv("Final CPTS/TotalMortality.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SM_t2.p<-as.matrix(SM_t2.p[,3:7])
colnames(SM_t2.p) <- NULL
SM_t2.pr<-array(t(SM_t2.p),dim=c(5,5,5),dimnames = list(SM_t2=SM_indir, SM_dir=SM_dir, SM_indir=SM_indir))

SM_t3.p<-read.csv("Final CPTS/SM_t3.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SM_t3.p<-as.matrix(SM_t3.p[,3:7])
colnames(SM_t3.p) <- NULL
SM_t3.pr<-array(t(SM_t3.p),dim=c(5,2,5),dimnames = list(SM_t3=SM_t3, Sediment_changes=SChan, SM_t2=SM_t2))

# Deep meiofauna

DM_dir.p<-read.csv("Final CPTS/SSI_dir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
DM_dir.p<-as.matrix(DM_dir.p[,3:7])
colnames(DM_dir.p) <- NULL
DM_dir.pr<-array(t(DM_dir.p),dim=c(5,3,2),dimnames = list(DM_dir=DM_dir, Mining_intensity=MI,Mined=Mined)) 

DM_indir.p<-read.csv("Final CPTS/DM_indir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
DM_indir.p<-as.matrix(DM_indir.p[,4:8])
colnames(DM_indir.p) <- NULL
DM_indir.pr<-array(t(DM_indir.p),dim=c(5,2,3,3),dimnames = list(DM_indir=DM_indir, Contaminant_release=CRel,Sediment_deposition=SDep, Suspended_sediment=SS)) 

DM_t2.p<-read.csv("Final CPTS/TotalMortality.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
DM_t2.p<-as.matrix(DM_t2.p[,3:7])
colnames(DM_t2.p) <- NULL
DM_t2.pr<-array(t(DM_t2.p),dim=c(5,5,5),dimnames = list(DM_t2=DM_indir, DM_dir=DM_dir, DM_indir=DM_indir))

DM_t3.p<-read.csv("Final CPTS/DM_t3.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
DM_t3.p<-as.matrix(DM_t3.p[,3:7])
colnames(DM_t3.p) <- NULL
DM_t3.pr<-array(t(DM_t3.p),dim=c(5,2,5),dimnames = list(DM_t3=DM_t3, Sediment_changes=SChan, DM_t2=DM_t2))

# Large mobile macrofauna
LMM_dir.p<-read.csv("Final CPTS/SSI_dir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
LMM_dir.p<-as.matrix(LMM_dir.p[,3:7])
colnames(LMM_dir.p) <- NULL
LMM_dir.pr<-array(t(LMM_dir.p),dim=c(5,3,2),dimnames = list(LMM_dir=LMM_dir, Mining_intensity=MI,Mined=Mined)) 

LMM_indir.p<-read.csv("Final CPTS/LMM_indir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
LMM_indir.p<-as.matrix(LMM_indir.p[,4:8])
colnames(LMM_indir.p) <- NULL
LMM_indir.pr<-array(t(LMM_indir.p),dim=c(5,2,3,3),dimnames = list(LMM_indir=LMM_indir, Contaminant_release=CRel,Sediment_deposition=SDep, Suspended_sediment=SS)) 

LMM_t2.p<-read.csv("Final CPTS/TotalMortality.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
LMM_t2.p<-as.matrix(LMM_t2.p[,3:7])
colnames(LMM_t2.p) <- NULL
LMM_t2.pr<-array(t(LMM_t2.p),dim=c(5,5,5),dimnames = list(LMM_t2=LMM_indir, LMM_dir=LMM_dir, LMM_indir=LMM_indir))

LMM_t3.p<-read.csv("Final CPTS/LMM_t3.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
LMM_t3.p<-as.matrix(LMM_t3.p[,3:7])
colnames(LMM_t3.p) <- NULL
LMM_t3.pr<-array(t(LMM_t3.p),dim=c(5,2,5),dimnames = list(LMM_t3=LMM_t3, Sediment_changes=SChan, LMM_t2=LMM_t2))

# Large sessile infauna

LSI_dir.p<-read.csv("Final CPTS/SSI_dir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
LSI_dir.p<-as.matrix(LSI_dir.p[,3:7])
colnames(LSI_dir.p) <- NULL
LSI_dir.pr<-array(t(LSI_dir.p),dim=c(5,3,2),dimnames = list(LSI_dir=LSI_dir, Mining_intensity=MI,Mined=Mined)) 

LSI_indir.p<-read.csv("Final CPTS/LSI_indir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
LSI_indir.p<-as.matrix(LSI_indir.p[,4:8])
colnames(LSI_indir.p) <- NULL
LSI_indir.pr<-array(t(LSI_indir.p),dim=c(5,2,3,3),dimnames = list(LSI_indir=LSI_indir, Contaminant_release=CRel,Sediment_deposition=SDep, Suspended_sediment=SS)) 

LSI_t2.p<-read.csv("Final CPTS/TotalMortality.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
LSI_t2.p<-as.matrix(LSI_t2.p[,3:7])
colnames(LSI_t2.p) <- NULL
LSI_t2.pr<-array(t(LSI_t2.p),dim=c(5,5,5),dimnames = list(LSI_t2=LSI_indir, LSI_dir=LSI_dir, LSI_indir=LSI_indir))

LSI_t3.p<-read.csv("Final CPTS/LSI_t3.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
LSI_t3.p<-as.matrix(LSI_t3.p[,3:7])
colnames(LSI_t3.p) <- NULL
LSI_t3.pr<-array(t(LSI_t3.p),dim=c(5,2,5),dimnames = list(LSI_t3=LSI_t3, Sediment_changes=SChan, LSI_t2=LSI_t2))

# Sessile soft-bodied megafauna

SSBM_dir.p<-read.csv("Final CPTS/SESF_dir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SSBM_dir.p<-as.matrix(SSBM_dir.p[,3:7])
colnames(SSBM_dir.p) <- NULL
SSBM_dir.pr<-array(t(SSBM_dir.p),dim=c(5,3,2),dimnames = list(SSBM_dir=SSBM_dir, Mining_intensity=MI,Mined=Mined)) 

SSBM_indir.p<-read.csv("Final CPTS/SSBM_indir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SSBM_indir.p<-as.matrix(SSBM_indir.p[,4:8])
colnames(SSBM_indir.p) <- NULL
SSBM_indir.pr<-array(t(SSBM_indir.p),dim=c(5,2,3,3),dimnames = list(SSBM_indir=SSBM_indir, Contaminant_release=CRel,Sediment_deposition=SDep, Suspended_sediment=SS)) 

SSBM_t2.p<-read.csv("Final CPTS/TotalMortality.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SSBM_t2.p<-as.matrix(SSBM_t2.p[,3:7])
colnames(SSBM_t2.p) <- NULL
SSBM_t2.pr<-array(t(SSBM_t2.p),dim=c(5,5,5),dimnames = list(SSBM_t2=SSBM_indir, SSBM_dir=SSBM_dir, SSBM_indir=SSBM_indir))

SSBM_t3.p<-read.csv("Final CPTS/SSBM_t3.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SSBM_t3.p<-as.matrix(SSBM_t3.p[,3:7])
colnames(SSBM_t3.p) <- NULL
SSBM_t3.pr<-array(t(SSBM_t3.p),dim=c(5,2,5),dimnames = list(SSBM_t3=SSBM_t3, Sediment_changes=SChan, SSBM_t2=SSBM_t2))


# Mobile predatory epifauna

MPE_dir.p<-read.csv("Final CPTS/GH_dir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
MPE_dir.p<-as.matrix(MPE_dir.p[,3:7])
colnames(MPE_dir.p) <- NULL
MPE_dir.pr<-array(t(MPE_dir.p),dim=c(5,3,2),dimnames = list(MPE_dir=MPE_dir, Mining_intensity=MI,Mined=Mined)) 

MPE_indir.p<-read.csv("Final CPTS/MPE_indir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
MPE_indir.p<-as.matrix(MPE_indir.p[,4:8])
colnames(MPE_indir.p) <- NULL
MPE_indir.pr<-array(t(MPE_indir.p),dim=c(5,2,3,3),dimnames = list(MPE_indir=MPE_indir, Contaminant_release=CRel,Sediment_deposition=SDep, Suspended_sediment=SS)) 

MPE_t2.p<-read.csv("Final CPTS/TotalMortality.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
MPE_t2.p<-as.matrix(MPE_t2.p[,3:7])
colnames(MPE_t2.p) <- NULL
MPE_t2.pr<-array(t(MPE_t2.p),dim=c(5,5,5),dimnames = list(MPE_t2=MPE_indir, MPE_dir=MPE_dir, MPE_indir=MPE_indir))

MPE_t3.p<-read.csv("Final CPTS/MPE_t3.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
MPE_t3.p<-as.matrix(MPE_t3.p[,3:7])
colnames(MPE_t3.p) <- NULL
MPE_t3.pr<-array(t(MPE_t3.p),dim=c(5,2,5),dimnames = list(MPE_t3=MPE_t3, Sediment_changes=SChan, MPE_t2=MPE_t2))

# Mobile grazing epifauna

MGE_dir.p<-read.csv("Final CPTS/GH_dir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
MGE_dir.p<-as.matrix(MGE_dir.p[,3:7])
colnames(MGE_dir.p) <- NULL
MGE_dir.pr<-array(t(MGE_dir.p),dim=c(5,3,2),dimnames = list(MGE_dir=MGE_dir, Mining_intensity=MI,Mined=Mined)) 

MGE_indir.p<-read.csv("Final CPTS/MGE_indir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
MGE_indir.p<-as.matrix(MGE_indir.p[,4:8])
colnames(MGE_indir.p) <- NULL
MGE_indir.pr<-array(t(MGE_indir.p),dim=c(5,2,3,3),dimnames = list(MGE_indir=MGE_indir, Contaminant_release=CRel,Sediment_deposition=SDep, Suspended_sediment=SS)) 

MGE_t2.p<-read.csv("Final CPTS/TotalMortality.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
MGE_t2.p<-as.matrix(MGE_t2.p[,3:7])
colnames(MGE_t2.p) <- NULL
MGE_t2.pr<-array(t(MGE_t2.p),dim=c(5,5,5),dimnames = list(MGE_t2=MGE_indir, MGE_dir=MGE_dir, MGE_indir=MGE_indir))

MGE_t3.p<-read.csv("Final CPTS/MGE_t3.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
MGE_t3.p<-as.matrix(MGE_t3.p[,3:7])
colnames(MGE_t3.p) <- NULL
MGE_t3.pr<-array(t(MGE_t3.p),dim=c(5,2,5),dimnames = list(MGE_t3=MGE_t3, Sediment_changes=SChan, MGE_t2=MGE_t2))

# Grazing hyperbenthos

GH_dir.p<-read.csv("Final CPTS/GH_dir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
GH_dir.p<-as.matrix(GH_dir.p[,3:7])
colnames(GH_dir.p) <- NULL
GH_dir.pr<-array(t(GH_dir.p),dim=c(5,3,2),dimnames = list(GH_dir=GH_dir, Mining_intensity=MI,Mined=Mined)) 

GH_indir.p<-read.csv("Final CPTS/GH_indir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
GH_indir.p<-as.matrix(GH_indir.p[,4:8])
colnames(GH_indir.p) <- NULL
GH_indir.pr<-array(t(GH_indir.p),dim=c(5,2,3,3),dimnames = list(GH_indir=GH_indir, Contaminant_release=CRel,Sediment_deposition=SDep, Suspended_sediment=SS)) 

GH_t2.p<-read.csv("Final CPTS/TotalMortality.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
GH_t2.p<-as.matrix(GH_t2.p[,3:7])
colnames(GH_t2.p) <- NULL
GH_t2.pr<-array(t(GH_t2.p),dim=c(5,5,5),dimnames = list(GH_t2=GH_indir, GH_dir=GH_dir, GH_indir=GH_indir))

GH_t3.p<-read.csv("Final CPTS/GH_t3.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
GH_t3.p<-as.matrix(GH_t3.p[,3:7])
colnames(GH_t3.p) <- NULL
GH_t3.pr<-array(t(GH_t3.p),dim=c(5,2,5),dimnames = list(GH_t3=GH_t3, Sediment_changes=SChan, GH_t2=GH_t2))

# Preadatory hyperbenthos

PH_dir.p<-read.csv("Final CPTS/PH_dir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
PH_dir.p<-as.matrix(PH_dir.p[,3:7])
colnames(PH_dir.p) <- NULL
PH_dir.pr<-array(t(PH_dir.p),dim=c(5,3,2),dimnames = list(PH_dir=PH_dir, Mining_intensity=MI,Mined=Mined)) 

PH_indir.p<-read.csv("Final CPTS/PH_indir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
PH_indir.p<-as.matrix(PH_indir.p[,4:8])
colnames(PH_indir.p) <- NULL
PH_indir.pr<-array(t(PH_indir.p),dim=c(5,2,3,3),dimnames = list(PH_indir=PH_indir, Contaminant_release=CRel,Sediment_deposition=SDep, Suspended_sediment=SS)) 

PH_t2.p<-read.csv("Final CPTS/TotalMortality.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
PH_t2.p<-as.matrix(PH_t2.p[,3:7])
colnames(PH_t2.p) <- NULL
PH_t2.pr<-array(t(PH_t2.p),dim=c(5,5,5),dimnames = list(PH_t2=PH_indir, PH_dir=PH_dir, PH_indir=PH_indir))

PH_t3.p<-read.csv("Final CPTS/PH_t3.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
PH_t3.p<-as.matrix(PH_t3.p[,3:7])
colnames(PH_t3.p) <- NULL
PH_t3.pr<-array(t(PH_t3.p),dim=c(5,2,5),dimnames = list(PH_t3=PH_t3, Sediment_changes=SChan, PH_t2=PH_t2))

# Sessile erect suspension feeders

SESF_dir.p<-read.csv("Final CPTS/SESF_dir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SESF_dir.p<-as.matrix(SESF_dir.p[,3:7])
colnames(SESF_dir.p) <- NULL
SESF_dir.pr<-array(t(SESF_dir.p),dim=c(5,3,2),dimnames = list(SESF_dir=SESF_dir, Mining_intensity=MI,Mined=Mined)) 

SESF_indir.p<-read.csv("Final CPTS/SESF_indir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SESF_indir.p<-as.matrix(SESF_indir.p[,4:8])
colnames(SESF_indir.p) <- NULL
SESF_indir.pr<-array(t(SESF_indir.p),dim=c(5,2,3,3),dimnames = list(SESF_indir=SESF_indir, Contaminant_release=CRel,Sediment_deposition=SDep, Suspended_sediment=SS)) 

SESF_t2.p<-read.csv("Final CPTS/TotalMortality.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SESF_t2.p<-as.matrix(SESF_t2.p[,3:7])
colnames(SESF_t2.p) <- NULL
SESF_t2.pr<-array(t(SESF_t2.p),dim=c(5,5,5),dimnames = list(SESF_t2=SESF_indir, SESF_dir=SESF_dir, SESF_indir=SESF_indir))

SESF_t3.p<-read.csv("Final CPTS/SESF_t3.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SESF_t3.p<-as.matrix(SESF_t3.p[,3:7])
colnames(SESF_t3.p) <- NULL
SESF_t3.pr<-array(t(SESF_t3.p),dim=c(5,2,5),dimnames = list(SESF_t3=SESF_t3, Sediment_changes=SChan, SESF_t2=SESF_t2))

# Sessile encrusting suspension feeders

SCSF_dir.p<-read.csv("Final CPTS/SESF_dir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SCSF_dir.p<-as.matrix(SCSF_dir.p[,3:7])
colnames(SCSF_dir.p) <- NULL
SCSF_dir.pr<-array(t(SCSF_dir.p),dim=c(5,3,2),dimnames = list(SCSF_dir=SCSF_dir, Mining_intensity=MI,Mined=Mined)) 

SCSF_indir.p<-read.csv("Final CPTS/SCSF_indir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SCSF_indir.p<-as.matrix(SCSF_indir.p[,4:8])
colnames(SCSF_indir.p) <- NULL
SCSF_indir.pr<-array(t(SCSF_indir.p),dim=c(5,2,3,3),dimnames = list(SCSF_indir=SCSF_indir, Contaminant_release=CRel,Sediment_deposition=SDep, Suspended_sediment=SS)) 

SCSF_t2.p<-read.csv("Final CPTS/TotalMortality.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SCSF_t2.p<-as.matrix(SCSF_t2.p[,3:7])
colnames(SCSF_t2.p) <- NULL
SCSF_t2.pr<-array(t(SCSF_t2.p),dim=c(5,5,5),dimnames = list(SCSF_t2=SCSF_indir, SCSF_dir=SCSF_dir, SCSF_indir=SCSF_indir))

SCSF_t3.p<-read.csv("Final CPTS/SCSF_t3.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SCSF_t3.p<-as.matrix(SCSF_t3.p[,3:7])
colnames(SCSF_t3.p) <- NULL
SCSF_t3.pr<-array(t(SCSF_t3.p),dim=c(5,2,5),dimnames = list(SCSF_t3=SCSF_t3, Sediment_changes=SChan, SCSF_t2=SCSF_t2))

# Sessile erect filter feeders

SEFF_dir.p<-read.csv("Final CPTS/SESF_dir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SEFF_dir.p<-as.matrix(SEFF_dir.p[,3:7])
colnames(SEFF_dir.p) <- NULL
SEFF_dir.pr<-array(t(SEFF_dir.p),dim=c(5,3,2),dimnames = list(SEFF_dir=SEFF_dir, Mining_intensity=MI,Mined=Mined)) 

SEFF_indir.p<-read.csv("Final CPTS/SEFF_indir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SEFF_indir.p<-as.matrix(SEFF_indir.p[,4:8])
colnames(SEFF_indir.p) <- NULL
SEFF_indir.pr<-array(t(SEFF_indir.p),dim=c(5,2,3,3),dimnames = list(SEFF_indir=SEFF_indir, Contaminant_release=CRel,Sediment_deposition=SDep, Suspended_sediment=SS)) 

SEFF_t2.p<-read.csv("Final CPTS/TotalMortality.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SEFF_t2.p<-as.matrix(SEFF_t2.p[,3:7])
colnames(SEFF_t2.p) <- NULL
SEFF_t2.pr<-array(t(SEFF_t2.p),dim=c(5,5,5),dimnames = list(SEFF_t2=SEFF_indir, SEFF_dir=SEFF_dir, SEFF_indir=SEFF_indir))

SEFF_t3.p<-read.csv("Final CPTS/SEFF_t3.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SEFF_t3.p<-as.matrix(SEFF_t3.p[,3:7])
colnames(SEFF_t3.p) <- NULL
SEFF_t3.pr<-array(t(SEFF_t3.p),dim=c(5,2,5),dimnames = list(SEFF_t3=SEFF_t3, Sediment_changes=SChan, SEFF_t2=SEFF_t2))

# Sessile encrusting filter feeders

SCFF_dir.p<-read.csv("Final CPTS/SESF_dir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SCFF_dir.p<-as.matrix(SCFF_dir.p[,3:7])
colnames(SCFF_dir.p) <- NULL
SCFF_dir.pr<-array(t(SCFF_dir.p),dim=c(5,3,2),dimnames = list(SCFF_dir=SCFF_dir, Mining_intensity=MI,Mined=Mined)) 

SCFF_indir.p<-read.csv("Final CPTS/SCFF_indir.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SCFF_indir.p<-as.matrix(SCFF_indir.p[,4:8])
colnames(SCFF_indir.p) <- NULL
SCFF_indir.pr<-array(t(SCFF_indir.p),dim=c(5,2,3,3),dimnames = list(SCFF_indir=SCFF_indir, Contaminant_release=CRel,Sediment_deposition=SDep, Suspended_sediment=SS)) 

SCFF_t2.p<-read.csv("Final CPTS/TotalMortality.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SCFF_t2.p<-as.matrix(SCFF_t2.p[,3:7])
colnames(SCFF_t2.p) <- NULL
SCFF_t2.pr<-array(t(SCFF_t2.p),dim=c(5,5,5),dimnames = list(SCFF_t2=SCFF_indir, SCFF_dir=SCFF_dir, SCFF_indir=SCFF_indir))

SCFF_t3.p<-read.csv("Final CPTS/SCFF_t3.csv", header=TRUE,fill=TRUE, sep=',', stringsAsFactors = FALSE)
SCFF_t3.p<-as.matrix(SCFF_t3.p[,3:7])
colnames(SCFF_t3.p) <- NULL
SCFF_t3.pr<-array(t(SCFF_t3.p),dim=c(5,2,5),dimnames = list(SCFF_t3=SCFF_t3, Sediment_changes=SChan, SCFF_t2=SCFF_t2))

#--------------------------------------------------------------------------#

# COMBINING CPTS AND MODEL STRUCTURE TO A BN

# Combine dag and the local distributions into a CPT

cpt2 <- list(Mined=Mp, 
            Mining_intensity=MIp,
            Sediment_type = STp,
            Sediment_contaminants =SC.pr,
            Plume_release = PRp,
            Nodule_removal=NRp,
            Depth_extraction=DEp,
            Volume_extraction = VE.pr,
            Distance_mining_site=DistMS.pr,
            Contaminant_release = CRel.pr,
            Sediment_changes=SCh.pr,
            Suspended_sediment= SS.pr,
            Sediment_deposition=SD.pr,
            SSI_indir=SSI_indir.pr,
            SSI_dir=SSI_dir.pr,
            SSI_t2=SSI_t2.pr,
            SSI_t3=SSI_t3.pr,
            SMI_indir=SMI_indir.pr,
            SMI_dir=SMI_dir.pr,
            SMI_t2=SMI_t2.pr,
            SMI_t3=SMI_t3.pr,
            SSBM_indir=SSBM_indir.pr,
            SSBM_dir=SSBM_dir.pr,
            SSBM_t2=SSBM_t2.pr,
            SSBM_t3=SSBM_t3.pr, 
            LMM_indir=LMM_indir.pr,
            LMM_dir=LMM_dir.pr,
            LMM_t2=LMM_t2.pr,
            LMM_t3=LMM_t3.pr,
            GH_indir=GH_indir.pr,
            GH_dir=GH_dir.pr,
            GH_t2=GH_t2.pr,
            GH_t3=GH_t3.pr,
            SESF_indir=SESF_indir.pr,
            SESF_dir=SESF_dir.pr,
            SESF_t2=SESF_t2.pr,
            SESF_t3=SESF_t3.pr,
            DM_indir=DM_indir.pr,
            DM_dir=DM_dir.pr,
            DM_t2=DM_t2.pr,
            DM_t3=DM_t3.pr,
            SM_indir=SM_indir.pr,
            SM_dir=SM_dir.pr,
            SM_t2=SM_t2.pr,
            SM_t3=SM_t3.pr,
            LSI_indir=LSI_indir.pr,
            LSI_dir=LSI_dir.pr,
            LSI_t2=LSI_t2.pr,
            LSI_t3=LSI_t3.pr,
            MPE_indir=MPE_indir.pr,
            MPE_dir=MPE_dir.pr,
            MPE_t2=MPE_t2.pr,
            MPE_t3=MPE_t3.pr,
            MGE_indir=MGE_indir.pr,
            MGE_dir=MGE_dir.pr,
            MGE_t2=MGE_t2.pr,
            MGE_t3=MGE_t3.pr,
            PH_indir=PH_indir.pr,
            PH_dir=PH_dir.pr,
            PH_t2=PH_t2.pr,
            PH_t3=PH_t3.pr,
            SEFF_indir=SEFF_indir.pr,
            SEFF_dir=SEFF_dir.pr,
            SEFF_t2=SEFF_t2.pr,
            SEFF_t3=SEFF_t3.pr,
            SCFF_indir=SCFF_indir.pr,
            SCFF_dir=SCFF_dir.pr,
            SCFF_t2=SCFF_t2.pr,
            SCFF_t3=SCFF_t3.pr,
            SCSF_indir=SCSF_indir.pr,
            SCSF_dir=SCSF_dir.pr,
            SCSF_t2=SCSF_t2.pr,
            SCSF_t3=SCSF_t3.pr
            )
CRbn <- custom.fit(dagCR, cpt2)

junction <- compile(as.grain(CRbn))

#--------------------------------------------------------------------------#

# USING THE BN FOR PROBABILITY QUERIES

# Probability of a node given evidence

###############################################
##### SCENARIO 1 HIGH DISTURBANCE #######
################################################

# Change mining intensity to illusatre changes in pressures
# Set evidence for each spatial domain separately to obtain results for all

Sce1in <- setEvidence(junction, nodes = c("Mining_intensity","Depth_extraction", "Distance_mining_site", "Plume_release"), states =c("100%", ">30cm", "Inside block", "At the seafloor"))
Sce1out <- setEvidence(junction, nodes = c("Mining_intensity","Depth_extraction", "Distance_mining_site", "Plume_release"), states =c("100%", ">30cm", "Outside block", "At the seafloor"))
Sce1near <- setEvidence(junction, nodes = c("Mining_intensity","Depth_extraction", "Distance_mining_site", "Plume_release"), states =c("100%", ">30cm", "Nearfield", "At the seafloor"))

# example of checking query results for a specific node
querygrain(Sce1out, nodes = "Sediment_deposition")$Sediment_deposition

#--------------------------------------------------------------------------#

# VISUALISING THE RESULTS FOR SCENARIO 1

# Physicochemical parameters

# Probabilities of pressures from mining
r11<-as.data.frame.table(querygrain(Sce1in, nodes = "Sediment_deposition")$Sediment_deposition)
r12<-as.data.frame.table(querygrain(Sce1out, nodes = "Sediment_deposition")$Sediment_deposition)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "Sediment_deposition")$Sediment_deposition)

SDdf<-cbind(r11,r12n[,2],r12[,2]) # combine spatial domains to one dataframe
names(SDdf)<-c("Sediment deposition", "In", "Nearfield","Farfield")
SDdf2<-melt(SDdf)
names(SDdf2)<-c("Sediment_deposition", "Spatial_domain", "Probability")

SDepp<-ggplot(SDdf2,aes(x = Sediment_deposition, y=Probability, fill = Spatial_domain)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_fill_manual(values=c("#FC4E07", "#E7B800", "#00AFBB"))+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")

#old colours
#SDepp<-ggplot(SDdf2,aes(x = Sediment_deposition, y=Probability, fill = Spatial_domain)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_fill_manual(values=c("#FC4E07", "#E7B800", "#00AFBB"))+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")

r13<-as.data.frame.table(querygrain(Sce1in, nodes = "Suspended_sediment")$Suspended_sediment)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "Suspended_sediment")$Suspended_sediment)
r14n<-as.data.frame.table(querygrain(Sce1near, nodes = "Suspended_sediment")$Suspended_sediment)

SSdf<-cbind(r13, r14n[,2],r14[,2]) # combine spatial domains to one dataframe
names(SSdf)<-c("Suspended sediment","In","Nearfield","Farfield")
SSdf2<-melt(SSdf)
names(SSdf2)<-c("Suspended_sediment", "Spatial_domain", "Probability")

SSp<-ggplot(SSdf2,aes(x = Suspended_sediment, y=Probability, fill = Spatial_domain)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_fill_manual(values=c("#FC4E07", "#E7B800", "#00AFBB"))+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")

# Contaminant release
r15<-as.data.frame.table(querygrain(Sce1in, nodes = "Contaminant_release")$Contaminant_release)
r16<-as.data.frame.table(querygrain(Sce1out, nodes = "Contaminant_release")$Contaminant_release)
r16n<-as.data.frame.table(querygrain(Sce1near, nodes = "Contaminant_release")$Contaminant_release)

CRdf<-cbind(r15, r16n[,2],r16[,2]) # combine spatial domains to one dataframe
names(CRdf)<-c("Contaminant_release","In","Nearfield","Farfield")
CRdf2<-melt(CRdf)
names(CRdf2)<-c("Contaminant_release", "Spatial_domain", "Probability")

CRp<-ggplot(CRdf2,aes(x =Contaminant_release, y=Probability, fill = Spatial_domain)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_fill_manual(values=c("#FC4E07", "#E7B800", "#00AFBB"))+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")

# Habitat change
r17<-as.data.frame.table(querygrain(Sce1in, nodes = "Sediment_changes")$Sediment_changes)
r18<-as.data.frame.table(querygrain(Sce1out, nodes = "Sediment_changes")$Sediment_changes)
r18n<-as.data.frame.table(querygrain(Sce1near, nodes = "Sediment_changes")$Sediment_changes)

SCdf<-cbind(r17, r18n[,2],r18[,2]) # combine spatial domains to one dataframe
names(SCdf)<-c("Sediment changes","In","Nearfield","Farfield")
SCdf2<-melt(SCdf)
names(SCdf2)<-c("Sediment_changes", "Spatial_domain", "Probability")

SCp<-ggplot(SCdf2,aes(x=Sediment_changes, y=Probability, fill = Spatial_domain)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_fill_manual(values=c("#FC4E07", "#E7B800", "#00AFBB"))+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")


pressuresS1<-ggarrange(SSp+theme(plot.title = element_text(hjust = 0.5))+labs(title="Suspended sediment")+theme(legend.position="none"),
          SDepp +theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank())+labs(y="SMI", title="Sediment deposition")+theme(legend.position="none"),
           CRp+theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank())+labs(y="SMI", title="Contaminant release")+theme(legend.position="none"),
          SCp+theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank())+labs(y="SMI", title="Sediment changes"),
          nrow=1,widths = c(1,1,1,1), common.legend=TRUE)

png(file="C:/Users/kaikkonenl/Documents/CR ERA model/Data/Plots/Pressures-S1.png",
    width=1200, height=600)
pressuresS1
dev.off()

# IMPACTS ON BENTHIC FAUNA
# SMI

r11<-as.data.frame.table(querygrain(Sce1in, nodes = "SMI_t2")$SMI_t2)
r12<-as.data.frame.table(querygrain(Sce1in, nodes = "SMI_t3")$SMI_t3)

  SMIdf<-cbind(r11, r12[,2]) # combine timesteps to one dataframe
  names(SMIdf)<-c("Mortality", "Immediately", "1 year after")
  SMIdf2<-melt(SMIdf)
  names(SMIdf2)<-c("Mortality", "Timestep", "Probability")
  
r13<-as.data.frame.table(querygrain(Sce1out, nodes = "SMI_t2")$SMI_t2)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "SMI_t3")$SMI_t3)

  SMIdf3<-cbind(r13, r14[,2]) # combine timesteps to one dataframe
  names(SMIdf3)<-c("Mortality", "Immediately", "1 year after")
  SMIdf4<-melt(SMIdf3)
  names(SMIdf4)<-c("Mortality", "Timestep", "Probability")
  
  r15<-as.data.frame.table(querygrain(Sce1near, nodes = "SMI_t2")$SMI_t2)
  r16<-as.data.frame.table(querygrain(Sce1near, nodes = "SMI_t3")$SMI_t3)
  
  SMIdf5<-cbind(r15, r16[,2]) # combine timesteps to one dataframe
  names(SMIdf5)<-c("Mortality", "Immediately", "1 year after")
  SMIdf6<-melt(SMIdf5)
  names(SMIdf6)<-c("Mortality", "Timestep", "Probability")
  
r11n<-as.data.frame.table(querygrain(Sce1near, nodes = "SMI_t2")$SMI_t2)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "SMI_t3")$SMI_t3)
  
  SMIdf7<-cbind(r11n, r12n[,2]) # combine timesteps to one dataframe
  names(SMIdf7)<-c("Mortality", "Immediately", "1 year after")
  SMIdf8<-melt(SMIdf7)
  names(SMIdf8)<-c("Mortality", "Timestep", "Probability")

S1_SMI_near<-ggplot(SMIdf8,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S1_SMI_in<-ggplot(SMIdf2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))

S1_SMI_out<-ggplot(SMIdf4,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# SSI

r11<-as.data.frame.table(querygrain(Sce1in, nodes = "SSI_t2")$SSI_t2)
r12<-as.data.frame.table(querygrain(Sce1in, nodes = "SSI_t3")$SSI_t3)

SSIdf<-cbind(r11, r12[,2]) # combine timesteps to one dataframe
names(SSIdf)<-c("Mortality", "Immediately", "1 year after")
SSIdf2<-melt(SSIdf)
names(SSIdf2)<-c("Mortality", "Timestep", "Probability")

r13<-as.data.frame.table(querygrain(Sce1out, nodes = "SSI_t2")$SSI_t2)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "SSI_t3")$SSI_t3)

SSIdf3<-cbind(r13, r14[,2]) # combine timesteps to one dataframe
names(SSIdf3)<-c("Mortality", "Immediately", "1 year after")
SSIdf4<-melt(SSIdf3)
names(SSIdf4)<-c("Mortality", "Timestep", "Probability")

r15<-as.data.frame.table(querygrain(Sce1near, nodes = "SSI_t2")$SSI_t2)
r16<-as.data.frame.table(querygrain(Sce1near, nodes = "SSI_t3")$SSI_t3)

SSIdf5<-cbind(r15, r16[,2]) # combine timesteps to one dataframe
names(SSIdf5)<-c("Mortality", "Immediately", "1 year after")
SSIdf6<-melt(SSIdf5)
names(SSIdf6)<-c("Mortality", "Timestep", "Probability")

r11n<-as.data.frame.table(querygrain(Sce1near, nodes = "SSI_t2")$SSI_t2)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "SSI_t3")$SSI_t3)

SSIdf7<-cbind(r11n, r12n[,2]) # combine timesteps to one dataframe
names(SSIdf7)<-c("Mortality", "Immediately", "1 year after")
SSIdf8<-melt(SSIdf7)
names(SSIdf8)<-c("Mortality", "Timestep", "Probability")

S1_SSI_near<-ggplot(SSIdf8,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S1_SSI_in<-ggplot(SSIdf2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))

S1_SSI_out<-ggplot(SSIdf4,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# Surface meiofauna

r11<-as.data.frame.table(querygrain(Sce1in, nodes = "SM_t2")$SM_t2)
r12<-as.data.frame.table(querygrain(Sce1in, nodes = "SM_t3")$SM_t3)

SMdf<-cbind(r11, r12[,2]) # combine timesteps to one dataframe
names(SMdf)<-c("Mortality", "Immediately", "1 year after")
SMdf2<-melt(SMdf)
names(SMdf2)<-c("Mortality", "Timestep", "Probability")

r13<-as.data.frame.table(querygrain(Sce1out, nodes = "SM_t2")$SM_t2)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "SM_t3")$SM_t3)

SMdf3<-cbind(r13, r14[,2]) # combine timesteps to one dataframe
names(SMdf3)<-c("Mortality", "Immediately", "1 year after")
SMdf4<-melt(SMdf3)
names(SMdf4)<-c("Mortality", "Timestep", "Probability")

r15<-as.data.frame.table(querygrain(Sce1near, nodes = "SM_t2")$SM_t2)
r16<-as.data.frame.table(querygrain(Sce1near, nodes = "SM_t3")$SM_t3)

SMdf5<-cbind(r15, r16[,2]) # combine timesteps to one dataframe
names(SMdf5)<-c("Mortality", "Immediately", "1 year after")
SMdf6<-melt(SMdf5)
names(SMdf6)<-c("Mortality", "Timestep", "Probability")

r11n<-as.data.frame.table(querygrain(Sce1near, nodes = "SM_t2")$SM_t2)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "SM_t3")$SM_t3)

SMdf7<-cbind(r11n, r12n[,2]) # combine timesteps to one dataframe
names(SMdf7)<-c("Mortality", "Immediately", "1 year after")
SMdf8<-melt(SMdf7)
names(SMdf8)<-c("Mortality", "Timestep", "Probability")

S1_SM_near<-ggplot(SMdf8,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S1_SM_in<-ggplot(SMdf2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))

S1_SM_out<-ggplot(SMdf4,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# Deep meiofauna

r11<-as.data.frame.table(querygrain(Sce1in, nodes = "DM_t2")$DM_t2)
r12<-as.data.frame.table(querygrain(Sce1in, nodes = "DM_t3")$DM_t3)

DMdf<-cbind(r11, r12[,2]) # combine timesteps to one dataframe
names(DMdf)<-c("Mortality", "Immediately", "1 year after")
DMdf2<-melt(DMdf)
names(DMdf2)<-c("Mortality", "Timestep", "Probability")

r13<-as.data.frame.table(querygrain(Sce1out, nodes = "DM_t2")$DM_t2)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "DM_t3")$DM_t3)

DMdf3<-cbind(r13, r14[,2]) # combine timesteps to one dataframe
names(DMdf3)<-c("Mortality", "Immediately", "1 year after")
DMdf4<-melt(DMdf3)
names(DMdf4)<-c("Mortality", "Timestep", "Probability")

r15<-as.data.frame.table(querygrain(Sce1near, nodes = "DM_t2")$DM_t2)
r16<-as.data.frame.table(querygrain(Sce1near, nodes = "DM_t3")$DM_t3)

DMdf5<-cbind(r15, r16[,2]) # combine timesteps to one dataframe
names(DMdf5)<-c("Mortality", "Immediately", "1 year after")
DMdf6<-melt(DMdf5)
names(DMdf6)<-c("Mortality", "Timestep", "Probability")

r11n<-as.data.frame.table(querygrain(Sce1near, nodes = "DM_t2")$DM_t2)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "DM_t3")$DM_t3)

DMdf7<-cbind(r11n, r12n[,2]) # combine timesteps to one dataframe
names(DMdf7)<-c("Mortality", "Immediately", "1 year after")
DMdf8<-melt(DMdf7)
names(DMdf8)<-c("Mortality", "Timestep", "Probability")

S1_DM_near<-ggplot(DMdf8,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S1_DM_in<-ggplot(DMdf2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))

S1_DM_out<-ggplot(DMdf4,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# Small soft-bodied megafauna
r11<-as.data.frame.table(querygrain(Sce1in, nodes = "SSBM_t2")$SSBM_t2)
r12<-as.data.frame.table(querygrain(Sce1in, nodes = "SSBM_t3")$SSBM_t3)

SSBMdf<-cbind(r11, r12[,2]) # combine timesteps to one dataframe
names(SSBMdf)<-c("Mortality", "Immediately", "1 year after")
SSBMdf2<-melt(SSBMdf)
names(SSBMdf2)<-c("Mortality", "Timestep", "Probability")

r13<-as.data.frame.table(querygrain(Sce1out, nodes = "SSBM_t2")$SSBM_t2)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "SSBM_t3")$SSBM_t3)

SSBMdf3<-cbind(r13, r14[,2]) # combine timesteps to one dataframe
names(SSBMdf3)<-c("Mortality", "Immediately", "1 year after")
SSBMdf4<-melt(SSBMdf3)
names(SSBMdf4)<-c("Mortality", "Timestep", "Probability")

r11n<-as.data.frame.table(querygrain(Sce1near, nodes = "SSBM_t2")$SSBM_t2)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "SSBM_t3")$SSBM_t3)

SSBMdf7<-cbind(r11n, r12n[,2]) # combine timesteps to one dataframe
names(SSBMdf7)<-c("Mortality", "Immediately", "1 year after")
SSBMdf8<-melt(SSBMdf7)
names(SSBMdf8)<-c("Mortality", "Timestep", "Probability")

S1_SSBM_near<-ggplot(SSBMdf8,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S1_SSBM_in<-ggplot(SSBMdf2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))

S1_SSBM_out<-ggplot(SSBMdf4,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# SCSF
r11<-as.data.frame.table(querygrain(Sce1in, nodes = "SCSF_t2")$SCSF_t2)
r12<-as.data.frame.table(querygrain(Sce1in, nodes = "SCSF_t3")$SCSF_t3)

SCSFdf<-cbind(r11, r12[,2]) # combine timesteps to one dataframe
names(SCSFdf)<-c("Mortality", "Immediately", "1 year after")
SCSFdf2<-melt(SCSFdf)
names(SCSFdf2)<-c("Mortality", "Timestep", "Probability")

r13<-as.data.frame.table(querygrain(Sce1out, nodes = "SCSF_t2")$SCSF_t2)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "SCSF_t3")$SCSF_t3)

SCSFdf3<-cbind(r13, r14[,2]) # combine timesteps to one dataframe
names(SCSFdf3)<-c("Mortality", "Immediately", "1 year after")
SCSFdf4<-melt(SCSFdf3)
names(SCSFdf4)<-c("Mortality", "Timestep", "Probability")

r11n<-as.data.frame.table(querygrain(Sce1near, nodes = "SCSF_t2")$SCSF_t2)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "SCSF_t3")$SCSF_t3)

SCSFdf7<-cbind(r11n, r12n[,2]) # combine timesteps to one dataframe
names(SCSFdf7)<-c("Mortality", "Immediately", "1 year after")
SCSFdf8<-melt(SCSFdf7)
names(SCSFdf8)<-c("Mortality", "Timestep", "Probability")

S1_SCSF_near<-ggplot(SCSFdf8,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S1_SCSF_in<-ggplot(SCSFdf2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))

S1_SCSF_out<-ggplot(SCSFdf4,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# SEFF
r11<-as.data.frame.table(querygrain(Sce1in, nodes = "SEFF_t2")$SEFF_t2)
r12<-as.data.frame.table(querygrain(Sce1in, nodes = "SEFF_t3")$SEFF_t3)

SEFFdf<-cbind(r11, r12[,2]) # combine timesteps to one dataframe
names(SEFFdf)<-c("Mortality", "Immediately", "1 year after")
SEFFdf2<-melt(SEFFdf)
names(SEFFdf2)<-c("Mortality", "Timestep", "Probability")

r13<-as.data.frame.table(querygrain(Sce1out, nodes = "SEFF_t2")$SEFF_t2)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "SEFF_t3")$SEFF_t3)

SEFFdf3<-cbind(r13, r14[,2]) # combine timesteps to one dataframe
names(SEFFdf3)<-c("Mortality", "Immediately", "1 year after")
SEFFdf4<-melt(SEFFdf3)
names(SEFFdf4)<-c("Mortality", "Timestep", "Probability")

r11n<-as.data.frame.table(querygrain(Sce1near, nodes = "SEFF_t2")$SEFF_t2)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "SEFF_t3")$SEFF_t3)

SEFFdf7<-cbind(r11n, r12n[,2]) # combine timesteps to one dataframe
names(SEFFdf7)<-c("Mortality", "Immediately", "1 year after")
SEFFdf8<-melt(SEFFdf7)
names(SEFFdf8)<-c("Mortality", "Timestep", "Probability")

S1_SEFF_near<-ggplot(SEFFdf8,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S1_SEFF_in<-ggplot(SEFFdf2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))

S1_SEFF_out<-ggplot(SEFFdf4,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# SCFF
r11<-as.data.frame.table(querygrain(Sce1in, nodes = "SCFF_t2")$SCFF_t2)
r12<-as.data.frame.table(querygrain(Sce1in, nodes = "SCFF_t3")$SCFF_t3)

SCFFdf<-cbind(r11, r12[,2]) # combine timesteps to one dataframe
names(SCFFdf)<-c("Mortality", "Immediately", "1 year after")
SCFFdf2<-melt(SCFFdf)
names(SCFFdf2)<-c("Mortality", "Timestep", "Probability")

r13<-as.data.frame.table(querygrain(Sce1out, nodes = "SCFF_t2")$SCFF_t2)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "SCFF_t3")$SCFF_t3)

SCFFdf3<-cbind(r13, r14[,2]) # combine timesteps to one dataframe
names(SCFFdf3)<-c("Mortality", "Immediately", "1 year after")
SCFFdf4<-melt(SCFFdf3)
names(SCFFdf4)<-c("Mortality", "Timestep", "Probability")

r11n<-as.data.frame.table(querygrain(Sce1near, nodes = "SCFF_t2")$SCFF_t2)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "SCFF_t3")$SCFF_t3)

SCFFdf7<-cbind(r11n, r12n[,2]) # combine timesteps to one dataframe
names(SCFFdf7)<-c("Mortality", "Immediately", "1 year after")
SCFFdf8<-melt(SCFFdf7)
names(SCFFdf8)<-c("Mortality", "Timestep", "Probability")

S1_SCFF_near<-ggplot(SCFFdf8,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S1_SCFF_in<-ggplot(SCFFdf2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))

S1_SCFF_out<-ggplot(SCFFdf4,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# SESF
r11<-as.data.frame.table(querygrain(Sce1in, nodes = "SESF_t2")$SESF_t2)
r12<-as.data.frame.table(querygrain(Sce1in, nodes = "SESF_t3")$SESF_t3)

SESFdf<-cbind(r11, r12[,2]) # combine timesteps to one dataframe
names(SESFdf)<-c("Mortality", "Immediately", "1 year after")
SESFdf2<-melt(SESFdf)
names(SESFdf2)<-c("Mortality", "Timestep", "Probability")

r13<-as.data.frame.table(querygrain(Sce1out, nodes = "SESF_t2")$SESF_t2)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "SESF_t3")$SESF_t3)

SESFdf3<-cbind(r13, r14[,2]) # combine timesteps to one dataframe
names(SESFdf3)<-c("Mortality", "Immediately", "1 year after")
SESFdf4<-melt(SESFdf3)
names(SESFdf4)<-c("Mortality", "Timestep", "Probability")


r11n<-as.data.frame.table(querygrain(Sce1near, nodes = "SESF_t2")$SESF_t2)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "SESF_t3")$SESF_t3)

SESFdf7<-cbind(r11n, r12n[,2]) # combine timesteps to one dataframe
names(SESFdf7)<-c("Mortality", "Immediately", "1 year after")
SESFdf8<-melt(SESFdf7)
names(SESFdf8)<-c("Mortality", "Timestep", "Probability")

S1_SESF_near<-ggplot(SESFdf8,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S1_SESF_in<-ggplot(SESFdf2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))
S1_SESF_out<-ggplot(SESFdf4,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# LMM

r11<-as.data.frame.table(querygrain(Sce1in, nodes = "LMM_t2")$LMM_t2)
r12<-as.data.frame.table(querygrain(Sce1in, nodes = "LMM_t3")$LMM_t3)

LMMdf<-cbind(r11, r12[,2]) # combine timesteps to one dataframe
names(LMMdf)<-c("Mortality", "Immediately", "1 year after")
LMMdf2<-melt(LMMdf)
names(LMMdf2)<-c("Mortality", "Timestep", "Probability")

r13<-as.data.frame.table(querygrain(Sce1out, nodes = "LMM_t2")$LMM_t2)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "LMM_t3")$LMM_t3)

LMMdf3<-cbind(r13, r14[,2]) # combine timesteps to one dataframe
names(LMMdf3)<-c("Mortality", "Immediately", "1 year after")
LMMdf4<-melt(LMMdf3)
names(LMMdf4)<-c("Mortality", "Timestep", "Probability")

r11n<-as.data.frame.table(querygrain(Sce1near, nodes = "LMM_t2")$LMM_t2)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "LMM_t3")$LMM_t3)

LMMdf7<-cbind(r11n, r12n[,2]) # combine timesteps to one dataframe
names(LMMdf7)<-c("Mortality", "Immediately", "1 year after")
LMMdf8<-melt(LMMdf7)
names(LMMdf8)<-c("Mortality", "Timestep", "Probability")

S1_LMM_near<-ggplot(LMMdf8,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S1_LMM_in<-ggplot(LMMdf2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))
S1_LMM_out<-ggplot(LMMdf4,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# LSI

r11<-as.data.frame.table(querygrain(Sce1in, nodes = "LSI_t2")$LSI_t2)
r12<-as.data.frame.table(querygrain(Sce1in, nodes = "LSI_t3")$LSI_t3)

LSIdf<-cbind(r11, r12[,2]) # combine timesteps to one dataframe
names(LSIdf)<-c("Mortality", "Immediately", "1 year after")
LSIdf2<-melt(LSIdf)
names(LSIdf2)<-c("Mortality", "Timestep", "Probability")

r13<-as.data.frame.table(querygrain(Sce1out, nodes = "LSI_t2")$LSI_t2)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "LSI_t3")$LSI_t3)

LSIdf3<-cbind(r13, r14[,2]) # combine timesteps to one dataframe
names(LSIdf3)<-c("Mortality", "Immediately", "1 year after")
LSIdf4<-melt(LSIdf3)
names(LSIdf4)<-c("Mortality", "Timestep", "Probability")

r15<-as.data.frame.table(querygrain(Sce1near, nodes = "LSI_t2")$LSI_t2)
r16<-as.data.frame.table(querygrain(Sce1near, nodes = "LSI_t3")$LSI_t3)

LSIdf5<-cbind(r15, r16[,2]) # combine timesteps to one dataframe
names(LSIdf5)<-c("Mortality", "Immediately", "1 year after")
LSIdf6<-melt(LSIdf5)
names(LSIdf6)<-c("Mortality", "Timestep", "Probability")

r11n<-as.data.frame.table(querygrain(Sce1near, nodes = "LSI_t2")$LSI_t2)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "LSI_t3")$LSI_t3)

LSIdf7<-cbind(r11n, r12n[,2]) # combine timesteps to one dataframe
names(LSIdf7)<-c("Mortality", "Immediately", "1 year after")
LSIdf8<-melt(LSIdf7)
names(LSIdf8)<-c("Mortality", "Timestep", "Probability")

S1_LSI_near<-ggplot(LSIdf8,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S1_LSI_in<-ggplot(LSIdf2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))

S1_LSI_out<-ggplot(LSIdf4,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# MPE

r11<-as.data.frame.table(querygrain(Sce1in, nodes = "MPE_t2")$MPE_t2)
r12<-as.data.frame.table(querygrain(Sce1in, nodes = "MPE_t3")$MPE_t3)

MPEdf<-cbind(r11, r12[,2]) # combine timesteps to one dataframe
names(MPEdf)<-c("Mortality", "Immediately", "1 year after")
MPEdf2<-melt(MPEdf)
names(MPEdf2)<-c("Mortality", "Timestep", "Probability")

r13<-as.data.frame.table(querygrain(Sce1out, nodes = "MPE_t2")$MPE_t2)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "MPE_t3")$MPE_t3)

MPEdf3<-cbind(r13, r14[,2]) # combine timesteps to one dataframe
names(MPEdf3)<-c("Mortality", "Immediately", "1 year after")
MPEdf4<-melt(MPEdf3)
names(MPEdf4)<-c("Mortality", "Timestep", "Probability")

r15<-as.data.frame.table(querygrain(Sce1near, nodes = "MPE_t2")$MPE_t2)
r16<-as.data.frame.table(querygrain(Sce1near, nodes = "MPE_t3")$MPE_t3)

MPEdf5<-cbind(r15, r16[,2]) # combine timesteps to one dataframe
names(MPEdf5)<-c("Mortality", "Immediately", "1 year after")
MPEdf6<-melt(MPEdf5)
names(MPEdf6)<-c("Mortality", "Timestep", "Probability")

r11n<-as.data.frame.table(querygrain(Sce1near, nodes = "MPE_t2")$MPE_t2)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "MPE_t3")$MPE_t3)

MPEdf7<-cbind(r11n, r12n[,2]) # combine timesteps to one dataframe
names(MPEdf7)<-c("Mortality", "Immediately", "1 year after")
MPEdf8<-melt(MPEdf7)
names(MPEdf8)<-c("Mortality", "Timestep", "Probability")

S1_MPE_near<-ggplot(MPEdf8,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S1_MPE_in<-ggplot(MPEdf2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))

S1_MPE_out<-ggplot(MPEdf4,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# MGE
r11<-as.data.frame.table(querygrain(Sce1in, nodes = "MGE_t2")$MGE_t2)
r12<-as.data.frame.table(querygrain(Sce1in, nodes = "MGE_t3")$MGE_t3)

MGEdf<-cbind(r11, r12[,2]) # combine timesteps to one dataframe
names(MGEdf)<-c("Mortality", "Immediately", "1 year after")
MGEdf2<-melt(MGEdf)
names(MGEdf2)<-c("Mortality", "Timestep", "Probability")

r13<-as.data.frame.table(querygrain(Sce1out, nodes = "MGE_t2")$MGE_t2)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "MGE_t3")$MGE_t3)

MGEdf3<-cbind(r13, r14[,2]) # combine timesteps to one dataframe
names(MGEdf3)<-c("Mortality", "Immediately", "1 year after")
MGEdf4<-melt(MGEdf3)
names(MGEdf4)<-c("Mortality", "Timestep", "Probability")

r15<-as.data.frame.table(querygrain(Sce1near, nodes = "MGE_t2")$MGE_t2)
r16<-as.data.frame.table(querygrain(Sce1near, nodes = "MGE_t3")$MGE_t3)

MGEdf5<-cbind(r15, r16[,2]) # combine timesteps to one dataframe
names(MGEdf5)<-c("Mortality", "Immediately", "1 year after")
MGEdf6<-melt(MGEdf5)
names(MGEdf6)<-c("Mortality", "Timestep", "Probability")

r11n<-as.data.frame.table(querygrain(Sce1near, nodes = "MGE_t2")$MGE_t2)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "MGE_t3")$MGE_t3)

MGEdf7<-cbind(r11n, r12n[,2]) # combine timesteps to one dataframe
names(MGEdf7)<-c("Mortality", "Immediately", "1 year after")
MGEdf8<-melt(MGEdf7)
names(MGEdf8)<-c("Mortality", "Timestep", "Probability")

S1_MGE_near<-ggplot(MGEdf8,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S1_MGE_in<-ggplot(MGEdf2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))

S1_MGE_out<-ggplot(MGEdf4,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# GH
r11<-as.data.frame.table(querygrain(Sce1in, nodes = "GH_t2")$GH_t2)
r12<-as.data.frame.table(querygrain(Sce1in, nodes = "GH_t3")$GH_t3)

GHdf<-cbind(r11, r12[,2]) # combine timesteps to one dataframe
names(GHdf)<-c("Mortality", "Immediately", "1 year after")
GHdf2<-melt(GHdf)
names(GHdf2)<-c("Mortality", "Timestep", "Probability")

r13<-as.data.frame.table(querygrain(Sce1out, nodes = "GH_t2")$GH_t2)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "GH_t3")$GH_t3)

GHdf3<-cbind(r13, r14[,2]) # combine timesteps to one dataframe
names(GHdf3)<-c("Mortality", "Immediately", "1 year after")
GHdf4<-melt(GHdf3)
names(GHdf4)<-c("Mortality", "Timestep", "Probability")

r11n<-as.data.frame.table(querygrain(Sce1near, nodes = "GH_t2")$GH_t2)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "GH_t3")$GH_t3)

GHdf7<-cbind(r11n, r12n[,2]) # combine timesteps to one dataframe
names(GHdf7)<-c("Mortality", "Immediately", "1 year after")
GHdf8<-melt(GHdf7)
names(GHdf8)<-c("Mortality", "Timestep", "Probability")

S1_GH_near<-ggplot(GHdf8,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S1_GH_in<-ggplot(GHdf2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))
S1_GH_out<-ggplot(GHdf4,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# PH

r11<-as.data.frame.table(querygrain(Sce1in, nodes = "PH_t2")$PH_t2)
r12<-as.data.frame.table(querygrain(Sce1in, nodes = "PH_t3")$PH_t3)

PHdf<-cbind(r11, r12[,2]) # combine timesteps to one dataframe
names(PHdf)<-c("Mortality", "Immediately", "1 year after")
PHdf2<-melt(PHdf)
names(PHdf2)<-c("Mortality", "Timestep", "Probability")

r13<-as.data.frame.table(querygrain(Sce1out, nodes = "PH_t2")$PH_t2)
r14<-as.data.frame.table(querygrain(Sce1out, nodes = "PH_t3")$PH_t3)

PHdf3<-cbind(r13, r14[,2]) # combine timesteps to one dataframe
names(PHdf3)<-c("Mortality", "Immediately", "1 year after")
PHdf4<-melt(PHdf3)
names(PHdf4)<-c("Mortality", "Timestep", "Probability")

r15<-as.data.frame.table(querygrain(Sce1near, nodes = "PH_t2")$PH_t2)
r16<-as.data.frame.table(querygrain(Sce1near, nodes = "PH_t3")$PH_t3)

PHdf5<-cbind(r15, r16[,2]) # combine timesteps to one dataframe
names(PHdf5)<-c("Mortality", "Immediately", "1 year after")
PHdf6<-melt(PHdf5)
names(PHdf6)<-c("Mortality", "Timestep", "Probability")

r11n<-as.data.frame.table(querygrain(Sce1near, nodes = "PH_t2")$PH_t2)
r12n<-as.data.frame.table(querygrain(Sce1near, nodes = "PH_t3")$PH_t3)

PHdf7<-cbind(r11n, r12n[,2]) # combine timesteps to one dataframe
names(PHdf7)<-c("Mortality", "Immediately", "1 year after")
PHdf8<-melt(PHdf7)
names(PHdf8)<-c("Mortality", "Timestep", "Probability")

S1_PH_near<-ggplot(PHdf8,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S1_PH_in<-ggplot(PHdf2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))

S1_PH_out<-ggplot(PHdf4,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))


##### SEPARATE PLOTS for ALL FUNCTIONAL GROUPS

## INFAUNA

Infaunaplot<-ggarrange(S1_SMI_in+theme(plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Inside block")+theme(legend.position="none"),
                       S1_SMI_near + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Near-field"),
                       S1_SMI_out + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Far-field"),
                       S1_SSI_in + labs(y="LMM")+
                         theme(legend.position="none"),
                       S1_SSI_near + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                       S1_SSI_out+ theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank(),legend.position="none"),
                       S1_LMM_in + 
                         theme(legend.position="none"),
                       S1_LMM_near+ theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank(),legend.position="none"),
                       S1_LMM_out + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                       S1_LSI_in + labs(y="GH")+
                         theme(legend.position="none"),
                       S1_LSI_near + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                       S1_LSI_out+ theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank(),legend.position="none"),
                       S1_DM_in + 
                         theme(legend.position="none"),
                       S1_DM_near + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="bottom",plot.title = element_text(hjust = 0.5)),
                       S1_DM_out+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),legend.position="none"),
                       S1_SM_in + 
                         theme(legend.position="none"),
                       S1_SM_near + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                       S1_SM_out+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),legend.position="none"),
                       ncol=3,nrow=6, widths = c(1,1,1), font.label = list(size = 10, color = "black", family = NULL),vjust = 0.2,hjust=-0.1,labels = c("Small mobile infauna","","","Small sessile infauna","", "","Large mobile macrofauna"," ","","Large sessile macrofauna","","","Deep meiofauna","","", "Surface meiofauna","",""), common.legend=TRUE)


### Mobile benthos

Mobileplot<-ggarrange(S1_MPE_in+theme(plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Inside block")+theme(legend.position="none"),
                       S1_MPE_near + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Near-field"),
                       S1_MPE_out + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Far-field"),
                       S1_MGE_in + labs(y="LMM")+
                         theme(legend.position="none"),
                       S1_MGE_near + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                       S1_MGE_out+ theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank(),legend.position="none"),
                       S1_PH_in + 
                         theme(legend.position="none"),
                       S1_PH_near+ theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank(),legend.position="none"),
                       S1_PH_out + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                       S1_GH_in + labs(y="GH")+
                         theme(legend.position="none"),
                       S1_GH_near + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                       S1_GH_out+ theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank(),legend.position="none"),
                       
                       ncol=3,nrow=4, widths = c(1,1,1), font.label = list(size = 10, color = "black", family = NULL),vjust = 0.2,hjust=-0.1,labels = c("Mobile predatory epifauna","","","Mobile grazing epifauna","", "","Predatory hyperbenthos"," ","","Grazing hyperbenthos","",""), common.legend=TRUE)


# Sessile fauna
Sessileplot<-ggarrange(S1_SESF_in+theme(plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Inside block")+theme(legend.position="none"),
                       S1_SESF_near + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Near-field"),
                       S1_SESF_out + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Far-field"),
                       S1_SCSF_in + labs(y="LMM")+
                         theme(legend.position="none"),
                       S1_SCSF_near + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                       S1_SCSF_out+ theme(axis.text.y = element_blank(),
                                          axis.ticks.y = element_blank(),
                                          axis.title.y = element_blank(),legend.position="none"),
                       S1_SEFF_in + 
                         theme(legend.position="none"),
                       S1_SEFF_near+ theme(axis.text.y = element_blank(),
                                           axis.ticks.y = element_blank(),
                                           axis.title.y = element_blank(),legend.position="none"),
                       S1_SEFF_out + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                       S1_SCFF_in + labs(y="GH")+
                         theme(legend.position="none"),
                       S1_SCFF_near + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                       S1_SCFF_out+ theme(axis.text.y = element_blank(),
                                          axis.ticks.y = element_blank(),
                                          axis.title.y = element_blank(),legend.position="none"),
                       S1_SSBM_in + labs(y="GH")+
                         theme(legend.position="none"),
                       S1_SSBM_near + 
                         theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                       S1_SSBM_out+ theme(axis.text.y = element_blank(),
                                          axis.ticks.y = element_blank(),
                                          axis.title.y = element_blank(),legend.position="none"),
                       
                       ncol=3,nrow=5, widths = c(1,1,1), font.label = list(size = 10, color = "black", family = NULL),vjust = 0.2,hjust=-0.1,labels = c("Sessile erect suspension feeders","","","Sessile encrusting suspension feeders","", "","Sessile erect filter feeders"," ","","Sessile encrusting filter feeders","","", "Small soft-bodied suspension/filter feeders", "", ""), common.legend=TRUE)

# Combination plot

compS1<-ggarrange(S1_SM_in+theme(plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Inside block")+theme(legend.position="none"),
                  S1_SM_near + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Near-field"),
                  S1_SM_out + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Far-field"),
                  S1_SMI_in + labs(y="LMM")+
                    theme(legend.position="none"),
                  S1_SMI_near + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                  S1_SMI_out+ theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank(),legend.position="none"),
                  S1_LSI_in + 
                    theme(legend.position="none"),
                  S1_LSI_near+ theme(axis.text.y = element_blank(),
                                     axis.ticks.y = element_blank(),
                                     axis.title.y = element_blank(),legend.position="none"),
                  S1_LSI_out + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                  S1_SESF_in + labs(y="GH")+
                    theme(legend.position="none"),
                  S1_SESF_near + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                  S1_SESF_out+ theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank(),legend.position="none"),
                  S1_SSBM_in + 
                    theme(legend.position="none"),
                  S1_SSBM_near + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+theme(legend.position="bottom",plot.title = element_text(hjust = 0.5)),
                  S1_SSBM_out+ theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank(),legend.position="none"),
                  S1_GH_in + 
                    theme(legend.position="none"),
                  S1_GH_near + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                  S1_GH_out+ theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank(),legend.position="none"),
                  ncol=3,nrow=6, widths = c(1,1,1), font.label = list(size = 10, color = "black", family = NULL),vjust = 0.2,hjust=-0.1,labels = c("Surface meiofauna","","","Small mobile infauna","", "","Large sessile macrofauna"," ","","Sessile erect suspension feeders","","","Small soft-bodied suspension /filter feeders","","", "Grazing hyperbenthos","",""), common.legend=TRUE)

png(file="Sessileplot.png",
    width=800, height=1200)
Sessileplot
dev.off()

png(file="C:/Users/kaikkonenl/Documents/CR ERA model/Data/Plots/Infaunaplot.png",
    width=800, height=1200)
Infaunaplot
dev.off()

png(file="C:/Users/kaikkonenl/Documents/CR ERA model/Data/Plots/Mobileplot.png",
    width=800, height=1200)
Mobileplot
dev.off()

png(file="C:/Users/kaikkonenl/Documents/CR ERA model/Data/Plots/Compilation-S1.png",
    width=800, height=1200)
compS1
dev.off()


###############################################
##### SCENARIO 2 INTERMEDIATE DISTURBANCE #######
################################################

# Probability of a node given evidence

# Scenarios 2: Surface collector mining or low-penetration bottom trawling
# Change mining intensity to illustrate changes in pressures

Sce2in <- setEvidence(junction, nodes = c("Mining_intensity","Depth_extraction", "Distance_mining_site", "Plume_release"), states =c("50%", "1-10cm", "Inside block", "At the seafloor"))
Sce2out <- setEvidence(junction, nodes = c("Mining_intensity","Depth_extraction", "Distance_mining_site", "Plume_release"), states =c("50%", "1-10cm", "Outside block", "At the seafloor"))
Sce2near <- setEvidence(junction, nodes = c("Mining_intensity","Depth_extraction", "Distance_mining_site", "Plume_release"), states =c("50%", "1-10cm", "Nearfield", "At the seafloor"))

# Physicochemical parameters

# Probabilities of pressures from mining
r21<-as.data.frame.table(querygrain(Sce2in, nodes = "Sediment_deposition")$Sediment_deposition)
r22<-as.data.frame.table(querygrain(Sce2out, nodes = "Sediment_deposition")$Sediment_deposition)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "Sediment_deposition")$Sediment_deposition)

SDdf<-cbind(r21,r22n[,2],r22[,2]) # combine spatial domains to one dataframe
names(SDdf)<-c("Sediment deposition", "In", "Nearfield","Farfield")
SDdf2<-melt(SDdf)
names(SDdf2)<-c("Sediment_deposition", "Spatial_domain", "Probability")

SDepp2<-ggplot(SDdf2,aes(x = Sediment_deposition, y=Probability, fill = Spatial_domain)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_fill_manual(values=c("#FC4E07", "#E7B800", "#00AFBB"))+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")

r23<-as.data.frame.table(querygrain(Sce2in, nodes = "Suspended_sediment")$Suspended_sediment)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "Suspended_sediment")$Suspended_sediment)
r24n<-as.data.frame.table(querygrain(Sce2near, nodes = "Suspended_sediment")$Suspended_sediment)

SSdf<-cbind(r23, r24n[,2],r24[,2]) # combine spatial domains to one dataframe
names(SSdf)<-c("Suspended sediment","In","Nearfield","Farfield")
SSdf2<-melt(SSdf)
names(SSdf2)<-c("Suspended_sediment", "Spatial_domain", "Probability")

  SSp2<-ggplot(SSdf2,aes(x = Suspended_sediment, y=Probability, fill = Spatial_domain)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_fill_manual(values=c("#FC4E07", "#E7B800", "#00AFBB"))+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")

# Contaminant release
r25<-as.data.frame.table(querygrain(Sce2in, nodes = "Contaminant_release")$Contaminant_release)
r26<-as.data.frame.table(querygrain(Sce2out, nodes = "Contaminant_release")$Contaminant_release)
r26n<-as.data.frame.table(querygrain(Sce2near, nodes = "Contaminant_release")$Contaminant_release)

CRdf<-cbind(r25, r26n[,2],r26[,2]) # combine spatial domains to one dataframe
names(CRdf)<-c("Contaminant_release","In","Nearfield","Farfield")
CRdf2<-melt(CRdf)
names(CRdf2)<-c("Contaminant_release", "Spatial_domain", "Probability")

  CRp2<-ggplot(CRdf2,aes(x =Contaminant_release, y=Probability, fill = Spatial_domain)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_fill_manual(values=c("#FC4E07", "#E7B800", "#00AFBB"))+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")

# Habitat change
r27<-as.data.frame.table(querygrain(Sce2in, nodes = "Sediment_changes")$Sediment_changes)
r18<-as.data.frame.table(querygrain(Sce2out, nodes = "Sediment_changes")$Sediment_changes)
r18n<-as.data.frame.table(querygrain(Sce2near, nodes = "Sediment_changes")$Sediment_changes)

SCdf<-cbind(r27, r18n[,2],r18[,2]) # combine spatial domains to one dataframe
names(SCdf)<-c("Sediment changes","In","Nearfield","Farfield")
SCdf2<-melt(SCdf)
names(SCdf2)<-c("Sediment_changes", "Spatial_domain", "Probability")

  SCp2<-ggplot(SCdf2,aes(x=Sediment_changes, y=Probability, fill = Spatial_domain)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_fill_manual(values=c("#FC4E07", "#E7B800", "#00AFBB"))+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")


pressuresS2<-ggarrange(SSp2+theme(plot.title = element_text(hjust = 0.5))+labs(title="Suspended sediment")+theme(legend.position="none"),
                       SDepp2 +theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank())+labs(y="SMI", title="Sediment deposition")+theme(legend.position="none"),
                       CRp2+theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank())+labs(y="SMI", title="Contaminant release")+theme(legend.position="none"),
                       SCp2+theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank())+labs(y="SMI", title="Sediment changes"),
                       nrow=1,widths = c(1,1,1,1), common.legend=TRUE)

# IMPACTS ON BENTHIC FAUNA

# SMI

r21<-as.data.frame.table(querygrain(Sce2in, nodes = "SMI_t2")$SMI_t2)
r22<-as.data.frame.table(querygrain(Sce2in, nodes = "SMI_t3")$SMI_t3)

SMIdfS2<-cbind(r21, r22[,2]) # combine timesteps to one dataframe
names(SMIdfS2)<-c("Mortality", "Immediately", "1 year after")
SMIdf2S2<-melt(SMIdfS2)
names(SMIdf2S2)<-c("Mortality", "Timestep", "Probability")

r23<-as.data.frame.table(querygrain(Sce2out, nodes = "SMI_t2")$SMI_t2)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "SMI_t3")$SMI_t3)

SMIdf3S2<-cbind(r23, r24[,2]) # combine timesteps to one dataframe
names(SMIdf3S2)<-c("Mortality", "Immediately", "1 year after")
SMIdf4S2<-melt(SMIdf3S2)
names(SMIdf4S2)<-c("Mortality", "Timestep", "Probability")

r25<-as.data.frame.table(querygrain(Sce2near, nodes = "SMI_t2")$SMI_t2)
r26<-as.data.frame.table(querygrain(Sce2near, nodes = "SMI_t3")$SMI_t3)

SMIdf5S2<-cbind(r25, r26[,2]) # combine timesteps to one dataframe
names(SMIdf5S2)<-c("Mortality", "Immediately", "1 year after")
SMIdf6S2<-melt(SMIdf5S2)
names(SMIdf6S2)<-c("Mortality", "Timestep", "Probability")

r21n<-as.data.frame.table(querygrain(Sce2near, nodes = "SMI_t2")$SMI_t2)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "SMI_t3")$SMI_t3)

SMIdf7S2<-cbind(r21n, r22n[,2]) # combine timesteps to one dataframe
names(SMIdf7S2)<-c("Mortality", "Immediately", "1 year after")
SMIdf8S2<-melt(SMIdf7S2)
names(SMIdf8S2)<-c("Mortality", "Timestep", "Probability")

S2_SMI_near<-ggplot(SMIdf8S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S2_SMI_in<-ggplot(SMIdf2S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))

S2_SMI_out<-ggplot(SMIdf4S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# SSI

r21<-as.data.frame.table(querygrain(Sce2in, nodes = "SSI_t2")$SSI_t2)
r22<-as.data.frame.table(querygrain(Sce2in, nodes = "SSI_t3")$SSI_t3)

SSIdfS2<-cbind(r21, r22[,2]) # combine timesteps to one dataframe
names(SSIdfS2)<-c("Mortality", "Immediately", "1 year after")
SSIdf2S2<-melt(SSIdfS2)
names(SSIdf2S2)<-c("Mortality", "Timestep", "Probability")

r23<-as.data.frame.table(querygrain(Sce2out, nodes = "SSI_t2")$SSI_t2)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "SSI_t3")$SSI_t3)

SSIdf3S2<-cbind(r23, r24[,2]) # combine timesteps to one dataframe
names(SSIdf3S2)<-c("Mortality", "Immediately", "1 year after")
SSIdf4S2<-melt(SSIdf3S2)
names(SSIdf4S2)<-c("Mortality", "Timestep", "Probability")

r25<-as.data.frame.table(querygrain(Sce2near, nodes = "SSI_t2")$SSI_t2)
r26<-as.data.frame.table(querygrain(Sce2near, nodes = "SSI_t3")$SSI_t3)

SSIdf5S2<-cbind(r25, r26[,2]) # combine timesteps to one dataframe
names(SSIdf5S2)<-c("Mortality", "Immediately", "1 year after")
SSIdf6S2<-melt(SSIdf5S2)
names(SSIdf6S2)<-c("Mortality", "Timestep", "Probability")

r21n<-as.data.frame.table(querygrain(Sce2near, nodes = "SSI_t2")$SSI_t2)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "SSI_t3")$SSI_t3)

SSIdf7S2<-cbind(r21n, r22n[,2]) # combine timesteps to one dataframe
names(SSIdf7S2)<-c("Mortality", "Immediately", "1 year after")
SSIdf8S2<-melt(SSIdf7S2)
names(SSIdf8S2)<-c("Mortality", "Timestep", "Probability")

S2_SSI_near<-ggplot(SSIdf8S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))
S2_SSI_in<-ggplot(SSIdf2S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))
S2_SSI_out<-ggplot(SSIdf4S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# Surface meiofauna

r21<-as.data.frame.table(querygrain(Sce2in, nodes = "SM_t2")$SM_t2)
r22<-as.data.frame.table(querygrain(Sce2in, nodes = "SM_t3")$SM_t3)

SMdfS2<-cbind(r21, r22[,2]) # combine timesteps to one dataframe
names(SMdfS2)<-c("Mortality", "Immediately", "1 year after")
SMdf2S2<-melt(SMdfS2)
names(SMdf2S2)<-c("Mortality", "Timestep", "Probability")

r23<-as.data.frame.table(querygrain(Sce2out, nodes = "SM_t2")$SM_t2)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "SM_t3")$SM_t3)

SMdf3S2<-cbind(r23, r24[,2]) # combine timesteps to one dataframe
names(SMdf3S2)<-c("Mortality", "Immediately", "1 year after")
SMdf4S2<-melt(SMdf3S2)
names(SMdf4S2)<-c("Mortality", "Timestep", "Probability")

r25<-as.data.frame.table(querygrain(Sce2near, nodes = "SM_t2")$SM_t2)
r26<-as.data.frame.table(querygrain(Sce2near, nodes = "SM_t3")$SM_t3)

SMdf5S2<-cbind(r25, r26[,2]) # combine timesteps to one dataframe
names(SMdf5S2)<-c("Mortality", "Immediately", "1 year after")
SMdf6S2<-melt(SMdf5)
names(SMdf6S2)<-c("Mortality", "Timestep", "Probability")

r21n<-as.data.frame.table(querygrain(Sce2near, nodes = "SM_t2")$SM_t2)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "SM_t3")$SM_t3)

SMdf7S2<-cbind(r21n, r22n[,2]) # combine timesteps to one dataframe
names(SMdf7S2)<-c("Mortality", "Immediately", "1 year after")
SMdf8S2<-melt(SMdf7S2)
names(SMdf8S2)<-c("Mortality", "Timestep", "Probability")

S2_SM_near<-ggplot(SMdf8S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S2_SM_in<-ggplot(SMdf2S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))

S2_SM_out<-ggplot(SMdf4S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# Deep meiofauna

r21<-as.data.frame.table(querygrain(Sce2in, nodes = "DM_t2")$DM_t2)
r22<-as.data.frame.table(querygrain(Sce2in, nodes = "DM_t3")$DM_t3)

DMdfS2<-cbind(r21, r22[,2]) # combine timesteps to one dataframe
names(DMdfS2)<-c("Mortality", "Immediately", "1 year after")
DMdf2S2<-melt(DMdfS2)
names(DMdf2S2)<-c("Mortality", "Timestep", "Probability")

r23<-as.data.frame.table(querygrain(Sce2out, nodes = "DM_t2")$DM_t2)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "DM_t3")$DM_t3)

DMdf3S2<-cbind(r23, r24[,2]) # combine timesteps to one dataframe
names(DMdf3S2)<-c("Mortality", "Immediately", "1 year after")
DMdf4S2<-melt(DMdf3S2)
names(DMdf4S2)<-c("Mortality", "Timestep", "Probability")

r25<-as.data.frame.table(querygrain(Sce2near, nodes = "DM_t2")$DM_t2)
r26<-as.data.frame.table(querygrain(Sce2near, nodes = "DM_t3")$DM_t3)

DMdf5S2<-cbind(r25, r26[,2]) # combine timesteps to one dataframe
names(DMdf5S2)<-c("Mortality", "Immediately", "1 year after")
DMdf6S2<-melt(DMdf5S2)
names(DMdf6S2)<-c("Mortality", "Timestep", "Probability")

r21n<-as.data.frame.table(querygrain(Sce2near, nodes = "DM_t2")$DM_t2)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "DM_t3")$DM_t3)

DMdf7S2<-cbind(r21n, r22n[,2]) # combine timesteps to one dataframe
names(DMdf7S2)<-c("Mortality", "Immediately", "1 year after")
DMdf8S2<-melt(DMdf7S2)
names(DMdf8S2)<-c("Mortality", "Timestep", "Probability")

S2_DM_near<-ggplot(DMdf8S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))
S2_DM_in<-ggplot(DMdf2S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))
S2_DM_out<-ggplot(DMdf4S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# SSBM
r21<-as.data.frame.table(querygrain(Sce2in, nodes = "SSBM_t2")$SSBM_t2)
r22<-as.data.frame.table(querygrain(Sce2in, nodes = "SSBM_t3")$SSBM_t3)

SSBMdfS2<-cbind(r21, r22[,2]) # combine timesteps to one dataframe
names(SSBMdfS2)<-c("Mortality", "Immediately", "1 year after")
SSBMdf2S2<-melt(SSBMdfS2)
names(SSBMdf2S2)<-c("Mortality", "Timestep", "Probability")

r23<-as.data.frame.table(querygrain(Sce2out, nodes = "SSBM_t2")$SSBM_t2)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "SSBM_t3")$SSBM_t3)

SSBMdf3S2<-cbind(r23, r24[,2]) # combine timesteps to one dataframe
names(SSBMdf3S2)<-c("Mortality", "Immediately", "1 year after")
SSBMdf4S2<-melt(SSBMdf3S2)
names(SSBMdf4S2)<-c("Mortality", "Timestep", "Probability")

r21n<-as.data.frame.table(querygrain(Sce2near, nodes = "SSBM_t2")$SSBM_t2)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "SSBM_t3")$SSBM_t3)

SSBMdf7S2<-cbind(r21n, r22n[,2]) # combine timesteps to one dataframe
names(SSBMdf7S2)<-c("Mortality", "Immediately", "1 year after")
SSBMdf8S2<-melt(SSBMdf7S2)
names(SSBMdf8S2)<-c("Mortality", "Timestep", "Probability")

S2_SSBM_near<-ggplot(SSBMdf8S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))

S2_SSBM_in<-ggplot(SSBMdf2S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))

S2_SSBM_out<-ggplot(SSBMdf4S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# SCSF
r21<-as.data.frame.table(querygrain(Sce2in, nodes = "SCSF_t2")$SCSF_t2)
r22<-as.data.frame.table(querygrain(Sce2in, nodes = "SCSF_t3")$SCSF_t3)

SCSFdfS2<-cbind(r21, r22[,2]) # combine timesteps to one dataframe
names(SCSFdfS2)<-c("Mortality", "Immediately", "1 year after")
SCSFdf2S2<-melt(SCSFdfS2)
names(SCSFdf2S2)<-c("Mortality", "Timestep", "Probability")

r23<-as.data.frame.table(querygrain(Sce2out, nodes = "SCSF_t2")$SCSF_t2)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "SCSF_t3")$SCSF_t3)

SCSFdf3S2<-cbind(r23, r24[,2]) # combine timesteps to one dataframe
names(SCSFdf3S2)<-c("Mortality", "Immediately", "1 year after")
SCSFdf4S2<-melt(SCSFdf3S2)
names(SCSFdf4S2)<-c("Mortality", "Timestep", "Probability")

r21n<-as.data.frame.table(querygrain(Sce2near, nodes = "SCSF_t2")$SCSF_t2)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "SCSF_t3")$SCSF_t3)

SCSFdf7S2<-cbind(r21n, r22n[,2]) # combine timesteps to one dataframe
names(SCSFdf7S2)<-c("Mortality", "Immediately", "1 year after")
SCSFdf8S2<-melt(SCSFdf7S2)
names(SCSFdf8S2)<-c("Mortality", "Timestep", "Probability")

S2_SCSF_near<-ggplot(SCSFdf8S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))
S2_SCSF_in<-ggplot(SCSFdf2S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))
S2_SCSF_out<-ggplot(SCSFdf4S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# SEFF
r21<-as.data.frame.table(querygrain(Sce2in, nodes = "SEFF_t2")$SEFF_t2)
r22<-as.data.frame.table(querygrain(Sce2in, nodes = "SEFF_t3")$SEFF_t3)

SEFFdfS2<-cbind(r21, r22[,2]) # combine timesteps to one dataframe
names(SEFFdfS2)<-c("Mortality", "Immediately", "1 year after")
SEFFdf2S2<-melt(SEFFdfS2)
names(SEFFdf2S2)<-c("Mortality", "Timestep", "Probability")

r23<-as.data.frame.table(querygrain(Sce2out, nodes = "SEFF_t2")$SEFF_t2)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "SEFF_t3")$SEFF_t3)

SEFFdf3S2<-cbind(r23, r24[,2]) # combine timesteps to one dataframe
names(SEFFdf3S2)<-c("Mortality", "Immediately", "1 year after")
SEFFdf4S2<-melt(SEFFdf3S2)
names(SEFFdf4S2)<-c("Mortality", "Timestep", "Probability")

r21n<-as.data.frame.table(querygrain(Sce2near, nodes = "SEFF_t2")$SEFF_t2)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "SEFF_t3")$SEFF_t3)

SEFFdf7S2<-cbind(r21n, r22n[,2]) # combine timesteps to one dataframe
names(SEFFdf7S2)<-c("Mortality", "Immediately", "1 year after")
SEFFdf8S2<-melt(SEFFdf7S2)
names(SEFFdf8S2)<-c("Mortality", "Timestep", "Probability")

S2_SEFF_near<-ggplot(SEFFdf8S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))
S2_SEFF_in<-ggplot(SEFFdf2S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))
S2_SEFF_out<-ggplot(SEFFdf4S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# SCFF
r21<-as.data.frame.table(querygrain(Sce2in, nodes = "SCFF_t2")$SCFF_t2)
r22<-as.data.frame.table(querygrain(Sce2in, nodes = "SCFF_t3")$SCFF_t3)

SCFFdfS2<-cbind(r21, r22[,2]) # combine timesteps to one dataframe
names(SCFFdfS2)<-c("Mortality", "Immediately", "1 year after")
SCFFdf2S2<-melt(SCFFdfS2)
names(SCFFdf2S2)<-c("Mortality", "Timestep", "Probability")

r23<-as.data.frame.table(querygrain(Sce2out, nodes = "SCFF_t2")$SCFF_t2)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "SCFF_t3")$SCFF_t3)

SCFFdf3S2<-cbind(r23, r24[,2]) # combine timesteps to one dataframe
names(SCFFdf3S2)<-c("Mortality", "Immediately", "1 year after")
SCFFdf4S2<-melt(SCFFdf3S2)
names(SCFFdf4S2)<-c("Mortality", "Timestep", "Probability")

r21n<-as.data.frame.table(querygrain(Sce2near, nodes = "SCFF_t2")$SCFF_t2)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "SCFF_t3")$SCFF_t3)

SCFFdf7S2<-cbind(r21n, r22n[,2]) # combine timesteps to one dataframe
names(SCFFdf7S2)<-c("Mortality", "Immediately", "1 year after")
SCFFdf8S2<-melt(SCFFdf7S2)
names(SCFFdf8S2)<-c("Mortality", "Timestep", "Probability")

S2_SCFF_near<-ggplot(SCFFdf8S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))
S2_SCFF_in<-ggplot(SCFFdf2S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))
S2_SCFF_out<-ggplot(SCFFdf4S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# SESF
r21<-as.data.frame.table(querygrain(Sce2in, nodes = "SESF_t2")$SESF_t2)
r22<-as.data.frame.table(querygrain(Sce2in, nodes = "SESF_t3")$SESF_t3)

SESFdfS2<-cbind(r21, r22[,2]) # combine timesteps to one dataframe
names(SESFdfS2)<-c("Mortality", "Immediately", "1 year after")
SESFdf2S2<-melt(SESFdfS2)
names(SESFdf2S2)<-c("Mortality", "Timestep", "Probability")

r23<-as.data.frame.table(querygrain(Sce2out, nodes = "SESF_t2")$SESF_t2)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "SESF_t3")$SESF_t3)

SESFdf3S2<-cbind(r23, r24[,2]) # combine timesteps to one dataframe
names(SESFdf3S2)<-c("Mortality", "Immediately", "1 year after")
SESFdf4S2<-melt(SESFdf3S2)
names(SESFdf4S2)<-c("Mortality", "Timestep", "Probability")

r21n<-as.data.frame.table(querygrain(Sce2near, nodes = "SESF_t2")$SESF_t2)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "SESF_t3")$SESF_t3)

SESFdf7S2<-cbind(r21n, r22n[,2]) # combine timesteps to one dataframe
names(SESFdf7S2)<-c("Mortality", "Immediately", "1 year after")
SESFdf8S2<-melt(SESFdf7S2)
names(SESFdf8S2)<-c("Mortality", "Timestep", "Probability")

S2_SESF_near<-ggplot(SESFdf8S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))
S2_SESF_in<-ggplot(SESFdf2S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))
S2_SESF_out<-ggplot(SESFdf4S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# LMM

r21<-as.data.frame.table(querygrain(Sce2in, nodes = "LMM_t2")$LMM_t2)
r22<-as.data.frame.table(querygrain(Sce2in, nodes = "LMM_t3")$LMM_t3)

LMMdfS2<-cbind(r21, r22[,2]) # combine timesteps to one dataframe
names(LMMdfS2)<-c("Mortality", "Immediately", "1 year after")
LMMdf2S2<-melt(LMMdfS2)
names(LMMdf2S2)<-c("Mortality", "Timestep", "Probability")

r23<-as.data.frame.table(querygrain(Sce2out, nodes = "LMM_t2")$LMM_t2)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "LMM_t3")$LMM_t3)

LMMdf3S2<-cbind(r23, r24[,2]) # combine timesteps to one dataframe
names(LMMdf3S2)<-c("Mortality", "Immediately", "1 year after")
LMMdf4S2<-melt(LMMdf3S2)
names(LMMdf4S2)<-c("Mortality", "Timestep", "Probability")

r21n<-as.data.frame.table(querygrain(Sce2near, nodes = "LMM_t2")$LMM_t2)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "LMM_t3")$LMM_t3)

LMMdf7S2<-cbind(r21n, r22n[,2]) # combine timesteps to one dataframe
names(LMMdf7S2)<-c("Mortality", "Immediately", "1 year after")
LMMdf8S2<-melt(LMMdf7S2)
names(LMMdf8S2)<-c("Mortality", "Timestep", "Probability")

S2_LMM_near<-ggplot(LMMdf8S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))
S2_LMM_in<-ggplot(LMMdf2S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))
S2_LMM_out<-ggplot(LMMdf4S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# LSI

r21<-as.data.frame.table(querygrain(Sce2in, nodes = "LSI_t2")$LSI_t2)
r22<-as.data.frame.table(querygrain(Sce2in, nodes = "LSI_t3")$LSI_t3)

LSIdfS2<-cbind(r21, r22[,2]) # combine timesteps to one dataframe
names(LSIdfS2)<-c("Mortality", "Immediately", "1 year after")
LSIdf2S2<-melt(LSIdfS2)
names(LSIdf2S2)<-c("Mortality", "Timestep", "Probability")

r23<-as.data.frame.table(querygrain(Sce2out, nodes = "LSI_t2")$LSI_t2)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "LSI_t3")$LSI_t3)

LSIdf3S2<-cbind(r23, r24[,2]) # combine timesteps to one dataframe
names(LSIdf3S2)<-c("Mortality", "Immediately", "1 year after")
LSIdf4S2<-melt(LSIdf3S2)
names(LSIdf4S2)<-c("Mortality", "Timestep", "Probability")

r25<-as.data.frame.table(querygrain(Sce2near, nodes = "LSI_t2")$LSI_t2)
r26<-as.data.frame.table(querygrain(Sce2near, nodes = "LSI_t3")$LSI_t3)

LSIdf5S2<-cbind(r25, r26[,2]) # combine timesteps to one dataframe
names(LSIdf5S2)<-c("Mortality", "Immediately", "1 year after")
LSIdf6S2<-melt(LSIdf5S2)
names(LSIdf6S2)<-c("Mortality", "Timestep", "Probability")

r21n<-as.data.frame.table(querygrain(Sce2near, nodes = "LSI_t2")$LSI_t2)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "LSI_t3")$LSI_t3)

LSIdf7S2<-cbind(r21n, r22n[,2]) # combine timesteps to one dataframe
names(LSIdf7S2)<-c("Mortality", "Immediately", "1 year after")
LSIdf8S2<-melt(LSIdf7S2)
names(LSIdf8S2)<-c("Mortality", "Timestep", "Probability")

S2_LSI_near<-ggplot(LSIdf8S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))
S2_LSI_in<-ggplot(LSIdf2S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))
S2_LSI_out<-ggplot(LSIdf4S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# MPE

r21<-as.data.frame.table(querygrain(Sce2in, nodes = "MPE_t2")$MPE_t2)
r22<-as.data.frame.table(querygrain(Sce2in, nodes = "MPE_t3")$MPE_t3)

MPEdfS2<-cbind(r21, r22[,2]) # combine timesteps to one dataframe
names(MPEdfS2)<-c("Mortality", "Immediately", "1 year after")
MPEdf2S2<-melt(MPEdfS2)
names(MPEdf2S2)<-c("Mortality", "Timestep", "Probability")

r23<-as.data.frame.table(querygrain(Sce2out, nodes = "MPE_t2")$MPE_t2)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "MPE_t3")$MPE_t3)

MPEdf3S2<-cbind(r23, r24[,2]) # combine timesteps to one dataframe
names(MPEdf3S2)<-c("Mortality", "Immediately", "1 year after")
MPEdf4S2<-melt(MPEdf3S2)
names(MPEdf4S2)<-c("Mortality", "Timestep", "Probability")

r25<-as.data.frame.table(querygrain(Sce2near, nodes = "MPE_t2")$MPE_t2)
r26<-as.data.frame.table(querygrain(Sce2near, nodes = "MPE_t3")$MPE_t3)

MPEdf5S2<-cbind(r25, r26[,2]) # combine timesteps to one dataframe
names(MPEdf5S2)<-c("Mortality", "Immediately", "1 year after")
MPEdf6S2<-melt(MPEdf5S2)
names(MPEdf6S2)<-c("Mortality", "Timestep", "Probability")

r21n<-as.data.frame.table(querygrain(Sce2near, nodes = "MPE_t2")$MPE_t2)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "MPE_t3")$MPE_t3)

MPEdf7S2<-cbind(r21n, r22n[,2]) # combine timesteps to one dataframe
names(MPEdf7S2)<-c("Mortality", "Immediately", "1 year after")
MPEdf8S2<-melt(MPEdf7S2)
names(MPEdf8S2)<-c("Mortality", "Timestep", "Probability")

S2_MPE_near<-ggplot(MPEdf8S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))
S2_MPE_in<-ggplot(MPEdf2S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))
S2_MPE_out<-ggplot(MPEdf4S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# MGE
r21<-as.data.frame.table(querygrain(Sce2in, nodes = "MGE_t2")$MGE_t2)
r22<-as.data.frame.table(querygrain(Sce2in, nodes = "MGE_t3")$MGE_t3)

MGEdfS2<-cbind(r21, r22[,2]) # combine timesteps to one dataframe
names(MGEdfS2)<-c("Mortality", "Immediately", "1 year after")
MGEdf2S2<-melt(MGEdfS2)
names(MGEdf2S2)<-c("Mortality", "Timestep", "Probability")

r23<-as.data.frame.table(querygrain(Sce2out, nodes = "MGE_t2")$MGE_t2)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "MGE_t3")$MGE_t3)

MGEdf3S2<-cbind(r23, r24[,2]) # combine timesteps to one dataframe
names(MGEdf3S2)<-c("Mortality", "Immediately", "1 year after")
MGEdf4S2<-melt(MGEdf3S2)
names(MGEdf4S2)<-c("Mortality", "Timestep", "Probability")

r25<-as.data.frame.table(querygrain(Sce2near, nodes = "MGE_t2")$MGE_t2)
r26<-as.data.frame.table(querygrain(Sce2near, nodes = "MGE_t3")$MGE_t3)

MGEdf5S2<-cbind(r25, r26[,2]) # combine timesteps to one dataframe
names(MGEdf5S2)<-c("Mortality", "Immediately", "1 year after")
MGEdf6S2<-melt(MGEdf5S2)
names(MGEdf6S2)<-c("Mortality", "Timestep", "Probability")

r21n<-as.data.frame.table(querygrain(Sce2near, nodes = "MGE_t2")$MGE_t2)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "MGE_t3")$MGE_t3)

MGEdf7S2<-cbind(r21n, r22n[,2]) # combine timesteps to one dataframe
names(MGEdf7S2)<-c("Mortality", "Immediately", "1 year after")
MGEdf8S2<-melt(MGEdf7S2)
names(MGEdf8S2)<-c("Mortality", "Timestep", "Probability")

S2_MGE_near<-ggplot(MGEdf8S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))
S2_MGE_in<-ggplot(MGEdf2S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))
S2_MGE_out<-ggplot(MGEdf4S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# GH
r21<-as.data.frame.table(querygrain(Sce2in, nodes = "GH_t2")$GH_t2)
r22<-as.data.frame.table(querygrain(Sce2in, nodes = "GH_t3")$GH_t3)

GHdfS2<-cbind(r21, r22[,2]) # combine timesteps to one dataframe
names(GHdfS2)<-c("Mortality", "Immediately", "1 year after")
GHdf2S2<-melt(GHdfS2)
names(GHdf2S2)<-c("Mortality", "Timestep", "Probability")

r23<-as.data.frame.table(querygrain(Sce2out, nodes = "GH_t2")$GH_t2)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "GH_t3")$GH_t3)

GHdf3S2<-cbind(r23, r24[,2]) # combine timesteps to one dataframe
names(GHdf3S2)<-c("Mortality", "Immediately", "1 year after")
GHdf4S2<-melt(GHdf3S2)
names(GHdf4S2)<-c("Mortality", "Timestep", "Probability")

r21n<-as.data.frame.table(querygrain(Sce2near, nodes = "GH_t2")$GH_t2)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "GH_t3")$GH_t3)

GHdf7S2<-cbind(r21n, r22n[,2]) # combine timesteps to one dataframe
names(GHdf7S2)<-c("Mortality", "Immediately", "1 year after")
GHdf8S2<-melt(GHdf7S2)
names(GHdf8S2)<-c("Mortality", "Timestep", "Probability")

S2_GH_near<-ggplot(GHdf8S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))
S2_GH_in<-ggplot(GHdf2S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))
S2_GH_out<-ggplot(GHdf4S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))

# PH

r21<-as.data.frame.table(querygrain(Sce2in, nodes = "PH_t2")$PH_t2)
r22<-as.data.frame.table(querygrain(Sce2in, nodes = "PH_t3")$PH_t3)

PHdfS2<-cbind(r21, r22[,2]) # combine timesteps to one dataframe
names(PHdfS2)<-c("Mortality", "Immediately", "1 year after")
PHdf2S2<-melt(PHdfS2)
names(PHdf2S2)<-c("Mortality", "Timestep", "Probability")

r23<-as.data.frame.table(querygrain(Sce2out, nodes = "PH_t2")$PH_t2)
r24<-as.data.frame.table(querygrain(Sce2out, nodes = "PH_t3")$PH_t3)

PHdf3S2<-cbind(r23, r24[,2]) # combine timesteps to one dataframe
names(PHdf3S2)<-c("Mortality", "Immediately", "1 year after")
PHdf4S2<-melt(PHdf3S2)
names(PHdf4S2)<-c("Mortality", "Timestep", "Probability")

r25<-as.data.frame.table(querygrain(Sce2near, nodes = "PH_t2")$PH_t2)
r26<-as.data.frame.table(querygrain(Sce2near, nodes = "PH_t3")$PH_t3)

PHdf5S2<-cbind(r25, r26[,2]) # combine timesteps to one dataframe
names(PHdf5S2)<-c("Mortality", "Immediately", "1 year after")
PHdf6S2<-melt(PHdf5S2)
names(PHdf6S2)<-c("Mortality", "Timestep", "Probability")

r21n<-as.data.frame.table(querygrain(Sce2near, nodes = "PH_t2")$PH_t2)
r22n<-as.data.frame.table(querygrain(Sce2near, nodes = "PH_t3")$PH_t3)

PHdf7S2<-cbind(r21n, r22n[,2]) # combine timesteps to one dataframe
names(PHdf7S2)<-c("Mortality", "Immediately", "1 year after")
PHdf8S2<-melt(PHdf7S2)
names(PHdf8S2)<-c("Mortality", "Timestep", "Probability")

S2_PH_near<-ggplot(PHdf8S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("Decrease in abundance")+scale_fill_manual(values = c("#E7B800", "#EBCC2A"))
S2_PH_in<-ggplot(PHdf2S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#FC4E07", "#FD9A72"))
S2_PH_out<-ggplot(PHdf4S2,aes(x = Mortality, y=Probability, fill = Timestep)) + geom_bar(stat="identity",position = "dodge")+theme_classic()+ scale_y_continuous(name="Probability", limits=c(0, 1))+xlab("")+scale_fill_manual(values = c("#00AFBB","#78B7C5"))


##### SEPARATE PLOTS for ALL FUNCTIONAL GROUPS

## INFAUNA

InfaunaplotS2<-ggarrange(S2_SMI_in+theme(plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Inside block")+theme(legend.position="none"),
                         S2_SMI_near + 
                           theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Near-field"),
                         S2_SMI_out + 
                           theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Far-field"),
                         S2_SSI_in + labs(y="LMM")+
                           theme(legend.position="none"),
                         S2_SSI_near + 
                           theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                         S2_SSI_out+ theme(axis.text.y = element_blank(),
                                           axis.ticks.y = element_blank(),
                                           axis.title.y = element_blank(),legend.position="none"),
                         S2_LMM_in + 
                           theme(legend.position="none"),
                         S2_LMM_near+ theme(axis.text.y = element_blank(),
                                            axis.ticks.y = element_blank(),
                                            axis.title.y = element_blank(),legend.position="none"),
                         S2_LMM_out + 
                           theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                         S2_LSI_in + labs(y="GH")+
                           theme(legend.position="none"),
                         S2_LSI_near + 
                           theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                         S2_LSI_out+ theme(axis.text.y = element_blank(),
                                           axis.ticks.y = element_blank(),
                                           axis.title.y = element_blank(),legend.position="none"),
                         S2_SM_in + 
                           theme(legend.position="none"),
                         S2_SM_near + 
                           theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                         S2_SM_out+ theme(axis.text.y = element_blank(),
                                          axis.ticks.y = element_blank(),
                                          axis.title.y = element_blank(),legend.position="none"),
                         S2_DM_in + 
                           theme(legend.position="none"),
                         S2_DM_near + 
                           theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())+theme(legend.position="bottom",plot.title = element_text(hjust = 0.5)),
                         S2_DM_out+ theme(axis.text.y = element_blank(),
                                          axis.ticks.y = element_blank(),
                                          axis.title.y = element_blank(),legend.position="none"),
                         ncol=3,nrow=6, widths = c(1,1,1), font.label = list(size = 10, color = "black", family = NULL),vjust = 0.2,hjust=-0.1,labels = c("Small mobile infauna","","","Small sessile infauna","", "","Large mobile macrofauna"," ","","Large sessile macrofauna","","","Deep meiofauna","","", "Surface meiofauna","",""), common.legend=TRUE)

### Mobile benthos

MobileplotS2<-ggarrange(S2_MPE_in+theme(plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Inside block")+theme(legend.position="none"),
                        S2_MPE_near + 
                          theme(axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Near-field"),
                        S2_MPE_out + 
                          theme(axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Far-field"),
                        S2_MGE_in + labs(y="LMM")+
                          theme(legend.position="none"),
                        S2_MGE_near + 
                          theme(axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                        S2_MGE_out+ theme(axis.text.y = element_blank(),
                                          axis.ticks.y = element_blank(),
                                          axis.title.y = element_blank(),legend.position="none"),
                        S2_PH_in + 
                          theme(legend.position="none"),
                        S2_PH_near+ theme(axis.text.y = element_blank(),
                                          axis.ticks.y = element_blank(),
                                          axis.title.y = element_blank(),legend.position="none"),
                        S2_PH_out + 
                          theme(axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                        S2_GH_in + labs(y="GH")+
                          theme(legend.position="none"),
                        S2_GH_near + 
                          theme(axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                        S2_GH_out+ theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank(),legend.position="none"),
                        
                        ncol=3,nrow=4, widths = c(1,1,1), font.label = list(size = 10, color = "black", family = NULL),vjust = 0.2,hjust=-0.1,labels = c("Mobile predatory epifauna","","","Mobile grazing epifauna","", "","Predatory hyperbenthos"," ","","Grazing hyperbenthos","",""), common.legend=TRUE)


# Sessile fauna
SessileplotS2<-ggarrange(S2_SESF_in+theme(plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Inside block")+theme(legend.position="none"),
                         S2_SESF_near + 
                           theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Near-field"),
                         S2_SESF_out + 
                           theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Far-field"),
                         S2_SCSF_in + labs(y="LMM")+
                           theme(legend.position="none"),
                         S2_SCSF_near + 
                           theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                         S2_SCSF_out+ theme(axis.text.y = element_blank(),
                                            axis.ticks.y = element_blank(),
                                            axis.title.y = element_blank(),legend.position="none"),
                         S2_SEFF_in + 
                           theme(legend.position="none"),
                         S2_SEFF_near+ theme(axis.text.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.title.y = element_blank(),legend.position="none"),
                         S2_SEFF_out + 
                           theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                         S2_SCFF_in + labs(y="GH")+
                           theme(legend.position="none"),
                         S2_SCFF_near + 
                           theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                         S2_SCFF_out+ theme(axis.text.y = element_blank(),
                                            axis.ticks.y = element_blank(),
                                            axis.title.y = element_blank(),legend.position="none"),
                         S2_SSBM_in + labs(y="GH")+
                           theme(legend.position="none"),
                         S2_SSBM_near + 
                           theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                         S2_SSBM_out+ theme(axis.text.y = element_blank(),
                                            axis.ticks.y = element_blank(),
                                            axis.title.y = element_blank(),legend.position="none"),
                         
                         ncol=3,nrow=5, widths = c(1,1,1), font.label = list(size = 10, color = "black", family = NULL),vjust = 0.2,hjust=-0.1,labels = c("Sessile erect suspension feeders","","","Sessile encrusting suspension feeders","", "","Sessile erect filter feeders"," ","","Sessile encrusting filter feeders","","", "Small soft-bodied suspension/filter feeders", "", ""), common.legend=TRUE)

# Combination plot

compS2<-ggarrange(S2_SM_in+theme(plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Inside block")+theme(legend.position="none"),
                  S2_SM_near + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Near-field"),
                  S2_SM_out + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5))+labs(y="SMI", title="Far-field"),
                  S2_SMI_in + labs(y="LMM")+
                    theme(legend.position="none"),
                  S2_SMI_near + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                  S2_SMI_out+ theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank(),legend.position="none"),
                  S2_LSI_in + 
                    theme(legend.position="none"),
                  S2_LSI_near+ theme(axis.text.y = element_blank(),
                                     axis.ticks.y = element_blank(),
                                     axis.title.y = element_blank(),legend.position="none"),
                  S2_LSI_out + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                  S2_SESF_in + labs(y="GH")+
                    theme(legend.position="none"),
                  S2_SESF_near + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                  S2_SESF_out+ theme(axis.text.y = element_blank(),
                                     axis.ticks.y = element_blank(),
                                     axis.title.y = element_blank(),legend.position="none"),
                  S2_SSBM_in + 
                    theme(legend.position="none"),
                  S2_SSBM_near + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+theme(legend.position="bottom",plot.title = element_text(hjust = 0.5)),
                  S2_SSBM_out+ theme(axis.text.y = element_blank(),
                                     axis.ticks.y = element_blank(),
                                     axis.title.y = element_blank(),legend.position="none"),
                  S2_GH_in + 
                    theme(legend.position="none"),
                  S2_GH_near + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank())+theme(legend.position="none",plot.title = element_text(hjust = 0.5)),
                  S2_GH_out+ theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank(),legend.position="none"),
                  ncol=3,nrow=6, widths = c(1,1,1), font.label = list(size = 10, color = "black", family = NULL),vjust = 0.2,hjust=-0.1,labels = c("Surface meiofauna","","","Small mobile infauna","", "","Large sessile macrofauna"," ","","Sessile erect suspension feeders","","","Small soft-bodied suspension /filter feeders","","", "Grazing hyperbenthos","",""), common.legend=TRUE)



