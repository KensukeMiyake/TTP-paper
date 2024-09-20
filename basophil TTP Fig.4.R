library(grandR)
library(ggplot2)
library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)
library(minpack.lm)
suppressPackageStartupMessages({
  library(grandR)
  library(ggplot2)
  library(patchwork)
})

set.seed(1234)

#================calculate the normalized count of labeled transcripts==================
#generate the grandR object
BMBA <- ReadGRAND("BMBA.tsv.gz",
                  design=c(Design$dur.4sU,Design$Condition,Design$Replicate)) #load GRAND-SLAM output table
Coldata(BMBA) #metadata
BMBA
BMBA <- FilterGenes(BMBA,minval = 20, mincol = 6) 
BMBA
PlotPCA(BMBA,aest=aes(color=duration.4sU.original,shape=Condition))

BMBA<-Normalize(BMBA)
BMBA

#the normalized count of labeled transcripts
norm = BMBA$data$norm
ntr = BMBA$data$ntr
label.norm = norm*ntr
rownames(label.norm) = BMBA$gene.info$Symbol
label.norm = as.data.frame(label.norm)
label.norm$symbol = rownames(label.norm)

View(label.norm)

#the average of normalized count of labeled transcripts at each time point in each genotype
label.norm$WT_0 = apply(label.norm[, c(1:3)], 1, mean)
label.norm$KO_0 = apply(label.norm[, c(4:6)], 1, mean)
label.norm$WT_1 = apply(label.norm[, c(7:9)], 1, mean)
label.norm$KO_1 = apply(label.norm[, c(10:12)], 1, mean)
label.norm$WT_2 = apply(label.norm[, c(13:15)], 1, mean)
label.norm$KO_2 = apply(label.norm[, c(16:18)], 1, mean)
label.norm$WT_4 = apply(label.norm[, c(19:21)], 1, mean)
label.norm$KO_4 = apply(label.norm[, c(22:24)], 1, mean)

View(label.norm)

#normalize to the average of the labeled count at 0h
label.norm$WT_0_A = label.norm[, 1]/label.norm$WT_0
label.norm$WT_0_B = label.norm[, 2]/label.norm$WT_0
label.norm$WT_0_C = label.norm[, 3]/label.norm$WT_0
label.norm$KO_0_A = label.norm[, 4]/label.norm$KO_0
label.norm$KO_0_B = label.norm[, 5]/label.norm$KO_0
label.norm$KO_0_C = label.norm[, 6]/label.norm$KO_0
label.norm$WT_1_A = label.norm[, 7]/label.norm$WT_0
label.norm$WT_1_B = label.norm[, 8]/label.norm$WT_0
label.norm$WT_1_C = label.norm[, 9]/label.norm$WT_0
label.norm$KO_1_A = label.norm[, 10]/label.norm$KO_0
label.norm$KO_1_B = label.norm[, 11]/label.norm$KO_0
label.norm$KO_1_C = label.norm[, 12]/label.norm$KO_0
label.norm$WT_2_A = label.norm[, 13]/label.norm$WT_0
label.norm$WT_2_B = label.norm[, 14]/label.norm$WT_0
label.norm$WT_2_C = label.norm[, 15]/label.norm$WT_0
label.norm$KO_2_A = label.norm[, 16]/label.norm$KO_0
label.norm$KO_2_B = label.norm[, 17]/label.norm$KO_0
label.norm$KO_2_C = label.norm[, 18]/label.norm$KO_0
label.norm$WT_4_A = label.norm[, 19]/label.norm$WT_0
label.norm$WT_4_B = label.norm[, 20]/label.norm$WT_0
label.norm$WT_4_C = label.norm[, 21]/label.norm$WT_0
label.norm$KO_4_A = label.norm[, 22]/label.norm$KO_0
label.norm$KO_4_B = label.norm[, 23]/label.norm$KO_0
label.norm$KO_4_C = label.norm[, 24]/label.norm$KO_0

View(label.norm)
names(label.norm)

#================calculate the half-lives of transcripts==================
#remove transcripts whose labeled count are 0 at 0h
dim(label.norm)
label.norm = filter(label.norm, WT_0 != 0)
label.norm = filter(label.norm, KO_0 != 0)
dim(label.norm)
View(label.norm)

label_WT = label.norm[, c(c(34:36, 40:42, 46:48, 52:54))]
label_KO = label.norm[, c(c(37:39, 43:45, 49:51, 55:57))]

#generate the function for calculating half-lives of each transcript
#REF:https://github.com/t-neumann/slamdunk/issues/70
half_calc <- function(y) {
  replicate = c("A","B","C","A","B","C","A","B","C","A","B","C")
  timepoints = c(0,0,0,1,1,1,2,2,2,4,4,4)
  df = data.frame(replicate, timepoints, y)
  model = nlsLM(y~Plat + (y0 - Plat) * exp(-k * (timepoints)),
                data = df,
                start=list(
                  y0= 1,
                  Plat = 0,
                  k= 0.5),
                
                upper = c(1,0,Inf),
                lower = c(1,0,0),
                control = nls.lm.control(maxiter = 1000),
                na.action = na.omit)
  a = summary(model)
  return(a$parameters[3,1])
}

#calculate the variable k of each transcript
#applying a single exponential decay model
#transcripts with errors when applying this model were removed
half_WT = c()
for(i in 1:nrow(label_WT)){
  tryCatch(
    {half_WT[i] = half_calc(as.numeric(label_WT[i,c(1:12)]))}
    , error = function(e){half_WT[i] = NA}
  ) 
}

length(half_WT)

half_KO = c()
for(i in 1:nrow(label_KO)){
  tryCatch(
    {half_KO[i] = half_calc(as.numeric(label_KO[i,c(1:12)]))}
    , error = function(e){half_KO[i] = NA}
  ) 
}

length(half_KO)

symbol = as.vector(rownames(label.norm))
length(symbol)

#calculate the half-lives from the variable k
#half-life = loge(2)/k
half_time = data.frame(symbol, half_WT, half_KO)
half_time$time_WT = log(2)/half_time$half_WT
half_time$time_KO = log(2)/half_time$half_KO
half_time$ratio_KO_WT = half_time$time_KO/half_time$time_WT
rownames(half_time) = half_time$symbol
dim(half_time)

View(half_time)

#remove transcripts with errors when applying a single exponential decay model 
half_time <- na.omit(half_time)
dim(half_time)

View(half_time)

#half-life (hour to minute)
half_time$min_WT = half_time$time_WT*60
half_time$min_KO = half_time$time_KO*60

#log-transformation of half-life (min)
half_time$logminWT = log2(half_time$min_WT)
half_time$logminKO = log2(half_time$min_KO)


#================compare half-lives of WT and TTP-KO basophils (Fig.4C)==================
interested = c("Il4", "Il13", "Cxcl2", "Ccl3", "Ptgs2")

#annotate genes whose stability was changed by TTP-deficiency
half_time$stability <- "NA"
half_time$stability[half_time$ratio_KO_WT > 1.1] <- "UP" #cut-off value = 1.1
half_time$stability[half_time$ratio_KO_WT < 1/1.1] <- "DOWN"
half_time$stability[half_time$symbol %in% interested] <- "INTEREST"
half_time$stability = factor(half_time$stability, levels = c("NA", "UP", "DOWN", "INTEREST"))

#plot (Fig.4C)
p=ggplot(data=half_time, aes(x=logminKO, y=logminWT, col=stability, size=stability))+
  geom_point() +
  scale_color_manual(values=c("gray", "#fed0e0", "gray40", "red")) +
  scale_size_manual(values=c(0.5,0.5,0.5,2))+
  geom_abline(intercept = -log2(1.1), slope = 1, linewidth = 0.5, linetype = 5)+
  geom_abline(intercept = log2(1.1), slope = 1, linewidth = 0.5, linetype = 5)+
  xlim(1.8,10)+
  ylim(1.8,10)+
  coord_fixed()+
  theme_classic()+
  theme(legend.position = "none")

p = Seurat::LabelPoints(plot = p, points = interested, repel = TRUE)
p

#================enrichment analysis (Fig.4C)==================
#pick up transcripts whose stability was upregulated by TTP-deficiency
upgenes <- half_time$symbol[half_time$ratio_KO_WT >1.1] #cut-off value = 1.1

#convert gene symbol to ENTREZID
upgenes <- bitr(upgenes, 
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)

#GO BP enrichment analysis
ego_up <- enrichGO(gene = upgenes$ENTREZID,
                   OrgDb         = org.Mm.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05, 
                   readable      = TRUE)
head(as.data.frame(ego_up))
View(as.data.frame(ego_up))

terms = c("regulation of inflammatory response", "production of molecular mediator of immune response", "cytokine production involved in immune response", "cytoplasmic translation", "ribosome biogenesis")

#plot (Fig.4D)
clusterProfiler::dotplot(ego_up, showCategory = terms, font.size=8)+ geom_count() + scale_size_area()