---
title: "Filtering_for_cell_types"
author: "Laura Wolbeck"
date: "2023-09-25"
output: html_document
---
 
The dataset is filtered for ambient RNA, low UMI count, high mitochondrial gene fraction and doublets.
Now follows filtering based on marker genes and cell types. 

# Setup
```{r}
#Load helper functions
source("/people/nrq364/embryonic_PolyIC/scRNA_helper_new.R")

library(CRMetrics)
library(magrittr)
library(dplyr)
library(conos)
library(pagoda2)
library(qs)
library(ggplot2)
library(ggrastr)
library(RColorBrewer)
library(Seurat)
```



#1.  Read data
read in count matrix from previous notebook
```{r}
cms <- qread("cms_filtered_cellbender.qs")
```

#2. Conos #1 with all cells
create initial embedding with all cells 

vector with sample names
```{r, eval=F}
names =c("control_1_",
           "control_2_",
           "control_3_",
         "control_4_",
         "control_5_",
         "control_6_",
          "treated_1_",
           "treated_2_",
           "treated_3_",
           "treated_4_",
           "treated_5_",
           "treated_6_")
```

using a helper function that combines normalization, integration of samples, embedding into graph and clustering
```{r, eval=FALSE}
con <- quickConos(cms,
                  names,
                  n.cores.p2=10,
                  n.cores.con=20, get.tsne = F, alignment.strength=0.2)

con <- con$con
```

## additional clusters
Rerun leiden clustering to find additional clusters
```{r, eval=FALSE}
con$findCommunities(method = leiden.community, resolution = 4, min.group.size = 50)
```

plot graph with new clusters
```{r, eval=FALSE}
con$plotGraph(groups = con$clusters$leiden$groups %>% factor)
```

table of number of nuclei per cluster
```{r, eval=FALSE}
table(con$clusters$leiden$groups %>% factor)
```

save factor of clusters and original/initial conos object
```{r}
leiden25 <- con$clusters$leiden$groups %>% factor
qsave(leiden25, "leiden25.qs")
qsave(con, "con_org.qs")
```

## annotation of initial dataset
see functions used for annotation of dataset in the notebook "Annotation.Rmd"
```{r, fig.height=10}
con_org$plotGraph(groups=anno)
```

### Figure 19A
export UMAP representation without labeling of clusters
```{r}
plot <- con_org$plotGraph(mark.groups=F) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
rasterize(plot, layers='Point', dpi=1000)
ggsave("UMAP_con_org_without_labels.pdf", width=14, height=8.6 )
```

##Figure 19B: plot marker genes on UMAP
only of the cell types that are not retained in final UMAP, i.e. non-telencephalic and non-neuronal cells

```{r}
genes <- c("Grhl2","Epcam", "Hba-x", "Hbb-y", "Runx1", "Fyb","En1", "Fgf8", "Alx4", "Pdgfra", "Sox10", "Col9a1", "Dkk3", "Cd34","Pecam1")
```

```{r}
plots <- lapply(genes, function(gene) {
  data <- con_org$plotGraph(gene=gene)
  data <- data$data
  ggplot(data, aes(x = x, y =y, color= Color))+ geom_point(size=0.1, alpha=0.5)+ theme_light()+ theme(panel.grid.major = element_blank(), panel.grid.minor =    element_blank())+ xlab(NULL)+ylab(NULL) +scale_color_gradient(low="#dbd9d9",high="#FF0000")+ theme(axis.ticks = element_blank(), axis.text =element_blank(),legend.position="none",plot.title = element_text(size = 8))+labs(color = gene)+ggtitle(gene)+theme( plot.margin = margin(t = -1,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))
})

```

plot expression of marker genes on UMAP and export
```{r, fig.height=8}
plots <- lapply(plots, function(x) {rasterize(x,layers='Point', dpi=1000)})
cowplot::plot_grid(plotlist = plots, ncol=3)
ggsave("UMAPs_con_org_marker_genes.pdf", width=5, height = 7)
```

#3. Filter for telencephalic, neurogenic and microglial genes
vector with names of telencephalic genes, we want to keep all cells expressing at least one of these genes
```{r}
genes <- c("Emx1", "Emx2", "Foxg1", "Six3", "Eomes", "Neurog1", "Neurog2","Slc17a6", "Slc17a7", "Gad1","Gad2", "Il6ra", "Spi1", "Mrc1", "Cd86")
```

check number of cells per sample before filtering
```{r, eval=FALSE}
sapply(cms, function(d) dim(d)[2])
sum(sapply(cms, function(d) dim(d)[2]))
```

Extract all cell names where at least one of the genes in the vector "genes" has a count > 0
```{r}
cell_names <- lapply(cms, function(cm) {
  gene_counts <- sapply(genes, function(gene) cm[gene, ])
  colnames(cm)[rowSums(gene_counts) > 0]
})

cell_names <- unlist(cell_names)
```

Subset the count matrices to include only the cells with counts > 0 for at least one of the genes in the vector
```{r}
cms_genes <- lapply(cms, function(cm) cm[, colnames(cm) %in% cell_names])
```

check number of cells per sample after filtering
```{r, eval=FALSE}
sapply(cms_genes, function(d) dim(d)[2])
sum(sapply(cms_genes, function(d) dim(d)[2]))
```
```{r}
qsave(cms_genes, "cms_genes.qs")
```

#4. Conos #2: selected genes
new embedding on remaining cells
```{r}
cms <- qs::qread("cms_genes.qs")
```

```{r, eval=F}
#vector with sample names 
names =c("control_1_",
           "control_2_",
           "control_3_",
         "control_4_",
         "control_5_",
         "control_6_",
          "treated_1_",
           "treated_2_",
           "treated_3_",
           "treated_4_",
           "treated_5_",
           "treated_6_")
```
```{r, eval=FALSE}
con <- quickConos(cms,
                  names,
                  n.cores.p2=10,
                  n.cores.con=20, get.tsne = F, alignment.strength=0.2)

con <- con$con

```
```{r}
qsave(con, "con_genes.qs")
con_genes <- qread("con_genes.qs")
```

inspect UMAP and clustering
```{r}
con_genes$plotGraph() 
```

Colour graph per sample
```{r, echo=F, fig.height=15}
con_genes$plotGraph(groups=con_genes$getDatasetPerCell(), size=1, mark.groups=F, legend.position="right", show.legend = T, title="") + dotSize(2)
``` 

annotate dataset with help of previous annotation of initial dataset
annotation is named "anno_genes"
dataset still includes cell types that need to be removed

```{r, echo=F}
con_genes$plotGraph(groups=anno_genes) 
```

##Figure 20: plot intermediate UMAP and export
```{r}
plot <- con_genes$plotGraph(mark.groups=F) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
rasterize(plot, layers='Point', dpi=1000)
ggsave("UMAP_con_genes_without_labels.pdf", width=14, height=8.6 )
```

## remove unwanted cell types
get cell names to be removed
```{r,eval=F}
remove_cluster <- names(anno_genes[anno_genes %in% c("Vascular_endothelial_cells", "Nasal_process_cells", "Epithelial_cells"  , "Optic_cup_and_retinal_progenitor_cells","Midbrain_hinbrain_progenitors", "Erythrocyte_progenitors")])
```

```{r}
length(remove_cluster)
```

```{r, eval=FALSE}
sapply(cms, function(d) dim(d)[2])
sum(sapply(cms, function(d) dim(d)[2]) %>% setNames(names))
```

```{r, eval=FALSE}
#remove "remove_cluster" from cms
cms %<>% lapply(function(s) s[,!colnames(s) %in% remove_cluster])
```
```{r, eval=FALSE}
sapply(cms, function(d) dim(d)[2])
sum(sapply(cms, function(d) dim(d)[2]) )
```
```{r, eval=FALSE, echo=FALSE}
qsave(cms, "cms_final.qs") 
```

#5. Conos #3: final embedding 
```{r}
cms <- qread("cms_final.qs")
```

```{r, eval=F}
#vector with sample names 
names =c("control_1_",
           "control_2_",
           "control_3_",
         "control_4_",
         "control_5_",
         "control_6_",
          "treated_1_",
           "treated_2_",
           "treated_3_",
           "treated_4_",
           "treated_5_",
           "treated_6_")
```
```{r, eval=FALSE}
con <- quickConos(cms,
                  names,
                  n.cores.p2=10,
                  n.cores.con=20, get.tsne = F, alignment.strength=0.2)

con <- con$con
```

find additional clusters and save conos object
```{r, eval=FALSE}
con$findCommunities(method = leiden.community, resolution = 3, min.group.size = 15)
qsave(con, "con.qs")
```

save clustering
```{r}
leiden <- con$clusters$leiden$groups %>% factor
qsave(leiden, "leiden.qs")
```
annotation is done using functions in  "Annotation.Rmd" and is called "anno"

## Figure 24: final UMAP
define colours per sample
```{r}
colours <- c( "#FF0000" ,"#FF6000", "#FFBF00", "#DFFF00" ,"#00FFFF" , "#20FF00", "#00FF40", "#DF00FF","#80FF00","#009FFF", "#0040FF","#FF0060" , "#8000FF", "#00FF9F", "#FF00BF","#2000FF" ) %>% setNames(c("Cajal_Retzius_cells", "CGE_IPC","Cortical_hem"  ,  "Dorsal_IPC" ,  "Immature_GABAergic_neurons","Highly_proliferative_RG_S" , "Immature_principal_neurons",    "Immature_principal_neurons_Tbr1","MGE_IPC","Microglia",
 "Dorsal_RG" , "Anteriomedial_RG", "Lateral_RG", "PIC_enriched_immature_principal_neurons", "Roof_plate","Highly_proliferative_RG_G2M"))
```

plot and export UMAP without labeling
```{r}
plot<- con$plotGraph(groups=anno, color.by="cluster",mark.groups=F, size=0.2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ scale_colour_manual(values=colours)
plot
rasterize(plot, layers='Point', dpi=1000)
ggsave("UMAP_con_without_labels.pdf", width=14, height=8.6 )
```

## Figure 21: UMAP coloured by samples, condition, sex, and Xist

UMAP coloured by sample
```{r}
data <- con$plotGraph(groups=con$getDatasetPerCell())
data <- data$data

plot1 <- ggplot(data, aes(x = x, y =y, color= Group))+ geom_point(size=0.2, alpha=0.5)+ theme_light()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ xlab("")+ylab("") + theme(axis.ticks = element_blank(), axis.text =element_blank())+labs(color = "Sample") +guides(color = guide_legend(override.aes = list(size=3, alpha=1)))+scale_color_manual(values=rainbow(12),labels = new_labels)
plot1
rasterize(plot1, layers='Point', dpi=1000)
```

UMAP coloured by condition
```{r}
data <- con$plotGraph(groups=condition)
data <- data$data

plot2<-ggplot(data, aes(x = x, y =y, color= Group))+ geom_point(size=0.2, alpha=0.5)+ theme_light()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ xlab("")+ylab("") + theme(axis.ticks = element_blank(), axis.text =element_blank())+labs(color = "Condition") +guides(color = guide_legend(override.aes = list(size=3, alpha=1)))+scale_color_manual(values=rainbow(2)) 


plot2
rasterize(plot2, layers='Point', dpi=1000)
```

UMAP coloured by gender
create factor with gender per cell
```{r, eval=T}
gender <- con$getDatasetPerCell() 
levels(gender)
levels(gender) <- c("male", "mix", "male", "mix", "mix", "mix", "mix", "female", "mix", "mix", "mix", "mix")
```
```{r}
data <- con$plotGraph(groups=gender)
data <- data$data

plot3 <- ggplot(data, aes(x = x, y =y, color= Group))+ geom_point(size=0.2, alpha=0.5)+ theme_light()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ xlab("")+ylab("") + theme(axis.ticks = element_blank(), axis.text =element_blank())+labs(color = "Sex") +guides(color = guide_legend(override.aes = list(size=3, alpha=1)))+scale_color_manual(values=rainbow(3))
plot1
rasterize(plot3, layers='Point', dpi=1000)
```

UMAP with Xist expression
```{r}
data <- con$plotGraph(gene="Xist")
data <- data$data

plot4 <- ggplot(data, aes(x = x, y =y, color= Color))+ geom_point(size=0.2, alpha=0.5)+ theme_light()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ xlab("")+ylab("") +scale_color_gradient(low="#dbd9d9",high="#FF0000")+ theme(axis.ticks = element_blank(), axis.text =element_blank())+labs(color = "Xist")
plot4
rasterize(plot4, layers='Point', dpi=1000)
```

figure with 4 UMAPS
```{r, fig.height=12}
gridplot <- cowplot::plot_grid(plot1, plot2, plot3, plot4, align = "v", axis = "lr")
gridplot
ggsave("4_UMAPs_coloured.png", width=14, height=8.6 )
```

## Figure 22: dotplot of sex genes
get count matrix
```{r}
cm <- con$getJointCountMatrix(raw=F) 
#%>% 
  #Matrix::t() %>% 
  #as.sparse()
```
```{r}
sex_genes <- c ("Xist", "Uty") 
sccore:::dotPlot(cm, markers=sex_genes, cell.groups=con$getDatasetPerCell(), dot.scale = 10, cols=c("white","#0E4B99"), gene.order = T)
ggsave("dotplot_sex.pdf", width=4, height=4)
```

## Figure 25: dotplot of cluster marker genes
```{r}
#reorder anno factor levels
anno_order <- factor(anno, levels=c('Cajal_Retzius_cells', 'Immature_PN_Tbr1' , 'Immature_PN','PIC_immature_PN', 'Dorsal_IPC','Immature_GABAergic_N','MGE_progenitors','CGE_progenitors', 'Roof_plate', 'Cortical_hem','Highly_prolif_RG_S','Highly_prolif_RG_G2M','Lateral_RG', 'Anteriomedial_RG', 'Dorsal_RG', 'Microglia'))
```

```{r}
markers <- c("Trp73", "Reln", 'Eomes','Tbr1','Slc17a6' ,'Dcx','Thsd7b','Neurod1' ,'Neurog2','Btg2','Gad2','Dlx1','Nkx2-1','Shh','Nr2f2','Ptprz1','Wnt3','Rspo3', 'Rspo2', 'Bmp6',"Pax6",'Mcm4','Sox2','Nes','Rplp1','Meis2','Mpped1','Zic1','Zic4','Emx2', 'Creb5','Csf1r', "Spi1")
```
```{r}
RColorBrewer::brewer.pal.info
RColorBrewer::brewer.pal(3,"Blues")
```

plot dotplot and export
```{r, fig.width=7}
sccore::dotPlot(markers, cm, cell.groups = anno_order, gene.order = T, xlab = "", ylab = "", cols=c("white","#0E4B99")) +theme_light()+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+ theme(text = element_text(size = 11)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("dotplot_marker_genes_2.png", width=9.5, height=4.0)
```

## Figure 26
### barplot number of cells per cluster
```{r}
anno_table <- as.data.frame(table(anno_order))
ggplot(anno_table, aes(anno_order, Freq)) + geom_bar(stat = "identity",fill="#2F5496", alpha=0.2, color="#2F5496", width=0.8) + xlab("") + ylab("number of nuclei") + theme_light()+ theme(axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.title.y = element_text(size = 9), axis.text.x =  element_text(size = 9),axis.text.y =  element_text(size = 9))
ggsave("barplot_nuclei_per_cluster.pdf", width=4, height=3 )  
```

### stacked barplot: samplecomposition per cluster
```{r}
con$clusters$leiden$groups <- anno_order
data <- plotClusterBarplots(con, legend.height = 0.1, show.entropy=F, show.composition = T, show.size = F)
data <- data$data

new_labels <- c("control_1", "control_2","control_3", "control_4","control_5","control_6","treated_1","treated_2","treated_3","treated_4","treated_5","treated_6")

p1 <- ggplot(data, aes(fill=sample, y=f, x=cluster))+ theme_light()+ geom_bar(position="fill", stat="identity", alpha=1)+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))  + xlab("")+ ylab("Fraction of nuclei") + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+theme(legend.position="bottom")+guides(
         fill = guide_legend(override.aes = list(size = 0.5),title="Sample", title.position="top", title.hjust = 0.5))+ theme(text = element_text(size = 14)) + scale_y_continuous(expand = c(0,0)) +
  scale_fill_discrete(labels = new_labels)
p1
```

### violin plots
using seurat object created further down
```{r}
 iqr = function(z, lower = 0.25, upper = 0.75) {
  data.frame(
    y = median(z),
    ymin = quantile(z, lower),
    ymax = quantile(z, upper)
  )
}
```
```{r}
p<- VlnPlot(seurat, features="nFeature_RNA", pt.size = 0)
p<- p$data
ggplot(p, aes(y = nFeature_RNA, x = ident, color = ident)) +
  geom_violin(scale = "width", aes(fill = ident), 
                color = NA, alpha = 0.2, 
              trim = F)+theme_light()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme( legend.position = "none")+  ggtitle(" ") + ylab("genes per nucleus") + xlab(" ")+ stat_summary(fun.data = iqr, geom = "errorbar", show.legend = F, color = "grey30", width = 0.25, stroke=0.5) +stat_summary(fun.data = iqr, geom = "point", show.legend = F, shape = 23, stroke = 0.5, color = "black", alpha = 0.9,aes(fill = ident), size=0.7)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values =colours)+ylim(0, 10000) +
  theme(axis.title.y = element_text(size = 9), axis.text.x =  element_text(size = 9),axis.text.y =  element_text(size = 9))
ggsave("vln_genes_per_cluster.pdf", width=4, height=3 ) 
```
```{r}
p<- VlnPlot(seurat, features="nCount_RNA", pt.size = 0)
p<- p$data
ggplot(p, aes(y = nCount_RNA, x = ident, color = ident)) +
  geom_violin(scale = "width", aes(fill = ident), 
                color = NA, alpha = 0.2, 
              trim = F)+theme_light()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme( legend.position = "none")+  ggtitle(" ") + ylab("UMIs per nucleus") + xlab(" ")+ stat_summary(fun.data = iqr, geom = "errorbar", show.legend = F, color = "grey30", width = 0.25) +stat_summary(fun.data = iqr, geom = "point", show.legend = F, shape = 23, stroke = 0.5, color = "black", alpha = 0.9,aes(fill = ident), size=0.7)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ scale_fill_manual(values =colours)+#+ylim(0, 100000)+
  scale_y_continuous(trans = "log10") +
  theme(axis.title.y = element_text(size = 9), axis.text.x =  element_text(size = 9),axis.text.y =  element_text(size = 9))
  
ggsave("vln_UMI_per_cluster.pdf", width=4, height=3 ) 
```

## Figure 38: Cytokine receptor expression in ctrl samples
```{r}
cm <- con$getJointCountMatrix(raw=F) 
#%>% 
  #Matrix::t() %>% 
  #as.sparse()
```
```{r}
cm_ctrl <- cm[grepl("control", rownames(cm)), ]
qsave(cm_ctrl, "cm_ctrl.qs")
head (anno_order)
anno_order_ctrl <- anno_order[grepl("control", names(anno_order))]
```

```{r}
receptors <- c(
"Il6ra",
"Il6st",
"Il10rb",
"Il11ra1",
"Il13ra1",
"Il15ra",
"Il17ra",
"Il17rb",
"Il17rc",
"Il17rd",
"Il21r",
"Il20rb",
"Il12rb2",
"Il4ra",
"Cx3cr1",
"Cxcr2",
"Cxcr4",
"Ccr1",
"Ccr9",
"Tnfrsf1a",
"Tnfrsf1b",
"Tnfrsf13b",
"Tnfrsf19",
"Tnfrsf21",
"Ifnar1",
"Ifnar2",
"Ifngr1",
"Ifngr2",
"Ifnlr1")
```

```{r, fig.height=12}
sccore::dotPlot(receptors, cm_ctrl, cell.groups = anno_order_ctrl, gene.order = T, xlab = "", ylab = "", cols=c("white","#0E4B99")) +theme_light()+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+ theme(text = element_text(size = 11)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("dotplot_cytokine_receptors.png", width=8.8, height=3.8)
```

# 6. Pagoda web app
create pagoda web app to interactively explore dataset
```{r, eval=FALSE, echo=FALSE}
# Merge count matrices
cm <- conos:::mergeCountMatrices(cms,F)
```

```{r, eval=FALSE, echo=FALSE}
# Pagoda2 processing on merged matrices
p2 <- basicP2proc(cm, n.cores = 20, n.odgenes = 3e3, nPcs = 100, k = 30, perplexity = 50, log.scale = TRUE, trim = 10, keep.genes = NULL, 
min.cells.per.gene = 3, min.transcripts.per.cell = 200, get.largevis = T, get.tsne = F, make.geneknn = TRUE)

p2 <- extendedP2proc(p2, organism="mm")
go.env2 <- p2$go.env
p2 <- p2$p2

#Make cell metadata
mit_frac <- mitoFraction(con, "mouse")
sample <- con$getDatasetPerCell()
mit_frac <- as.factor(round(mit_frac, digits=3))
metadata.listfactors <- list(sample = sample,
                             mit_frac = mit_frac,
                             annotation = anno,
                             leiden= con$clusters$leiden$groups)
metadata.forweb <- factorListToMetadata(metadata.listfactors)
```

Create web object
```{r, eval=FALSE}
p2.web <- webP2proc(p2, additionalMetadata = metadata.forweb, title = 'PolyIC E10.5 2024', go.env = go.env2, make.go.sets = F)

#Add UMAP embedding
p2.web <- addEmbeddingP2Web(p2.web, con)

#Create .bin-file
p2.web$serializeToStaticFast('p2web_final.bin',verbose=T)
```

# 7. Cell cycle seurat
## generate seurat object
get just raw counts,
```{r}
cm_raw <- con$getJointCountMatrix(raw=T) %>% 
  Matrix::t() %>% 
  as.sparse()
```

Create Seurat object from raw cm
```{r}
seurat <- CreateSeuratObject(counts = cm_raw, project = "Embryonic_PolyIC", min.cells = 0, min.features = 0)
```

add normalized cm
```{r}
cm <- con$getJointCountMatrix(raw=F) %>% 
  Matrix::t() %>% 
  as.sparse()
seurat <- SetAssayData(object= seurat, slot= "data", cm )
```

add annotation as metadata
```{r}
seurat$annotation <- anno
seurat$orig.ident <- anno
```
```{r}
level_order <- c('Cajal_Retzius_cells', 'Immature_principal_neurons_Tbr1' , 'Immature_principal_neurons','PIC_enriched_immature_principal_neurons', 'Dorsal_IPC','Immature_GABAergic_neurons','MGE_IPC','CGE_IPC', 'Roof_plate', 'Cortical_hem','Highly_proliferative_RG_S','Highly_proliferative_RG_G2M','Lateral_RG', 'Anteriomedial_RG', 'Dorsal_RG', 'Microglia')
```

reorder annotation
```{r}
Idents(seurat) <- "annotation" 
seurat$annotation <- factor(x = seurat$annotation, levels = level_order)
```

add conditions meta data
```{r, eval=T}
condition <- con$getDatasetPerCell() %>% substr(1, 7) %>% 
  setNames(names(con$getDatasetPerCell()))
condition <- as.factor(condition)
```
```{r}
seurat$condition <- condition
```
```{r}
Idents(seurat) <- "annotation"
```
```{r}
qsave(seurat, "seurat.qs")
```


## perform cell cycle scoring
```{r}
# A list of cell cycle markers, from Tirosh et al, 2015 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4944528/), is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
```
```{r}
seurat <- CellCycleScoring(seurat, g2m.features = g2m.genes, s.features = s.genes)
```

### Figure 27A: cellcycle composition
```{r}
data <- plotClusterBarplots(con, groups = anno_order , show.entropy = F, show.size = F, sample.factor = seurat$Phase)
data <- data$data
p3 <- ggplot(data, aes(fill=sample, y=f, x=cluster)) + geom_bar(position="fill", stat="identity", alpha=1)+ theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +theme(text=element_text(size=14)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+ xlab("") + ylab("Fraction of nuclei") +theme(legend.position="bottom")+guides(
         fill = guide_legend(title="Cell cycle phase",title.position="top", title.hjust = 0.5)) + scale_y_continuous(expand = c(0,0))
p3
```

### Figure 27B: UMAP coloured by cell cycle
```{r}
data <- con$plotGraph(groups=seurat$Phase)
data <- data$data

plot3 <- ggplot(data, aes(x = x, y =y, color= Group))+ geom_point(size=0.2, alpha=0.5)+ theme_light()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ xlab("")+ylab("") + theme(axis.ticks = element_blank(), axis.text =element_blank())+labs(color = "Cell cycle") +guides(color = guide_legend(override.aes = list(size=3,alpha=1)))+scale_color_manual(values=rainbow(3), labels = c("G0", "G2M", "S"))
plot3
rasterize(plot3, layers='Point', dpi=1000)
ggsave("UMAP_cell_cycle.pdf", width=14, height=8.6 )
```

# 8. Prepare cell names with annotation, embedding coordinates, condition and sample metadata for python (scVelo)

##annotation
```{r}
head(anno)
```
```{r}
annoscv <- anno
names(annoscv) <- gsub("!!", ":", names(annoscv))
names(annoscv) <- gsub("-1", "x", names(annoscv))
head(annoscv)
```
```{r}
write.csv(annoscv, 'anno.csv')
```

## cell names for scVelo
```{r}
cells <- names(anno)
head(cells)
```

```{r}
cells <- gsub("!!", ":", cells)
cells <- gsub("-1", "x", cells)
head(cells)
```
```{r}
write.csv(cells, 'cells.csv')
```

## embedding for scVelo
```{r}
con_embedding <- con$embedding
head(con_embedding)
```

```{r}
# order names of embedding  as in "anno"
con_embedding %<>% as.data.frame()
con_embedding <- con_embedding[match(names(anno), rownames(con_embedding)),] 
head(con_embedding)
```
```{r}
rownames(con_embedding) <- gsub("!!", ":", rownames(con_embedding))
rownames(con_embedding) <- gsub("-1", "x", rownames(con_embedding))
head(con_embedding)
```

```{r}
write.csv(con_embedding, 'embedding.csv')
```

## sample and condition per cell for scvelo
```{r}
sample <- con$getDatasetPerCell()
head(sample)
```

```{r}
# order cells in "sample" as in anno
sample <- sample[match(names(anno), names(sample))] 
head(sample)
```
```{r}
condition <- con$getDatasetPerCell() %>% substr(1, 7) %>% 
  setNames(names(con$getDatasetPerCell()))
head(condition)
```

```{r}
# order cells in "condition" as in "anno"
condition %<>% .[match(names(anno), names(.))] 
head(condition)
```

create one dataframe with sample info, condition and anotation per cell
```{r}
metadata <- data.frame( sample=as.character(sample),condition=condition, annotation=anno, row.names = names(anno))
head(metadata)

```
```{r}
rownames(metadata) <- gsub("!!", ":", rownames(metadata))
rownames(metadata) <- gsub("-1", "x", rownames(metadata))
head(metadata)
```
```{r}
write.csv(metadata, 'metadatascv.csv')
```

the saved files are used for RNA velocity analysis with scvelo in python: "scVelo.ipynb"


#9. Cellchat DB
```{r}
cellchat_ctrls <-qread("/people/nrq364/embryonic_PolyIC_2/cellchat_ctrls2.qs")
cellchat_treated <- qread("/people/nrq364/embryonic_PolyIC_2/cellchat_treated2.qs")
```
```{r}
library(CellChat) #version 1.6.1
library(patchwork)
options(stringsAsFactors = FALSE)
library(magrittr)
library(ComplexHeatmap)
```

## prepare ctrl and treated cm
```{r}
cm <- con$getJointCountMatrix(raw=F)%>% 
  Matrix::t() %>% 
  as.sparse() 
saveRDS(cm, "cm2_normalized.rds")
```

```{r}
anno_ctrls <- anno[grep("control", names(anno))]
ctrl_cells <- names(anno_ctrls)
anno_treated <- anno[grep("treated", names(anno))]
treated_cells <- names(anno_treated)
```
```{r}
length(ctrl_cells)
length(treated_cells)
```

```{r, eval=FALSE}
#keep "keep_cluster" from cms
cm_ctrl <- cm[,colnames(cm) %in% ctrl_cells]
cm_treated <- cm[,colnames(cm) %in% treated_cells]
```

```{r}
qsave(cm_ctrl, "cm_ctrl.qs")
qsave(cm_treated, "cm_treated.qs")
```

## create one seurat object per condition
create seurat object with cm normalized in pagoda/conos
```{r}
seurat_ctrls <- CreateSeuratObject(counts = cm_ctrl, project = "Embryonic_PolyIC", min.cells = 0, min.features = 0)
seurat_treated <- CreateSeuratObject(counts = cm_treated, project = "Embryonic_PolyIC", min.cells = 0, min.features = 0)
```

annotation as meta data
```{r}
seurat_ctrls$annotation <- anno_ctrls
seurat_treated$annotation <- anno_treated
```

```{r}
seurat_ctrls$orig.ident <- anno_ctrls
seurat_treated$orig.ident <- anno_treated
```
```{r}
qsave(seurat_ctrls, "seurat_ctrls.qs")
qsave(seurat_treated, "seurat_treated.qs")
```
## create cell chat objects
```{r}
cellchat_ctrls <- createCellChat(object = seurat_ctrls, group.by = "annotation")
cellchat_treated <- createCellChat(object = seurat_treated, group.by = "annotation")
```
## set database of ligand-receptor pairs
Our database CellChatDB is a manually curated database of literature-supported ligand-receptor interactions in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60% of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions and 19% of cell-cell contact interactions
```{r}
CellChatDB <- CellChatDB.mouse
```
```{r}
# set the used database in the object 
cellchat_ctrls@DB <- CellChatDB
cellchat_treated@DB <- CellChatDB
```
```{r}
# subset the expression data of signaling genes for saving computation cost
cellchat_ctrls %<>% subsetData() 
cellchat_treated %<>% subsetData()
```
```{r}
future::plan("multiprocess", workers = 4) # do parallel
cellchat_ctrls %<>% identifyOverExpressedGenes()
cellchat_ctrls %<>% identifyOverExpressedInteractions()
cellchat_treated %<>% identifyOverExpressedGenes()
cellchat_treated %<>% identifyOverExpressedInteractions()
```

## Compute the communication probability
By default, CellChat uses a statistically robust mean method called ‘trimean’, which produces fewer interactions than other methods. However, we find that CellChat performs well at predicting stronger interactions, which is very helpful for narrowing down on interactions for further experimental validations. Of note, ‘trimean’ approximates 25% truncated mean, implying that the average gene expression is zero if the percent of expressed cells in one group is less than 25%.
When setting population.size = TRUE, abundant cell types tend to have strong interactions than rare cells.
```{r}
cellchat_ctrls %<>% computeCommunProb(population.size=F, raw.use= TRUE) 
cellchat_treated %<>% computeCommunProb(population.size=F, raw.use= TRUE)
```
##Infer the cell-cell communication at a signaling pathway level
CellChat computes the communication probability on signaling pathway level by summarizing the communication probabilities of all ligands-receptors interactions associated with each signaling pathway.
```{r}
cellchat_ctrls %<>% computeCommunProbPathway()
cellchat_treated %<>% computeCommunProbPathway()
```

##Calculate the aggregated cell-cell communication network
We can calculate the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability. 
```{r}
cellchat_ctrls %<>% aggregateNet()
cellchat_treated %<>% aggregateNet()
```

```{r}
# Compute the network centrality scores
cellchat_ctrls <- netAnalysis_computeCentrality(cellchat_ctrls, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_treated <- netAnalysis_computeCentrality(cellchat_treated, slot.name = "netP")
```
```{r}
qsave(cellchat_ctrls, "cellchat_ctrls.qs")
qsave(cellchat_treated, "cellchat_treated.qs")
```
## Comparison analysis
### merge ctrl and treated
```{r}
cellchat <- qread("cellchat.qs")
```
```{r}
object.list <- list(control = cellchat_ctrls, treated = cellchat_treated)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
```

### Figure 32
```{r}
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
g <- gg1 + gg2
g
ggsave("CellChat.pdf")
```

```{r}
qsave(cellchat, "cellchat.qs")
```

