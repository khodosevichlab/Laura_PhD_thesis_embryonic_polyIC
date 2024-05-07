---
title: "CRMetrics_preprocessing"
author: "Laura Wolbeck"
date: "2023-09-15"
output: html_document
---

This script performs quality control steps on the count matrices retrieved from Cell Ranger pipeline.

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
```

# 1. Create CRMetrics object
```{r, eval=F}
crm <- CRMetrics$new(data.path = "/data/counts/embryonic_PolyIC_counts_v7.1.0/", 
                     metadata = "metadata.csv", 
                     n.cores = 60)
```

total number of sequenced nuclei
```{r}
sum(sapply(crm$cms, function(d) dim(d)[2]))
```
```{r}
metrics <- crm$summary.metrics
metrics %<>% as.data.frame()
```

median genes and UMIs across all samples
```{r}
#median genes per cell across all samples
median(metrics[metrics$metric=="median genes per cell",]$value)
#median UMIs per cell across all samples
median(metrics[metrics$metric=="median umi counts per cell",]$value)
```

# 2. cellbender v0.2.0
```{r, eval=F, fig.height=8}
crm$prepareCellbender(shrinkage = 100, 
                      show.expected.cells = TRUE, 
                      show.total.droplets = TRUE, n.cores = 60)
```

prepare vector of total droplets included per sample, we use three times the cell count estimated from Cell Ranger
```{r}
droplets <- crm$getTotalDroplets()
droplets
```

save cellbender script and then run it in the terminal
```{r, eval=F}
crm$saveCellbenderScript(file = "cellbender_script.sh", 
                         fpr = 0.01, 
                         epochs = 150, 
                         use.gpu = TRUE, total.droplets = droplets)
```

inspect cellbenders results
## Figure 16
```{r}
crm$plotCbCells() +
  scale_fill_manual(values = c("#1B9E77","#1B9E77","#1B9E77","#1B9E77","#1B9E77","#1B9E77", "#9E0142","#9E0142","#9E0142","#9E0142","#9E0142","#9E0142"))
```

Plot cellbenders training performance
```{r, fig.width = 12, fig.height = 10}
crm$plotCbTraining()# + theme(legend.text = element_text(size=20), legend.key.size = unit(2, 'cm'))
```
```{r}
crm$plotCbAmbExp()
```
```{r}
crm$plotCbAmbGenes()
```

from crm object remove current count matrices and add count matrices after cellbender correction
```{r, eval=F}
crm$cms = NULL
crm$detailed.metrics = NULL
crm$addDetailedMetrics(cellbender = T, n.cores = 60)
qsave(crm, "crm_cellbender.qs", nthreads = 10)
```
```{r}
crm <- qread("crm_cellbender.qs", nthreads = 10)
```


# 3. Conos 1 (before filtering)
in this section a first embedding is estimated to visualize nuclei that need to be removed on the UMAP embedding

perform standard preprocessing using pagoda, i.e. normalization
```{r, eval=F}
crm$doPreprocessing(preprocess = "pagoda")
```

```{r, eval=F}
crm$createEmbedding()
```

```{r, eval=F}
crm$plotEmbedding()
```

# 4. Cell filtering
## scrublet
estimate doublets using scrublet
```{r, eval=F}
crm$detectDoublets(env = "r-reticulate", conda.path = "/opt/software/miniconda/4.12.0/condabin/conda", method = "scrublet")
```

### Figure 18B
```{r}
plot <- crm$plotEmbedding(doublet.method = "scrublet")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
rasterize(plot, layers='Point', dpi=1000)
ggsave("UMAP_scrublet.pdf", width=7, height=4.3 )
```

## depth (UMI counts per nucleus)
inspect density histogram of UMI counts per sample
```{r, fig.heigth=12, fig.width=12, warning=F}
crm$plotDepth(cutoff = 1000)
```

### Figure 18C
plot doublets on UMAP 
```{r}
plot <- crm$plotEmbedding(depth = TRUE, 
             depth.cutoff = 1000) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot
```

```{r}
rasterize(plot, layers='Point', dpi=1000)
ggsave("UMAP_depth.pdf", width=7, height=4.3 )
```

## Mitochondrial gene fraction
plot nuclei >5% mt genes on UMAP
```{r}
crm$plotEmbedding(mito.frac = TRUE, 
             mito.cutoff = 0.05, 
             species = "mouse")
```

## Plot filtered cells
And we can plot the cells to be filtered per sample where `combination` means a cell that has at least two filter labels, e.g. `mito` and `depth`.
```{r}
crm$plotFilteredCells(type = "bar", 
                      doublet.method = "scrublet", 
                      depth = TRUE, 
                      depth.cutoff = 1000,  mito.frac = TRUE, mito.cutoff = 0.05, species="mouse")
```

# 5. Save filtered cms
```{r eval = FALSE}
crm$filterCms(depth.cutoff = 1000, 
              doublets = "scrublet",
              samples.to.exclude = NULL,
              species = "mouse", mito.cutoff = 0.05)
```

The filtered list of count matrices is stored in $cms.filtered which can be saved on disk afterwards.
```{r, eval=F}
qsave(crm$cms.filtered, "cms_filtered_cellbender.qs", 
      nthreads = 10)
```

number of cells after filtering
```{r}
sum(sapply(cms, function(d) dim(d)[2])) 
```

# 6. Plot genes/nucleus and UMIs/nucleus after all filtering steps
```{r}
cms <- qread("cms_filtered_cellbender.qs")
```

create a new crm object just for plotting the VlnPlots
```{r, eval=F}
crm <- CRMetrics$new(cms =  cms, 
                     n.cores = 60, metadata="metadata.csv")
crm$addSummaryFromCms()
crm$addDetailedMetrics()
```

number of nuclei after filtering 
```{r}
sum(sapply(crm$cms, function(d) dim(d)[2]))
```
```{r}
metrics <- crm$summary.metrics
metrics %<>% as.data.frame()
```

```{r}
#median genes per cell across all samples
median(metrics[metrics$metric=="median genes per cell",]$value)
#median UMIs per cell across all samples
median(metrics[metrics$metric=="median umi counts per cell",]$value)
```

```{r}
 iqr = function(z, lower = 0.25, upper = 0.75) {
  data.frame(
    y = median(z),
    ymin = quantile(z, lower),
    ymax = quantile(z, upper)
  )
}
```

## plot genes/nucleus
```{r}
p <- crm$plotDetailedMetrics(comp.group = "condition",
                       metrics = "gene_count", 
                        plot.geom = "violin" , hline = F)
p <- p$data
```

```{r}
p1 <- p %>% ggplot(aes( y = value, x = sample,fill=condition )) + geom_violin(alpha=0.5, color = NA)+ ylab("genes per nucleus")+xlab("") + scale_fill_manual(values = c("#1B9E77","#9E0142"))+ theme_light()   + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun.data = iqr, geom = "errorbar", show.legend = F, color = "grey30", width = 0.25) +stat_summary(fun.data = iqr, geom = "point", show.legend = F,aes(fill = condition), shape = 21, stroke = 0.5, color = "black", alpha = 1, size = 0.5)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ scale_x_discrete(labels = c("control_1", "control_2","control_3", "control_4","control_5","control_6","treated_1","treated_2","treated_3","treated_4","treated_5","treated_6")) +guides(fill = FALSE)
p1
```

## plot UMIs/nucleus
```{r, fig.height=3}
plot <- crm$plotDetailedMetrics(comp.group = "condition",
                        metrics = "UMI_count", 
                        plot.geom = "violin" , hline = F)

plot <- plot$data
```

```{r}

p2 <- plot %>% ggplot(aes( y = value, x = sample,fill=condition )) + geom_violin(alpha=0.5, color = NA)+ ylab("UMIs per nucleus")+xlab("") + scale_fill_manual(values = c("#1B9E77","#9E0142"))+ theme_light()   + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun.data = iqr, geom = "errorbar", show.legend = F, color = "grey30", width = 0.25) +stat_summary(fun.data = iqr, geom = "point", show.legend = F,aes(fill = condition), shape = 21, stroke = 0.5, color = "black", alpha = 1, size = 0.5)+coord_cartesian(ylim = c(1000,100000))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ scale_y_log10()+ scale_x_discrete(labels = c("control_1", "control_2","control_3", "control_4","control_5","control_6","treated_1","treated_2","treated_3","treated_4","treated_5","treated_6"))
```

### Figure 17
plot both plots 
```{r}
gridplot <- cowplot::plot_grid(p1, p2, align = "vh", axis = "lr", label_colour ="#2F5496", label_y = 1)
gridplot
ggsave("Vln_UMI_genes_per_nucleus.pdf", width=20, height=6, units = "cm" )
```


Now the dataset is filtered for ambient RNA, low UMI count, high mitochondrial gene fraction and doublets.
Further filtering based on cell types/marker genes follows in the notebook "Notebook_2_Filtering_for_cell_types"