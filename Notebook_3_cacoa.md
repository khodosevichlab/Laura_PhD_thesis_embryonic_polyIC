---
title: "Cacoa"
author: "Laura Wolbeck"
date: "2024-01-08"
output: html_document
---
This script performs the differential analysis using the R package cacoa.

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
library(openxlsx)
library(writexl)
library(EnhancedVolcano)
library(ggrepel)
library(cacoa)
```

#1.  Read data

```{r}
con <- qread("/people/nrq364/embryonic_PolyIC_2/con.qs")
anno <- qread("anno.qs")
```

#2. Generate cacoa object
```{r,eval=FALSE }
condition <- setNames(c("control", "control", "control", "control","control", "control", "treated", "treated", "treated","treated", "treated", "treated"), names(con$samples))
```

```{r,eval=FALSE}
cao <- Cacoa$new(con, cell.groups = anno, sample.groups=condition, n.cores = 10, ref.level = "control", target.level = "treated")
```

define colors
```{r}
cao$cell.groups.palette <- c( "#FF0000" ,"#FF6000", "#FFBF00", "#DFFF00" ,"#00FFFF" , "#20FF00", "#00FF40", "#DF00FF","#80FF00","#009FFF", "#0040FF","#FF0060" , "#8000FF", "#00FF9F", "#FF00BF","#2000FF" ) %>% setNames(c("Cajal_Retzius_cells", "CGE_IPC","Cortical_hem"  ,  "Dorsal_IPC" ,  "Immature_GABAergic_neurons","Highly_proliferative_RG_S" , "Immature_principal_neurons",    "Immature_principal_neurons_Tbr1","MGE_IPC","Microglia",
 "Dorsal_RG" , "Anteriomedial_RG", "Lateral_RG", "PIC_enriched_immature_principal_neurons", "Roof_plate","Highly_proliferative_RG_G2M"))
cao$sample.groups.palette <- c("#9E0142", "#1B9E77") %>% setNames(c("treated", "control"))
qsave(cao, "cao_without_treated5_and_ctrl3.qs")
```

#3. Cluster-based expression shifts
```{r}
cao$estimateExpressionShiftMagnitudes(top.n.genes = 300, n.permutations = 1000) #n.permutations 1000 is the default, for top.n.genes set a lower number because neurons are very similar, if you normalize to a big chunk of similar genes it will not give you anything
```

## Figure 30
```{r}
p3 <- cao$plotExpressionShiftMagnitudes(type = "box")

 
 iqr = function(z, lower = 0.25, upper = 0.75) {
  data.frame(
    y = median(z),
    ymin = quantile(z, lower),
    ymax = quantile(z, upper)
  )
}
```

```{r}
p <- p3$data %>% ggplot(aes( y = value, x = Type, color = Type)) +
  geom_violin(scale = "width", 
                color = NA, aes(fill = Type), alpha = 0.2, 
              trim = F) + theme_light() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + scale_color_manual(values = cao$cell.groups.palette)+ scale_fill_manual(values = cao$cell.groups.palette) +
  geom_hline(yintercept = 0, linewidth = 0.1) +
  stat_summary(fun.data = iqr, geom = "errorbar", show.legend = F, aes(color = Type),stroke = 0.5, alpha = 0.9,color = "grey30", width = 0.25) + 
  stat_summary(fun.data = iqr, geom = "point", show.legend = F,aes(fill = Type), shape = 23, color = "grey30", alpha = 0.9)  + ylab(("normalized expression\ndistance")) + xlab("") + theme(plot.margin = margin(0.5,0,0,2, "cm"))+scale_x_discrete(labels= c("Immature_GABAergic_N", "Highly_prolif_RG_S", "PIC_immature_PN","Dorsal_RG","Highly_prolif_RG_G2/M","Lateral_RG", "Cortical_hem","MGE_progenitors", "Dorsal_IPC", "Microglia", "CGE_progenitors", "Roof_plate", "Anteriomedial_RG", "Immature_PN","Immature_PN_Tbr1", "Cajal_Retzius_cells"))+theme(axis.text.x = element_text( vjust = 0.5)) +theme(axis.text.x =  element_text(size = 10))
p
ggsave("Expressionshifts_Vln.png",width=7, height=4.5)
# the more to the right the more affected, below 0 does not mean the controls are more different than treated, just means there is nothing

#with raw p.value the 2 most right cell types are significant
```

#4. Pairwise expression distances between samples 

```{r}
new_labels <- c("control_1", "control_2","control_3", "control_4","control_5","control_6","treated_1","treated_2","treated_3","treated_4","treated_5","treated_6")
```

## Figure 23
colouring by 10x version, sequencing batch, and sex was done in adobe illustrator
```{r}
cao$plotSampleDistances(space='expression.shifts', font.size=4, show.sample.size=T, method="UMAP") + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +theme( legend.text = element_text(size = 12))
ggsave("UMAP_expression_shifts.pdf",width=7, height=4.3 )
```

#5. Cluster-based compositional changes (cell loadings)

```{r}
cao$estimateCellLoadings()
```

```{r, fig.height=6}
 cao$plotCellLoadings(signif.threshold=0.05, show.pvals = T) 
```

## Figure 29
```{r}
p2 <- cao$plotCellLoadings(show.pvals = F)
p2 <- p2$data

library(data.table)
p2 <- data.table(p2)
lvl <- p2[, abs(median(values)), by = ind][order(-V1),ind] %>% as.character()
levels(p2$ind) <- rev(lvl)

p2$ind <- factor(p2$ind, levels = rev(lvl))

palcomp <- setNames(ifelse(p2[, (median(values)), by = ind][order(-V1),V1] > 0, "#9E0142", "#1B9E77"), p2[, (median(values)), by = ind][order(-V1),ind])
#palcomp2 <- setNames(ifelse(p2[, (median(values)), by = ind][order(-V1),V1] > 0, "darkblue", "grey50"), p2[, (median(values)), by = ind][order(-V1),ind])

ord <- p2[, abs(median(values)), by = ind][order(V1)]$ind
```

```{r}
p <- p2 %>% ggplot( aes( x = values, y = factor(ind, levels = ord), color = ind)) +
  geom_violin(scale = "width", 
                #fill = "#ededed", 
                color = NA, aes(fill = ind), alpha = 0.2) +
  theme_light() + theme(legend.position = "none") + scale_fill_manual(values = palcomp) +
  geom_vline(xintercept = 0, linewidth = 0.1) +
  #new_scale_fill() + 
  stat_summary(fun.data = iqr, geom = "errorbar", show.legend = F, aes(color = ind), stroke = 0.6, alpha = 0.9, color = "grey30", width = 0.25) + 
  stat_summary(fun.data = iqr, geom = "point", show.legend = F,aes(fill = ind), shape = 23, stroke = 1, color = "grey30", alpha = 0.9, size = 2) +  scale_fill_manual(values = palcomp) + xlab("separation coefficient") + ylab("") + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

p+scale_y_discrete(labels= rev(c("Immature_PN_Tbr1", "PIC_immature_PN","Highly_prolif_RG_S", "CGE_progenitors", "Cajal_Retzius_cells", "Roof_plate", "Immature_PN","MGE_progenitors", "Dorsal_IPC","Microglia","Dorsal_RG","Highly_prolif_RG_G2/M","Anteriomedial_RG","Lateral_RG","Cortical_hem","Immature_GABAergic_N")))+theme(axis.text.y = element_text( vjust = 0.5)) +theme(axis.text.y =  element_text(size = 10))
ggsave("Compositionshifts.pdf",width=7, height=4.5)
```

#6. DEG analysis
DEG per cluster, number of cells per cluster not fixed so number of DEGs are size dependent
if you would do with resampling this is rather a measure of expression shift and you would not look at the specific DEGs that pop up 
```{r,eval=FALSE}
cao$estimateDEPerCellType(n.cores=50, test = "DESeq2.Wald",  n.resampling=0) 
```

explanation of output column names:

    baseMean: mean of normalized counts for all samples
    log2FoldChange: log2 fold change
    lfcSE: standard error
    stat: Wald statistic
    pvalue: Wald test p-value
    padj: BH adjusted p-values


save all DE genes, not only with padj < 0.05
```{r}
DEG <- cao$test.results$de %>% lapply("[[", "res")
names(DEG)[names(DEG) == 'PIC_enriched_immature_principal_neurons'] <- 'PIC_immature_principal_neu'
#DEG_sig <- DEG[lapply(DEG, function(x) {x$padj < 0.05})]
write.xlsx(x = DEG, file= "DEG_all_genes_embryonic_PIC.xlsx")
```

## Figure 33

transform DEG list into one big dataframe

```{r}
big_df <- bind_rows(DEG, .id = "celltype")
#rename the celltypes
big_df %<>% 
  mutate(celltype = case_when(
    celltype == "CGE_IPC" ~ "CGE_progenitors",
    celltype ==  "Immature_GABAergic_neurons" ~ "Immature_GABAergic_N",
    celltype == "Highly_proliferative_RG_S"  ~ "Highly_prolif_RG_S",
    celltype == "Immature_principal_neurons" ~ "Immature_PN",
    celltype == "Immature_principal_neurons_Tbr1" ~ "Immature_PN_Tbr1",
    celltype == "MGE_IPC"  ~ "MGE_progenitors",
    celltype == "PIC_enriched_immature_principal_neurons"  ~ "PIC_immature_PN",
    celltype == "Highly_proliferative_RG_G2M"  ~ "Highly_prolif_RG_G2M",
    TRUE ~ celltype  # Keep other values unchanged
  ))
```

count significant DEGs
```{r}
#upregulated
sum(big_df$padj < 0.05 & big_df$log2FoldChange > 0)
#downregulated
sum(big_df$padj < 0.05 & big_df$log2FoldChange < 0)
```

```{r}
# genes to label
deg <- c("Bcl11b", "Bcl11a", "Negr1", "Pcdh10", "Gria", "Tlr1", "Il1rapl1", "Sp8", "Nfib", "Pcdh17", "Cux2", "Nfia", "Kcnh7", "Dscaml1", "Myt1l", "Gpc6", "Ncam2", "Dpyd", "Mctp1", "Cacna1e", "Irs1", "Trpm3", "Wnt7b", "Robo1", "Epha5","Syt7", "Ptn", "Nckap5", "Entpd1", "Kcnk10", "Plxna2", "Nrxn3", "Sall1","Sema3e", "Ttyh1", "Slc2a3", "Nfix", "Siah3", "Dnmt3a", "Csmd1","Creb5", "Ntrk2")
```


```{r, fig.height=9, fig.width=7}
p<- ggplot(data=big_df, aes(x = Gene, y=log2FoldChange)) + geom_point(
    position = position_jitter(width = 0.01),
    size = 0.7,
    alpha = 0.5, color="grey")+geom_point(
    data = big_df %>% filter(padj < 0.05 & log2FoldChange > 1),
    position = position_jitter(width = 0.01),
    size = 0.7,
    alpha = 0.7, color= "red"
  ) + geom_point(
    data = big_df %>% filter(padj < 0.05 & log2FoldChange < (-1)),
    position = position_jitter(width = 0.01),
    size = 0.7,
    alpha = 0.7, color="blue"
  ) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x.bottom = element_blank() ) +theme_minimal()+
  theme(axis.ticks.y = element_line())+
  geom_hline(yintercept = c(-1, 1), linetype = "dashed",linewidth=0.3) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
```


```{r, fig.height=9, fig.width=6}
d<- p+  geom_label_repel(
    data = big_df[big_df$Gene %in% deg & big_df$padj < 0.05 & abs(big_df$log2FoldChange) > 1, ],aes(label = Gene),max.overlaps = 100,
    box.padding = 0.5,label.padding = 0.2,min.segment.length = 0, size=2)+ facet_wrap(.~celltype, nrow = 2,strip.position = "bottom") +theme(legend.position = "none")+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x.bottom = element_blank() )+
  theme(strip.text.x = element_text(angle = 90, hjust=1), panel.spacing = unit(0.1, "lines"))
d
rasterize(d, layers='Point', dpi=1000)
ggsave("DEG_manhatten.pdf",width=5.0, height = 9)
```

#7. GSEA

```{r}
cao$estimateOntology(type="GSEA", org.db=org.Mm.eg.db::org.Mm.eg.db)
```
```{r}
cao$saveOntologyAsTable("GSEA.tsv", name="GSEA", p.adj=0.05)#saves only sig terms
GSEA <- read.table("GSEA.tsv", sep="\t") 
```

save GSEA as excel file
```{r}
cell_types <- levels(anno)
GSEA_list <- lapply(cell_types, {function(x) GSEA[GSEA$V1 == x, ]})
names(GSEA_list) <- cell_types
names(GSEA_list)[names(GSEA_list) == 'PIC_enriched_immature_principal_neurons'] <- 'PIC_immature_principal_neu'
write.xlsx(x = GSEA_list, file= "GSEA_embryonic_PIC.xlsx")#headers missing
```

```{r}
celltypes <- c("CGE_IPC","Cortical_hem","Dorsal_IPC", "Highly_proliferative_RG_S","MGE_IPC","Dorsal_RG",                      "Anteriomedial_RG","Lateral_RG","Roof_plate","Highly_proliferative_RG_G2M","Cajal_Retzius_cells" , "Immature_GABAergic_neurons" ,"Immature_principal_neurons"   ,           "Immature_principal_neurons_Tbr1", "PIC_enriched_immature_principal_neurons")
```

## down 
### Figure 34A: top 10
```{r}
cao$plotOntologyHeatmap( name="GSEA", genes="down",  subtype = "BP", cluster=T, top.n = 10)+theme(axis.text.x.top = element_text( vjust = 0.5))

ggsave("GSEA_top10_down.pdf",width=6, height = 2.8)
```

### Figure 34B: other interesting terms
```{r, fig.height=6}
cao$plotOntologyHeatmap( name="GSEA", genes="down",  subtype = "BP", cluster=T, top.n = 25, description.regex="interleukin|Wnt| migration|neuro|neural| differentiation| pallium| development|brain|development|social |behavior|proliferation|immune|telencephalon|inflammat|cytokine|response|plexin", description.exclude.regex="muscle|lung|skin|cartilage|vasclua|eye|mesenchym|hair|peripheral|connective|skeletal|epidermal|epidermis|sensory|astrocyte|epithelial|glial")+theme(axis.text.x.top = element_text( vjust = 0.5))
ggsave("GSEA_25_down.pdf",width=6.8, height = 5.0)
```

## up
### Figure 35: top 20
```{r}
cao$plotOntologyHeatmap( name="GSEA", genes="up",  subtype = "BP", cluster=T, top.n = 20,readjust.p = T)+theme(axis.text.x.top = element_text( vjust = 0.5))
ggsave("GSEA_top20_up.pdf",width=6.5, height = 4.3)
```