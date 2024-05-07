---
title: "Annotation"
output: html_document
date: "2024-05-07"
---

This notebook summarizes the strategy and fucntions used for annotation of clusters.
# Setup
```{r}
#Load helper functions
source("/people/nrq364/embryonic_PolyIC/scRNA_helper_new.R")

library(magrittr)
library(dplyr)
library(conos)
library(ggplot2)
library(shiny)
library(plotly)
```

# Plot marker genes
generally known marker genes can be plotted on UMAP or dotplots to identify the cell types
```{r}
#plot one gene
con$plotGraph(gene="Tcf4")

#plot many genes
genes <- c("Zic1", "Zic4", "Fgf8", "Fgf17", "Il17rd", "Ptprm")
plots <- lapply(genes, function(gene) {
  con$plotGraph(gene = gene, title = gene, size = 0.2)
})
cowplot::plot_grid(plotlist = plots)
```

# Estimate upregulated genes per cluster
estimate upregulated DEGs per cluster (=marker genes per cluster) and then identify the cell type; consider only using the control cells to estimate marker genes per cluster 
```{r}
marker <- con$getDifferentialGenes(groups=anno, append.auc=TRUE, upregulated.only=TRUE)
```

# Annotation strategy
To refine clustering first...
...start with using the clusters estimated by conos... 
```{r}
#save original clusters estimated by leiden algorithm
leiden22 <- con$clusters$leiden$groups %>% factor
qsave(leiden22, "leiden22.qs")

#rename to anno
anno <- leiden22
```

...potentially rerun the clustering for the whole dataset (change resolution and min.group.size)...
```{r, eval=FALSE}
#find new clusters
con$findCommunities(method = leiden.community, resolution = 3, min.group.size = 15)

#plot new clusters in UMAP and as table
con$plotGraph(groups = con$clusters$leiden$groups %>% factor, title="")
table(con$clusters$leiden$groups %>% factor)
```

...or subcluster specific clusters.
```{r,eval=FALSE, echo=F}
#subcluster cluster "18"
anno <- findSubcommunities(con, "18", groups=anno, resolution=1)
anno %<>% as.factor()
levels(anno)
#plot new clusters in UMAP
con$plotGraph(groups = anno, title="", size=0.2)

#plot only specific cluster in UMAP
con$plotGraph(groups=anno, subgroups="18_3", size=0.2)
```


 use helper functions from scRNA_helper.R to rename or collapse clusters
```{r,eval=FALSE, echo=F}
# rename "18_3" to "immune_cells"
anno %<>% renameAnnotation("18_3", "immune_cells")
anno %<>% renameAnnotation("18", "blood_cells")
anno %<>% renameAnnotation("20", "neural_crest")

#collapse all clusters with "18" in name
anno %<>% collapseAnnotation("18")
```

If the upper strategy does not give the desired clusters based on marker genes, manual selection of cells can be done:
to select specific cells and name them the following code can be used, shiny and plotly need to be installed
a new window will open in which desired cells will be manually selected
the cells will be saved as a vector/factor

```{r}
x <- con$plotGraph() 
sd <- x$data


library(shiny)
library(plotly)

ui <- fluidPage(
    plotlyOutput("plot"),
    verbatimTextOutput("click"),
    verbatimTextOutput("brush")
)

server <- function(input, output, session) {
    frame1 <- data.frame()
    frame2 <- data.frame()
    nms <- row.names(mtcars)

    output$plot <- renderPlotly({
        p <- ggplot(sd, aes(x = x, y = y, key = CellName)) + geom_point()
        ggplotly(p) %>% layout(dragmode = "lasso")
    })

    output$click <- renderPrint({
        d <- event_data("plotly_click")
        if (!is.null(d)) {
            frame1 <<- frame1[is.null(frame1$pointNumber), ] # Optional line to remove the previous selections
            frame1 <<- rbind(frame1, d) 
        }
            frame1
        })

    output$brush <- renderPrint({
        d <- event_data("plotly_selected")
        if (!is.null(d)) {
            frame2 <<- frame2[is.null(frame2$pointNumber), ] # Optional line to remove the previous selections 
            frame2 <<- rbind(frame2, d)
        }
            frame2
saveRDS(frame2, file = "immature_principal_neurons.rds") #remember to change this dataframe name once you do a new cell selection
    })

}

shinyApp(ui, server)
```

add selected cells to the annotation factor:
```{r}
#define function to add cell selection from shiney app to anno
add_selection_to_anno <- function(x, name) { #name should be in "
  x <- readRDS(paste0(name, ".rds"))
  levels(anno) <- c(levels(anno), name) # to add a new level to anno
  anno[x$key] <- name
 return(anno)
}

#use like this: 
#first argument is the name of the file you defined in the code above with the shiny app, second argument is the name of the cluster in your annotation
add_selection_to_anno(immature_principal_neurons, "immature_principal_neurons")
```