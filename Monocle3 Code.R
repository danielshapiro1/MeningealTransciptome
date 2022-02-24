#### TRAJECTORY ####
#### Library and data loading ----
library(Seurat)
library(patchwork)
library(readr)
# library(scCATCH)
library(SingleR)
library(tidyverse)
library(monocle3)
# library(SeuratData)
library(magrittr)
# library(colourpicker)
# library(rDGIdb)

ashleytbi <- readRDS(file = "/Users/danielshapiro/Desktop/meninges.all.names.080320.rds")
ashleytbi <- UpdateSeuratObject(ashleytbi)

#### Prelim Analysis of populations in each conditoon ####
prop.table(table(Idents(ashleytbi))) # show cell pop by type
table(Idents(ashleytbi), ashleytbi$OriginalClusters) # show cell pop in each cluster
table(Idents(ashleytbi), ashleytbi@meta.data[["orig.ident"]]) # show cell pop in each treatment

# lets just pull out the treatment into its own slot
ashleytbi$treatment <- ashleytbi@meta.data[["orig.ident"]]

#### SEPARATE TREATMENT GROUPS ####
tbi.SHAM <- subset(ashleytbi, subset = treatment == c("sham.meninges"))
tbi.REAL <- subset(ashleytbi, subset = treatment == c("TBI.meninges"))

prop.table(table(Idents(tbi.SHAM), tbi.SHAM$treatment)) # show cell pop in each treatment
prop.table(table(Idents(tbi.REAL), tbi.REAL$treatment)) # show cell pop in each treatment

#### TOTAL Trajectory ----
#### Seurat -> Monocle ----
# in previous versions we tried the seurat wrapper it just didnt work
# below we manually wrap the data ourselves

# convert to monocle cds object 
# Extract data, phenotype data, and feature data from the SeuratObject
expressiondata <- ashleytbi@assays[["RNA"]]@data

cellmd <- ashleytbi@meta.data

genemd <- data.frame(gene_short_name = row.names(expressiondata), 
                     row.names = row.names(expressiondata))

# Construct monocle cds
ashleytbi.cds <- new_cell_data_set(expression_data = expressiondata,
                                  cell_metadata = cellmd,
                                  gene_metadata = genemd)
ashleytbi.cds <- preprocess_cds(ashleytbi.cds, num_dim = 30) # we used 30 in earlier seurat scripts

# 
# run clustering again (didnt transfer from seurat)
ashleytbi.cds <- reduce_dimension(ashleytbi.cds, max_components = 3, reduction_method = "UMAP")
ashleytbi.cds <- cluster_cells(ashleytbi.cds, max_components = 3, reduction_method = "UMAP")

#### optional: Transfer Seurat embeddings and plot #####
# Note that these may be calculated on the Integrated object, not the counts
# and thus will involve fewer genes
temp.cds <- ProjectDim(ashleytbi, reduction = "pca") # this will be removed
reducedDim(ashleytbi.cds, type = "PCA") <- temp.cds@reductions$pca@cell.embeddings 
ashleytbi.cds@preprocess_aux$prop_var_expl <- temp.cds@reductions$pca@stdev
plot_pc_variance_explained(ashleytbi.cds)

# Transfer Seurat UMAP embeddings
ashleytbi.cds@int_colData@listData$reducedDims$UMAP <- temp.cds@reductions$umap@cell.embeddings

#### Learning trajectory ----
# now learn the PATH (trajectory)
ashleytbi.cds <- learn_graph(ashleytbi.cds)

# this calls up a shiny app, choose the ROOT NODE
ashleytbi.cds <- order_cells(ashleytbi.cds, reduction_method = "UMAP")

## transfer singleR labels to moncle3 object
colData(ashleytbi.cds)$singleRCalls <- ashleytbi@meta.data[["SingleR.calls"]] # call this by opening the object
colData(ashleytbi.cds)$assigned_cell_type <- ashleytbi@active.ident # call this by opening the object

# finally, you can visualize the learned path
pdf("monocle3UMAP_projectedgraph", width=6, height=6)
plot_cells(ashleytbi.cds,
           color_cells_by = "assigned_cell_type",
           label_groups_by_cluster=F,
           show_trajectory_graph = T,
           trajectory_graph_segment_size = 1,
           label_leaves=F, # this gives a little node label (outcome)
           label_roots = T,
           label_branch_points = F,
           graph_label_size = 1, # size of # in circle
           group_label_size = 3,
           cell_size = 1,
           alpha = 0.7,
           scale_to_range = T) # sync color scheme
dev.off()

# now you can show pseudotime
pdf("Figure_images/monocle3_pseudotime_seuratpartition.pdf", width=7, height=6)
plot_cells(ashleytbi.cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = F,
           trajectory_graph_segment_size = 1,
           label_leaves=F, # this gives a little node label (outcome)
           label_roots = T,
           label_branch_points = F,
           graph_label_size = 1, # size of # in circle
           group_label_size = 3,
           cell_size = 1,
           alpha = 1,
           scale_to_range = T) 
dev.off()


#### Subset Trajectory & analysis of SMC----
ashleytbi_subset.cds <- choose_cells(ashleytbi.cds) # calls up shiny app

plot_cells(ashleytbi_subset.cds,
           color_cells_by = "pseudotime", # cell type or none
           show_trajectory_graph = F,
           trajectory_graph_segment_size = 1,
           label_leaves=F, # this gives a little node label (outcome)
           label_roots = T,
           label_branch_points = F,
           graph_label_size = 1, # size of # in circle
           group_label_size = 3,
           cell_size = 1,
           alpha = 0.7,
           scale_to_range = T) 


### MORAN's I Test of Autocorrelation ###
# now we can extrapolate genes that are differentially expressed in this region
# Moran’s I is a measure of multi-directional and multi-dimensional spatial autocorrelation. 
# the statistic tells you whether cells at nearby positions on a 
# trajectory will have similar (or dissimilar) +
# expression levels for the gene being tested.
## first lets do the whole dataset
# a special gene module score heatmap (for the whole dataset)
pr_graph_test_res <- graph_test(ashleytbi.cds, neighbor_graph="principal_graph", cores=2)
write.csv(pr_graph_test_res, file = "ashleytbi_moransI_in_SMC_subcluster.csv")
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.00000001)) # you can adjust the p-value here
head(pr_deg_ids)
gene_module_df <- find_gene_modules(ashleytbi.cds[pr_deg_ids,], resolution=1e-3)
cell_group_df <- tibble::tibble(cell=row.names(colData(ashleytbi.cds)), 
                                cell_group=colData(ashleytbi.cds)$assigned_cell_type)
agg_mat <- aggregate_gene_expression(ashleytbi.cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

# which then can be visualized like so;
# this can show you the different gene modules that can are responsible for changes over pseudotime
plot_cells(ashleytbi.cds,
           genes=gene_module_df %>% filter(module %in% c(2,3,7)), # specify the module you want to examine
           label_cell_groups=T,
           show_trajectory_graph=F)

subset(gene_module_df, module == 2)

## now lets do the subsets
pr_graph_test_res.sub <- graph_test(ashleytbi.cds_subset, neighbor_graph="principal_graph", cores=2)

pr_deg_ids.sub <- row.names(subset(pr_graph_test_res.sub, q_value < 0.00000001))
head(pr_deg_ids.sub)
# collect the trajectory-variable genes into modules
gene_module_df.sub <- find_gene_modules(ashleytbi.cds_subset[pr_deg_ids.sub,], resolution=1e-3)
# visualize these genes
# here I am just pulling out genes that have high moran's i and might be helpful in the paper
pdf("monocle3_genesoverpseudotime_seuratpartition.pdf", width=7, height=6)
plot_cells(ashleytbi.cds_subset, 
           genes=c("SERPINF1", "FBLN1", "PPP1R14A","MYH11", "RAMP1", "IGFBP2"), # this is faceting by the genes that are DE
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)

dev.off()

pdf("monocle3_genesoverpseudotime_seuratpartition_extended.pdf", width=7, height=6)
plot_cells(ashleytbi.cds_subset, 
           genes=c("MYH11", "RAMP1", 'IGFBP2',
                   "LUM", "TCF21", "C7",
                   "SERPINF1",  "FBLN1", "PPP1R14A",
                   "FN1", "CNN1", "TNFRSF11B"), # this is faceting by the genes that are DE
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)

dev.off()
# recluster at higher definition
ashleytbi.cds_subset = cluster_cells(ashleytbi.cds_subset, resolution=1e-2)

pdf("monocle3_RNAvelocitySUBSET_seuratpartition.pdf", width=6, height=6)
plot_cells(ashleytbi.cds_subset, 
           color_cells_by="cluster",
           label_groups_by_cluster=F,
           show_trajectory_graph = T,
           trajectory_graph_segment_size = 1,
           label_leaves=F, # this gives a little node label (outcome)
           label_roots = F,
           label_branch_points = F,
           graph_label_size = 1, # size of # in circle
           group_label_size = 4,
           cell_size = 1,
           alpha = 0.5,
           scale_to_range = T)
dev.off()

pdf("monocle3_RNAvelocitySUBSET_MODULES_seuratpartition.pdf", width=6, height=6)
plot_cells(ashleytbi.cds_subset,
           genes=gene_module_df %>% filter(module %in% c(2,3)), # specify the module you want to examine
           label_cell_groups=T,
           graph_label_size = 1, # size of # in circle
           group_label_size = 4,
           alpha = 0.5,
           show_trajectory_graph=F)
dev.off()

#### additonal search for cluster markers ----
# first set active idents
Idents(ashleytbi) <- ashleytbi@meta.data[["SingleR.calls"]]
markers_SMC <- FindMarkers(ashleytbi, ident.1 = "SMC") # since we didnt specific ident.2, this shows relative to all other cell types

write.csv(markers_SMC, file = "SMC_specific_DEGs.csv")

features <- c("MYH11", "FN1",  "COL6A1", 'COL6A2', "PPP1R14A",
              "TNFRSF11B", #TOP HEAVY
              "FBLN1","LUM", "TCF21", "C7", "C6", # BOTTOM HEAVY
              "SERPINF1" )

RidgePlot(ashleytbi, features = features, ncol = 1, group.by = "SingleR.calls", log = TRUE, y.max = 20)
VlnPlot(ashleytbi, features = features, ncol = 1, group.by = "SingleR.calls", log = TRUE)
FeaturePlot(ashleytbi, features = features)
DotPlot(ashleytbi, features = features, group.by = "SingleR.calls")

pdf("Figure_images/monocle3_genesoverpseudotime_seuratpartition_extended_v2.pdf", width=7, height=6)
plot_cells(ashleytbi.cds_subset, 
           genes=features, # this is faceting by the genes that are DE
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)

dev.off()


#### output seurat object for future use ####
saveRDS(ashleytbi, file = "final_ashleytbi_labeled.rds")





#### SHAM TRAJECTORY ####

#### TBI TRAJECTORY ####







