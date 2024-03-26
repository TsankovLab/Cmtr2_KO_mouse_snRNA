library (Seurat)
library (ggplot2)
library (patchwork)
library (tidyr)
library (plyr)
library (tidyverse)
library (rstatix)

# Load functions
# Generate barplots or boxplots of meta groups proportions specified
# meta_groups: vector or meta group names:
# 1) first meta_group is the group on which proportions are calculated
# 2) second meta_group split first meta_group on x axes  
# 3) third meta_group will group barplots separately
# if splits include only one value runs barplot instead
cellComp = function (
  seurat_obj = NULL, 
  metaGroups = NULL, # vector of at least 3 metaGroups e.g. c('orig.ident','celltypes','celltypes'),
  plot_as = 'box', # box or bar 
  pal = NULL,
  prop = TRUE,
  ptable_factor = 1, # specify which column of the data.frame or seurat object metadata should be used to compute proportions
  facet_ncol = 20,
  facet_scales = 'free',
  subset_prop = NULL, # subset prop table by any group in any column
  removeNA = TRUE,
  returnDF = FALSE
  )
  {
  require (ggplot2) 
  require (ggpubr)  
  if (is.data.frame (seurat_obj))
    {
    meta_groups_df = seurat_obj[,metaGroups]  
    } else {
    meta_groups_df = seurat_obj@meta.data[,metaGroups]
    }
  # Refactor to remove 0 groups
  #meta_groups_df =  as.data.frame(lapply(unclass(meta_groups_df),as.character),stringsAsFactors=T)
  if(is.null(pal)) pal = rainbow (length(unique(meta_groups_df[,2])))
  #if(is.null(pal) & plot_as == 'box') pal = rainbow (length(unique(meta_groups_df[,3])))
  if (prop)
    {
    ccomp_df = as.data.frame (prop.table (table (meta_groups_df),ptable_factor))
    ccomp_df = na.omit (ccomp_df) # this is to remove NaN somehow produced from the line above 
    } else {
    ccomp_df = as.data.frame (table (meta_groups_df)) 
    }

  if(removeNA) ccomp_df = ccomp_df[ccomp_df$Freq != 0, ] # remove 0s from proportions
  if (!is.null (subset_prop)) 
    {
    subset_col = unlist(sapply (seq(ncol(ccomp_df)), function(x) if(any(ccomp_df[,x] %in% subset_prop)) colnames(ccomp_df)[x]))
    ccomp_df = ccomp_df[ccomp_df[,subset_col] %in% subset_prop,]
    }
  #colnames (ccomp_df) = c(paste0('Var_',seq_along(metaGroups)), 'proportion')  
  if (plot_as == 'box')
    {
    p = ggplot (ccomp_df, aes_string (x= metaGroups[2], y= 'Freq')) +
        theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_fill_manual (values= pal) + xlab (metaGroups[2]) + ylab (ifelse (prop, 'proportion','counts'))
    if (length(metaGroups > 2)) p = p + geom_boxplot(aes_string (fill= metaGroups[3]), outlier.size=.2, alpha = 0.7, lwd=.2) 
    else p = p + geom_boxplot(aes_string (fill= metaGroups[2]), outlier.size=.2, alpha = 0.7, lwd=.2)   
    if (length(metaGroups) > 3) p = p + facet_wrap (as.formula(paste("~", metaGroups[4])), scales=facet_scales, ncol=facet_ncol)
    } 
  if (plot_as == 'bar')
    {
    p = ggplot (ccomp_df, aes_string (x= metaGroups[1], y= 'Freq')) +
        geom_bar(position="stack", stat="identity", aes_string(fill= metaGroups[2])) +
        #geom_bar(position="dodge", stat="identity", aes_string(fill= metaGroups[2])) +
        theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        scale_fill_manual (values= pal) + xlab (metaGroups[2]) + ylab (ifelse (prop, 'proportion','counts'))
    if (length(metaGroups) == 3) p = p + facet_wrap (as.formula(paste("~", metaGroups[3])), scales=facet_scales, ncol=facet_ncol)
    }
  if (returnDF) return(ccomp_df) else 
  return (p)
  }
  
  
### START ANALYSIS ####


## specify folder, metadata project name and other variables
proj_name = 'Alena_snRNA'
projdir = paste0("/ahg/regevdata/projects/ICA_Lung/Bruno/Meyerson_scRNA/",proj_name,"_analysis/")
meta = read.csv (paste0(projdir,'metadata.csv'), row.names=NULL)
projdir = paste0(projdir,'test')
dir.create (paste0(projdir,'/Plots'), recursive=T)

#meta = meta[]
samples_path = sapply (meta$sampleID, function(x)
	paste0('/ahg/regevdata/projects/ICA_Lung/Bruno/Meyerson_scRNA/raw_data/',x,'_raw_feature_bc_matrix'))
meta$sample_path = samples_path
# Set project folder

org = 'mouse'
samples_path = samples_path

nFeat = 400 # Number of features per cells. default 400
nCounts = 800 # Number of UMI per cell. Default 800
pchM = 25 # Percent mitochondrial genes. Default 25 

variablefeatures = 'seurat'
nfeat = 2000 # number of variable genes to consider for dimentionality reduction
sigPCs = 15
vars_to_regress = NULL
ccRegress = FALSE # Regress cell cycle gene expression 
metaGroupNames = c('Genotype','sampleID','Litter.ID','Sex')
res = c(0.2, 0.8, 2, 3, 5) # denovo cluster resolutions 

# Import count matrices
scrna_mats = sapply (meta$sample_path, function (dir) 
    {
    if (grepl('.h5',dir)) 
      {
      Seurat::Read10X_h5 (dir, use.names = TRUE, unique.features = TRUE)    
      } else {
      Seurat::Read10X (data.dir = dir)
      }
    })

srt = sapply (seq_along (scrna_mats), function(i) CreateSeuratObject (counts = scrna_mats[[i]], min.cells = 0, min.features = nFeat, project = proj_name))  
names (srt) = meta$sampleID

###--- Hard filtering based on nfeat and UMI and p.mito ---###')
if (org == 'mouse') for (i in seq_along (srt)) srt[[i]]$percent.mt = PercentageFeatureSet (srt[[i]], pattern = "^mt-")
  
# Filter data
cell_idx = sapply (lapply (srt, function(x) x$nFeature_RNA > nFeat #& nFeature_RNA < 6000 
  & x$percent.mt < pchM & x$nCount_RNA > nCounts), sum)
#write.table (cell_idx, paste0(projdirF0,'cell_idx.csv'))
srt = lapply (srt[cell_idx > 0], function (x) x[,x$nFeature_RNA > nFeat
  & x$percent.mt < pchM & x$nCount_RNA > nCounts])

set.seed (1234)
#*- Data Processing step -*


message ('Merge samples in one seurat object and map metadata')
srtM = merge (srt[[1]], y = srt[2:length(srt)],
              add.cell.ids = names (srt), project = proj_name)        
srtM$sampleID = rep (meta$sampleID, sapply (srt, ncol)) # attach sampleIDs to merged object
srtM$orig.ident = rep (meta$sampleID, sapply (srt, ncol)) # attach sampleIDs to merged object
srt = srtM

  # Add additional metadata to merged obj ####
  if (exists('meta'))
    {
    if (ncol (meta) > 1)
      {
      for (x in colnames (meta[,colnames(meta) != 'sampleID', drop=F])) 
        { 
        metaGroup = meta[,x]
        names(metaGroup) = meta$sampleID
        srt@meta.data[,x] = metaGroup[srt$sampleID]
        }
      }
    } 

# Process merged data
srt = NormalizeData (object = srt, normalization.method = "LogNormalize", scale.factor = 10000)

if (variablefeatures == 'seurat') srt = FindVariableFeatures (srt, selection.method = "vst", nfeat = nfeat)
# Find variable features using scran (works bettern than FindVariableFeatures with vst as this finds very lowly expressed vf)
  
###--- Cell cycle scoring ---###
if (org == 'mouse')
  {
  s.genes = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/cellcycle_mouse.s.genes.rds')
  g2m.genes = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/cellcycle_mouse.g2m.genes.rds')
  }
# if (org == 'human')
#   {
#   cc.genes <- readLines('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/regev_lab_cell_cycle_genes.txt')
#   s.genes <- cc.genes[1:43]
#   g2m.genes <- cc.genes[44:97]
#   }

srt = CellCycleScoring (object = srt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
srt$cc = srt$S.Score + srt$G2M.Score


if (!is.null(vars_to_regress)) srt = ScaleData (srt, features = VariableFeatures (object=srt), vars.to.regress=vars_to_regress) else 
srt = ScaleData (srt, features = VariableFeatures (object=srt))  
srt = RunPCA (srt, features = VariableFeatures (object = srt), npcs = ifelse(ncol(srt) <= 30,ncol(srt)-1,30), ndims.print = 1:5, nfeat.print = 5, verbose = FALSE)

# Plot distribution of nFeature and nCounts
# metaGroups = srt@meta.data[,metaGroupNames, drop=F]
# if (!any (colnames(metaGroups) == 'sampleID')) cbind (data.frame (sampleID = srt$sampleID), metaGroups)

ccomp_df = as.data.frame (table(srt$sampleID))
cc_p2 = ggplot (ccomp_df, aes (x= Var1, y= Freq)) +
        geom_bar (position="stack", stat="identity") +
        theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=1)) + 
        ggtitle (paste('Tot cells',ncol(srt))) + theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
srt$nFeature_RNAL = log10 (srt$nFeature_RNA)
srt$nCount_RNAL = log10 (srt$nCount_RNA)
vln_p = VlnPlot (srt, features = c("nFeature_RNAL", "nCount_RNAL", "percent.mt"), combine=T,group.by = 'sampleID',pt.size = 0, ncol = 3)
png (paste0(projdir, "/Plots/QC_nFeat_nCount_m.percent_vlnPlot.png"), 2800, 1000, res=300)
print (cc_p2 | vln_p[[1]] | vln_p[[2]] | vln_p[[3]]) + plot_layout (widths=c(1,2,2,2))
dev.off()


# Set variables for UMAP
reductionName = 'umap'
reductionSave = 'pca'
reductionGraphKnn = 'RNA_knn'
reductionGraphSnn = 'RNA_snn' 

srt = RunUMAP (object = srt, reduction = reductionSave, dims = 1:sigPCs)

# UMAP of non corrected and harmony corrected clusters
humap_p = lapply (metaGroupNames, function(y) DimPlot (object = srt, reduction = reductionName, pt.size = .01, group.by = y) + theme_classic())

png (paste0(projdir,'/Plots/',paste(metaGroupNames, collapse='_'),'_umap.png'), width = 2500, height = 1200, pointsize=10, res = 300, type="cairo")
print (wrap_plots (humap_p), ncol=3)
dev.off()
  
# Color UMAP by number of detected genes and color from 5% to 95% of counts
qc_metrics = c('nCount_RNA', 'nFeature_RNA', 'percent.mt')
umap_df = data.frame (srt[[reductionName]]@cell.embeddings, nCount_RNA = log10(srt$nCount_RNA+1), nFeature_RNA = log10(srt$nFeature_RNA+1), percent.mt = srt$percent.mt)
umap_p = lapply (qc_metrics, function(x) ggplot(data = umap_df) + 
geom_point (mapping = aes_string (x = colnames(umap_df)[1], y=colnames(umap_df)[2], color = x), size = .01) + 
scale_colour_gradientn (colours = rainbow(7)) +
theme_classic() + 
theme(
   plot.background = element_blank()
  ))

png (paste0(projdir,'/Plots/QC_umap.png'), width = 3500, height = 1200, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap_p), ncol=3)
dev.off()
  
# Run denovo clustering on non-adjusted reductions
srt = FindNeighbors (object = srt, reduction = reductionSave, dims = 1:sigPCs, k.param = 30,
                              verbose = TRUE, force.recalc = T, graph.name=c(reductionGraphKnn,reductionGraphSnn))
for (i in seq_along(res)) srt = FindClusters (srt, resolution = res[i], verbose = T, n.start = 100, graph.name=reductionGraphSnn)
  #srt = AddMetaData (srt, metadata= srt$seurat_clusters, col.name = paste0('seurat_clusters_',res[i]))

clust_p1 = list()
for (i in seq_along(res)) clust_p1[[i]] = DimPlot (srt, pt.size = .01, label = T, group.by= paste0(reductionGraphSnn,'_res.',res[i]), reduction = reductionName) + NoLegend()

  
png (paste0(projdir,'/Plots/denovo_clusters_',reductionSave,'_umaps.png'), 3000, 3000, pointsize=10, res = 300, type="cairo")
print (wrap_plots (clust_p1, ncol=2))
dev.off() 
  

# set additional column metadata
srt$Genotype2 = ifelse (srt$Genotype == 'Cmtr2 WT','WT','KO')
ccomp_df = as.data.frame (table(srt$sampleID))
ccomp_df$sampleID2 = srt$Genotype2 [match (ccomp_df$Var1, srt$sampleID)]
ccomp_df = ccomp_df[order (ccomp_df$sampleID2),]
ccomp_df$sampleID2 = paste0(ccomp_df$sampleID2, c(1,2,3,4,5,1,2,3,4))
srt$sampleID2 = ccomp_df$sampleID2[factor(srt$sampleID, levels = unique(ccomp_df$Var1))]

genotype_pal = c(WT = 'lightyellow',KO='orchid')
library (circlize)
library (paletteer)
celltype_pal_fun = colorRampPalette (paletteer_d("ggthemes::Summer"))
celltype_pal = celltype_pal_fun(26)



###--- Use JShendure mouse dev data to annotate cell types ---# 
### Annotate cell clusters using label lransfer from JShendure data and Marioni after subsampling ###
JshendureDir = '/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/JShendure_mouse_dev/'
srt_js = readRDS (paste0(JshendureDir,'srt_JShendure_mouse_dev_E9.5.rds'))
srt_js$Main_cell_type[srt_js$Main_cell_type == 'Chondroctye progenitors'] = 'Chondrocyte progenitors'
# marioniDir = '/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/'
# srt_mg = readRDS (paste0(marioniDir,'srt_marioni_gastrulation.rds'))
srt_refs = list (srt_js)#, srt_mg)  
names (srt_refs) = c('JShendure')#,'Marioni')
metaGroupNames = c('Main_cell_type')#,'celltype')
subsample = Inf
set.seed (123)
small_ann = 10
for (i in seq_along (srt_refs))
  {
  ref = srt_refs[[i]]
  metaGroup = ref@meta.data[,metaGroupNames[i]]
  metaGroup[is.na(metaGroup)] = 'not_assigned' # rename NA  
  sub_cells = lapply (unique(metaGroup), function(x) {
    sub_celltype = colnames(ref)[metaGroup == x]
    if (length(sub_celltype) < subsample) {
      sub_celltype } else {
    sub_celltype = sample (sub_celltype, subsample)
    }})
  sub_cells = unlist (sub_cells)
  ref = ref[,sub_cells[sub_cells != 'not_assigned']] # remove NA

  metaGroup = ref@meta.data[,metaGroupNames[i]]
  query = srt
  DefaultAssay (query) = 'RNA'
  anchors = FindTransferAnchors (reference = ref, query = query, dims = 1:30)
  predictions = TransferData (anchorset = anchors, refdata = metaGroup,
  dims = 1:30)
  small_ann2 = table (predictions$predicted.id) < small_ann
  predictions$predicted.id[predictions$predicted.id %in% names(which(small_ann2))] = 'Others'
  write.csv (predictions, paste0(projdir,'/prediction_scores_',names (srt_refs)[i],'_subsampled_',subsample,'.csv'))
  srt@meta.data [,paste0('predicted_',names (srt_refs)[i], '_',subsample)] = predictions$predicted.id
  }

#saveRDS (srt, paste0(projdir, 'srt.rds'))

for (i in seq_along (srt_refs))
  {
  predictions = read.csv (paste0 (projdir,'/prediction_scores_',names (srt_refs)[i],'_subsampled_',subsample,'.csv'))
  umap_df = data.frame (srt[[reductionName]]@cell.embeddings, prediction = srt@meta.data[,paste0('predicted_',names (srt_refs)[i], '_',subsample)])
    direction_col = umap_df %>% group_by(prediction) %>% summarise (mean (UMAP_1))
    umap_df$prediction = factor (umap_df$prediction, levels = direction_col$prediction[order(-direction_col[,2])])
    labeltransfer_p1 = ggplot(data = umap_df) + 
    geom_point (mapping = aes_string (x = colnames(umap_df)[1], y= colnames(umap_df)[2], color = 'prediction'), size = .1) + 
    #scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    ggtitle (names(srt_refs)[1]) + 
    theme_void() + 
    scale_color_manual (values = setNames(celltype_pal, levels(umap_df$prediction))) +
    theme(legend.key.size = unit(0.2, "cm")) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme(legend.position="bottom")
    umap_df = data.frame (srt[[reductionName]]@cell.embeddings, Genotype = srt@meta.data[,'Genotype2'])
  labeltransfer_p2 = ggplot(data = umap_df) + 
    geom_point (mapping = aes_string (x = colnames(umap_df)[1], y= colnames(umap_df)[2], fill = 'Genotype'),color = 'black',  size = 2, shape=21, stroke = 0.1) + 
    #scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    ggtitle (names(srt_refs)[1]) + 
    theme_void() +
    scale_fill_manual (values= genotype_pal)
    
  #labeltransfer_p1 = DimPlot (srt, pt.size = .05, label = F, group.by= paste0('predicted_',names (srt_refs)[i],'_',subsample), reduction = reductionName) 
  
  pred_mat = data.frame (predicted.id = predictions$predicted.id, score = predictions$prediction.score.max)
  predictions_p = ggplot (pred_mat, aes (x= predicted.id, y = score)) +
  geom_boxplot () + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  png (paste0(projdir,'/Plots/Label_tranfer_',names(srt_refs)[i],'_subsampled_',subsample,'UMAP.png'), width=7500, height=3000,res=500)
  print (wrap_plots (labeltransfer_p2))#, predictions_p, ncol=3))
  dev.off() 
  }

# Plot both label transfers UMAPs
metaGroupNames = c('predicted_JShendure_Inf','RNA_snn_res.0.8', 'Genotype')
labeltransfer_p1 = lapply (metaGroupNames, function(x) DimPlot (srt, pt.size = .05, label = T, group.by= x, reduction = reductionName))

png (paste0(projdir,'/Plots/Cell_annotations_',paste(metaGroupNames, collapse='-'),'_umaps.png'), width=3500, height=2000,res=100)
print (wrap_plots (labeltransfer_p1, ncol=2))
dev.off()

### Cell composition analysis
metaGroupName = 'predicted_JShendure_Inf' 
metaGroupName2 = 'Genotype2'
# metaGroupName3 = 'predicted_Marioni'
cc_df = srt@meta.data[srt@meta.data[,metaGroupName] != 'Others',]
cc_df[,metaGroupName] = factor (cc_df[,metaGroupName], levels =names(table (cc_df[,metaGroupName])[order (-table (cc_df[,metaGroupName]))]))
cc_box1 = cellComp (
  seurat_obj = cc_df, 
  metaGroups = c('sampleID', metaGroupName, 'Genotype2'),
  plot_as = 'box',
  ptable_factor = c(1,3),
  pal = genotype_pal,
  ) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# cc_box2 = cellComp (
#   seurat_obj = srt, 
#   metaGroups = c(metaGroupName, metaGroupName1),
#   plot_as = 'bar',
#   ptable_factor = c(1),
#   #pal = viridis::mako (length(unique(srt@meta.data[,'Genotype2']))),
#   ) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   theme(legend.key.size = unit(.2, 'cm'))

# cc_box3 = cellComp (
#   seurat_obj = srt, 
#   metaGroups = c(metaGroupName, metaGroupName3),
#   plot_as = 'bar',
#   ptable_factor = c(1),
#   #pal = viridis::mako (length(unique(srt@meta.data[,'Genotype2']))),
#   ) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   theme(legend.key.size = unit(.2, 'cm'))


cc_df2 = cc_box1$data
cc_df2 = cc_df2[!cc_df2$predicted_JShendure_Inf %in% c('Premature oligodendrocyte','Inhibitory neuron progenitors'),]
stat.test = cc_df2 %>%
  group_by (predicted_JShendure_Inf) %>%
  t_test(Freq ~ Genotype2) %>%
  adjust_pvalue(method = "none") %>%
  add_significance()
stat.test = stat.test %>% add_xy_position (x = 'predicted_JShendure_Inf', step.increase=0.1)

cc_box1 = cc_box1 + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
   bracket.nudge.y = 0, hide.ns = TRUE,
    label = "p.adj.signif")

png (paste0(projdir,'/Plots/cell_composition_',metaGroupName,'2.png'), width=3500, height=2000, res=500)
cc_box1
dev.off()

# use speckle package to test for cell composition differences ####
library (speckle)
library (limma)

cell_info <- speckle_example_data()
head(cell_info)

metaGroupNames = c('sampleID','predicted_JShendure_Inf','Genotype2')

ccomp_df = srt@meta.data [,metaGroupNames]
cc_res = propeller (clusters = ccomp_df[,metaGroupNames[2]], sample = ccomp_df[,metaGroupNames[1]], group = ccomp_df[,metaGroupNames[3]])
write.csv (cc_res, 'cell_composition_analysis.csv')




### Compute proliferating cells fraction between genotypes
metaGroupName1 = 'predicted_JShendure_Inf'
metaGroupName2 = 'sampleID2'
metaGroupName3 = 'Genotype2'

pdf (paste0('Plots/cycling_distribution.pdf'))
hist (srt$S.Score)
hist (srt$G2M.Score)
dev.off()
srt$S.Score2 = ifelse (srt$S.Score > .2, 'cycling','G1')
srt$G2M.Score2 = ifelse (srt$G2M.Score > .2, 'cycling','G1')
srt$Phase2 = ifelse (srt$S.Score2 != 'G1' & srt$G2M.Score2 != 'G1', 'cycling','G1')
metaGroupNames = c(metaGroupName2, 'Phase2', metaGroupName3, metaGroupName1)
ccc_box1 = cellComp (
  seurat_obj = srt, 
  metaGroups = metaGroupNames,
  plot_as = 'box',
#  pal = pal1,
  prop = TRUE,
  ptable_factor = c(1,4),
  subset_prop = 'cycling',
  facet_ncol = 6,
  facet_scales = 'free'
  ) + theme_classic()#+ NoLegend() + ggtitle (i)

metaGroupNames = c(metaGroupName2, 'Phase2', metaGroupName1)
ccc_bar1 = cellComp (
  seurat_obj = srt, 
  metaGroups = metaGroupNames,
  plot_as = 'bar',
  pal = c(cycling = 'blue', G1= 'grey'),
  prop = TRUE,
  ptable_factor = c(1,3),
  #subset_prop = 'cycling',
  facet_ncol = 6,
  facet_scales = 'free'
  ) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+ NoLegend() + ggtitle (i)


metaGroupNames = c(metaGroupName2, 'Phase2', metaGroupName3)
metaGroups = metaGroupNames
meta_groups_df = srt@meta.data[,metaGroupNames]
pal = genotype_pal
ccomp_df = as.data.frame (prop.table (table (meta_groups_df),1))
ccomp_df = na.omit (ccomp_df) # this is to remove NaN somehow produced from the line above 
ccomp_df = ccomp_df[ccomp_df$Freq != 0, ] # remove 0s from proportions
ccomp_df = ccomp_df[ccomp_df$Phase2 == 'cycling',]
ccomp_df$Phase2 = as.factor (as.character(ccomp_df$Phase2))
#comp_df$Var_2 = as.factor (as.character(comp_df$Phase2))
stat.test2 = ccomp_df %>%
  group_by (Phase2) %>%
  t_test(Freq ~ Genotype2) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
ccc_box2 = ggplot (ccomp_df, aes (x= Phase2, y= Freq)) +
        geom_boxplot(aes (fill= Genotype2), outlier.size=.2)  + 
        theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_fill_manual (values= pal) + xlab (metaGroups[2]) + ylab ('proportion') + 
        theme_classic() + 
        geom_point(aes_string(fill = 'Genotype2'), size = 2, shape = 21, position = position_jitterdodge()) +
        ylim (c(0, 0.008))

stat.test2 =  stat.test2 %>% add_xy_position (x = "Phase2", step.increase=0.001,dodge=1)
ccc_box2 = ccc_box2 + stat_pvalue_manual (stat.test2, remove.bracket=F,
   bracket.nudge.y = .0008, hide.ns = F,
    label = "p.adj") 

png(paste0 ('Plots/cc_genotype_composition2.png'), width=1500, height=1300, res=500)
ccc_box2
dev.off()

