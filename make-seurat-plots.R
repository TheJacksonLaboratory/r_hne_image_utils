# Pipeline to integrate H&E, CNN imaging features, and HoVer-Net cell counts, 
# using Seurat tools.

#install.packages('Seurat')
suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(Seurat))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(patchwork))
suppressPackageStartupMessages(p_load(data.table))  
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(reticulate))
suppressPackageStartupMessages(library("optparse"))

# Rscript ./seurat-plotting-utils.R --wsi-file=/Users/whitebr/work/jax/pdxnet/r_hne_image_utils/test/Images/9468.svs --hovernet-thumbnail-file=/Users/whitebr/work/jax/pdxnet/r_hne_image_utils/test/Hovernet/9468_wsi/thumb/9468.png --inception-feature-file=/Users/whitebr/work/jax/pdxnet/r_hne_image_utils/test/Inception/Features/features_9468.svs.csv --hovernet-counts-file=/Users/whitebr/work/jax/pdxnet/r_hne_image_utils/test/Hovernet/tile_512_mask/9468_cellcount.csv 

option_list <- list(
    make_option(c("--wsi-file"), action="store",
                default=NULL,
                help="Path to whole slide image (WSI) in format compatible with python openslide (e.g., SVS, TIF, or NDPI)"),
    make_option(c("--hovernet-thumbnail-file"), action="store",
                default=NULL,
                help="Path to thumbnail file output by Hover-Net (in the thumb/ output directory)."),
    make_option(c("--hovernet-counts-file"), action="store",
                default=NULL,
                help="Path to file containing post-processed, per-tile Hover-Net results (e.g., created by https://cgc.sbgenomics.com/u/pdxnet/pdxnet-image/analysis/cruncher/hovernet-wsi-out#hovernet_tile_tissue_cell_count_persample_workflow.py)"),
    make_option(c("--inception-feature-file"), action="store",
                default=NULL,
                help="Path to CSV file containing per-tile Inception features (e.g., created by run-inception-v3-on-wsi.py)")
#    make_option(c("--hovernet-mask-file"), action="store",
#                default=NULL,
#                help="Path to mask file output by Hover-Net (in the mask/ output directory).")
)

descr <- "\
  Plot Inception-V3 features and HoVer-Net cell counts overlayed on H&E image using Seurat tools.
"

parser <- OptionParser(usage = "%prog [options]", option_list=option_list, description=descr)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

wsi.file <- opt$`wsi-file`
thumbnail.file <- opt$`hovernet-thumbnail-file`
mask.file <- opt$`hovernet-mask-file`
hovernet.cnts.file <- opt$`hovernet-counts-file`
inception.feature.file <- opt$`inception-feature-file`

if(is.null(wsi.file)) { print_help(parser); stop("Must pass WSI file with --wsi-file") }
#if(is.null(mask.file)) { print_help(parser); stop("Must pass Hover-Net mask file with --hovernet-mask-file") }
if(is.null(thumbnail.file)) { print_help(parser); stop("Must pass Hover-Net thumbnail file with --hovernet-thumbail-file") }
if(is.null(hovernet.cnts.file)) { print_help(parser); stop("Must pass Hover-Net post-processed output file with --hovernet-counts-file") }
if(is.null(inception.feature.file)) { print_help(parser); stop("Must pass inception feature file with --inception-feature-file") }

# Read in the hovernet data

# Expected columns in hovernet post-processed output holding the tile coordinates and per-cell type counts for that tile
hovernet.cell.type.cols <- c("nolabel", "neoplastic", "inflammatory", "connective", "necrosis", "non.neoplastic")
hovernet.cols <- c("X", "Y", "XX", "YY", hovernet.cell.type.cols)
# Expected column in hovernet post-processed output holding the fraction of the tile that overlaps the tissue mask
hovernet.mask.fraction.col <- "mask_fraction"

hovernet.tbl <- read.csv(hovernet.cnts.file, header = T, as.is = T)
missing.cols <- hovernet.cols[!(hovernet.cols %in% colnames(hovernet.tbl))]
if(length(missing.cols) > 0) {
  stop(paste0("Missing the following expected columns from hovernet counts file: ", paste0(missing.cols, collapse = ", "), "\n"))
}

if(!(hovernet.mask.fraction.col %in% colnames(hovernet.tbl))) {
  stop(paste0("Was expecting ", hovernet.mask.fraction.col, " in hovernet counts file\n"))
}

hovernet.data <- hovernet.tbl[, hovernet.cols]
hovernet.mask.data <- hovernet.tbl[, hovernet.mask.fraction.col, drop = TRUE]

# Subset the hovernet data to that passing the mask threshold
mask_thresh=1.0
hovernet.data <- hovernet.data[hovernet.mask.data >= mask_thresh,]

# Read in the inception features
inception.features <- as.data.frame(fread(inception.feature.file))

#' Extract dimensions from whole slide image (WSI) file
#'
#' \code{get.wsi.dimensions} returns the x and y dimensions of a WSI.
#'
#' This function uses python's openslide (via reticulate) to
#' extract the dimensions. Hence, the WSI should be in a format
#' compatible with openslide (e.g., SVS, TIF, or NDPI).
#'
#' @param wsi.file A string. _Absolute_ path to the WSI image compatible with
#'   openslide (e.g., SVS, TIF, or NDPI).
#' @return vector with x and y coordinates of WSI
get.wsi.dimensions <- function(wsi.file) {
  openslide <- import("openslide")
  slide <- openslide$OpenSlide(wsi.file)
  dims <- unlist(slide$dimensions)
  return(dims)
}

dims <- get.wsi.dimensions(wsi.file)
wsi_x <- dims[1]
wsi_y <- dims[2]

sample <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(thumbnail.file))

# Load a thumbnail of the H&E-stained image file
hne.image <- png::readPNG(thumbnail.file)

mag_ratio=1

# Name tiles by number: tile_x_y
# Translate from pixels to tile number
tileNames <- paste("tile_", hovernet.data$X/(hovernet.data$XX-hovernet.data$X), "_", hovernet.data$Y/(hovernet.data$YY-hovernet.data$Y), sep="") 
xCoords <- ((hovernet.data$X+hovernet.data$XX)/2)*(ncol(hne.image)/wsi_x)
yCoords <- ((hovernet.data$Y+hovernet.data$YY)/2)*(nrow(hne.image)/wsi_y)
coordinates <- data.frame(cbind(xCoords,yCoords))
names(coordinates) <- c("imagecol","imagerow")
row.names(coordinates) <- tileNames
row.names(hovernet.data) <- tileNames

# Remove unused hovNet columns, add total cell count, Tumor/Stroma Fraction
hovernet.data <- hovernet.data[, hovernet.cell.type.cols]
hovernet.data$Total <- apply(hovernet.data,1,sum)
hovernet.data <- hovernet.data %>% rowwise() %>% mutate(TumorFraction = (neoplastic)/(Total+0.01)) 
hovernet.data <- hovernet.data %>% rowwise() %>% mutate(StromaFraction = (connective)/(Total+0.01))
hovernet.data <- hovernet.data %>% rowwise() %>% mutate(NecroticFraction = (necrosis)/(Total+0.01))
hovernet.data <- as.data.frame(hovernet.data)
row.names(hovernet.data) <- tileNames

# Remove pixel coordinates from inception.features
# Translate from pixels to tile number
tileNames <- paste("tile_",inception.features$minx/(inception.features$maxx-inception.features$minx), "_", inception.features$miny/(inception.features$maxy-inception.features$miny), sep="") 
row.names(inception.features) <- tileNames
feature.cols <- colnames(inception.features)[grepl(colnames(inception.features), pattern="feat")]
inception.features <- inception.features[,feature.cols]

# Remove tiles with tissue threshold (per hovernet cell counts)
inception.features <- inception.features[row.names(hovernet.data), ]
coordinates <- coordinates[row.names(hovernet.data), ]
hovernet.data <- as.data.frame(t(hovernet.data))
inception.features <- as.data.frame(t(inception.features))

# Create Seurat image object
radius=1/60 #just did trial and error here for dot size on plots
img<-new(
  Class='VisiumV1',
  image=hne.image, # an array with the image data
  #assay="Spatial",
  key="image_",
  scale.factors=scalefactors(spot=1,fiducial=1,hires=1,lowres=1),
  coordinates=coordinates, #data.frame with rows as spots (tiles) and axes in the columns
  spot.radius=radius # radius size, set by trial and error above
)

# Make Seurat object: create spatial object, normalize, scale, find features, dimensionality reduction, clustering
seuratWSI1 <- CreateSeuratObject(counts = inception.features, assay = "Spatial", project = "PDXNetPDMR")
seuratWSI1@images$image <- img #add the H&E image

seuratWSI1 <- NormalizeData(seuratWSI1, normalization.method = "LogNormalize")
seuratWSI1 <- FindVariableFeatures(seuratWSI1, selection.method = "vst", nfeatures = 500)
seuratWSI1 <- ScaleData(seuratWSI1, features = rownames(seuratWSI1))
seuratWSI1 <- RunPCA(seuratWSI1, features = VariableFeatures(seuratWSI1))
seuratWSI1 <- FindNeighbors(seuratWSI1, dims = 1:15)
seuratWSI1 <- FindClusters(seuratWSI1, resolution = 0.5, graph.name="Spatial_snn")
seuratWSI1 <- RunUMAP(seuratWSI1, dims=1:15)

# Add HoVer-Net cell count data to Seurat object
# Create a new assay to store HoVerNet counts in 
hovernet_assay <- CreateAssayObject(counts = hovernet.data)
# Add this assay to the previously created Seurat object
seuratWSI1[["Hovernet"]] <- hovernet_assay 

# Cluster by inception features 
wsi_ratio <- wsi_y/wsi_x
pt_ratio <- (30000/wsi_x)*mag_ratio
s1 <- SpatialDimPlot(seuratWSI1,images="image",label=F,label.size=3,alpha=c(1,1),crop=F, pt.size.factor=pt_ratio, stroke=NA)+plot_annotation(title=sample)+theme(aspect.ratio = wsi_ratio)
ggsave(paste(sample,"_inceptionfeatures_seurat.pdf",sep=""), s1, width=8, height=8)

# Plots of HoVer-Net counts
DefaultAssay(seuratWSI1) <- "Hovernet"
#SpatialFeaturePlot(seuratWSI1, features=c("neoplastic"), alpha=c(0.1,1)) + theme(legend.position="right") #tile-level HoverNet counts
#SpatialFeaturePlot(seuratWSI1, features=c("connective"), alpha=c(0.1,1)) + theme(legend.position="right") #tile-level HoverNet counts
#SpatialFeaturePlot(seuratWSI1, features=c("inflammatory"), alpha=c(0.1,1)) + theme(legend.position="right") #tile-level HoverNet counts
#SpatialFeaturePlot(seuratWSI1, features=c("necrosis"), alpha=c(0.1,1)) + theme(legend.position="right") #tile-level HoverNet counts

# Creates all in 1 figure
#SpatialFeaturePlot(seuratWSI1,crop=F, pt.size.factor=1.8, features=c("Total","TumorFraction","StromaFraction","NecroticFraction","neoplastic","connective","necrosis","inflammatory","non.neoplastic","nolabel"),alpha=c(0.1,1))+plot_annotation(title=sample)+theme(aspect.ratio = 5)
g1 <- SpatialFeaturePlot(seuratWSI1,crop=F, pt.size.factor=pt_ratio,stroke=NA, features=c("Total"),alpha=c(0.1,1))+plot_annotation(title=sample)+theme(aspect.ratio = wsi_ratio)
g2 <- SpatialFeaturePlot(seuratWSI1,crop=F, pt.size.factor=pt_ratio,stroke=NA, features=c("TumorFraction"),alpha=c(0.1,1))+plot_annotation(title=sample)+theme(aspect.ratio = wsi_ratio)
g3 <- SpatialFeaturePlot(seuratWSI1,crop=F, pt.size.factor=pt_ratio,stroke=NA, features=c("StromaFraction"),alpha=c(0.1,1))+plot_annotation(title=sample)+theme(aspect.ratio = wsi_ratio)
g4 <- SpatialFeaturePlot(seuratWSI1,crop=F, pt.size.factor=pt_ratio,stroke=NA, features=c("NecroticFraction"),alpha=c(0.1,1))+plot_annotation(title=sample)+theme(aspect.ratio = wsi_ratio)
g5 <- SpatialFeaturePlot(seuratWSI1,crop=F, pt.size.factor=pt_ratio,stroke=NA, features=c("neoplastic"),alpha=c(0.1,1))+plot_annotation(title=sample)+theme(aspect.ratio = wsi_ratio)
g6 <- SpatialFeaturePlot(seuratWSI1,crop=F, pt.size.factor=pt_ratio,stroke=NA, features=c("connective"),alpha=c(0.1,1))+plot_annotation(title=sample)+theme(aspect.ratio = wsi_ratio)
g7 <- SpatialFeaturePlot(seuratWSI1,crop=F, pt.size.factor=pt_ratio,stroke=NA, features=c("necrosis"),alpha=c(0.1,1))+plot_annotation(title=sample)+theme(aspect.ratio = wsi_ratio)
g8 <- SpatialFeaturePlot(seuratWSI1,crop=F, pt.size.factor=pt_ratio,stroke=NA, features=c("inflammatory"),alpha=c(0.1,1))+plot_annotation(title=sample)+theme(aspect.ratio = wsi_ratio)
g9 <- SpatialFeaturePlot(seuratWSI1,crop=F, pt.size.factor=pt_ratio,stroke=NA, features=c("non.neoplastic"),alpha=c(0.1,1))+plot_annotation(title=sample)+theme(aspect.ratio = wsi_ratio)
g10 <- SpatialFeaturePlot(seuratWSI1,crop=F, pt.size.factor=pt_ratio,stroke=NA, features=c("nolabel"),alpha=c(0.1,1))+plot_annotation(title=sample)+theme(aspect.ratio = wsi_ratio)

figure <- grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,nrow = 3, ncol=4, top=textGrob(sample))
ggsave(paste(sample,"_hovernet_cellcount_seurat.pdf",sep=""), figure, width=12, height=12)

# Plot umap with Inception features, colored by Hovernet cell counts
FeaturePlot(seuratWSI1,c("Total","neoplastic","inflammatory","necrosis","connective","TumorFraction","StromaFraction","NecroticFraction"))+plot_annotation(title=sample)
ggsave(paste(sample,"_inceptionfeatures_seurat_umap.pdf",sep=""), width=12, height=10)
#FeaturePlot(seuratWSI1,c("Total","TumorFraction","StromaFraction","NecroticFraction"))+plot_annotation(title=sample)
#ggsave(paste(sample,"_inceptionfeatures_seurat_umap_fraction.pdf",sep=""), width=10, height=10)

#function to assign predominant HoVerNet cell type by tile, with prioritized tie breaker
getTileCellType<-function(x,typeNames){
  x<-as.numeric(x)
  ret<-""
  tmp1<-typeNames[which(x==max(x))]
  ret<-tmp1[1]
  if(length(tmp1)>1){
    if("inflammatory" %in% tmp1){ret<-"inflammatory"}else{
      if("connective" %in% tmp1){ret<-"connective"}else{
        if("neoplastic" %in% tmp1){ret<-"neoplastic"}else{
          if("necrosis" %in% tmp1){ret<-"necrosis"}else{
            if("non.neoplastic" %in% tmp1){ret<-"non.neoplastic"}
          }}}}}
  return(ret)
}

#function to assign predominant HoVerNet cell type by tile if >70%, else it is called "mix"
getTileCellType1<-function(x,typeNames){
  x<-as.numeric(x)
  y=x/sum(x)
  ret<-"NA"
  ret<-typeNames[which(y>0.70)]
  if(length(ret)==0){
    ret<-"mix"
  }
  return(ret)
}


#plot umap with Inception features, colored by tile-level majority Hovernet cell count
#need to update this so color is consistent by category, i.e. neoplastic always red, etc.
tileCellType <- data.frame(rbind(apply(hovernet.data[1:6,],2,getTileCellType1,row.names(hovernet.data)[1:6])))
Idents(seuratWSI1) <- tileCellType
DimPlot(seuratWSI1,pt.size=2)+plot_annotation(title=sample)
ggsave(paste(sample,"_inceptionfeatures_seurat_umap_cell.pdf",sep=""), width=8, height=8)

