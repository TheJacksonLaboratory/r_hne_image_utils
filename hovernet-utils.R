#' Tile (post-processed) Hover-Net output.
#'
#' @param hovernet_centroid_file A string holding the post-processed Hover-Net results from an image.
#'                               This will have been created by the python function
#'                               extract_cell_info_from_hovernet_output defined in hovernet_utils.py in the
#'                               library hne_image_utils. See the repo: 
#'                               https://github.com/TheJacksonLaboratory/hne_image_utils
#'                               It is expected to have one row per nucleus with have columns centroid_x and centroid_y giving the
#'                               x and y coordinates of the centroid of that nucleus.
#' @param json_cell_info_file A string holding the path to Hover-Net's type_info.json file describing each of the cell types.
#'                            e.g., from cell type id to cell type name.
#' @param target.image.resolution An int giving the resolution at which we would like to tile the image.
#' @param hovernet.resolution An int giving the resolution at which Hover-Net was run (likely, 40 as this is the
#'                            resolution at which Hover-Net was trained on the PanNuke dataset)
#' @param tile_size An int giving the length and width of the tile in the _target.image.resolution_
#' @return a list with three entries:
#'         hovernet_tiles is a data.frame, where each row corresponds to a particular cell type in a particular tile
#'                         and the columns are:
#'                         image_id: the id of the image
#'                         tile_x, tile_y: the coordinates of the upper left hand corner of the tile_size x tile_size tile in the target.image.resolution.
#'                         cell_type_label: the label of the cell type,
#'                         count: the number of cells of that type in the _tile_,
#'                         area: the total area of cells of that type in the _tile_ _scaled_ to be the area in the target.image.resolution,
#'                         median_area: the median area of cells of that type in the _tile_ _scaled_ to be the area in the target.image.resolution,
#'         hovernet_slides is a data.frame, where each row corresponds to a cell type
#'                         and the columns are:
#'                         image_id: the id of the image
#'                         cell_type_label: the label of the cell type,
#'                         count: the number of cells of that type in the _image_,
#'                         area: the total area of cells of that type in the _image_ _scaled_ to be the area in the target.image.resolution,
#'                         median_area: the median area of cells of that type in the _image_ _scaled_ to be the area in the target.image.resolution,
#'         hovernet_nuclei is a data.frame, where each row corresponds to a nucleus
#'                         and the columns are:
#'                         image_id: the id of the image
#'                         cell_type_label: the label of the cell,
#'                         centroid_x, centroid_y: the x, y coordinates of the cell
#'                         tile_x, tile_y: the top-left coordinates of the tile to which the cell belongs
#'                         cell_area: the area of the cell
process.hovernet.centroids <- function(hovernet_centroid_file, json_cell_info_file, target.image.resolution = 20, hovernet.resolution = 40, tile_size = 512) {
  # hovernet_centroids <- read.table(hovernet_centroid_file, sep=",", header=TRUE)
  hovernet_centroids <- fread(hovernet_centroid_file)
  
  # Convert from 40x hovernet.resolution at which HoverNet was run to possibly lower resolution of original image
  scale.factor = hovernet.resolution / target.image.resolution
  # Convert the centroids to coordinates
  # The tile_size / 2 puts them in the middle of the tile ... that's OK for seurat, but not here.
  # Instead, assign to the tile "upper left corner"
  # hovernet_centroids$tile_x <- ( tile_size * floor(hovernet_centroids$centroid_x / scale.factor / tile_size) ) + ( tile_size / 2 )
  # hovernet_centroids$tile_y <- ( tile_size * floor(hovernet_centroids$centroid_y / scale.factor / tile_size) ) + ( tile_size / 2 )
  hovernet_centroids$tile_x <- ( tile_size * floor(hovernet_centroids$centroid_x / scale.factor / tile_size) )
  hovernet_centroids$tile_y <- ( tile_size * floor(hovernet_centroids$centroid_y / scale.factor / tile_size) )
  
  # Read in the labels for the cell types
  cell_info_tbl = fromJSON(file = json_cell_info_file)
  cell_type = as.numeric(names(cell_info_tbl))
  cell_type_label = as.vector(unlist(lapply(cell_info_tbl, function(x) x[[1]])))
  cell_info_tbl <- data.frame(cell_type_label = cell_type_label, cell_type = cell_type)
  hovernet_centroids <- merge(hovernet_centroids, cell_info_tbl)
  
  # Count the number of cells of each type and the total area _at the tile level_
  # NB: need to scale the entire area
  hovernet_tiles <-
    ddply(hovernet_centroids[, c("image_id", "tile_x", "tile_y", "cell_type_label", "cell_type_area")],
          .variables = c("image_id", "tile_x", "tile_y", "cell_type_label"),
          .parallel = FALSE,
          .fun = function(df) {
            cell_type_label = df[1, "cell_type_label"]
            tile_count = nrow(df)
            tile_areas = df$cell_type_area / ( scale.factor * scale.factor )
            tile_area = sum(tile_areas)
            median_tile_area = median(tile_areas)
            data.frame(cell_type_label = cell_type_label, count = tile_count, area = tile_area, median_area = median_tile_area)
          })
  
  # Collapse results to the slide / image level
  hovernet_slides <-
    ddply(hovernet_tiles,
          .variables = c("image_id", "cell_type_label"),
          .parallel = FALSE,
          .fun = function(df) {
            cell_type_label = df[1, "cell_type_label"]
            slide_count = sum(df$count)
            slide_area = sum(df$area)
            data.frame(cell_type_label = cell_type_label, count = slide_count, area = slide_area)
          })

  hovernet_centroids <- hovernet_centroids[, c("image_id", "centroid_x", "centroid_y", "tile_x", "tile_y", "cell_type_label", "cell_type_area")]
  colnames(hovernet_centroids) <- c("image_id", "centroid_x", "centroid_y", "tile_x", "tile_y", "cell_type_label", "cell_area")
  ret.list <- list("hovernet_tiles" = hovernet_tiles, "hovernet_slides" = hovernet_slides, "hovernet_nuclei" = hovernet_centroids)
  return(ret.list)
}

#' Extract cell type info from Hover-Net type_info.json file
#'
#' @param json_path A string holding the path to the type_info.json file input to Hover-Net
#' @return A data.frame whose rows correspond to a cell type and with columns:
#            id: the cell type id (e.g., as included in the Hover-Net JSON output)
#            label: the label of the cell type (e.g., neopla, etc)
#            r_color: the red channel of the RGB color specification used by Hover-Net for this cell type
#            g_color: the green channel of the RGB color specification used by Hover-Net for this cell type
#            b_color: the blue channel of the RGB color specification used by Hover-Net for this cell type
extract.hovernet.cell.info <- function(json.type.info.path) {

  cell.dict <- fromJSON(file=json.type.info.path)
  df <- 
    ldply(cell.dict,
          .fun = function(elem) {
            data.frame(label = elem[[1]], 
                       r_color = elem[[2]][1], g_color = elem[[2]][2], b_color = elem[[2]][3])
          })
  colnames(df)[1] <- "id"
  df
}

#' Extract cell info (e.g., centroid location and cell type) from Hover-Net JSON output.
#'
#' @param json_path A string holding the path to the json file output by Hover-Net
#' @return A list containing
#'           mag: the magnification at which Hover-Net was run
#'           df: A data.frame whose rows correspond to cells and with columns:
#'             id: the id of the cell
#'             centroid_x, centroid_y: the x and y coordinates of the cell's centroid
#'             type: cell type identifier to be cross referenced with HoverNet type_info.json file
#'             type_prob: the probability of cell type
extract.cell.info.from.hovernet.output <- function(json.path) {
  js <- fromJSON(file=json.path)
  mag <- js["mag"]
  # length(js[["nuc"]])
  df <- 
    ldply(js[["nuc"]],
          .fun = function(elem) {
            data.frame(centroid.x = elem$centroid[1], centroid.y = elem$centroid[2],
                       type = elem$type, type_prob = elem$type_prob)
          })
  colnames(df)[1] <- "id"
  list("mag" = mag, "df" = df)
}
