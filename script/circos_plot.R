library(circlize)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
circos.genomicHeatmap <- function (bed, col, numeric.column = NULL, border = NA, border_lwd = par("lwd"), border_lty = par("lty"), 
                                   connection_height = convert_height(5, "mm"), line_col = par("col"), line_lwd = par("lwd"), 
                                   line_lty = par("lty"), heatmap_height = 0.15, side = c("inside", "outside"), 
                                   track.margin = circos.par("track.margin")) {
  mat <- bed[, -(1:3), drop = FALSE]
  if (is.null(numeric.column)) {
    numeric.column <- which(apply(mat, 2, "is.numeric"))
    if (length(numeric.column) == 0) {
      stop("You don't have numeric columns in `bed`.")
    }
  }
  else {
    if (is.numeric(numeric.column)) {
      numeric.column <- numeric.column - 3
    }
  }
  mat <- mat[, numeric.column, drop = FALSE]
  mat <- as.matrix(mat)
  bed <- cbind(bed[, 1:3], mat)
  side <- match.arg(side)
  if (missing(col)) {
    col <- colorRamp2(seq(min(na.rm = TRUE), max(mat, na.rm = TRUE), 
                         length = 3), c("blue", "#EEEEEE", "red"))
  }
  if (is.function(col)) {
    col <- col(mat)
  }
  if (length(border) == 1) {
    if (is.na(border)) {
      border <- col
    }
  }
  if (length(border) == 1) {
    border <- rep(border, length(col))
    dim(border) <- dim(col)
  }
  if (side == "inside") {
    circos.genomicPosTransformLines(bed, posTransform = posTransform.default, 
                                    horizontalLine = "top", track.height = connection_height, 
                                    track.margin = c(convert_height(1, "mm"), track.margin[2]), 
                                    cell.padding = c(0, 0, 0, 0), col = line_col, lwd = line_lwd, 
                                    lty = line_lty)
    circos.genomicTrackPlotRegion(bed, stack = TRUE, panel.fun = function(region, value, ...) {
      l = bed[, 1] == CELL_META$sector.index
      i = getI(...)
      circos.genomicRect(region, value, col = col[l, i], 
                         border = border[l, i], lwd = border_lwd, lty = border_lty, 
                         posTransform = posTransform.default, ...)
    }, bg.border = NA, track.height = heatmap_height, track.margin = c(track.margin[1], 
                                                                       0), cell.padding = c(0, 0, 0, 0))
  }
  else {
    circos.genomicTrackPlotRegion(bed, stack = TRUE, panel.fun = function(region, value, ...) {
      l = bed[, 1] == CELL_META$sector.index
      i = getI(...)
      circos.genomicRect(region, value, col = col[l, i], 
                         border = border[l, i], lwd = border_lwd, lty = border_lty, 
                         posTransform = posTransform.default, ...)
    }, bg.border = NA, track.height = heatmap_height, track.margin = c(0, track.margin[2]), cell.padding = c(0, 0, 0, 0))
  }
}

argv <- commandArgs(T)
input <- argv[1]
output <- argv[2]
# input <- "circos/input_H3K4me3_cov_1000000.tsv"
# output <- "circos/input_H3K4me3_cov_1000000.pdf"

coverage <- read.table(input)
coverage_norm <- apply(coverage[,5:9], 2, function(x){(x-coverage[,4])/sum(x-coverage[,4])*1000})
cov_max <- ceiling(max(coverage_norm))
coverage_chrom <- apply(coverage[,2:3], 2, function(x){x/10000000})
coverage_norm <- data.frame(coverage$V1, coverage_chrom, coverage_norm)
col_value <- c(min(coverage_norm[,4:8]), median(as.matrix(coverage_norm[,4:8])), max(coverage_norm[,4:8]))
legend_col_value <- c(min(coverage_norm[,4:8]), median(as.matrix(coverage_norm[,4:8]))-min(coverage_norm[,4:8]), median(as.matrix(coverage_norm[,4:8])), max(coverage_norm[,4:8])-median(as.matrix(coverage_norm[,4:8])), max(coverage_norm[,4:8]))
col_fun <- colorRamp2(col_value, c("green", "yellow", "red"))

pdf(output)
col_set <- brewer.pal(12, 'Set3')
chrom_size <- read.table("Bdi.chromSize.for.circos")
circos.clear()
circos.par("track.height" = 0.05, start.degree = 80, gap.degree = c(1,1,1,1,20))
circos.initialize(factors = chrom_size$V1, x = chrom_size$V2/10000000)
circos.track(factors = chrom_size$V1, y = chrom_size$V2,bg.col = col_set[11],
              panel.fun = function(x,y){
                n=floor(CELL_META$xlim[2])
                circos.axis(major.at=seq(0,CELL_META$xlim[2]+1),labels = seq(0,n))
                circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(7, "mm"), 
                            CELL_META$sector.index)
                }
             )
circos.genomicHeatmap(coverage_norm[,c(1,2,3,8)], col = col_fun, side = "outside", heatmap_height = 0.1, line_col = "white", connection_height = 0.005)
circos.genomicHeatmap(coverage_norm[,c(1,2,3,7)], col = col_fun, side = "outside", heatmap_height = 0.1, line_col = "white", connection_height = 0.005)
circos.genomicHeatmap(coverage_norm[,c(1,2,3,6)], col = col_fun, side = "outside", heatmap_height = 0.1, line_col = "white", connection_height = 0.005)
circos.genomicHeatmap(coverage_norm[,c(1,2,3,5)], col = col_fun, side = "outside", heatmap_height = 0.1, line_col = "white", connection_height = 0.005)
circos.genomicHeatmap(coverage_norm[,1:4], col = col_fun, side = "outside", heatmap_height = 0.1, line_col = "white", connection_height = 0.005)

lgd_links <- Legend(at = c(0:cov_max), col_fun = col_fun, title_position = "topleft", title = "")
lgd_list_vertical <- packLegend(lgd_links)
pushViewport(viewport(x = unit(2, "mm"), y = unit(4, "mm"), 
                      width = grobWidth(lgd_list_vertical), 
                      height = grobHeight(lgd_list_vertical), 
                      just = c('left', 'bottom')))
grid.draw(lgd_list_vertical)
p <- seq(0.42, 0.85, length=5)
text(0, p, LETTERS[1:5], cex = 1.2, font = 2)
dev.off()


