#' Title plot_scRNA_reducion
#'
#' @param scRNA A seurat object
#' @param reduction umap tsne
#' @param plot.column Which column in metadata to plot, default is "seurat_clusters"
#' @param label Whether to plot cell labels, default is TRUE
#' @param label.text.size The cell label size, default is 6
#' @param cellcolors A color vector. cellcolors=c("#8B0000","#E18727FF","#0072B5FF","#800080","#BC3C29FF","#20854EFF")
#'
#' @return
#' @export
#'
#' @examples plot_umap(scRNA = scRNA,plot.column = "group",label = F)
plot_scRNA_reducion <- function(scRNA,reduction='umap',plot.column="seurat_clusters",
                      label=T,label.text.size=6,
                      cellcolors=c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF", "#008280FF", "#BB0021FF", "#5F559BFF", "#A20056FF", "#808180FF", "#1B1919FF")) {
  # https://mp.weixin.qq.com/s/5euLsq08ckP_-WDE819K1Q
  library(tidyverse)
  library(cowplot)
  library(ggplot2)
  scRNA <-  scRNA
  #提取细胞类型信息 “labels”
  data1 <- scRNA@meta.data
  table(data1$labels)
  #提取x、y轴需要的数据 tSNE_1、tSNE_2
  data2 <- scRNA@reductions[[reduction]]@cell.embeddings%>%as.data.frame()
  #绘图数据合并
  mydata <- merge(data2,data1,by=0,all.x = T)%>%
    column_to_rownames("Row.names")
  mydata$labels <- mydata[,plot.column]
  # set colors
  sort(unique(mydata$labels))
  x.axis = names(mydata)[1]
  y.axis = names(mydata)[2]

  p <- ggplot(data = mydata,aes_string(x.axis,y.axis,fill='labels',colour='labels')) +
    geom_point(shape=21,size=1.5,alpha=0.5)+ #点形状、大小、透明度
    scale_fill_manual(values = cellcolors)+
    scale_colour_manual(values = cellcolors)+
    theme_bw(base_rect_size = 1)+ #主题设置可根据自己的喜好设置
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.position = "right",
          legend.key.height = unit(1,'cm'),
          legend.key.width = unit(0.5,'cm'))+
    guides(fill=guide_legend(override.aes = list(size=3,alpha=1))) #为了美观，需要修改图例点的大小和透明度

  if (label) {
    # 设置文本位置，这里设置中位数为例
    median_df <- mydata %>% group_by(labels) %>%
      summarise(median.1 = median(UMAP_1),
                median.2 = median(UMAP_2))
    head(median_df)
    p <- p + ggrepel::geom_text_repel(data = median_df,
                                      aes(median.1,median.2, label =labels), #文字注释：位置+内容
                                      size=label.text.size,
                                      color= "gray20",
                                      min.segment.length = 0)
  }
  return(p)
}



