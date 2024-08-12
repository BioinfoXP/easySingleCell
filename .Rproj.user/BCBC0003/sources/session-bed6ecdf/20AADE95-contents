do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          meta.cluster = cellInfo.tb$meta.cluster,
                          colname.patient = "patient",
                          loc = cellInfo.tb$loc,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){
  ##input data
  library(data.table)
  dir.create(dirname(out.prefix),F,T)

  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)

  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}

  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)

  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }

  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx))
  }else{
    return(OR.dist.mtx)
  }
}

test.dist.table <- function(count.dist,min.rowSum=0)
{
  count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb,keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid","cid","count")
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col]-this.c
    this.m <- matrix(c(this.c,
                       sum.row[this.row]-this.c,
                       other.col.c,
                       sum(sum.col)-sum.row[this.row]-other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)
    data.frame(rid=this.row,
               cid=this.col,
               p.value=res.test$p.value,
               OR=res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
  count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
  return(count.dist.melt.ext.tb)
}


# 加载必要的包
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)

# 定义主分析函数
analyze_tissue_dist <- function(meta_data, output_prefix, pdf_width = 8, pdf_height = 4, verbose = 1) {

  # 调用 do.tissueDist 函数进行主要分析
  OR_immune_list <- do.tissueDist(cellInfo.tb = meta_data,
                                  out.prefix = sprintf("%s.Immune_cell", output_prefix),
                                  pdf.width = pdf_width, pdf.height = pdf_height, verbose = verbose)

  # 返回分析结果
  return(OR_immune_list)
}

# 定义绘图函数
plot_heatmap <- function(OR_list) {
  # 提取 OR 值结果
  a <- OR_list[["OR.dist.tb"]] %>%
    as.data.frame() %>%
    column_to_rownames(var = "rid") %>%
    na.omit()

  # 提取 P 值结果
  b <- OR_list$count.dist.melt.ext.tb[, c(1, 2, 6)] %>%
    spread(key = "cid", value = "adj.p.value") %>%
    column_to_rownames(var = "rid")

  # 只选择在a中的行
  b <- b[rownames(a),]

  # 调整 P 值符号表示
  col <- viridis(11, option = "D")
  b <- ifelse(b >= 0.05 & (a > 1.5 | a < 0.5), "",
              ifelse(b < 0.0001 & (a > 1.5 | a < 0.5), "****",
                     ifelse(b < 0.001 & (a > 1.5 | a < 0.5), "***",
                            ifelse(b < 0.01 & (a > 1.5 | a < 0.5), "**",
                                   ifelse(b < 0.05 & (a > 1.5 | a < 0.5), "*", "")))))

  bk <- c(seq(0, 0.99, by = 0.01), seq(1, 2, by = 0.01))

  # 绘制热图
  pheatmap(a, border_color = NA, fontsize = 9, cellheight = 12, cellwidth = 20,
           clustering_distance_rows = "correlation", display_numbers = b,
           number_color = "black", fontsize_number = 10, cluster_col = FALSE,
           cluster_rows = TRUE, breaks = bk, treeheight_row = 20, treeheight_col = 20,
           color = c(colorRampPalette(colors = col[1:6])(length(bk) / 2),
                     colorRampPalette(colors = col[6:11])(length(bk) / 2)))
}


