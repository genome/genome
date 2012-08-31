args = commandArgs(TRUE)
library(stats)
library(ggplot2)
f <- function(x, height) {
    ans <- median(x)
        data.frame(ymin=ans, ymax=ans, y=ans)
}
genes = read.table(args[1], header=FALSE)
fpkm_matrix = read.table(args[2], header=TRUE)
fpkm.subset = fpkm_matrix[fpkm_matrix$mapped_gene_name %in% genes$V1,]
tumor_cols = grep("tumor", names(fpkm_matrix))
normal_cols = grep("normal", names(fpkm_matrix))
tumor_names <- names(fpkm.subset[,tumor_cols])
normal_names <- names(fpkm.subset[,normal_cols])
tumor_data <- fpkm.subset[,!(names(fpkm.subset) %in% normal_names)]
tumor.melt <- melt(tumor_data, id.vars=c("tracking_id", "ensg_name", "mapped_gene_name", "locus"))
tumor.melt$TISSUE <- "TUMOR"
normal_data <- fpkm.subset[,!(names(fpkm.subset) %in% tumor_names)]
normal.melt <- melt(normal_data, id.vars=c("tracking_id", "ensg_name", "mapped_gene_name", "locus"))
normal.melt$TISSUE <- "NORMAL"
patty.melt <- rbind(tumor.melt, normal.melt)
names(patty.melt)[names(patty.melt)=="variable"]<-"SAMPLE NAME"
names(patty.melt)[names(patty.melt)=="value"]<-"FPKM"
p <- ggplot(patty.melt, aes(x=TISSUE,y=FPKM))
col_num = ceiling(nrow(genes)/2)
preplot <- p + geom_jitter(position=position_jitter(width=.01,height=.01)) + facet_wrap(facets="mapped_gene_name", scales="free_y", ncol=col_num ) + opts(strip.text.x = theme_text(size = 20), axis.text.x=theme_text(size=15), axis.text.y=theme_text(size=13)) 
plot <- preplot + stat_summary(fun.data = f, geom = "crossbar", colour = "red", width=0.3) 
height =12
if(nrow(genes) == 1) { height = 6 }
pdf(args[3], height=height, width=(col_num * 5))
print(plot)
dev.off()

