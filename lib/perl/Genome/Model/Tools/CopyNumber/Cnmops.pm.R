#!/usr/bin/env Rscript

usage <- function()
{
  writeLines("Usage:\n\tRscript Cnmops.pm.R tumor.bam normal.bam capture.bed outdir --test[or --notest]")
}

convert_start_end_character <-function(t) {
  t$start <- sprintf("%d", t$start)
  t$end <- sprintf("%d", t$end)
  return(t)
}

get_capture_regions <- function()
{
  segments <- read.table(capture_bed,sep="\t",as.is=TRUE)
  gr <- GRanges(segments[,1],IRanges(segments[,2],segments[,3]))
  if(test) {
    gr <- gr[1:15000]
  }
  return(gr)
}

get_tumor_normal_logrr <- function(ref_analysis_norm)
{
  R  <- cn.mops:::.makeLogRatios(ref_analysis_norm, "CN2")
  regions <- data.frame(seqnames=as.vector(seqnames(R)),
      starts=start(R)-1, ends=end(R))
  lrr <- attr(ref_analysis_norm, "individualCall")
  bed_logrr <- cbind(regions, lrr)
  colnames(bed_logrr) <- c("chr", "start", "end", "tumor_normal_log2r");
  bed_logrr<-convert_start_end_character(bed_logrr)
  op_bed_logrr = paste(out_dir, "/cnmops.bed_logrr.txt", sep = "");
  write.table(bed_logrr, file = op_bed_logrr, col.names = T,
      row.names = F, sep = "\t", quote = F)
  return(bed_logrr)
}

get_read_counts <- function()
{
  X.nor <- getSegmentReadCountsFromBAM(normal_bam, GR=gr, mode="paired")
  X.tum <- getSegmentReadCountsFromBAM(tumor_bam, GR=gr, mode="paired")
  X <- X.tum
  values(X) <- cbind(values(X.tum), values(X.nor))
  X <- normalizeGenome(X)
  save(X, file = paste(out_dir, "/cnmops.RD.Robject", sep = ""))
  seqlevels(X, force = TRUE)<-c(1:22, "X", "Y")
  return(X)
}

plot_segments <- function()
{
  segment_file=paste(out_dir, "/cnmops.segplot_WG.pdf", sep = "")
  pdf(segment_file, width = 14, height = 14)
  segplot(ref_analysis_norm, segStat = "median")
  dev.off()
  for (i in 1:22) {
    if(i %in% as.vector(seqnames(localAssessments(ref_analysis_norm)))) {
      file = paste(out_dir, "/cnmops.segplot_chr", i, ".pdf", sep = "");
      pdf(file, width = 10, height = 10);
      segplot(ref_analysis_norm, seqnames=i, segStat = "median");
      title(ylab = "log2 RD(tumor/normal)")
      abline(h=log2(1/2), col = "green")
      abline(h=log2(3/2), col = "blue")
      abline(h=log2(4/2), col = "yellow")
      legend("bottomright", inset=.05, title="Copy Number difference (tumor - normal)", c("1","2", "CNV segment", "-1"), horiz=TRUE, col = c("blue", "yellow", "red", "green"), lty=1)
      dev.off()
    }
  }
}

save_ref_analysis_norm <- function(ref_analysis_norm)
{
  ref_analysis_norm_file = paste(out_dir, "/cnmops.ref_analysis_norm.Robject", sep = "")
  if(length(cnvr(ref_analysis_norm)) != 0) {
    ref_analysis_norm <- cn.mops:::.replaceNames(ref_analysis_norm,
        colnames(ref_analysis_norm@normalizedData),"Tumor_Normal_Log2Ratio")
  }
  save(ref_analysis_norm, file=ref_analysis_norm_file)
}

write_cnvr<-function()
{
  cnv_r<-cnvs(ref_analysis_norm)
  df <- data.frame(seqnames=as.vector(seqnames(cnv_r)), starts=start(cnv_r)-1, ends=end(cnv_r))
  colnames(df) <- c("chr", "start", "end")
  convert_start_end_character(df)
  df2 <- data.frame(seqnames=as.vector(seqnames(cnv_r)), starts=start(cnv_r)-1, ends=end(cnv_r), medians=values(cnv_r)$median, means=values(cnv_r)$mean, CN=values(cnv_r)$CN)
  convert_start_end_character(df2)
  colnames(df2) <- c("chr", "start", "end", "median", "mean", "CN")
  op_cnv_bed = paste(out_dir, "/cnmops.cnv.bed", sep = "")
  write.table(df, file=op_cnv_bed, quote=F, sep="\t", row.names=F, col.names=T)
  op_cnv_txt = paste(out_dir, "/cnmops.cnvs.txt", sep = "")
  write.table(df2, file=op_cnv_txt, quote=F, sep="\t", row.names=F, col.names=T)
}

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 5) {
  usage()
  quit()
}
library(cn.mops)

tumor_bam = args[1]
normal_bam = args[2]
capture_bed = args[3]
out_dir = args[4]

if(args[5] == "--test") {
  test = TRUE
} else {
  test = FALSE
}

set.seed(42)#use a deterministic seed, this enables testing and doesn't really affect the results.
gr<-get_capture_regions()
RD<-get_read_counts()
ref_analysis_norm <- referencecn.mops(RD[,1], RD[,2], returnPosterior = TRUE, segAlgorithm = "DNAcopy", norm = 0, lowerThreshold = -0.415, upperThreshold = 0.321, priorImpact = 20)
ref_analysis_norm <- calcIntegerCopyNumbers(ref_analysis_norm)
tumor_normal_logrr <- get_tumor_normal_logrr(ref_analysis_norm)
save_ref_analysis_norm(ref_analysis_norm)
plot_segments()
write_cnvr()
