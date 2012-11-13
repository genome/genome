### Survival analysis for mutation data ###

### original location of code: /gscuser/qzhang/gstat/survival/survival.R
### example input file: /gscuser/qzhang/gstat/survival/tcga.tsv

### Run it on command line like below
### for example,   R --no-save --args < survival.R vital_status.input mut_matrix.input legend.placement output_dir &

### clinical data /vital status input file, first three columns are sample_ID, survival_time, vital_status (0=living, 1=deceased)

######################## read input arguments

clinical.survival.data=commandArgs()[4];
mut.data=commandArgs()[5];
legend.placement=commandArgs()[6];
out.dir=commandArgs()[7];

######################## read and prepare data

vitals = read.table(clinical.survival.data,header=T);
mut_matrix = read.table(mut.data,header=T);
x = merge(vitals,mut_matrix,by.x=1,by.y=1);
write.table(x,file=paste(out.dir,"survival_analysis_data_matrix.tsv",sep="/"),quote=F,append=F,row.names=F,sep="\t")
colnames(x)[-c(1:3)]->phenos
if (class(x[,phenos])=="integer" & length(unique(x[,phenos]))<6) x[,phenos] [x[,phenos]>1]=1
# Make a list of distinctive colors using afriggeri.github.com/RYB
distinctColors = c(rgb(255,63,0,maxColorValue=255),rgb(51,23,0,maxColorValue=255),rgb(0,168,51,maxColorValue=255),rgb(41,95,153,maxColorValue=255),rgb(255,127,0,maxColorValue=255),rgb(255,127,127,maxColorValue=255),rgb(127,0,127,maxColorValue=255),rgb(127,211,25,maxColorValue=255),rgb(148,175,204,maxColorValue=255),rgb(153,75,0,maxColorValue=255),rgb(255,255,127,maxColorValue=255),rgb(89,11,63,maxColorValue=255),rgb(20,131,102,maxColorValue=255),rgb(191,0,63,maxColorValue=255),rgb(127,211,25,maxColorValue=255),rgb(255,255,0,maxColorValue=255),rgb(25,96,25,maxColorValue=255));

######################### survival analysis

library(survival)
logr=NULL

for (phenotype in phenos)
{
    #clean data
    loopdata <- x;
    loopdata <- loopdata[!is.na(loopdata[,phenotype]),];
    loopdata <- loopdata[!is.na(loopdata[,3]),];
    loopdata <- loopdata[!is.na(loopdata[,2]),];
    status=loopdata[,3];
    time=loopdata[,2];
    x1=loopdata[,phenotype];
    base.class = as.vector(sort(unique(x1)))[1];

    coxph(Surv(time, status) ~ x1, loopdata) -> co;
    summary(co)->co; co$conf->cox; co$logtest[3]->p; co$coef[5]->indv.p;
    rownames(cox) = sub("x1","",rownames(cox));
    if (length(rownames(cox))==1 && rownames(cox)[1]=="") { rownames(cox)[1] = "1"; }
    logr=rbind(logr,cbind(base.class,rownames(cox),phenotype,cox,indv.p,p))

    mfit.by <- survfit(Surv(time, status == 1) ~ x1, data = loopdata)
    ## file name for plot
    bitmap(file=paste(out.dir,"/",phenotype,"_survival_plot.png",sep=""))
    ## create survival plot
    plot(mfit.by,lty=1:10,ylab="Survival Probability",xlab="Time",col=distinctColors)
    if (dim(table(x1))>1) {
        title(paste(phenotype,", P=",signif(p,3),sep=""));
    }
    else {
        title(paste(phenotype));
    }
    legend(x=legend.placement, legend=names(table(x1)), lty = 1:10, col=distinctColors)
    dev.off()
}

########################## calculate fdr

logr=logr[,-5];
if (length(phenos) < 2) { logr=(t(logr)); }
fdr=p.adjust(as.numeric(logr[,"p"]),"fdr")
logr=cbind(logr,fdr)

######################### print output

colnames(logr)[1:9]=c("base.class","comparison.class","phenotype","hazard.ratio","lower.95","upper.95","2-class-p-value","p-value","fdr")
logr=logr[order(logr[,"p-value"]),]
write.table(logr,file=paste(out.dir,"survival_analysis_test_results.tsv",sep="/"),quote=F,append=F,row.names=F,sep="\t")
