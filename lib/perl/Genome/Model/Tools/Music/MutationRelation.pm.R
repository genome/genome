#get command line arguments
mutation_file = as.character(commandArgs()[4]);
permutations = as.numeric(commandArgs()[5]);
output_file = as.character(commandArgs()[6]);

# FUNCTION prepare mutation matrix data for permutation
seq.vjoint=function(x)
{
    tt=NULL;
    for (i in 1:length(x)) 
    {
        tt=c(tt,rep(names(x[i]),x[i]));
    }
    tt;
}
# END seq.vjoint

# FUNCTION randomly permute the mutation matrix, keeping mutation-number/sample unchanged
matrix.sample.v1=function(x,keep.freq=T)
{
    colSums(x)->pool;
    rowSums(x)->freq;
    freq=freq[order(freq,decreasing=T)];
    tt=NULL;
    for (i in c(1:length(freq)))
    {
        pooli=seq.vjoint(pool);
        xi=NULL;
        for (j in 1:freq[i])
        {
            if (length(pooli)>1) sample(pooli,1)->temp;
            if (length(pooli)==1) pooli->temp;
            xi=c(xi,temp);
            pooli=pooli[pooli!=temp];
        }
        tt=rbind(tt,cbind(i,xi));
        pool[xi]=pool[xi]-1;
        pool=pool[pool>0];
    }
    tt=table(tt[,1],tt[,2]);
    tt[tt>1]=1;
    tt;
}
# END matrix.sample.v1

# FUNCTION calculate the probability of having no correlation between any two genes
prob.concur=function(x)
{
    x>0->x;
    cols=colnames(x);
    cols=cols[order(cols)]; 
    n=length(cols);

    tt=NULL;
    for (i in 1:(n-1))
    {
        for (j in (i+1):n)
        {
            vi=cols[i];
            vj=cols[j];
            sum(x[,vi] & x[,vj])->nand;
            sum(x[,vi] != x[,vj])->nexc;
            ni=sum(x[,vi]);
            nj=sum(x[,vj]);
            temp=cbind(vi,vj,ni,nj,nand,nexc);
            tt=rbind(tt,temp);
        }
    }
    tt=as.data.frame(tt,stringsAsFactors=F);
    tt$nand=as.numeric(tt[,"nand"]);
    tt$nexc=as.numeric(tt[,"nexc"]);
    tt;
}
# END prob.concur

# FUNCTION performing mutation correlation
mut_cor_permu_test=function(x=NULL,n.permu=100,seed=NULL,mut.file=NULL,out.file=NULL,out.rdata=NULL)
{

    if (!is.null(mut.file)) read.table(mut.file,header=T)->x

    x=x[,-1];x[x>0]=1
    x=x[rowSums(x)>0,colSums(x)>0]
    prob.concur(x)->pc0

    if (!is.null(seed)) set.seed(seed)

    pp=0;en=0
    for (i in 1:n.permu)
    {

        doit=1; 
        while(doit==1) 
        {
            xi=matrix.sample.v1(x);
            if (max(xi)==1) doit=0;
        }

        prob.concur(xi)->pci;
        en=en+pci[,c("nand","nexc")];
        as.numeric(pci$nand>=pc0$nand)->pand;
        as.numeric(pci$nexc>=pc0$nexc)->pexc;
        pp=pp+cbind(pand,pexc);
    }

    en=en/n.permu;pp=pp/n.permu;
    pci[,c("nand","nexc")]=en;
    pci=cbind(pci,pp);
    tt=merge(pc0,pci,by.x=c("vi","vj"),by.y=c("vi","vj"));

    if (!is.null(out.rdata)) save(tt,file=out.rdata,compress=T);
    names(tt)= c("Gene1","Gene2","CntGene1","CntGene2","AndCnt","XorCnt","Perm_CntGene1","Perm_CntGene2","Perm_AndCnt","Perm_XorCnt","Pvalue_And","Pvalue_Xor")
    if (!is.null(out.file)) write.table(tt,file=out.file,quote=F,sep="\t",row.names=F);

    invisible (tt);
}
#END mut_cor_permu_test

#run test using mut_cor_permu_test
mut_cor_permu_test(mut.file=mutation_file,n.permu=permutations,out.file=output_file);
