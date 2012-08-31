#determine which test method to use
method = as.character(commandArgs()[4]);

if (method == "glm") {

    ### This program is for Generalized Linear Model (GLM) analysis
    ### Y = X + covar1 + covar2 ......
    ### Usually (not limited),
    ### Y is clinical trait (quantitative or binary)
    ### X is variant/gene/mutation etc.
    ### covar means covariate
    ### by Qunyuan Zhang (qunyuan@wustl.edu), 02/16/2012 updated 

    ### input options
    model.file = as.character(commandArgs()[5]);
    y.file = as.character(commandArgs()[6]);
    x.file = as.character(commandArgs()[7]);
    out.file = as.character(commandArgs()[8]);
    x.names="*";

    ### to run it on command line, type
    # R --no-save < glm.R model.file y.file x.file * out.csv
    # or (if you want to define x variables)
    # R --no-save < glm.R model.csv y.file x.file x1|x2|x3 out.csv
    # or (if x.file has been merged into y.file)
    # R --no-save < glm.R model.csv y.file * x1|x2|x3 out.csv
    # or (if you have defined x variable in model.file)
    # R --no-save < glm.R model.csv y.file * * out.csv

    ### about model.file 
    # column "type":  Q=quantitative trait, B=binary trait
    # column "y":     trait name
    # colum  "x":     variant/gene name; if x=NA or blank, it will determined by x.file and x.names
    # column "cvar":  covariate(s)
    # tab delimited 

    ### about y.file 
    # trait data file, column 1 must be sample id 
    # tab delimited 

    ### about x.file
    # usually mutation/variant/geene data file 
    # the first column must be sample id (the same as in y.file, ordered or not)
    # tab delimited  
    # x.file="*" if x.file already merged into y.file 

    ### about x.names
    # x.names="*" will use all column names in x.file as x variable names
    # or you can define it in the format x.names="gene1|gene2|gene3"
    # self-defined x.names have to be found in column names of y.file and/or x.file  

    #################### myglm fuction ##############
    myglm=function(z,trait,variant,covar=NA,ytype) {
        if (nchar(covar)==0 | is.na(covar) | is.null(covar)) { 
            model=formula(paste(trait,"~",variant)) 
        } else {
            model=formula(paste(trait,"~",variant,"+",covar))
        }
        if (ytype=="B") fit=glm(formula=model,data=z,family=binomial(link = "logit"))
        if (ytype=="Q") fit=glm(formula=model,data=z,family=gaussian(link = "identity"))
        fit
    }
    #################################################


    ### data input #####
    read.table(model.file,colClasses="character",na.strings = c("","NA"),sep="\t",header=T)->md
    read.table(y.file,na.strings = c("","NA"),sep="\t",header=T)->y
    if (x.names!="*") x.names=strsplit(x.names,split="[|]")[[1]]

    if (x.file!="*")
    {
        read.table(x.file,na.strings = c("","NA"),sep="\t",header=T)->x
        xid=colnames(x)[1]
        xs=colnames(x)[-1]
        if (x.names!="*") {x=x[,c(xid,x.names)];xs=colnames(x)[-1]}
        x.names=xs
        yid=colnames(y)[1]
        ysid = ! (colnames(y) %in% xs)
        y=y[,ysid]
        if (sum(ysid)==1) {y=data.frame(id=y);colnames(y)[1]=yid}
        #y=merge(y,x,by.x = xid, by.y = yid)
        y=merge(x,y,by.x = xid, by.y = yid)
    }

    ######### analysis ##########
    tt=NULL
    for (i in c(1:nrow(md)))
    {
        ytype=md[i,1];yi=md[i,2];xs=md[i,3];covi=md[i,4];memo=md[i,5]
        if (!is.na(xs) & nchar(xs)>0) xs=strsplit(xs,split="[|]")[[1]]
        if (is.na(xs)[1]|nchar(xs)[1]==0) xs=x.names 
        if (length(covi)==0) covi=NA
        for (xi in xs)
        {
            print(yi); print(xi); print(covi); print("******")
            if (ytype=="Q") try(anova(myglm(y,yi,xi,covi,ytype),test="F"))->fit
            if (ytype=="B") try(anova(myglm(y,yi,xi,covi,ytype),test="Chisq"))->fit
            if (class(fit)[1]!="try-error")
            {
                fit=as.matrix(fit)
                if (xi %in% rownames(fit)) tt=rbind(tt, cbind(yi,ytype,xi,as.data.frame(t(fit[xi,])),covi,memo))
            }
        }
    }
    #"yi","ytype","xi","Df","Deviance","Resid. Df","Resid. Dev","F","Pr(>F)","covi","memo"
    if (ytype=="Q") colnames(tt) = c("y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance","F_statistic","p-value","covariants","memo");
    if (ytype=="B") colnames(tt) = c("y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance","p-value","covariants","memo");
    write.table(tt,file=out.file,quote=F,sep="\t",row.names=F);

} else {

    # else, we process numerical or categorical clinical data correlation

    clinical_data = as.character(commandArgs()[5]);
    mutation_matrix = as.character(commandArgs()[6]);
    output_file = as.character(commandArgs()[7]);

    # FUNCTION finds the correlation between two variables
    cor2=function(ty,tx,method)
    {

        id=intersect(!is.na(ty),!is.na(tx));
        ty=ty[id];
        tx=tx[id];

        if(method=="cor")
        {
            tst=cor.test(tx,ty);
            s=tst$est;
            p=tst$p.value;
        }

        if(method=="wilcox")  #x must be (0,1) mutation data
        {
            tst=wilcox.test(x=ty[tx==0],y=ty[tx>=1])
            s=tst$stat
            p=tst$p.value
        }

        if(method=="chisq")
        {
            tst=chisq.test(tx,ty);
            s=tst$stat;
            p=tst$p.value;
        }

        if(method=="fisher")
        {
            tst=fisher.test(tx,ty)
            s=tst$p.value
            p=tst$p.value
        }

        if(method=="anova")
        {
            tst=summary(aov(ty~tx,as.data.frame(cbind(tx,ty))))
            s=tst[[1]]$F[1]
            p=tst[[1]]$Pr[1]
        }

        tt=c(p,s);
        tt;
    }
    # END cor2

    # FUNCTION runs correlation test on matrixes of data
    cor2test =function(y,x=NULL,method="cor",cutoff=1,sep="\t",outf=NULL)
    {

        if (!is.null(x))
        {

            if (length(x)==1) {read.table(x,header=T,sep=sep)->x;}
            if (length(y)==1) {read.table(y,header=T,sep=sep)->y;}
            colnames(y)[1]="id";
            colnames(x)[1]="id";
            tt=character(0);
            for (vi in colnames(x)[-1])
            {
                for (vj in colnames(y)[-1])
                {
                    tx=x[,c("id",vi)];
                    tx=tx[!is.na(tx[,vi]),];
                    tx=tx[!duplicated(tx[,"id"]),];
                    ty=y[,c("id",vj)];
                    ty=ty[!is.na(ty[,vj]),];
                    ty=ty[!duplicated(ty[,"id"]),];
                    xy=merge(tx,ty,by.x="id",by.y="id");
                    tx=xy[,2];
                    ty=xy[,3];
                    n=length(xy[,"id"]);
                    rst=try(cor2(ty,tx,method));
                    if (class(rst)=="try-error") {p=NA;s=NA;} else {p=rst[1];s=rst[2];}
                    t=c(vi,vj,method,n,s,p)

                    tt=rbind(tt,t);
                } #end vj
            } #end vi

            rownames(tt)=NULL;
            colnames(tt)=c("x","y","method","n","s","p");
            tt=as.data.frame(tt);
            tt[,"s"]=as.character(tt[,"s"]);
            tt[,"s"]=as.numeric(tt[,"s"]);
            tt[,"p"]=as.character(tt[,"p"]);
            tt[,"p"]=as.numeric(tt[,"p"]);
            fdr=p.adjust(tt[,"p"],method="fdr");
            bon=p.adjust(tt[,"p"],method="bon");
            tt=cbind(tt,fdr,bon);
            tt=tt[order(tt[,"p"]),];
        }

        if (is.null(x))
        {

            if (length(y)==1) {read.table(y,header=T,sep=sep)->y;}
            x=y;
            nxy=ncol(y)-1;
            colnames(y)[1]="id";
            colnames(x)[1]="id"
            tt=character(0);
            for (i in c(1:(nxy-1)))
            {
                for (j in c((i+1):nxy))
                {

                    vi=colnames(x)[-1][i];
                    vj=colnames(y)[-1][j];

                    tx=x[,c("id",vi)];
                    tx=tx[!is.na(tx[,vi]),];
                    tx=tx[!duplicated(tx[,"id"]),]
                    ty=y[,c("id",vj)];
                    ty=ty[!is.na(ty[,vj]),];
                    ty=ty[!duplicated(ty[,"id"]),];
                    xy=merge(tx,ty,by.x="id",by.y="id");
                    tx=xy[,2];
                    ty=xy[,3];
                    n=length(xy[,"id"]);
                    rst=try(cor2(ty,tx,method));
                    if (class(rst)=="try-error") {p=NA;s=NA;} else {p=rst[1];s=rst[2];}
                    t=c(vi,vj,method,n,s,p);
                    tt=rbind(tt,t);
                } #end vj
            } #end vi

            rownames(tt)=NULL;
            colnames(tt)=c("x","y","method","n","s","p");
            tt=as.data.frame(tt);
            tt[,"s"]=as.character(tt[,"s"]);
            tt[,"s"]=as.numeric(tt[,"s"]);
            tt[,"p"]=as.character(tt[,"p"]);
            tt[,"p"]=as.numeric(tt[,"p"]);
            fdr=p.adjust(tt[,"p"],method="fdr");
            bon=p.adjust(tt[,"p"],method="bon");
            tt=cbind(tt,fdr,bon);
            tt=tt[order(tt[,"p"]),];
        }


        if (!is.null(outf))
        {
            colnames(tt)=c("x","y","method","n","s","p","fdr","bon");

            #The amount of precision that R prints with is somehow machine dependent (or the R version?)
            tt[,"s"] = sapply(tt[,"s"], sprintf, fmt="%.4E");
            tt[,"p"] = sapply(tt[,"p"], sprintf, fmt="%.4E");
            tt[,"fdr"] = sapply(tt[,"fdr"], sprintf, fmt="%.2E");
            tt[,"bon"] = sapply(tt[,"bon"], sprintf, fmt="%.2E");

            #The ordering should be done after reformatting the precision (duh)
            tt=tt[order(tt[,"x"]),];
            tt=tt[order(tt[,"p"]),];

            write.table(tt,file=outf,quote=FALSE,row.names=FALSE,sep=",");
        }
        invisible(tt);
    }
    #END cor2test

    #run correlation test using function
    cor2test(y = clinical_data, x = mutation_matrix, method = method, outf = output_file);
}
