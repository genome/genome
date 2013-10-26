1. 
prepare a option file, mine is option_file_asms, most are self-explained

2. 
run it like below

bsub -e err 'R --no-save < burdentest.R option_file_asms Q trigRES ABCA1 10000:0:1'

or if you want to fix permutations for all analysis, use a non-zero seed (e.g. 123),like 10000:123:1

to repeat the analysis to get minimum permutation p value dsitribution for multiple trait, use -10000:123:100
which will repeat the analysis 100 times, thus the total permutation number is 10000*100=1M, very time consuming 
you may try -1000:123:100 to save 10 times of computer time

or 

bsub -e err 'R --no-save < burdentest.R option_file_asms Q trigRES ABCA1 10000:123 additional_covariates'

Input info:

option_file_asms (example option file)
Q (Q:quantitative trait, for binary trait, code it to 0,1, use "B" here)
trigRES (exampletrait name, must be in phenotype file)
ABCA1 (example gene name, must be in annotation file)
10000 (permutation number, increase it will incraese computer time, usually 10000, at least 100, can be more)

 
3.
Results will be generated in the output dir defined in the option file
As an example, mine is
out.dir="/gscmnt/sata424/info/medseq/Freimer-Boehnke/burdentest20120205/results"

*.single.csv  is single variant test result
*.burden.csv  is burden test result
*.null means no data 
*.error means unexpected data error

In burden.csv file, the columns below are included 

Trait
Gene
N : sample size
V : variant number for collapsing
MAF : Freq. cutoff
CMC : Li B, Leal SM. 2008. Methods for detecting associations with rare
variants for common diseases: application to analysis of sequence
data. Am J Hum Genet 83:311¨C321.

pCMC: Li and Leal's method with permuation

WSS: Madsen BE, Browning SR. 2009. A groupwise association test for rare
mutations using a weighted sum statistic. PLoS Genet 5:e1000384

aSum: Han F, Pan W. 2010. A data-adaptive sum test for disease association
with multiple common or rare variants. Hum Hered 70:42¨C54.

PWST : Qunyuan Zhang,Marguerite R. Irvin,2 Donna K. Arnett,2 Michael A. Province,1 and Ingrid Borecki1
Genetic Epidemiology (2011) A Data-Driven Method for Identifying Rare Variants with
Heterogeneous Trait Effects, also icluding the three modified versions below

SPWST
SPWST.up
SPWST.down   

4.
Individual result files need to be summarized. 

