#CLINSEQ CODE PATHS
/gscmnt/sata206/techd/git/genome/lib/perl/Genome/Model/ClinSeq/Commands/*.pm *.t
/gscmnt/sata206/techd/git/genome/lib/perl/Genome/Model/ClinSeq/Commands/Converge/*.pm *.t
/gscmnt/sata206/techd/git/genome/lib/perl/Genome/ProcessingProfile/ClinSeq.pm
/gscmnt/sata206/techd/git/genome/lib/perl/Genome/ProcessingProfile/ClinSeq.t
/gscmnt/sata206/techd/git/genome/lib/perl/Genome/Model/ClinSeq.pm
/gscmnt/sata206/techd/git/genome/lib/perl/Genome/Model/ClinSeq.t

#CLINSEQ DATA PATHS
#(data files used by the pipeline - not the input data, but annotation data, etc.)
/gscmnt/sata132/techd/mgriffit/reference_annotations/
/gscmnt/sata132/techd/mgriffit/reference_annotations/EnsemblGene/
/gscmnt/sata132/techd/mgriffit/reference_annotations/EntrezGene/
/gscmnt/sata132/techd/mgriffit/reference_annotations/GeneSymbolLists/
/gscmnt/sata132/techd/mgriffit/reference_annotations/hg18/
/gscmnt/sata132/techd/mgriffit/reference_annotations/hg18/ideogram/
/gscmnt/sata132/techd/mgriffit/reference_annotations/hg18/transcript_to_gene/
/gscmnt/sata132/techd/mgriffit/reference_annotations/hg19/
/gscmnt/sata132/techd/mgriffit/reference_annotations/hg19/ideogram/
/gscmnt/sata132/techd/mgriffit/reference_annotations/hg19/transcript_to_gene/
/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/DrugBank/query_files/

#COMPLETED TASKS
#As code writing tasks are completed, summarize them here
#Each of the following has been integrated into a single 'ClinSeq pipeline' - i.e. everything here gets run by a single command initiated for a single patient
#- Create methods to obtain and perform basic sanity checks on somatic variation and RNA-seq genome models
#- Summarize SNV lists.  Identify the subsets of SNVs that are protein coding relevant - produce a compact summary of each
#- Summarize InDel lists.  Similar to the summary for SNVs
#- Annotation of SNV and InDel lists at the gene level to identify potentially clinically relevant gene categories (e.g. kinases, druggable, etc.)
#- Create a 'database' of useful gene categories (kinases, phosphatases, etc.) using GO, DrugBank, etc.
#- Create method to intersect arbitrary lists of candidate genes with gene categories of interest
#- Create methods to parse and import DrugBank data and report on the drug interactions for candidate mutated, amplified, over-expressed genes
#- Create filtered versions of DrugBank that remove spurious or non-cancer relevant drug-gene interactions
#- Create method to identify copy number amplified and deleted regions on a gene-by-gene basis
#- Create visualizations of copy number results that are genome-wide, chromosome-by-chromosome, and specific genes or regions of the genome
#- Get read counts supporting SNVs called by WGS or Exome and intersect with read counts from RNA-Seq data
#- Summarize results from an RNA-seq model run for a single patient (organize into structured gene and isoform expression matrices)
#- Perform 'absolute outlier' gene/isoform expression analysis for a single patient and identify the most highly expressed kinases, druggable genes, etc.
#- Create visualizations of the most highly expressed genes as well as a predefined set of genes of interest
#- Create code to visualize chromosome ideograms and identify the cytoband of any gene or set of coordinates (in the format recognized by pathologists)
#- Automate creation of clonality plots
#- Automate identification of WGS, Exome and RNAseq datasets that have completed sequencing and are ready for ClinSeq analysis
#- Create methods to extract and parse the complete ClinicalTrials.gov database
#- Create methods to extract fastq files from BAM and trim reads for certain downstream analyses (defuse, ALEXA-seq, etc.)
#- Summary of TopHat RNA-seq alignment results
#  - Summarize read mapping, known vs. novel junctions, proportion of reads mapping across junctions, proportion of MT mapping reads
#- RNA-seq variant validation
#  - BAM read counts ... WGS, Exome, RNA
#- Added a detailed summary of builds going into a ClinSeq model, including extensive QC values
#- Added automatic creation of IGV XML sessions for each ClinSeq model
#- Refactored BamRead counts commands to use method code instead of shelling out to a system call
#- Run pairoscope on all predicted gene fusions from SV files
#- Scatter plot of Variant Allele Frequency in WGS Tumor vs Exome Tumor - report the R^2 and n
#- Scatter plot of Variant Allele Frequency in WGS Normal vs Exome Normal - report the R^2 and n
#- Experiment with some ways to visualize mutation variant allele frequency in combination with RNA expression
#  - WGS VAF (total RC) + Exome VAF (total RC) + RNA-seq VAF (total RC) + RNA-seq Gene FPKM...
#  - scatter plot of Exome VAF vs. RNAseq VAF and then color points according to gene FPKM


#TODO / FEATURE WISH LIST
#- Cosmetic.  Stop using Ansicolor warning and status messages in ClinSeq code.  Use $self->status_message and $self->warning_message instead

#- Summary stats.
#  - What is the non-synonymous to synonymous mutation rate in somatically mutated gene (i.e. tier1 mutations)
#    - i.e. what is the mutation rate per megabase: total tier 1, non-synonymous tier1, synonymous tier 1, total all tiers, etc

#- Generate lolliplot images for all SNVs to visualize the position of each variant
#  - Create a paired plot that shows the loliplot for data for patients that have been aggregated across projects
#- Create vcf formatted version of all SNVs combined, then tier 1,2,3 individually
#- Get all BAM read counts for all SNVs for convenience in downstream analysis (WGS, Exome & RNAseq) (Tumor and Normal)
#  - Create a grand list of all SNV coords from WGS and Exome, note whether each came from WGS, Exome or BOTH, Tier, read counts from all BAMs
#- CNV amplified / deleted + over-expressed / under-expressed
#- WGS vs. Exome Venn Diagrams : For SNVs & Indels
#- Mutated and expressed vs. not expressed summary
#  - Cufflinks expression.  FPKM + percentile
#  - Assess RNA-seq BAM read counts for Tier 2,3 SNVs.  High values might indicate an erroneous classification as Tier 2,3.  In theory RNA-seq should only cover Tier 1 ... 
#- RNA-seq outlier analysis
#  - Differential / relative comparisons to other samples / tumors of the same type
#- SVs.  Annotation strategies, validation, filtering

#- RNA-seq gene fusions
#  - Tophat fusion, Chimera scan
#- Previously discovered variants
#- Germline variants

#- Improve CNV analysis
#  - Find partial gene amplifications / deletions
#    - Use 'gmt copy-number cbs' and 'gmt copy-number cna-seg' to find segments that are copy-number amplified/deleted
#    - Find the genes that overlap these regions
#    - Make a new CNV gene script that produces one CNV value per gene & per transcript.  Use Ensembl annotations only.  Get annotations from annotator!

#- Enhanced SNV annotation
#  - Sift, PolyPHEN, or similar
#  - SNPeff
#  - Overlap of observed mutations with Cosmic / OMIM sites (use MUSIC?)
#  - Overlap of observed mutations with dbSNP
#  - Overlap of observed mutations with the NHLBI Exome Sequencing Project - Exome variant server
#  - Overlap of observed mutations with TGI recurrent sites
#  - Run snpEff on all variants
#  - Run the Ensembl variant effect predictor (VEP) on all variants
#  - Run regulatory annotation on Tier2 variants using: gmt annotate regulatory-features
#- Summarize gene annotation results (i.e. for each gene list, each gene is marked as Kinase, RTK, trancription factor, etc.)
#  - Summarize for each event list, how many genes belong to each catergory
#  - Converge results to gene-level for SNVs, InDels, CNVs (amps and dels), SVs, RNA-seq
#    - Similar code to what is already being done for converging of drug gene results across multiple patients


#- Druggable genes analysis
#  - 
#  - Are there any mutations that are in a kinase where that kinase itself is not druggable but:
#    - the region mutated is homologous to a mutation site in another kinase
#    - e.g. A HER2 mutation in a kinase domain and the corresponding region of that kinase domain is mutated in EGFR in another cancer type etc.
#    - Suggests that we could use an EGFR inhibitor to treat patients with this HER2 mutation.
#    - In general establish the cross-reactivity of kinase inhibitors or links that can be made between kinases by recurrence of mutation sites.

#SUMMARY STATS
#A descriptive stats, question/answer table
#Every record in the table will have the following:
#Question | Answer | Data type | analysis type | statistic type | Extra description
#What is the median coverage of mutated positions in the exome data? | 10 | Exome | SNV | Median | NA

#Question list:
#- Median coverage of SNV positions for WGS & Exome seperately, for both the Tumor and the Normal sample
#- Overall coverage of the genome or exome for both Tumor and Normal
#- Number and proportion of TopHat mapped reads for RNA-seq library
#- Number and Proportion of TopHat mapped reads for RNA-seq library that mapped to splice junctions
#- Number and Proportion of TopHat mapped reads for RNA-seq library that mapped to known splice junctions


#CLINSEQ - CLINICAL SEQUENCING REPORTS
#Executive summary
#- The most clinically relevant observations
#- The highest priority druggable events and diagnostic / prognostic observations
#- Summarize sample and data QC only if there are serious concerns that effect interpretation of the results

#Detailed decision tree summary
#- The positive/negative status of all steps in the decision tree
#- Where positive, summarize the event and the confidence we can assign to it and the resulting treatment recommendation
#- Where negative, summarize the power we had to assess that branch of the tree (e.g. no kinase mutations were detected but WGS coverage was low)
#- How deeply into the tree did we have to traverse to achieve a high scoring positive observation
#- The structure of the tree may have cancer site (or even subtype) specific properties
#- The structure of the tree will evolve as we go through the ClinSeq process with increasing numbers of patients
#  - For example, for cancer type X we may learn over time that mutations are the highest prioroty types of events but in cancer type Y it may be focal amplifications
#  - The initial structure of the tree will be based on existing knowledge from the literature and consulting with medical experts for each cancer type
#  - As ClinSeq data becomes available, we can informatically refactor the decision trees (e.g. by machine learning approaches)
#- The major branches of the tree might be: druggable, prognostic, diagnostic, treatment response, discovery
#- Each branch has sub-branches that effectively ask questions that drill down to progressively more detailed inquiries
#  - e.g. Is there a somatic mutation? Affecting protein sequence? Does the gene have an approved drug? At a known pathological site that responds to the drug?
#  - For each question, if the answer is NO, further questions may be asked that may still lead to a clinical action (but likely a lower priority one) 

#Clinical data summary
#- What do we need to know about the patient tumor sample(s) to help us interpret and maximize the value of the genome sequencing and analysis?
#- Some info from the can be directly assessed/corroborated by the WGS, Exome, RNAseq analysis...
#  - e.g. HER2 protein expression status from clinic -> HER2 genome amplification status from WGS, RNA expression status from RNA-seq
#  - e.g. BCR-ABL translocation status from clinic -> exact structure of translocation in both the genome and transcriptome from WGS and RNA-seq respectively
#  - and so on, there are many, many things being done in the clinic that we can verify and in some cases we can get a higher resolution assessment...
#- Pathology - show FFPE H&E image for the tumor.  If LCM was used, show the slide before/after LCM
#- Other clinical data on the sample. e.g. For breast cancer, results from ER, PR, HER2 staining
#- If sorting was used on the tumor cells (e.g. in ALL) - summarize the sort efficiency
#- Comments on expected purity and heterogeneity from the pathologist, surgeon, oncologist perspective
#- Comments on prior treatment status of the patient (what chemotherapy have they received before some procurement) from oncologist
#- Cytogenetics. Overall karyotype status.  Specific events assayed for the tumor type.  
#- e.g. overall genome stability, deletions, amplifications, translocations, microsatellite instability, etc.

#Production sequencing summary
#- Sample QC metrics and reports - Agilent traces etc.
#- Library types, coverage level, library complexity, etc.
#- Raw QC results -> sequencing error rates, read mapping efficiency
#- dbSNP concordance
#- Sequencing vs. microarray concordance -> sample swaps?, purity?
#- Etc.

#Production analysis summary
#- Analysis methods used -> processing profiles and parameters

#Reference alignment/assembly summary
#- How well have we interogated the genome, exome and transcriptome
#- i.e. what is our expectation of sensitivity...

#Somatic variation summary
#- Overview of SNV, InDel, CNV, SV, 'aberrant' gene expression

#Clinical sequencing summary
#- A synthesis of the medical/clinical implications of all analysis conducted above
#- The executive summary and decision tree summary will be compact representations of this analysis
#- This section will have the complete results from which these were generated
#- Candidate molecular event and gene lists from which clinically actionable events will be drawn
#- Where 'clinically actionable' implies some suggestion or modification of patient treatment
#- Many possible forms:
#  - Therapy that targets an activating mutation, amplified gene, over-expressed gene, activated pathway, etc.
#  - Observations of prognostic value -> e.g. low, medium, high risk of progression/metastasis
#  - Observations of diagnostic value -> subtyping of tumors
#  - Observations that influence dosing of therapies that are already part of standard of care
#- Observations relating to existing clinical trials

#Discovery summary
#- Observations that are not immediately clinically actionable but have a potential value in the development of novel therapies
#- Novel recurrent events
#- Novel private events in high impact genes
#- Previously known events that are novel to the current cancer type
#- Etc.

