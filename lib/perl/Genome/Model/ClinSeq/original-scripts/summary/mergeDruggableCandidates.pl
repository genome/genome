#!/usr/bin/env genome-perl -w
#Written by Malachi Griffith

#Join mutation, CNV and RNS-seq DE results into a single gene candidate list
#Map gene IDs to Entrez gene IDs where possible

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;

#Required at least one
my $mutant_file = '';
my $cnv_file = '';
my $rnaseq_file = '';

#Optional
my $cnv_upper_cutoff;
my $cnv_lower_cutoff;
my $rna_upper_cutoff;
my $rna_lower_cutoff;
my $verbose = '';
my $max_names = '';


GetOptions ('mutant_file=s'=>\$mutant_file, 'cnv_file=s'=>\$cnv_file, 'rnaseq_file=s'=>\$rnaseq_file,
            'cnv_upper_cutoff=f'=>\$cnv_upper_cutoff, 'cnv_lower_cutoff=f'=>\$cnv_lower_cutoff,
            'rna_upper_cutoff=f'=>\$rna_upper_cutoff, 'rna_lower_cutoff=f'=>\$rna_lower_cutoff,
            'verbose=s'=>\$verbose, 'max_names=i'=>\$max_names);

my $usage=<<INFO;

  Example usage: 
  
  mergeDruggableCandidates.pl  --mutant_file=/gscmnt/sata132/techd/mgriffit/hg1/snvs/SomaticMutations.txt  --cnv_file=/gscmnt/sata132/techd/mgriffit/hg1/cnvs/CNView_Ensembl_v58/CNView_Ensembl_v58.tsv  --rnaseq_file=/gscmnt/sata132/techd/mgriffit/hg1/rna_seq/hg1_vs_brcs/genes_fpkm_DE.tsv
  
  Notes:
  Must specify at least one of the following input files:

  Inputs:
  --mutant_file       Summarized mutations from summarizeSnvIndelFile.pl
  --cnv_file          Summarized CNV genes from CNView.pl
  --rnaseq_file       Summarized DE RNA-seq genes from buildCufflinksFpkmMatrix.pl and cufflinksDE.R

  --cnv_upper_cutoff  CNVs with a difference greater than this value will be allowed [optional, e.g. +2.0]
  --cnv_lower_cutoff  CNVs with a difference smaller than this value will be allowed [optional, e.g. -0.5]
                      To allow gains only, simply define an upper cutoff only

  --rna_upper_cutoff  RNA DE with a difference greater than this value will be allowed [optional, e.g. +2.0]
  --rna_lower_cutoff  RNA DE with a difference smaller than this value will be allowed [optional, e.g. -0.5]
                      To allow gains only, simply define an upper cutoff only

  --max_names         To allow multiple names per gene set this to greater than 1 (default is 1)

  --verbose=1         Verbose output

INFO

unless ($mutant_file || $cnv_file || $rnaseq_file){
  print GREEN, "$usage", RESET;
  exit();
}
unless($max_names){
  $max_names = 1;
}

#Parse Entrez flatfiles and Ensembl files from BioMart
#ftp://ftp.ncbi.nih.gov/gene/DATA/gene2accession.gz
#ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz
my $clinseq_annotations_dir = "/gscmnt/sata132/techd/mgriffit/reference_annotations/";
my $entrez_dir = "EntrezGene/";
my $ensembl_dir = "EnsemblGene/";
my $entrez_ensembl_data = &loadEntrezEnsemblData('-entrez_dir'=>$entrez_dir, '-ensembl_dir'=>$ensembl_dir);

#Create one merged list of gene candidates
my %gene_list;

#Import the mutations, get list of genes.  Watch out for cases where one position is associated with multiple genes.
#For mutations affecting multiple genes, divide into seperate gene records (could be mutations on opposite strands)
my $mutant_gene_count = 0;
my %gc_mutant;
if ($mutant_file){
  open (MUTANT, "$mutant_file") || die "\n\nMutant file could not be opened: $mutant_file\n\n";
  while(<MUTANT>){
    chomp($_);
    my @line = split("\t", $_);
    my $gene_string = $line[1];
    my @genes = split(",", $gene_string);
    foreach my $gene (@genes){
      if ($gene =~ /(\S+)/){
        $gene = $1;
        my $original_name = $gene_string;

        $gene = &fixGeneName('-entrez_ensembl_data'=>$entrez_ensembl_data, '-name'=>$gene);
        $gene = &fixGeneName('-entrez_ensembl_data'=>$entrez_ensembl_data, '-name'=>$gene); #Run again to give a chance for name -> Ensembl -> Entrez mappings...

        $gc_mutant{$gene}{mutant} = 1;
        $gc_mutant{$gene}{mutant_pos} = $line[0];
        $gc_mutant{$gene}{mutant_effect} = $line[2];
        $gc_mutant{$gene}{original_name} = $original_name;
        $mutant_gene_count++;
        $gene_list{$gene}=1;
      }
    }
  }
  close(MUTANT);
}

#Import CNV genes
my $cnv_gene_count = 0;
my %gc_cnv;
if ($cnv_file){
  open (CNV, "$cnv_file") || die "\n\nCNV file could not be opened: $cnv_file\n\n";
  my $header = 1;
  while(<CNV>){
    if ($header){$header = 0; next();}
    chomp($_);
    my @line = split("\t", $_);
    my $gene = $line[0];
    my $original_name = $gene;

    $gene = &fixGeneName('-entrez_ensembl_data'=>$entrez_ensembl_data, '-name'=>$gene);
    $gene = &fixGeneName('-entrez_ensembl_data'=>$entrez_ensembl_data, '-name'=>$gene); #Run again to give a chance for name -> Ensembl -> Entrez mappings...

    my $cnv_diff = $line[5];
    my $pass_filter = 0;
    if (defined($cnv_upper_cutoff)){
      if ($cnv_diff >= $cnv_upper_cutoff){
        $pass_filter = 1;
      }
    }
    if (defined($cnv_lower_cutoff)){
      if ($cnv_diff <= $cnv_lower_cutoff){
        $pass_filter = 1;
      }
    }
    #If neither filter was defined, allow everything to pass
    unless (defined($cnv_upper_cutoff) || defined($cnv_lower_cutoff)){
      $pass_filter = 1;
    }

    #Store all data
    $gc_cnv{$gene}{cnv} = 0;
    $gc_cnv{$gene}{cnv_diff} = $cnv_diff;
    $gc_cnv{$gene}{original_name} = $gene;
    
    $gene_list{$gene}=1;

    if ($pass_filter){
      $gc_cnv{$gene}{cnv} = 1;
      $cnv_gene_count++;
    }

  }
  close(CNV);
}


#Import the RNA-seq genes
my $rna_gene_count = 0;
my %gc_rna;
if ($rnaseq_file){
  open (RNA, "$rnaseq_file") || die "\n\nRNA-seq file could not be opened: $rnaseq_file\n\n";
  my $header = 1;
  my %columns;
  while(<RNA>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header){
      $header = 0; 
      my $cp = 0;
      foreach my $col_name (@line){
        $columns{$col_name}{col_pos} = $cp;
        $cp++;
      }
      next();
    }
    my $gene_id = $line[$columns{'gene_id'}{col_pos}];
    my $gene = $line[$columns{'gene_short_name'}{col_pos}];
    my $original_name = $gene;

    $gene = &fixGeneName('-entrez_ensembl_data'=>$entrez_ensembl_data, '-name'=>$gene, '-ensembl_id'=>$gene_id);
    $gene = &fixGeneName('-entrez_ensembl_data'=>$entrez_ensembl_data, '-name'=>$gene); #Run again to give a chance for name -> Ensembl -> Entrez mappings...

    my $rna_diff = $line[$columns{'Log2Diff_Medians'}{col_pos}];
    my $pass_filter = 0;
    if (defined($rna_upper_cutoff)){
      if ($rna_diff >= $rna_upper_cutoff){
        $pass_filter = 1;
      }
    }
    if (defined($rna_lower_cutoff)){
      if ($rna_diff <= $rna_lower_cutoff){
        $pass_filter = 1;
      }
    }
    #If neither filter was defined, allow everything to pass
    unless (defined($rna_upper_cutoff) || defined($rna_lower_cutoff)){
      $pass_filter = 1;
    }

    #Store all data
    $gc_rna{$gene}{rna} = 0;
    $gc_rna{$gene}{rna_diff} = $rna_diff;
    $gc_rna{$gene}{original_name} = $original_name;

    $gene_list{$gene}=1;

    if ($pass_filter){
      $gc_rna{$gene}{rna} = 1;
      $rna_gene_count++;
    }
  }
  close(RNA);
}


#Create a non-redundant list of genes that are mutated OR CNV gained OR over-expressed (Or any combination of the three)
my %gc;

foreach my $gene (keys %gene_list){

  my $candidate = 0;
  if ($gc_mutant{$gene}){
    if ($gc_mutant{$gene}{mutant}){
      $candidate = 1;
    }
  }
  if ($gc_cnv{$gene}){
    if ($gc_cnv{$gene}{cnv}){
      $candidate = 1;
    }
  }
  if ($gc_rna{$gene}){
    if ($gc_rna{$gene}{rna}){
      $candidate = 1;
    }
  }
  unless ($candidate){
    next();
  }

  $gc{$gene}{original_name} = "na";
  $gc{$gene}{mutant} = 0;
  $gc{$gene}{mutant_pos} = "na";
  $gc{$gene}{mutant_effect} = "na";
  $gc{$gene}{cnv} = 0;
  $gc{$gene}{cnv_diff} = "na";
  $gc{$gene}{rna} = 0;
  $gc{$gene}{rna_diff} = "na";

  if ($mutant_file){
    if ($gc_mutant{$gene}){
      $gc{$gene}{mutant} = $gc_mutant{$gene}{mutant};
      $gc{$gene}{mutant_pos} = $gc_mutant{$gene}{mutant_pos};
      $gc{$gene}{mutant_effect} = $gc_mutant{$gene}{mutant_effect};
      $gc{$gene}{original_name} = $gc_mutant{$gene}{original_name};
    }
  }

  if ($cnv_file){
    if ($gc_cnv{$gene}){
      $gc{$gene}{cnv} = $gc_cnv{$gene}{cnv};
      $gc{$gene}{cnv_diff} = $gc_cnv{$gene}{cnv_diff};
      $gc{$gene}{original_name} = $gc_cnv{$gene}{original_name};
    }
  }

  if ($rnaseq_file){
    if ($gc_rna{$gene}){
      $gc{$gene}{rna} = $gc_rna{$gene}{rna};
      $gc{$gene}{rna_diff} = $gc_rna{$gene}{rna_diff};
      $gc{$gene}{original_name} = $gc_rna{$gene}{original_name};
    }
  }

  $gc{$gene}{effect_count} = $gc{$gene}{mutant} + $gc{$gene}{cnv} + $gc{$gene}{rna};
}


#Write an output file, one gene per line, keyed on gene name with the mutation, CNV, expression status of each.
#print Dumper %gc;
print "gene\toriginal_name\tmutant?\tmutant_pos\tmutant_effect\tcnv?\tcnv_diff\trna_DE?\trna_diff\teffect_count\n";
foreach my $gene (sort keys %gc){

  #If the gene has more than max_names names, skip it
  my @g = split(",", $gene);
  my $g_count = @g;
  if ($g_count <= $max_names){
    print "$gene\t$gc{$gene}{original_name}\t$gc{$gene}{mutant}\t$gc{$gene}{mutant_pos}\t$gc{$gene}{mutant_effect}\t$gc{$gene}{cnv}\t$gc{$gene}{cnv_diff}\t$gc{$gene}{rna}\t$gc{$gene}{rna_diff}\t$gc{$gene}{effect_count}\n";
  }
}

if ($verbose){
  print BLUE, "\n\nSummary\n\tMutant gene count = $mutant_gene_count\n\tCNV gene count = $cnv_gene_count\n\tRNA gene count = $rna_gene_count\n\n", RESET;
}

exit();




#######################################################################################################################################################################
#Load Entrez Data from flatfiles                                                                                                                                      #
#######################################################################################################################################################################
sub loadEntrezEnsemblData{
  my %args = @_;
  my $entrez_dir = $args{'-entrez_dir'};
  my $ensembl_dir = $args{'-ensembl_dir'};
  my %edata;

  #Check input dirs
  unless (-e $entrez_dir && -d $entrez_dir){
    print RED, "\n\nEntrez dir not valid: $entrez_dir\n\n", RESET;
    exit();
  }
  unless ($entrez_dir =~ /\/$/){
    $entrez_dir .= "/";
  }
  unless (-e $ensembl_dir && -d $ensembl_dir){
    print RED, "\n\nEnsembl dir not valid: $ensembl_dir\n\n", RESET;
    exit();
  }
  unless ($ensembl_dir =~ /\/$/){
    $ensembl_dir .= "/";
  }

  #Load data from Ensembl files
  my %entrez_map;      #Entrez_id          -> symbol, synonyms
  my %ensembl_map;     #Ensembl_id         -> entrez_id(s) - from Entrez
  my %ensembl_map2;    #Ensembl_id         -> symbol(s) - from Ensembl
  my %symbols_map;     #Symbols            -> entrez_id(s)
  my %synonyms_map;    #Synonyms           -> entrez_id(s)
  my %p_acc_map;       #Protein accessions -> entrez_id(s)
  my %g_acc_map;       #Genomic accessions -> entrez_id(s)

  my $gene2accession_file = "$entrez_dir"."gene2accession.human";
  my $gene_info_file = "$entrez_dir"."gene_info.human";
  open (GENE, "$gene_info_file") || die "\n\nCould not open gene_info file: $gene_info_file\n\n";
  while(<GENE>){
    chomp($_);
    if ($_ =~ /^\#/){
      next();
    }
    my @line = split("\t", $_);
    my $tax_id = $line[0];
    #Skip all non-human records
    unless ($tax_id eq "9606"){
      next();
    }
    my $entrez_id = $line[1];
    my $symbol = $line[2];
    my $synonyms = $line[4];
    my $ext_ids = $line[5];

    #Get synonyms for each gene and divide each into a unique hash
    if ($synonyms eq "-"){
      $synonyms = "na";
    }
    my @synonyms_array = split("\\|", $synonyms);
    my %synonyms_hash;   
    foreach my $syn (@synonyms_array){
      $synonyms_hash{$syn} = 1;
    }

    #Parse the external IDs field for Ensembl gene IDs (Other possibilites include HGNC, MIM, HPRD)
    my %ensembl_hash;
    my @ext_ids_array = split("\\|", $ext_ids);
    $entrez_map{$entrez_id}{ensembl_id} = "na";
    foreach my $ext_string (@ext_ids_array){
      if ($ext_string =~ /Ensembl/i){
        if ($ext_string =~ /Ensembl\:(\w+)/){
          $entrez_map{$entrez_id}{ensembl_id} = $1;
          $ensembl_hash{$1} = 1;
        }else{
          print RED, "\n\nFormat of Ensembl field not understood: $ext_string\n\n", RESET;
          exit();
        }   
      }else{
        next();
      }
    }

    #Store entrez info keyed on entrez id
    #print "\n$entrez_id\t$symbol\t@synonyms_array";
    $entrez_map{$entrez_id}{symbol} = $symbol;
    $entrez_map{$entrez_id}{synonyms_string} = $synonyms;
    $entrez_map{$entrez_id}{synonyms_array} = \@synonyms_array;
    $entrez_map{$entrez_id}{synonyms_hash} = \%synonyms_hash;

    #Store entrez info keyed on symbol
    #print "\n$symbol\t$entrez_id";
    if ($symbols_map{$symbol}){
      my $ids = $symbols_map{$symbol}{entrez_ids};
      $ids->{$entrez_id} = 1;
    }else{
      my %tmp;
      $tmp{$entrez_id} = 1;
      $symbols_map{$symbol}{entrez_ids} = \%tmp;
    }

    #Store synonym to entrez_id mappings
    foreach my $syn (@synonyms_array){
      if ($synonyms_map{$syn}){
        my $ids = $synonyms_map{$syn}{entrez_ids};
        $ids->{$entrez_id} = 1;
      }else{
        my %tmp;
        $tmp{$entrez_id} = 1;
        $synonyms_map{$syn}{entrez_ids} = \%tmp;
      }
    }

    #Store ensembl to entrez_id mappings
    foreach my $ens (sort keys %ensembl_hash){
      if ($ensembl_map{$ens}){
        my $ids = $ensembl_map{$ens}{entrez_ids};
        $ids->{$entrez_id} = 1;
      }else{
        my %tmp;
        $tmp{$entrez_id} = 1;
        $ensembl_map{$ens}{entrez_ids} = \%tmp;
      }
    }
  }
  close (GENE);

  open (ACC, "$gene2accession_file") || die "\n\nCould not open gene2accession file: $gene2accession_file\n\n";
  while(<ACC>){
    chomp($_);
    if ($_ =~ /^\#/){
      next();
    }
    my @line = split("\t", $_);
    my $tax_id = $line[0];
    #Skip all non-human records
    unless ($tax_id eq "9606"){
      next();
    }
    my $entrez_id = $line[1];
    my $prot_id = $line[5];
    my $genome_id = $line[7];

    #Protein accession
    unless ($prot_id eq "-"){
      #If the prot is not defined, skip
      #Clip the version number
      if ($prot_id =~ /(\w+)\.\d+/){
        $prot_id = $1;
      }
      #print "\n$entrez_id\t$prot_id";
      if ($p_acc_map{$prot_id}){
        my $ids = $p_acc_map{$prot_id}{entrez_ids};
        $ids->{$entrez_id} = 1;
      }else{
        my %tmp;
        $tmp{$entrez_id} = 1;
        $p_acc_map{$prot_id}{entrez_ids} = \%tmp;
      }
    }

    #Genomic accession
    unless ($genome_id eq "-"){
      #If the genome accession is not defined, skip
      #Clip the version number
      if ($genome_id =~ /(\w+)\.\d+/){
        $genome_id = $1;
      }
      if ($g_acc_map{$genome_id}){
        my $ids = $g_acc_map{$genome_id}{entrez_ids};
        $ids->{$entrez_id} = 1;
      }else{
        my %tmp;
        $tmp{$entrez_id} = 1;
        $g_acc_map{$genome_id}{entrez_ids} = \%tmp;
      }
    }
  }
  close (ACC);

  #print Dumper %entrez_map;
  #print Dumper %symbols_map;
  #print Dumper %synonyms_map;
  #print Dumper %p_acc_map;

  #Now load ensembl gene id to gene name mappings from a series of legacy ensembl versions
  #Give preference to latest build
  my @files = qw (Ensembl_Genes_Human_v63.txt Ensembl_Genes_Human_v62.txt Ensembl_Genes_Human_v61.txt Ensembl_Genes_Human_v60.txt Ensembl_Genes_Human_v59.txt Ensembl_Genes_Human_v58.txt Ensembl_Genes_Human_v56.txt Ensembl_Genes_Human_v55.txt Ensembl_Genes_Human_v54.txt Ensembl_Genes_Human_v53.txt Ensembl_Genes_Human_v52.txt Ensembl_Genes_Human_v51.txt);

  foreach my $file (@files){
    my $path = "$ensembl_dir"."$file";
    open (ENSG, "$path") || die "\n\nCould not open file: $path\n\n";
    while(<ENSG>){
      chomp($_);
      my @line = split("\t", $_);
      my $ensg_id = $line[0];
      my $ensg_name = $line[1];
      if ($ensg_name =~ /(.*)\.\d+$/){
        $ensg_name = $1;
      }

      unless($ensembl_map2{$ensg_id}){
        $ensembl_map2{$ensg_id}{name}=$ensg_name;
        $ensembl_map2{$ensg_id}{source}=$file;
      }
    }
    close(ENSG);
  }

  $edata{'entrez_ids'} = \%entrez_map;
  $edata{'ensembl_ids'} = \%ensembl_map;
  $edata{'ensembl_ids2'} = \%ensembl_map2;
  $edata{'symbols'} = \%symbols_map;
  $edata{'synonyms'} = \%synonyms_map;
  $edata{'protein_accessions'} = \%p_acc_map;
  $edata{'genome_accessions'} = \%g_acc_map;

  return(\%edata);
}


#######################################################################################################################################################################
#If possible translate the current gene name or ID into an official gene name from Entrez                                                                             #
#######################################################################################################################################################################
sub fixGeneName{
  my %args = @_;
  my $edata = $args{'-entrez_ensembl_data'};
  my $original_name = $args{'-name'};
  my $ensembl_id;
  if (defined($args{'-ensembl_id'})){
     $ensembl_id = $args{'-ensembl_id'};
  }
  
  #Unless a better match is found, the original name will be returned
  my $corrected_name = $original_name; 

  #If the incoming gene name has a trailing version number, strip it off before comparison
  if ($original_name =~ /(.*)\.\d+$/){
    $original_name = $1;
  }

  #Load the mapping hashes
  my $entrez_map = $edata->{'entrez_ids'};
  my $ensembl_map = $edata->{'ensembl_ids'};
  my $ensembl_map2 = $edata->{'ensembl_ids2'};
  my $symbols_map = $edata->{'symbols'};
  my $synonyms_map = $edata->{'synonyms'};
  my $prot_acc_map = $edata->{'protein_accessions'};
  my $genome_acc_map = $edata->{'genome_accessions'};

  my $any_match = 0;
  my @entrez_symbols;
  my $entrez_name_string = '';

  #Try mapping directly to the entrez symbols
  my $entrez_match = 0;
  if ($symbols_map->{$original_name}){
    $entrez_match = 1;
    $any_match = 1;
    my $entrez_ids = $symbols_map->{$original_name}->{entrez_ids};
    foreach my $entrez_id (keys %{$entrez_ids}){
      my $entrez_symbol = $entrez_map->{$entrez_id}->{symbol};
      push (@entrez_symbols, $entrez_symbol);
    }
  }
  if ($entrez_match){
    $entrez_name_string = join(",", @entrez_symbols);
    $corrected_name = $entrez_name_string;
  }

  #Unless a match was already found, try mapping to ensembl IDs and then to entrez symbols
  #This assumes that the 'name' reported is actually an ensembl ID, something that happens routinely in the somatic variation pipeline...
  my $ensembl_match = 0;
  unless ($any_match){
    if ($ensembl_map->{$original_name}){
      $ensembl_match = 1;
      $any_match = 1;
      my $entrez_ids = $ensembl_map->{$original_name}->{entrez_ids};
      foreach my $entrez_id (keys %{$entrez_ids}){
        my $entrez_symbol = $entrez_map->{$entrez_id}->{symbol};
        push (@entrez_symbols, $entrez_symbol);
      }
    }
    if ($ensembl_match){
      $entrez_name_string = join(",", @entrez_symbols);
      $corrected_name = $entrez_name_string;
    }
  }


  #Unless a match was already found, try mapping to ensembl IDs (from Ensembl) and then to Ensembl symbols
  unless ($any_match){
    if ($ensembl_map2->{$original_name}){
      $ensembl_match = 1;
      $any_match = 1;
      $corrected_name = $ensembl_map2->{$original_name}->{name};
    }
  }

  #Unless a match was already found, try mapping to protein accession IDs, and then to Entrez symbols
  unless ($any_match){
    my $protein_acc_match = 0;
    if ($prot_acc_map->{$original_name}){
      $protein_acc_match = 1;
      $any_match = 1;
      my $entrez_ids = $prot_acc_map->{$original_name}->{entrez_ids};
      foreach my $entrez_id (keys %{$entrez_ids}){
        my $entrez_symbol = $entrez_map->{$entrez_id}->{symbol};
        push (@entrez_symbols, $entrez_symbol);
      }
    }
    if ($protein_acc_match){
      $entrez_name_string = join(",", @entrez_symbols);
      $corrected_name = $entrez_name_string;
    }
  }

  #Unless a match was already found, try mapping to genome IDs, and then to Entrez symbols
  unless ($any_match){
    my $genome_acc_match = 0;
    if ($genome_acc_map->{$original_name}){
      $genome_acc_match = 1;
      $any_match = 1;
      my $entrez_ids = $genome_acc_map->{$original_name}->{entrez_ids};
      foreach my $entrez_id (keys %{$entrez_ids}){
        my $entrez_symbol = $entrez_map->{$entrez_id}->{symbol};
        push (@entrez_symbols, $entrez_symbol);
      }
    }
    if ($genome_acc_match){
      $entrez_name_string = join(",", @entrez_symbols);
      $corrected_name = $entrez_name_string;
    }
  }



  #Unless a match was already found, try mapping to Entrez synonyms, and then to Entrez symbols
  #Only allow 1-to-1 matches for synonyms...
  unless ($any_match){
    my $synonyms_match = 0;
    if ($synonyms_map->{$original_name}){
      my $entrez_ids = $synonyms_map->{$original_name}->{entrez_ids};
      my $match_count = keys %{$entrez_ids};
      if ($match_count == 1){
        $synonyms_match = 1;
        $any_match = 1;
        foreach my $entrez_id (keys %{$entrez_ids}){
          my $entrez_symbol = $entrez_map->{$entrez_id}->{symbol};
          push (@entrez_symbols, $entrez_symbol);
        }
      }
    }
    if ($synonyms_match){
      $entrez_name_string = join(",", @entrez_symbols);
      $corrected_name = $entrez_name_string;
    }
  }

  #Unless a match was already found, try mapping to ensembl IDs (from Entrez) and then to entrez symbols - starting with an actual ensembl ID supplied separately
  if ($ensembl_id){
    unless ($any_match){
      if ($ensembl_map->{$ensembl_id}){
        $ensembl_match = 1;
        $any_match = 1;
        my $entrez_ids = $ensembl_map->{$ensembl_id}->{entrez_ids};
        foreach my $entrez_id (keys %{$entrez_ids}){
          my $entrez_symbol = $entrez_map->{$entrez_id}->{symbol};
          push (@entrez_symbols, $entrez_symbol);
        }
      }
      if ($ensembl_match){
        $entrez_name_string = join(",", @entrez_symbols);
        $corrected_name = $entrez_name_string;
      }
    }
  }

  #Unless a match was already found, try mapping to ensembl IDs (from Ensembl) and then to Ensembl symbols - starting with an actual ensembl ID supplied separately
  if ($ensembl_id){
    unless ($any_match){
      if ($ensembl_map2->{$ensembl_id}){
        $ensembl_match = 1;
        $any_match = 1;
        $corrected_name = $ensembl_map2->{$ensembl_id}->{name};
      }
    }
  }


  if ($verbose){
    if ($entrez_name_string eq $original_name){
      print BLUE, "\nSimple Entrez match: $original_name -> $corrected_name", RESET;
    }elsif($corrected_name eq $original_name){
      print YELLOW, "\nNo matches: $original_name -> $corrected_name", RESET;
    }else{
      print GREEN, "\nFixed name: $original_name -> $corrected_name", RESET;
    }
  }
  return($corrected_name);
}





