#!/usr/bin/env genome-perl
#Written by Malachi Griffith

#Before using this script, make sure the version of Ensembl you wish to use has been imported (both the mysql database and the corresponding mysql API)

#The user will supply:
#An ensembl version 
#An annotation build directory (e.g. /gscmnt/ams1102/info/model_data/2772828715/build106409296/annotation_data - NCBI-human.ensembl/58_37c_v2)
#A reference assembly build directory (e.g. /gscmnt/ams1102/info/model_data/2869585698/build106942997/ - GRCh37-lite-build37)

#1. Get chromosomes allowed for the specified build - get this from the actual bowtie index
#2. Get the master GTF for this version of Ensembl (either as a flat file from Ensembl FTP, or using the one already in the annotation build dir)
#3. Get all transcript IDs from Ensembl and divide them into categories by Biotype
#4. Produce filtered versions of this GTF file.  MT.gtf, rRNA.gtf, pseudo.gtf, rRNA_MT.gtf, rRNA_MT_pseudo.gtf

#Using hs_60_37e
#Possible ensembl GENE biotypes (for the chromosomes used in GRCh37-lite)
#IG_C_gene, IG_D_gene, IG_J_gene, IG_J_pseudogene, IG_V_gene, IG_V_pseudogene, IG_pseudogene, Mt_rRNA, Mt_tRNA, Mt_tRNA_pseudogene, TEC, lincRNA, miRNA, miRNA_pseudogene, misc_RNA, misc_RNA_pseudogene, non_coding, polymorphic_pseudogene, processed_transcript, protein_coding, pseudogene, rRNA, rRNA_pseudogene, scRNA_pseudogene, snRNA, snRNA_pseudogene, snoRNA, snoRNA_pseudogene, tRNA_pseudogene

#Possible ensembl TRANSCRIPT biotypes (for genes on the chromosomes used in GRCh37-lite)
#IG_C_gene, IG_D_gene, IG_J_gene, IG_J_pseudogene, IG_V_gene, IG_V_pseudogene, IG_pseudogene, Mt_rRNA, Mt_tRNA, Mt_tRNA_pseudogene, TEC, TR_gene, TR_pseudogene, ambiguous_orf, antisense, lincRNA, miRNA, miRNA_pseudogene, misc_RNA, misc_RNA_pseudogene, non_coding, nonsense_mediated_decay, polymorphic_pseudogene, processed_pseudogene, processed_transcript, protein_coding, pseudogene, rRNA, rRNA_pseudogene, retained_intron, retrotransposed, scRNA_pseudogene, snRNA, snRNA_pseudogene, snoRNA, snoRNA_pseudogene, tRNA_pseudogene, transcribed_processed_pseudogene, transcribed_unprocessed_pseudogene, unitary_pseudogene, unprocessed_pseudogene

#Define classes 'Mitochondrial', 'Ribosomal', 'Pseudogene'
#Mitochondrial genes:  On MT chromosome or of type: Mt_* (gene or transcript biotype)
#Ribosomal genes: rRNA 
#Pseudogenes: pseudogene and *_pseudogene

use DBI;
use strict;
use Data::Dumper;
use Getopt::Long;
use Benchmark;
use Term::ANSIColor qw(:constants);
use File::Basename;
use above 'Genome'; # remove 'above' when this is turned into a module

my $ensembl_version = '';
my $reference_build_dir = '';
my $genes_gtf = '';
my $species = '';
my $ribosomal_gene_ids = '';

GetOptions ('ensembl_version=s'=>\$ensembl_version, 'reference_build_dir=s'=>\$reference_build_dir, 'genes_gtf=s'=>\$genes_gtf, 'species=s'=>\$species, 'ribosomal_gene_ids=s'=>\$ribosomal_gene_ids);

unless ($ensembl_version && $reference_build_dir && $genes_gtf && $species && $ribosomal_gene_ids){
  print RED, "\n\nParameter missing", RESET;
  print GREEN, "\n\nUsage\n\ngetEnsemblRiboMtTranscripts.pl  --ensembl_version=58_37c  --reference_build_dir=/gscmnt/ams1102/info/model_data/2869585698/build106942997/  --genes_gtf=/gscmnt/ams1102/info/model_data/2772828715/build106409296/annotation_data/all_sequences.gtf  --species=homo_sapiens  --ribosomal_gene_ids=/gscmnt/ams1102/info/model_data/2772828715/build106409296/annotation_data/HumanRibosomalGeneNames.txt\n\n", RESET;
  exit();
}

#Get the ensembl API version from the user specified ensembl database version
my $ensembl_api_version;
if ($ensembl_version =~ /^(\d+)\_.*/){
  $ensembl_api_version = $1;
}else{
  print RED, "\n\nCould not determine Ensembl API version from Ensembl database version: $ensembl_version\n\n", RESET;
  exit();
}
  
#Check dirs
unless ($reference_build_dir =~ /\/$/){
  $reference_build_dir .= "/";
}
unless (-e $genes_gtf){
  print RED, "\n\nGenes gtf file not found: $genes_gtf\n\n", RESET;
  exit();
}
unless (-e $ribosomal_gene_ids){
  print RED, "\n\nNeed a list of ribosomal gene names for this species (where each line is: gene_name\tensg_id)\n\n", RESET;
  exit();
}

#Load Ensembl API for this ensembl version
my $api_path = "/gsc/scripts/share/ensembl-"."$ensembl_api_version/ensembl/modules";
my $bioperl_path = "/gsc/lib/perl5/bioperl/1.6.0/lib/perl5/";
unshift(@INC, "$api_path");
unshift(@INC, "$bioperl_path");
require Bio::EnsEMBL::DBSQL::DBAdaptor; #Used for local connections to EnsEMBL core databases

#1. Get chromosomes allowed for the specified build - get this from the actual bowtie index
my %bowtie_chrs;
my $bowtie_index = "$reference_build_dir";
my $cmd = "bowtie-inspect --names $reference_build_dir"."all_sequences.bowtie";
my @result = `$cmd`;
chomp(@result);
my $chr_count = scalar(@result);
foreach my $chr_line (@result){
  if ($chr_line =~ /^(\S+)\s+(.*)/){
    $bowtie_chrs{$1}{desc} = $2;
  }
}

#Load a list of ribosomal gene names and ids
my %ribo_names;
my %ribo_ids;
open (RIBO, "$ribosomal_gene_ids") || die "\n\nCould not load ribosomal gene ids: $ribosomal_gene_ids\n\n";
while(<RIBO>){
  chomp($_);
  my @line=split("\t", $_);
  $ribo_names{$line[0]}=1;
  $ribo_ids{$line[1]}=1;
}
close(RIBO);


#2. Get the master GTF for this version of Ensembl (either as a flat file from Ensembl FTP, or using the one already in the annotation build dir)
#Supplied by the user:
#Getting one directly from ensembl for testing purposes only:
#wget ftp://ftp.ensembl.org/pub/release-58/gtf/homo_sapiens/Homo_sapiens.GRCh37.58.gtf.gz
#Of from the annotation build dir:
#/gscmnt/ams1102/info/model_data/2772828715/build106409296/annotation_data/all_sequences.gtf


#3. Get all transcript IDs from Ensembl and divide them into categories by Biotype
my %ensg_ids;

#Get a normal DBI connection to the ensembl database
my $ensembl_database="$species"."_core_"."$ensembl_version";
my $ensembl_server="mysql1.gsc.wustl.edu";
my $ensembl_user="mse";
my $ensembl_password="";
my $ensembl_dbh = &connectDB('-database'=>$ensembl_database, '-server'=>$ensembl_server, '-user'=>$ensembl_user, '-password'=>$ensembl_password);

#Now use Ensembl API methods to get all transcripts for all genes and the biotype for each individual transcript
my $registry = "Bio::EnsEMBL::Registry";
$registry->load_registry_from_db(-host =>$ensembl_server, -user =>$ensembl_user, -pass=>$ensembl_password, -db_version =>$ensembl_api_version);

#Get 'slices' one at a time by chromosome, or failing that, by supercontig.  Using the list of chromosomes from the bowtie index as targets
#Make sure you can get a slice for every chromosome before proceeding
my $slice_adaptor = $registry->get_adaptor($species, 'Core', 'Slice');
foreach my $chr (sort keys %bowtie_chrs){
  my $slice;
  $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr );

  if ($slice_adaptor->fetch_by_region( 'chromosome', $chr )){
    #print BLUE, "\n\tChromosome slice found for: $chr", RESET;
  }elsif($slice_adaptor->fetch_by_region( 'supercontig', $chr )){
    #print BLUE, "\n\tSupercontig slice found for: $chr", RESET;
  }else{
    print RED, "\n\tSlice not found for: $chr", RESET;
    exit();
  }
}


#Store info from Ensembl API in a file.  Check for local copy of that before querying the database (slow...)
#Storing this info also makes it quicker to generate the ribo, mito, pseudo files for multiple input genes.gtf files...
my $transcript_info_file = "$species"."_"."$ensembl_version".".transcript_info.tsv";
if (-e $transcript_info_file){
  print YELLOW, "\n\nEnsembl transcript info file already exists - using that instead of importing from the API...\n\n", RESET;
}else{
  open(INFO, ">$transcript_info_file") || die "\n\nCould not open ensembl transcript info file for writing: $transcript_info_file\n\n";
  my $g = 0;
  my $t = 0;
  foreach my $chr (sort keys %bowtie_chrs){
    my $slice;
    if ($slice_adaptor->fetch_by_region( 'chromosome', $chr )){
      $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr );
    }elsif($slice_adaptor->fetch_by_region( 'supercontig', $chr )){
      $slice = $slice_adaptor->fetch_by_region( 'supercontig', $chr );
    }else{
      print RED, "\n\tSlice not found for: $chr", RESET;
      exit();
    }
    my @genes = @{ $slice->get_all_Genes() };
    foreach my $gene (@genes){
      $g++;
      #Get info on the gene object
      my $ensembl_g_id = $gene->stable_id();
      my $g_biotype = $gene->biotype();
      my $coord_sys  = $slice->coord_system()->name();
      my $seq_region = $slice->seq_region_name();
      my $gene_name = $gene->external_name();
      unless ($gene->is_known()){
        $gene_name = "Unknown";
      }

      #Get the transcripts associated with this gene object
      my %transcripts;
      my @trans_list = @{$gene->get_all_Transcripts()};
      my $trans_count = scalar(@trans_list);
      #Get the biotype for each transcript
      foreach my $transcript (@trans_list){
        $t++;
        my $ensembl_t_id = $transcript->stable_id();
        my $t_biotype = $transcript->biotype();
        print BLUE, "$g\t$ensembl_g_id\t$coord_sys\t$seq_region\t$g_biotype\t$t\t$ensembl_t_id\t$t_biotype\t$gene_name\n", RESET;
        print INFO "$g\t$ensembl_g_id\t$coord_sys\t$seq_region\t$g_biotype\t$t\t$ensembl_t_id\t$t_biotype\t$gene_name\n";
      }
    }
  }
  close (INFO);
}

#Load the results from the file generated above
my %transcripts;
my %ribo_transcripts;
my %mt_transcripts;
my %pseudo_transcripts;
open (TRANS, "$transcript_info_file") || die "\n\nCould not open ensembl transcript info file for reading\n\n";
while(<TRANS>){
  chomp($_);
  my @line = split("\t", $_);
  my $g_id = $line[1];
  my $chr = $line[3];
  my $g_biotype = $line[4];
  my $t_id = $line[6];
  my $t_biotype = $line[7];
  my $gene_name = $line[8];
  $transcripts{$t_id}{g_biotype} = $g_biotype;
  $transcripts{$t_id}{t_biotype} = $t_biotype;
  #Mitochondrial genes:  On MT chromosome or of type: Mt_* (gene or transcript biotype)
  if (($chr eq "MT") || ($g_biotype =~ /Mt/) || ($t_biotype =~ /MT/)){
    $mt_transcripts{$t_id}=1;
  }
  #Ribosomal genes: rRNA 
  if (($g_biotype =~ /rRNA/) || ($t_biotype =~ /rRNA/) || $ribo_names{$gene_name} || $ribo_ids{$g_id}){
    $ribo_transcripts{$t_id}=1;
  }
  #Pseudogenes: pseudogene and *_pseudogene
  if (($g_biotype =~ /pseudogene/) || ($t_biotype =~ /pseudogene/)){
    $pseudo_transcripts{$t_id}=1;
  }
}
close(TRANS);
my $trans_count = keys %transcripts;
my $ribo_count = keys %ribo_transcripts;
my $mt_count = keys %mt_transcripts;
my $pseudo_count = keys %pseudo_transcripts;

print BLUE, "\n\nFound the following:", RESET;
print BLUE, "\nTranscripts = $trans_count\nRibosomal transcripts = $ribo_count\nMitochondrial transcripts = $mt_count\nPseudogene transcripts = $pseudo_count\n\n", RESET;


#4. Produce filtered versions of this GTF file.  MT.gtf, rRNA.gtf, pseudo.gtf, rRNA_MT.gtf, rRNA_MT_pseudo.gtf
my $new_genes_gtf_file = "all_sequences.gtf";
my $ribo_file = "rRNA.gtf";
my $mt_file = "MT.gtf";
my $pseudo_file = "pseudogene.gtf";
my $ribo_mt_file = "rRNA_MT.gtf";
my $ribo_mt_pseudo_file = "rRNA_MT_pseudogene.gtf";
open(R, ">$ribo_file") || die "\n\nCould not open file for writing: $ribo_file\n\n";
open(M, ">$mt_file") || die "\n\nCould not open file for writing: $ribo_file\n\n";
open(P, ">$pseudo_file") || die "\n\nCould not open file for writing: $ribo_file\n\n";
open(RM, ">$ribo_mt_file") || die "\n\nCould not open file for writing: $ribo_file\n\n";
open(RMP, ">$ribo_mt_pseudo_file") || die "\n\nCould not open file for writing: $ribo_file\n\n";
open(GTF, "$genes_gtf") || die "\n\nCould not open GTF file: $genes_gtf\n\n";
open(NEWGTF, ">$new_genes_gtf_file") || die "\n\nCould not open file for writing: $new_genes_gtf_file\n\n";
while(<GTF>){
  my $t_id;
  if ($_ =~ /transcript\_id\s+\"(\S+)\"/){
    $t_id = $1;
  }else{
    print RED, "\n\nCould not identify ensembl transcript id in GTF line: $_\n\n", RESET;
    exit();
  }

  #Print out all 'CDS' and 'exon' records to the GTF file
  my @line = split("\t", $_);
  if ($line[2] =~ /exon|CDS/){

    #Skip any entry that does not correpond to a target chromosome
    if ($bowtie_chrs{$line[0]}){
      print NEWGTF "$_";
    }
  }

  if ($ribo_transcripts{$t_id}){
    print R "$_";
    print RM "$_";
    print RMP "$_";
  }
  if ($mt_transcripts{$t_id}){
    print M "$_";
    print RM "$_";
    print RMP "$_";
  }
  if ($pseudo_transcripts{$t_id}){
    print P "$_";
    print RMP "$_";
  }
}
close(R);
close(M);
close(P);
close(RM);
close(RMP);
close(GTF);
close(NEWGTF);

print YELLOW, "\n\nCheck resulting gtf file with: 'gffread -E $new_genes_gtf_file'\n\nRepair any errors or warnings and store a corrected version\n\n", RESET;


exit();


###############################################################################################################
#Create mysql database connection                                                                             #
###############################################################################################################
sub connectDB {
  my %args = @_;
  my $database_name = $args{'-database'};
  my $database_host = $args{'-server'};
  my $user_name = $args{'-user'};
  my $user_pw = $args{'-password'};

  my $dbh = DBI->connect( "dbi:mysql:database=$database_name;host=$database_host", $user_name, $user_pw, { PrintError => 1 } );
  my $connection_attempts = 1;

  #Check for failures to reconnect... if the connection failed, try again a few times
  my $connection_attempt_limit = 3;
  my $sleep_time = 15;
  while(!(defined($dbh)) && ($connection_attempts <= $connection_attempt_limit)){
    $connection_attempts++;
    print YELLOW, "\n\nDBH connect failed ... sleeping for $sleep_time and try attempt: $connection_attempts", RESET;
    sleep $sleep_time;
    $dbh = DBI->connect( "dbi:mysql:database=$database_name;host=$database_host", $user_name, $user_pw, { PrintError => 1 } );
  }

  return $dbh;
}


