#!/usr/bin/env genome-perl
#Written by Obi L. Griffith

#Purpose:
#This script queries the InterPro biomart website using a list InterPro accessions
#Sample perl snippet was obtained from Biomart website
#The result will be a list of UniProtKB protein accessions corresponding to each InterPro accession
#The UniProt ID will then be used to look up Entrez Gene details.

use strict;
use warnings;
use lib '/gscuser/ogriffit/lib/biomart-perl/lib';

use Term::ANSIColor qw(:constants);
use Getopt::Long;
use Data::Dumper;
use LWP;
use above 'Genome';
use Genome::Model::ClinSeq::Util qw(:all);

#Biomart specific stuff:
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;
my $confFile = "/gscuser/ogriffit/lib/biomart-perl/conf/biomart_Interpro_registry.xml"; #From: http://www.biomart.org/biomart/martservice?type=registry";
# Note: change action to 'clean' if you wish to start a fresh configuration  
# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry
my $action='cached';
my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
my $registry = $initializer->getRegistry;

binmode(STDOUT, ":utf8");
use Cwd 'abs_path';
my $lib_dir;
BEGIN{
  if (abs_path($0) =~ /(.*\/).*\/.*\.pl/){
    $lib_dir = $1."/druggenes";
  }
}
use lib $lib_dir;
use utility qw(:all);

my $infile = '';
my $outdir = '';
my $genes_file = '';
my $terms_file = '';

GetOptions ('infile=s'=>\$infile, 'outdir=s'=>\$outdir, 'genes_file=s'=>\$genes_file, 'terms_file=s'=>\$terms_file);

my $usage=<<INFO;

  Example usage: 
  
  getInterProMembers.pl  --infile=/gscuser/ogriffit/Projects/DruggableGenes/PotentiallyDruggable/Hopkins_and_Groom_2002/DruggableProteinFamiliesCurated.tsv --outdir=/gscuser/ogriffit/Projects/DruggableGenes/PotentiallyDruggable/Hopkins_and_Groom_2002/ --genes_file=HopkinsGroomGenes.tsv --terms_file=HopkinsGroomTerms.tsv 

  Details:
  --infile            FILE.  File containing list of InterPro protein families for which to look up member proteins
  --outdir                    PATH.  Directory to store results
  --genes_file		      FILE.  File to store resulting gene (uniprot) list from all interpro queries
  --terms_file          FILE.  File to store resulting terms (InterPro Names) from all interpro queries
INFO

unless ($infile && $outdir && $genes_file && $terms_file){
  print GREEN, "\n\n$usage\n\n", RESET;
  exit();
}

$outdir = &checkDir('-dir'=>$outdir, '-clear'=>"no");
my $genes_path="$outdir"."$genes_file";
my $terms_path="$outdir"."$terms_file";
my $tempfile1="biomart_temp.txt";
my $tempfile2="biomart_temp2.txt";

#Get InterPro ids to query
print BLUE, "\n\nLoading Interpro ids from: $infile\n\n", RESET;
my $interpro_data = &loadInterProData('-infile'=>$infile);
my @queryterms=@{$interpro_data};

#Query InterPro/Uniprot Biomart with each InterPro ID
#my @queryterms=("IPR000022","IPR002117"); #Used for testing purposes
my %InterPro2UniProt;
my %InterProDetails;
my %UniProtDetails;
foreach my $queryterm(@queryterms){
  $InterPro2UniProt{$queryterm}=&queryUniprotBiomart('-queryterm'=>$queryterm);
  my @Details=@{&queryInterProBiomart('-queryterm'=>$queryterm)};
  $InterProDetails{$queryterm}{Accession}=$Details[0];
  $InterProDetails{$queryterm}{Name}=$Details[1];
  $InterProDetails{$queryterm}{ShortName}=$Details[2];
  $InterProDetails{$queryterm}{Type}=$Details[3];
  $InterProDetails{$queryterm}{Count}=@{$InterPro2UniProt{$queryterm}};
}

#Get mapping of Uniprot Accessions to Entrez IDs, etc
#These can be obtained from here:
#ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/ 
#ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz
print "\nAttempting download of UniProt mapping file\n";
my $mapping_file_url="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz";
my $mapping_file_name="HUMAN_9606_idmapping_selected.tab.gz";
my $mapping_file_path = &download_file($mapping_file_url,$mapping_file_name);

print "\nParsing Uniprot mapping file\n";
my %UniProtMapping;
open (MAPPING, $mapping_file_path) or die "can't open $mapping_file_path\n";
while (<MAPPING>){
  my @data=split("\t",$_);
  my $uniprot_acc=$data[0];
  my $uniprot_id=$data[1];
  my $entrez_id=$data[2];
  unless ($entrez_id){$entrez_id="N/A";}
  my $ensembl_id=$data[19];
  unless ($ensembl_id){$ensembl_id="N/A";}
  $UniProtMapping{$uniprot_acc}{uniprot_acc}=$uniprot_acc;
  $UniProtMapping{$uniprot_acc}{entrez_id}=$entrez_id;
  $UniProtMapping{$uniprot_acc}{ensembl_id}=$ensembl_id;
}
close MAPPING;

#Write all results to file
open(GENES, ">$genes_path") or die "can't open $genes_path for write\n";
binmode(GENES, ":utf8");
print GENES "Uniprot_acc\tUniprot_id\tUniprot_protein_name\tUniprot_gene_name\tUniprot_evidence\tUniprot_status\tEntrez_id\tEnsembl_Id\tInterpro_acc\n";
foreach my $interpro (sort keys %InterPro2UniProt){
  my @uniprots=@{$InterPro2UniProt{$interpro}};
  foreach my $uniprot (@uniprots){
    print GENES "$uniprot\t$UniProtDetails{$uniprot}{Uniprot_id}\t$UniProtDetails{$uniprot}{Uniprot_protein_name}\t$UniProtDetails{$uniprot}{Uniprot_gene_name}\t$UniProtDetails{$uniprot}{Uniprot_evidence}\t$UniProtDetails{$uniprot}{Uniprot_status}\t$UniProtMapping{$uniprot}{entrez_id}\t$UniProtMapping{$uniprot}{ensembl_id}\t$interpro\n";
  }
}
close GENES;

open(TERMS, ">$terms_path") or die "can't open $terms_path for write\n";
binmode(TERMS, ":utf8");
print TERMS "Accession\tName\tShortName\tType\tCount\n";
foreach my $interpro (sort keys %InterProDetails){
  print TERMS "$InterProDetails{$interpro}{Accession}\t$InterProDetails{$interpro}{Name}\t$InterProDetails{$interpro}{ShortName}\t$InterProDetails{$interpro}{Type}\t$InterProDetails{$interpro}{Count}\n";
}
close TERMS;

#Clean up temp files
if (unlink($tempfile1) == 0) {print "$tempfile1 deleted successfully.";} else {print "$tempfile1 was not deleted.\n";}
if (unlink($tempfile2) == 0) {print "$tempfile2 deleted successfully.";} else {print "$tempfile2 was not deleted.\n";}
exit();

###################################################################################################
#Load InterPro IDs to query                                                                       #
###################################################################################################
sub loadInterProData{
  my %args = @_;
  my $infile = $args{'-infile'};
  my @accessions;
  open(INFILE, $infile) or die "can't open $infile\n";
  my $header=<INFILE>;
  while(<INFILE>){
    my @data=split("\t", $_);
    if ($data[1]=~/^IPR\d+/){
      push(@accessions,$data[1]);
    }else{
      print "IPR accession not recognized\n"; exit;
    }
  }
return(\@accessions);
}

###################################################################################################
#Query Uniprot Biomart with query term                                                            #
###################################################################################################
sub queryUniprotBiomart{
  my %args = @_;
  my $queryterm = $args{'-queryterm'};
  print "\nAttempting UniProt list query for $queryterm\n";

  my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');
		
	$query->setDataset("uniprot");
	$query->addFilter("interpro_id", [$queryterm]);
  $query->addFilter("proteome_name", ["Homo sapiens"]); #Note, this requires inclusion in "The complete human proteome", see: http://www.uniprot.org/faq/48
	$query->addFilter("entry_type", ["Swiss-Prot"]); #Only allow proteins with Swiss-prot (Reviewed) status and disqualify TrEMBL (Unreviewed) entries, see http://www.uniprot.org/faq/7
  #$query->addFilter("protein_evidence", ["1: Evidence at protein level"]); #Require evidence at protein level, see http://www.uniprot.org/docs/pe_criteria
  #See below for code to filter to multiple evidence levels after query result returned
  $query->addAttribute("accession");
  $query->addAttribute("name");
	$query->addAttribute("protein_name");
	$query->addAttribute("gene_name");
	$query->addAttribute("protein_evidence");
	$query->addAttribute("entry_type");
  my $query_runner = BioMart::QueryRunner->new();

  #Get count of expected results - use to make sure results are complete 
  $query->count(1);
  my $query_count;
  unless($query_count){
    $query_runner->execute($query);
    $query_count=$query_runner->getCount();
  }
  print "$query_count results expected for query\n";
  $query->count(0); #turn off counting so that full results can be obtained below 

  #Get results and store in temporary file
  open (BIOMART_OUT, ">$tempfile1") or die "Can't open $tempfile1 file for write\n";
  $query_runner->uniqueRowsOnly(1); #to obtain unique rows only
  $query_runner->execute($query);
  #$query_runner->printHeader(\*BIOMART_OUT);
  $query_runner->printResults(\*BIOMART_OUT);
  #$query_runner->printFooter(\*BIOMART_OUT);
  close BIOMART_OUT;

  #Read in results and check expected results against count above
  open (BIOMART_IN, "$tempfile1") or die "Can't open $tempfile1\n";
  my @results=<BIOMART_IN>;
  close BIOMART_IN;
  my $result_count=@results;
  print "$result_count results returned for query\n";
  unless ($result_count==$query_count){die "missing UniProt results for $queryterm\n";}
  chomp (@results);

  #Parse results
  my @accessions;
  foreach my $result (@results){
    my @data=split("\t", $result);
    my $Uniprot_acc=$data[0];
    my $Uniprot_id=$data[1];
    my $Uniprot_protein_name=$data[2]; unless($Uniprot_protein_name){$Uniprot_protein_name="N/A";}
    my $Uniprot_gene_name=$data[3]; unless($Uniprot_gene_name){$Uniprot_gene_name="N/A";}
    my $Uniprot_evidence=$data[4]; unless($Uniprot_evidence){$Uniprot_evidence="N/A";}
    my $Uniprot_status=$data[5]; unless($Uniprot_status){$Uniprot_status="N/A";}
    unless ($Uniprot_evidence=~/protein level/ || $Uniprot_evidence=~/transcript level/){next;} #Can't filter down to multiple evidence levels above, do here if desired
    push(@accessions,$Uniprot_acc);
    $UniProtDetails{$Uniprot_acc}{Uniprot_id}=$Uniprot_id;
    $UniProtDetails{$Uniprot_acc}{Uniprot_protein_name}=$Uniprot_protein_name;
    $UniProtDetails{$Uniprot_acc}{Uniprot_gene_name}=$Uniprot_gene_name;
    $UniProtDetails{$Uniprot_acc}{Uniprot_evidence}=$Uniprot_evidence;
    $UniProtDetails{$Uniprot_acc}{Uniprot_status}=$Uniprot_status;
  }
 return(\@accessions);
}


###################################################################################################
#Query InterPro Biomart with query term                                                           #
###################################################################################################
sub queryInterProBiomart{
  my %args = @_;
  my $queryterm = $args{'-queryterm'};
  print "Attempting InterPro details query for $queryterm\n";
  my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');
  $query->setDataset("entry");
  $query->addFilter("entry_id", [$queryterm]);
  $query->addAttribute("entry_id");
  $query->addAttribute("entry_name");
  $query->addAttribute("entry_short_name");
  $query->addAttribute("entry_type");
  my $query_runner = BioMart::QueryRunner->new();

  #Get results and store in temporary file
  open (BIOMART_OUT2, ">$tempfile2") or die "Can't open $tempfile2 file for write\n";
  $query_runner->uniqueRowsOnly(1); #to obtain unique rows only
  $query_runner->execute($query);
  #$query_runner->printHeader(\*BIOMART_OUT2);
  $query_runner->printResults(\*BIOMART_OUT2);
  #$query_runner->printFooter(\*BIOMART_OUT2);
  close BIOMART_OUT2;

  #Read in results
  open (BIOMART_IN2, "$tempfile2") or die "Can't open $tempfile2\n";
  my$line=<BIOMART_IN2>;
  chomp $line;
  my @details=split("\t",$line);
  close BIOMART_IN2;
  my $details_count=@details;
  print "Found $details_count details\n";
  unless ($details_count==4){die "missing details in $line\n";}
  return(\@details);
}

sub download_file {
    my $url = shift;
    my $targetfilename = shift;
    my $targetfilepath="$outdir"."$targetfilename";
    my $wget_cmd = "wget $url -O $targetfilepath";
    my $retval = Genome::Sys->shellcmd(cmd=>$wget_cmd);
    unless ($retval == 1){
      self->error_message('Failed to wget the specified URL');
      return;
    }
    #unzip if necessary
    if ($targetfilepath=~/\.gz$/){
      my $gunzip_cmd = "gunzip -f $targetfilepath";
      my $retval2 = Genome::Sys->shellcmd(cmd=>$gunzip_cmd);
      unless ($retval2 == 1){
        self->error_message('Failed to gunzip the specified file');
        return;
      }
      $targetfilepath=~s/\.gz$//;
    }
    print "Downloaded $targetfilepath\n";
    return $targetfilepath;
}

