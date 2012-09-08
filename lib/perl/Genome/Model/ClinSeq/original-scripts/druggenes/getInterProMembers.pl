#!/usr/bin/env genome-perl
#Written by Obi L. Griffith & Malachi Griffith

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
    $lib_dir = $1."/drug-genes";
  }
}
use lib $lib_dir;
use utility qw(:all);

my $infile = '';
my $outdir = '';
my $result_summary_file = '';


GetOptions ('infile=s'=>\$infile, 'outdir=s'=>\$outdir, 'result_file=s'=>\$result_summary_file);

my $usage=<<INFO;

  Example usage: 
  
  getInterProMembers.pl  --infile=/gscuser/ogriffit/Projects/DruggableGenes/PotentiallyDruggable/Hopkins_and_Groom_2002/DruggableProteinFamiliesCurated.tsv --outdir=/gscuser/ogriffit/Projects/DruggableGenes/PotentiallyDruggable/Hopkins_and_Groom_2002/ --result_file=DruggableProteins.tsv 

  Details:
  --infile            FILE.  File containing list of InterPro protein families for which to look up member proteins
  --outdir                    PATH.  Directory to store results
  --result_file		      FILE.  File to store results from all queries
INFO

unless ($infile && $outdir && $result_summary_file){
  print GREEN, "\n\n$usage\n\n", RESET;
  exit();
}

$outdir = &checkDir('-dir'=>$outdir, '-clear'=>"no");
my $result_path="$outdir"."$result_summary_file";


#Get InterPro ids to query
print BLUE, "\n\nLoading Interpro ids from: $infile\n", RESET;
#my $interpro_data = &loadInterProData($infile);


#Query InterPro Biomart with each InterPro ID
my $queryterm="IPR000022";

my $queryresult = &queryInterProBiomart('-queryterm'=>$queryterm);



#Write all results to file
open(RESULTS, ">$result_path") or die "can't open $result_path for write\n";
binmode(RESULTS, ":utf8");



print "\n\n";
close RESULTS;

exit();


###################################################################################################
#Load InterPro IDs to query                                                         #
###################################################################################################
sub loadInterProData{



}


###################################################################################################
#Query InterPro Biomart with query term                                                         #
###################################################################################################
sub queryInterProBiomart{
  my %args = @_;
  my $queryterm = $args{'-queryterm'};
  print "Attempting query for $queryterm\n";

  my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');
		
	$query->setDataset("uniprot");
	$query->addFilter("proteome_name", ["Homo sapiens"]);
	$query->addFilter("interpro_id", ["IPR002173"]);
	$query->addAttribute("accession");

  my $query_runner = BioMart::QueryRunner->new();

  #Get count of expected results - use to make sure results are complete 
  $query->count(1);
  $query_runner->execute($query);
  my $query_count=$query_runner->getCount();
  print "$query_count results found for query\n";
  $query->count(0);

  ############################## GET RESULTS ##########################
  # to obtain unique rows only
  # $query_runner->uniqueRowsOnly(1);

  $query_runner->execute($query);
#  $query_runner->printHeader();
  $query_runner->printResults();
#  $query_runner->printFooter();
  #####################################################################

 return($query_runner);
}

