#!/usr/bin/env genome-perl -w
#Written by Malachi Griffith

#Parse an XML database file from DrugBank
#Get drug, interaction and gene info for each drug-gene interaction in the database
#Get gene names and uniprot IDs from Entrez gene
#Add official 'EntrezGene' name to each gene record

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use XML::Simple;

my $script_dir;
use Cwd 'abs_path';
BEGIN{
  if (abs_path($0) =~ /(.*\/).*\.pl/){
    $script_dir = $1;
  }
}
use lib $script_dir;
use utility qw(:all);

binmode(STDOUT, ":utf8");

my $infile = '';
my $entrez_dir = '';
my $verbose = '';

GetOptions ('infile=s'=>\$infile, 'verbose=s'=>\$verbose, 'entrez_dir=s'=>\$entrez_dir);

my $usage=<<INFO;

  Example usage: 
  
  parseDrugBankXml.pl  --infile=drugbank.xml  --entrez_dir=./EntrezGene/

  
  Details:
  --infile                    PATH.  XML data file dowloaded from www.drugbank.ca
  --entrez_dir                PATH.  Directory containing gene info for gene name mapping ('gene2accession' and 'gene_info')
  --verbose=1                 Print more output while running

INFO

unless ($infile && $entrez_dir){
  print GREEN, "\n\n$usage\n\n", RESET;
  exit();
}

#Parse Entrez flatfiles
#ftp://ftp.ncbi.nih.gov/gene/DATA/gene2accession.gz
#ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz
my $entrez_data = &loadEntrezData('-entrez_dir'=>$entrez_dir);

#Instantiate an XML simple object
my $xs1 = XML::Simple->new();

#Create data object for the entire XML file, forcing array creation even if there is only one element
#my $xml = $xs1->XMLin($infile, ForceArray => 1);

#Create data object for the entire XML file, allowing mixed array/hash structure for the same field depending on whether there is only one element or more than one.
#my $xml = $xs1->XMLin($infile);

#Create data object for the entire XML file, allowing mixed array/hash structure for the same field depending of whether is only one element or more than one.
#Also actually specify the primary IDs you would like the resultng data structures to be keyed on ('drugbank-id' for drug records, 'id' for gene partners)
my $xml = $xs1->XMLin($infile, KeyAttr => ['nct_id', 'id'] );


#Get the 'drug' tree
my $info_source = "DrugBank";
my $drugs = $xml->{'drug'};

#Get the partners tree
#Unfortunately, the primary key of this tree is the gene name but drug intereactions seem to be linked by 'partner id'
#This means that we can not directly look up partner records for each drug.  We must traverse the partner hash for each drug (slow) or create a new data structure that is keyed in partner ID
my $partners = $xml->{'partners'};

#Build a new simpler partners object for convenience - since it is keyed on partner ID we could access directly as well...
my $partners_lite = &organizePartners('-partner_ref'=>$partners);

#Count distinct drug-gene interactions
my $ic = 0;

#Create three output files: Drugs, Targets, Interactions
my $drugs_outfile = "DrugBank_WashU_DRUGS.tsv";
my $targets_outfile = "DrugBank_WashU_TARGETS.tsv";
my $interactions_outfile = "DrugBank_WashU_INTERACTIONS.tsv";
open(DRUGS, ">$drugs_outfile") || die "\n\nCould not open outfile: $drugs_outfile\n\n";
binmode(DRUGS, ":utf8");
open(TARGETS, ">$targets_outfile") || die "\n\nCould not open outfile: $targets_outfile\n\n";
binmode(TARGETS, ":utf8");
open(INTERACTIONS, ">$interactions_outfile") || die "\n\nCould not open outfile: $interactions_outfile\n\n";
binmode(INTERACTIONS, ":utf8");

#Print out a header line foreach output file
my $interactions_header = "interaction_count\tdrug_id\tdrug_name\tdrug_synonyms\tdrug_brands\tdrug_type\tdrug_groups\tdrug_categories\tpartner_id\tknown_action?\ttarget_actions\tgene_symbol\tgene_name\tuniprot_id\tentrez_gene_symbols\tentrez_gene_synonyms\tlink";
print INTERACTIONS "$interactions_header\n";

my $drugs_header = "drug_id\tdrug_name\tdrug_synonyms\tdrug_brands\tdrug_type\tdrug_groups\tdrug_categories\ttarget_count\tlink";
print DRUGS "$drugs_header\n";

my $targets_header = "partner_id\tgene_symbol\tgene_name\tuniprot_id\tentrez_gene_symbols\tentrez_gene_synonyms";
print TARGETS "$targets_header\n";

foreach my $drug_id (sort {$a cmp $b} keys %{$drugs}){
  my $drug_name = $drugs->{$drug_id}->{'name'};
  my $drug_type = $drugs->{$drug_id}->{'type'};
  if ($verbose){
    print BLUE, "\n\n$drug_id\t$drug_type\t$drug_name", RESET;
  }

  #Get the drug groups
  my $group_list = $drugs->{$drug_id}->{'groups'};
  my @groups = @{&parseXmlTree('-ref'=>$group_list, '-value_name'=>'group')};
  if ($verbose){
    print BLUE, "\nGroups: @groups", RESET;
  }
  my $drug_groups_string = join(",", @groups);

  #Get the drug synonyms
  my $synonym_list = $drugs->{$drug_id}->{'synonyms'};
  my @synonyms = @{&parseXmlTree('-ref'=>$synonym_list, '-value_name'=>'synonym')};
  if ($verbose){
    print BLUE, "\nSynonyms: @synonyms", RESET;
  }
  my $drug_synonyms_string = join(",", @synonyms);

  #Get the drug brands
  my $brand_list = $drugs->{$drug_id}->{'brands'};
  my @brands = @{&parseXmlTree('-ref'=>$brand_list, '-value_name'=>'brand')};
  if ($verbose){
    print BLUE, "\nBrands: @brands", RESET;
  }
  my $drug_brands_string = join(",", @brands);

  #Get the drug categories
  my $category_list = $drugs->{$drug_id}->{'categories'};
  my @categories = @{&parseXmlTree('-ref'=>$category_list, '-value_name'=>'category')};
  if ($verbose){
    print BLUE, "\nCategories: @categories", RESET;
  }
  my $drug_categories_string = join(",", @categories);

  #Get the targets.  i.e. the gene parter ids for this drug
  my $targets = $drugs->{$drug_id}->{'targets'};
  my $target_list = $targets->{'target'};
  my @target_partners = @{&parseXmlTree('-ref'=>$target_list, '-value_name'=>'partner')};
  my $t_count = scalar(@target_partners);

  #Get the known action status (yes|unknown|no?) for this drug
  my @target_known_actions = @{&parseXmlTree('-ref'=>$target_list, '-value_name'=>'known-action')};

  #Get the actions associated with each drug->gene target interaction.  Note that there can be more than one action for each interaction - only defined if the action status was "yes"
  my @action_lists =  @{&parseXmlTree('-ref'=>$target_list, '-value_name'=>'actions')};
  my @target_actions_joined;
  foreach my $action_list (@action_lists){
    my @actions =  @{&parseXmlTree('-ref'=>$action_list, '-value_name'=>'action')};
    my $action_string = join(",", @actions);
    push(@target_actions_joined, $action_string);
  }
  if ($verbose){
    print "\n\tTarget Partners: @target_partners";
    print "\n\tTarget Known Actions: @target_known_actions";
    print "\n\tTarget Actions: @target_actions_joined";
  }

  #Strip <tabs> from string variables
  $drug_name =~ s/\t/ /g;
  $drug_synonyms_string =~ s/\t/ /g;
  $drug_brands_string =~ s/\t/ /g;
  $drug_type =~ s/\t/ /g;
  $drug_groups_string =~ s/\t/ /g;
  $drug_categories_string =~ s/\t/ /g;

  #Generate a URL link back to DrugBank using the DrugBank id
  my $link = "http://www.drugbank.ca/drugs/$drug_id";

  #Print out each drug
  my $drugs_line = "$drug_id\t$drug_name\t$drug_synonyms_string\t$drug_brands_string\t$drug_type\t$drug_groups_string\t$drug_categories_string\t$t_count\t$link";
  print DRUGS "$drugs_line\n";

  #Print out each drug-gene interaction...
  for (my $i = 0; $i < $t_count; $i++){
    my $target_pid = $target_partners[$i];
    my $target_known_action = $target_known_actions[$i];
    my $target_actions = $target_actions_joined[$i];

    unless ($partners_lite->{$target_pid}){
      print RED, "\n\nTarget PID: $target_pid is not defined in the partners hash\n\n", RESET;
      exit();
    }
    my $gene_symbol = $partners_lite->{$target_pid}->{gene_symbol};
    my $gene_name = $partners_lite->{$target_pid}->{gene_name};
    my $uniprotkb = $partners_lite->{$target_pid}->{uniprotkb};
    my $entrez_gene_symbols = $partners_lite->{$target_pid}->{entrez_gene_symbols};
    my $entrez_gene_synonyms = $partners_lite->{$target_pid}->{entrez_gene_synonyms};

    #Strip <tabs> from string variables
    $target_known_action =~ s/\t/ /g;
    $target_actions =~ s/\t/ /g;

    $ic++;
    my $interactions_line = "$ic\t$drug_id\t$drug_name\t$drug_synonyms_string\t$drug_brands_string\t$drug_type\t$drug_groups_string\t$drug_categories_string\t$target_pid\t$target_known_action\t$target_actions\t$gene_symbol\t$gene_name\t$uniprotkb\t$entrez_gene_symbols\t$entrez_gene_synonyms\t$link";
    print INTERACTIONS "$interactions_line\n";

  }
}

foreach my $pid (sort {$a <=> $b} keys %{$partners_lite}){
  my $gene_symbol = $partners_lite->{$pid}->{gene_symbol};
  my $gene_name = $partners_lite->{$pid}->{gene_name};
  my $uniprot_id = $partners_lite->{$pid}->{uniprotkb};
  my $entrez_gene_symbols = $partners_lite->{$pid}->{entrez_gene_symbols};
  my $entrez_gene_synonyms = $partners_lite->{$pid}->{entrez_gene_synonyms};

  my $targets_line = "$pid\t$gene_symbol\t$gene_name\t$uniprot_id\t$entrez_gene_symbols\t$entrez_gene_synonyms";
  print TARGETS "$targets_line\n";
}

close(DRUGS);
close(TARGETS);
close(INTERACTIONS);

print "\n\n";

exit();


#######################################################################################################################################################################
#Load Entrez Data from flatfiles                                                                                                                                      #
#######################################################################################################################################################################
sub loadEntrezData{
  my %args = @_;
  my $entrez_dir = $args{'-entrez_dir'};
  my %edata;

  #Check input dir
  unless (-e $entrez_dir && -d $entrez_dir){
    print RED, "\n\nEntrez dir not valid: $entrez_dir\n\n", RESET;
    exit();
  }
  unless ($entrez_dir =~ /\/$/){
    $entrez_dir .= "/";
  }

  my %entrez_map;      #Entrez_id -> symbol, synonyms
  my %symbols_map;     #Symbols   -> entrez_id(s)
  my %synonyms_map;    #Synonyms  -> entrez_id(s)
  my %p_acc_map;       #Protein accessions -> entrez_id

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
    if ($synonyms eq "-"){
      $synonyms = "na";
    }
    my @synonyms_array = split("\\|", $synonyms);
    my %synonyms_hash;   
    foreach my $syn (@synonyms_array){
      $synonyms_hash{$syn} = 1;
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
    #If the prot is not defined, skip
    if ($prot_id eq "-"){next();}
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
  close (ACC);

  #print Dumper %entrez_map;
  #print Dumper %symbols_map;
  #print Dumper %synonyms_map;
  #print Dumper %p_acc_map;

  $edata{'entrez_ids'} = \%entrez_map;
  $edata{'symbols'} = \%symbols_map;
  $edata{'synonyms'} = \%synonyms_map;
  $edata{'protein_accessions'} = \%p_acc_map;

  return(\%edata);
}


#######################################################################################################################################################################
#Attempt to map the gene_symbol / uniprotkb_id to Entrez Gene Symbol
#######################################################################################################################################################################
sub mapGeneToEntrez{
  my %args = @_;
  my $entrez_data = $args{'-entrez_data'};
  my $gene_symbol = $args{'-gene_symbol'};
  my $uniprot_id = $args{'-uniprot_id'};
  my $entrez_gene_symbol = "na";
  my $entrez_gene_synonyms = "na";

  my $entrez_map = $entrez_data->{'entrez_ids'};
  my $symbols_map = $entrez_data->{'symbols'};
  my $prot_acc_map = $entrez_data->{'protein_accessions'};

  if ($gene_symbol eq "na" && $uniprot_id eq "na"){
    #If gene_symbol and uniprot_id are 'na' ... nothing can be done.  Return 'na' for both symbol and synonyms

  }else{
    #If gene_symbol is not 'na', see if it matches an entrez gene symbol perfectly, return it and associated synonyms
    my $perfect_match = 0;
    my $uniprot_match = 0;
    unless ($gene_symbol eq "na"){
      if ($symbols_map->{$gene_symbol}){
        $perfect_match = 1;
        $entrez_gene_symbol = $gene_symbol;
        my $ids = $symbols_map->{$gene_symbol}->{entrez_ids};
        my %symbol_list;
        my %synonym_list;
        foreach my $entrez_id (keys %{$ids}){
          my $symbol = $entrez_map->{$entrez_id}->{symbol};
          $symbol_list{$symbol}=1;
          my @synonyms_array = @{$entrez_map->{$entrez_id}->{synonyms_array}};
          foreach my $syn (@synonyms_array){
            $synonym_list{$syn}=1;
          }
        }
        my @symbols1;
        my @synonyms1;
        foreach my $symbol (keys %symbol_list){
          push(@symbols1, $symbol);
        }
        foreach my $synonym (keys %synonym_list){
          push(@synonyms1, $synonym);
        }
        $entrez_gene_symbol = join(",", @symbols1);
        $entrez_gene_synonyms = join(",", @synonyms1);

      }else{
        #If gene_symbol does not match or is 'na', attempt to get a matching symbol using the uniprot ID
        unless ($uniprot_id eq "na"){
          if ($prot_acc_map->{$uniprot_id}){
            $uniprot_match = 1;
            my $ids = $prot_acc_map->{$uniprot_id}->{entrez_ids};
            my %symbol_list;
            my %synonym_list;
            foreach my $entrez_id (keys %{$ids}){
              my $symbol = $entrez_map->{$entrez_id}->{symbol};
              $symbol_list{$symbol}=1;
              my @synonyms_array = @{$entrez_map->{$entrez_id}->{synonyms_array}};
              foreach my $syn (@synonyms_array){
                $synonym_list{$syn}=1;
              }
            }
            my @symbols1;
            my @synonyms1;
            foreach my $symbol (keys %symbol_list){
              push(@symbols1, $symbol);
            }
            foreach my $synonym (keys %synonym_list){
              push(@synonyms1, $synonym);
            }
            $entrez_gene_symbol = join(",", @symbols1);
            $entrez_gene_synonyms = join(",", @synonyms1);

          }

        }
      }
    }

    #If the uniprot mapping fails, and the gene_symbol is not 'na' attempt a match to a synonym...
    #Too risky?? - could lead to false mappings?


  }

  my %result;
  $result{'entrez_gene_symbols'} = $entrez_gene_symbol;
  $result{'entrez_gene_synonyms'} = $entrez_gene_synonyms;
  
  #print YELLOW, "\n$gene_symbol\t$uniprot_id\t$entrez_gene_symbol\t$entrez_gene_synonyms", RESET;

  return(\%result);
}


#######################################################################################################################################################################
#Build a new partners object keyed on partner ID                                                                                                                      #
#######################################################################################################################################################################
sub organizePartners{
  my %args = @_;
  my $p_ref = $args{'-partner_ref'};
  my %p_lite;

  if ($verbose){
    print BLUE, "\n\nBuilding new 'partner lite' object\n\n", RESET;
  }
  my $p = $p_ref->{'partner'};

  foreach my $pid (keys %{$p}){

    #Get the gene symbol
    my $gene_symbol;
    my $gene_symbol_r = $p->{$pid}->{'gene-name'};
    if (ref($gene_symbol_r) eq "HASH"){
      $gene_symbol = %{$gene_symbol_r};
    }else{
      $gene_symbol = $gene_symbol_r;
    }
    unless ($gene_symbol){
      $gene_symbol = "na";
    }

    #Get the gene name
    my $gene_name;
    my $gene_name_r = $p->{$pid}->{'name'};
    if (ref($gene_name_r) eq "HASH"){
      $gene_name = %{$gene_name_r};
    }else{
      $gene_name = $gene_name_r;
    }
    unless ($gene_name){
      $gene_name = "na";
    }

    if ($verbose){
      print YELLOW, "\n\t$pid\t$gene_symbol\t$gene_name", RESET;
    }

    #Get the external identifiers - for now, just store the UniProt ID
    my %external_ids;
    my $ext_id_ref = $p->{$pid}->{'external-identifiers'};
    my $ext_id_list = $ext_id_ref->{'external-identifier'};
    if (ref($ext_id_list) eq "ARRAY"){
      foreach my $eid (@{$ext_id_list}){
        my $resource = $eid->{'resource'};
        my $identifier = $eid->{'identifier'};
        $external_ids{$resource} = $identifier;
        if ($verbose){
          print YELLOW, "\n\t\t$resource\t$identifier", RESET;
        }
      }
    }else{
        if ($ext_id_list->{'resource'}){
          my $resource = $ext_id_list->{'resource'};
          my $identifier = $ext_id_list->{'identifier'};
          if ($verbose){
            print YELLOW, "\n\t\t$resource\t$identifier", RESET;
          }
          $external_ids{$resource} = $identifier;
        }else{
          if ($verbose){
            print RED, "\n\t\tFound no external IDs at all", RESET;
          }
        }
    }
    my $uniprotkb = "na";
    if ($external_ids{'UniProtKB'}){
      $uniprotkb = $external_ids{'UniProtKB'};
    }

    #Store the info gathered thus far in the simplifed data structure keyed on partner ID
    #Remove <tabs> just in case
    $gene_symbol =~ s/\t/ /g;
    $gene_name =~ s/\t/ /g;
    $uniprotkb =~ s/\t/ /g;
    $p_lite{$pid}{gene_symbol} = $gene_symbol;
    $p_lite{$pid}{gene_name} = $gene_name;
    $p_lite{$pid}{uniprotkb} = $uniprotkb;

    #Attempt to map the gene_symbol / uniprotkb_id to Entrez Gene Symbol
    #Get the entrez symbol AND synonyms
    my $entrez_gene_info = &mapGeneToEntrez('-entrez_data'=>$entrez_data, '-gene_symbol'=>$gene_symbol, '-uniprot_id'=>$uniprotkb);
    my $entrez_gene_symbols = $entrez_gene_info->{'entrez_gene_symbols'};
    my $entrez_gene_synonyms = $entrez_gene_info->{'entrez_gene_synonyms'};
    $p_lite{$pid}{entrez_gene_symbols} = $entrez_gene_symbols;
    $p_lite{$pid}{entrez_gene_synonyms} = $entrez_gene_synonyms;




  }

  #print Dumper %p_lite;
  #foreach my $pid (sort {$a <=> $b} keys %p_lite){
  #  print CYAN, "\n$pid\t$p_lite{$pid}{gene_name}\t$p_lite{$pid}{drug_name}\t$p_lite{$pid}{uniprotkb}", RESET;
  #}


  return(\%p_lite);
}




