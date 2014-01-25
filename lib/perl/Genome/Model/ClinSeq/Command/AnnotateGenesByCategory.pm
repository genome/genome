package Genome::Model::ClinSeq::Command::AnnotateGenesByCategory;

use strict;
use warnings;
use Genome;
use Genome::Model::ClinSeq::Util qw(:all);

class Genome::Model::ClinSeq::Command::AnnotateGenesByCategory {
    is => 'Command::V2',
    has_input => [
        infile => {
            is => 'FilesystemPath',
            doc => 'tsv formatted file to added gene category annotations to',
        },
        cancer_annotation_db => {
            is => 'Genome::Db::Tgi::CancerAnnotation',
            doc => 'db of cancer annotation (see \'genome db list\' for latest version of desired database)', 
        },
        gene_name_columns => {
            is => 'Text',
            is_many => 1,
            doc => 'name(s) of column(s) containing gene names/symbols'
        },
        gene_groups => {
            is => 'Text',
            doc => 'name of list of gene groups to use',
            default => "Default",
        },
    ],
    has_optional_output => [
        category_outfile => {
            is => 'FilesystemPath',            
            doc => 'result file consisting of input file appended with gene category values',
        }
    ],
    doc => 'take a tsv file as input, indentify gene names, append matches between each of these genes and various gene categories (e.g. kinase, etc.)',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq annotate-genes-by-category --infile=/gscmnt/ams1108/info/model_data/2888708572/build134369422/AML103/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.readcounts.tsv --cancer_annotation_db='tgi/cancer-annotation/human/build37-20130711.1'  --gene-name-column='mapped_gene_name'

EOS
}

sub help_detail {
    return <<EOS
create an updated version of an input tsv file with gene names appending columns to indicate with gene categories each gene belongs to

EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);

  unless (-e $self->infile) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['infile'],
	                                          desc => "infile: " . $self->infile . " was not found",
                                          );
  }

  return @errors;
}


sub execute {
  my $self = shift;
  my $infile = $self->infile;
  my $cancer_annotation_db = $self->cancer_annotation_db;
  my $gene_symbol_lists_dir = $cancer_annotation_db->data_set_path("GeneSymbolLists");
  my @gene_name_columns = $self->gene_name_columns;


  #Get Entrez and Ensembl data for gene name mappings
  my $entrez_ensembl_data = $self->loadEntrezEnsemblData(-cancer_db => $cancer_annotation_db);

  my $symbol_list_names = $self->importSymbolListNames('-gene_symbol_lists_dir'=>$gene_symbol_lists_dir, '-verbose'=>0);
  my $master_list = $symbol_list_names->{master_list};
  my @symbol_list_names = sort {$master_list->{$a}->{order} <=> $master_list->{$b}->{order}} keys %{$master_list};
  my $gene_symbol_lists = $self->importGeneSymbolLists('-gene_symbol_lists_dir'=>$gene_symbol_lists_dir, '-symbol_list_names'=>\@symbol_list_names, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);

  my $sublists = $symbol_list_names->{sublists};
  unless (defined($sublists->{$self->gene_groups})) {
    $self->error_message("Could not find a gene sublist with the name ".$self->gene_groups);
    return;
  }

  my $target_groups = $sublists->{$self->gene_groups}->{groups};
  
  open (INFILE, "$infile") || die "\n\nCould not open input infile: $infile\n\n";
  my %data;
  my %cols;
  my $header = 1;
  my $header_line = '';
  my $l = 0;
  while(<INFILE>){
    $l++;    
    chomp($_);    
    my $record = $_;
    my @line = split("\t", $_);
    if ($header == 1){
      my $c = 0;
      $header_line = $_;
      foreach my $colname (@line){
        $cols{$colname}{position} = $c;
        $c++;
      }
      $header = 0;
      for my $gene_name_column (@gene_name_columns) {
          unless ($cols{$gene_name_column}){
              die $self->error_message("File has no $gene_name_column column: $infile");
          }
      }
      next();
    }
    $data{$l}{record} = $record;
    for my $gene_name_column (@gene_name_columns) {
        $data{$l}{$gene_name_column}{gene_name} = $line[$cols{$gene_name_column}{position}];
        $data{$l}{$gene_name_column}{sum} = 0;
    }
  }
  close(INFILE);

  #Figure out the gene matches to the gene symbol lists
  #Test each gene name in this column against those in the list and add a column with the match status (i.e. is it a kinase, cancer gene, etc.)
  my %member_list;
  foreach my $l (keys %data){
      for my $gene_name_column (@gene_name_columns) {
          my $gene_name = $data{$l}{$gene_name_column}{gene_name};
          for my $group_name (sort {$target_groups->{$a}->{order} <=> $target_groups->{$b}->{order}} keys %{$target_groups}) {
              my @group_members = @{$symbol_list_names->{master_group_list}->{$group_name}->{members}};
              for my $member (@group_members) {
                  $member_list{$member} = 1;
                  my $gene_symbols = $gene_symbol_lists->{$member}->{symbols};
                  if ($gene_symbols->{$gene_name}){
                      unless ($data{$l}{$gene_name_column}{$member}) {
                          $data{$l}{$gene_name_column}{$member} = 1;
                          $data{$l}{$gene_name_column}{sum}++;
                      }
                  }else{
                      $data{$l}{$gene_name_column}{$member} = 0;
                  }
              }
          }
      }
  }

  #Create a new file that has the same path and name as the input file except with '.catanno' injected into the end of the name
  my $base;
  my $extension;
  if ($infile =~ /(.*)(\.\w+)$/){
    $base = $1;
    $extension = $2;
  }else{
    die $self->error_message("Could not resolve base name path from infile: $infile");
  }
  my $category_outfile = $base . ".catanno" . $extension;

  #Print out a new file contain the extra columns
  open (OUT, ">$category_outfile") || die "\n\nCould not open output datafile: $category_outfile\n\n";
  my @gene_symbol_list_names = sort {$gene_symbol_lists->{$a}->{order} <=> $gene_symbol_lists->{$b}->{order}} keys %member_list;
  my @headers;
  if (@gene_name_columns > 1) {
      for my $gene_name_column (@gene_name_columns) {
          push @headers, map {$gene_name_column."-".$_} @gene_symbol_list_names;
      }
  }
  else {
      @headers = @gene_symbol_list_names;
  }
  my $gene_symbol_list_name_string = join("\t", @headers);
  print OUT "$header_line\t$gene_symbol_list_name_string\tgene_category_sum\n";
  foreach my $l (sort {$a <=> $b} keys %data){
    my @tmp;
    my $sum;
    for my $gene_name_column (@gene_name_columns) {
        foreach my $gene_symbol_list_name (@gene_symbol_list_names){
            push (@tmp, $data{$l}{$gene_name_column}{$gene_symbol_list_name});
        }
        $sum += $data{$l}{$gene_name_column}{sum};
    }
    my $new_cols_string = join("\t", @tmp);
    print OUT "$data{$l}{record}\t$new_cols_string\t$sum\n";
  }
  close(OUT);

  #Set output value
  $self->category_outfile($category_outfile);

  return 1;
}


1;

