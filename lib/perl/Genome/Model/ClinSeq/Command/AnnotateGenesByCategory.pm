package Genome::Model::ClinSeq::Command::AnnotateGenesByCategory;

use strict;
use warnings;
use Genome;
use Data::Dumper;
use Genome::Model::ClinSeq::Util qw(:all);

class Genome::Model::ClinSeq::Command::AnnotateGenesByCategory {
    is => 'Command::V2',
    has_input => [
        infile => {
            is => 'FilesystemPath',
            doc => 'tsv formatted file to added gene category annotations to',
        },
        gene_symbol_lists_dir => {
            is => 'FilesystemPath',
            doc => 'directory containing list of gene names by category',
        },
        gene_name_column => {
            is => 'Text',
            doc => 'name of column containing gene names/symbols'
        }
    ],
    has_output => [
        category_outfile => {
            is => 'FilesystemPath',            
            doc => 'result file consisting of input file appended with gene category values',
        }
    ],
    doc => 'take a tsv file as input, indentify gene names, append matches between each of these genes and various gene categories (e.g. kinase, etc.)',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq annotate-genes-by-category --infile=/gscmnt/ams1108/info/model_data/2888708572/build134369422/AML103/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.readcounts.tsv --gene-symbol-lists-dir=/gscmnt/sata132/techd/mgriffit/reference_annotations/GeneSymbolLists/  --gene-name-column='mapped_gene_name'

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

  unless (-e $self->gene_symbol_lists_dir && -d $self->gene_symbol_lists_dir) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['gene_symbol_lists_dir'],
	                                          desc => "symbol list dir: " . $self->gene_symbol_lists_dir . " not a valid directory",
                                          );
  }

  return @errors;
}


sub execute {
  my $self = shift;
  my $infile = $self->infile;
  my $gene_symbol_lists_dir = $self->gene_symbol_lists_dir;
  my $gene_name_column = $self->gene_name_column;

  $self->status_message("Adding gene category annotations to $infile using genes in this gene name column: $gene_name_column");

  #Get Entrez and Ensembl data for gene name mappings
  my $entrez_ensembl_data = &loadEntrezEnsemblData();

  my $symbol_list_names = &importSymbolListNames('-gene_symbol_lists_dir'=>$gene_symbol_lists_dir, '-verbose'=>0);
  my $master_list = $symbol_list_names->{master_list};
  my @symbol_list_names = sort {$master_list->{$a}->{order} <=> $master_list->{$b}->{order}} keys %{$master_list};
  my $gene_symbol_lists = &importGeneSymbolLists('-gene_symbol_lists_dir'=>$gene_symbol_lists_dir, '-symbol_list_names'=>\@symbol_list_names, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);


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
      unless ($cols{$gene_name_column}){
        die $self->error_message("File has no $gene_name_column column: $infile");
      }
      next();
    }
    $data{$l}{record} = $record;
    $data{$l}{gene_name} = $line[$cols{'mapped_gene_name'}{position}];
    $data{$l}{sum} = 0;
  }
  close(INFILE);

  #Figure out the gene matches to the gene symbol lists
  #Test each gene name in this column against those in the list and add a column with the match status (i.e. is it a kinase, cancer gene, etc.)
  foreach my $l (keys %data){
    my $gene_name = $data{$l}{gene_name};
    foreach my $gene_symbol_type (keys %{$gene_symbol_lists}){
      my $gene_symbols = $gene_symbol_lists->{$gene_symbol_type}->{symbols};
      if ($gene_symbols->{$gene_name}){
        $data{$l}{$gene_symbol_type} = 1;
        $data{$l}{sum}++;
      }else{
        $data{$l}{$gene_symbol_type} = 0;
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
  my @gene_symbol_list_names = sort {$gene_symbol_lists->{$a}->{order} <=> $gene_symbol_lists->{$b}->{order}} keys %{$gene_symbol_lists};
  my $gene_symbol_list_name_string = join("\t", @gene_symbol_list_names);
  print OUT "$header_line\t$gene_symbol_list_name_string\tgene_category_sum\n";
  foreach my $l (sort {$a <=> $b} keys %data){
    my @tmp;
    foreach my $gene_symbol_list_name (@gene_symbol_list_names){
      push (@tmp, $data{$l}{$gene_symbol_list_name});
    }
    my $new_cols_string = join("\t", @tmp);
    print OUT "$data{$l}{record}\t$new_cols_string\t$data{$l}{sum}\n";
  }
  close(OUT);

  #Set output value
  $self->category_outfile($category_outfile);

  return 1;
}


1;

