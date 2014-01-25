package Genome::Model::ClinSeq::RnaSeqUtil;

#Written by Malachi Griffith

require Exporter;

@ISA = qw( Exporter );
@EXPORT = qw();

@EXPORT_OK = qw(
                &parseFpkmFile &mergeIsoformsFile &calculateCufflinksStats
               );

%EXPORT_TAGS = (
                all => [qw(&parseFpkmFile &mergeIsoformsFile &calculateCufflinksStats)]
               );

use strict;
use warnings;
use Data::Dumper;
use Genome;
use Genome::Model::ClinSeq::Util qw(:all);
#TODO - Make a class.

#################################################################################################################
#Parse the genes.fpkm_trackin file                                                                              #
#################################################################################################################
sub parseFpkmFile{
  my $self = shift;
  my %args = @_;
  my $infile = $args{'-infile'};
  my $entrez_ensembl_data = $args{'-entrez_ensembl_data'};
  my $ensembl_tmap = $args{'-ensembl_map'};
  my $verbose = $args{'-verbose'};
  my $outfile;
  if (defined($args{'-outfile'})){
    $outfile = $args{'-outfile'};
  }

  if ($verbose){
    $self->status_message("\n\nParsing: $infile");
  }

  #Note that we could be parsing genes or transcripts here
  #Will therefore need ensembl maps keyed on both genes and transcripts
  my %ensembl_gmap;
  my $ensembl_gmap = \%ensembl_gmap;
  foreach my $tid (keys %{$ensembl_tmap}){
    my $gid = $ensembl_tmap->{$tid}->{ensg_id};
    my $gene_biotype = $ensembl_tmap->{$tid}->{gene_biotype};
    $ensembl_gmap->{$gid}->{gene_biotype} = $gene_biotype;
  }
  
  my %fpkm;
  my $header = 1;
  my $rc = 0;     #record count
  my %columns;
  open (FPKM, "$infile") || die "\n\nCould not open gene/isoform file: $infile\n\n";
  while(<FPKM>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header == 1){
      $header = 0;
      my $p = 0;
      foreach my $head (@line){
        $columns{$head}{position} = $p;
        $p++;
      }
      next();
    }
    $rc++;
    
    my $tracking_id = $line[$columns{'tracking_id'}{position}];
    my $gene_id = $line[$columns{'gene_id'}{position}];
    my $locus = $line[$columns{'locus'}{position}];
    my $length = $line[$columns{'length'}{position}];
    my $coverage = $line[$columns{'coverage'}{position}];
    my $FPKM = $line[$columns{'FPKM'}{position}];
    my $FPKM_conf_lo = $line[$columns{'FPKM_conf_lo'}{position}];
    my $FPKM_conf_hi = $line[$columns{'FPKM_conf_hi'}{position}];
    my $FPKM_status;
    if ($columns{'FPKM_status'}){
      $FPKM_status = $line[$columns{'FPKM_status'}{position}];
    }elsif($columns{'status'}){
      $FPKM_status = $line[$columns{'status'}{position}];
    }else{
      die $self->error_message("\n\nRequired column not found: 'FPKM_status' or 'status'");
    }

    #Fix gene name and create a new column for this name
    my $fixed_gene_name = $self->fixGeneName('-gene'=>$gene_id, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);

    #Key on tracking id AND locus coordinates
    my $key = "$tracking_id"."|"."$locus";

    #Get biotype.  Check for a match at gene or transcript level
    my $biotype = "na";
    if ($ensembl_gmap->{$tracking_id}){$biotype = $ensembl_gmap->{$tracking_id}->{gene_biotype};}
    if ($ensembl_tmap->{$tracking_id}){$biotype = $ensembl_tmap->{$tracking_id}->{transcript_biotype};}
    $fpkm{$key}{record_count} = $rc;
    $fpkm{$key}{tracking_id} = $tracking_id;
    $fpkm{$key}{mapped_gene_name} = $fixed_gene_name;
    $fpkm{$key}{gene_id} = $gene_id;
    $fpkm{$key}{biotype} = $biotype;
    $fpkm{$key}{locus} = $locus;
    $fpkm{$key}{length} = $length;
    $fpkm{$key}{coverage} = $coverage;
    $fpkm{$key}{FPKM} = $FPKM;
    $fpkm{$key}{FPKM_conf_lo} = $FPKM_conf_lo;
    $fpkm{$key}{FPKM_conf_hi} = $FPKM_conf_hi;
    $fpkm{$key}{FPKM_status} = $FPKM_status;
  }
  close(FPKM);

  my $gc = keys %fpkm;
  unless ($gc == $rc){
    die $self->error_message("\n\nFound $gc distinct gene|coord entries but $rc data lines - not good...\n\n");
  }

  #Print an outfile sorted on the key
  if ($outfile){
    open (OUT, ">$outfile") || die "\n\nCould not open gene file: $infile\n\n";
    print OUT "tracking_id\tmapped_gene_name\tgene_id\tbiotype\tlocus\tlength\tcoverage\tFPKM\tFPKM_conf_lo\tFPKM_conf_hi\tFPKM_status\n";
    foreach my $key (sort {$a cmp $b} keys %fpkm){
      print OUT "$fpkm{$key}{tracking_id}\t$fpkm{$key}{mapped_gene_name}\t$fpkm{$key}{gene_id}\t$fpkm{$key}{biotype}\t$fpkm{$key}{locus}\t$fpkm{$key}{length}\t$fpkm{$key}{coverage}\t$fpkm{$key}{FPKM}\t$fpkm{$key}{FPKM_conf_lo}\t$fpkm{$key}{FPKM_conf_hi}\t$fpkm{$key}{FPKM_status}\n";
    }
    close(OUT);
  }
  return(\%fpkm);
}


#################################################################################################################
#Merge the isoforms.fpkm_tracking file to the gene level                                                        #
#################################################################################################################
sub mergeIsoformsFile{
  my $self = shift;
  my %args = @_;
  my $infile = $args{'-infile'};
  my $entrez_ensembl_data = $args{'-entrez_ensembl_data'};
  my $ensembl_map = $args{'-ensembl_map'};
  my $verbose = $args{'-verbose'};

  my $status_file;
  $status_file = $args{'-status_file'} if (defined($args{'-status_file'}));
  
  my $outfile;
  $outfile = $args{'-outfile'} if (defined($args{'-outfile'}));

  if ($verbose){
    $self->status_message("\n\nParsing and merging to gene level: $infile");
  }

  my %trans;
  my %genes;
  my $header = 1;
  my $rc = 0;     #record count
  my %columns;
  open (TRANS, "$infile") || die "\n\nCould not open gene file: $infile\n\n";
  while(<TRANS>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header == 1){
      $header = 0;
      my $p = 0;
      foreach my $head (@line){
        $columns{$head}{position} = $p;
        $p++;
      }
      next();
    }
    $rc++;
   
    #Note.  The tracking ID in this file should be an Ensembl transcript id.  Use this ID and the specified ensembl version to look up the ENSG ID
    my $tracking_id = $line[$columns{'tracking_id'}{position}];
    my $original_gene_id = $line[$columns{'gene_id'}{position}];

    unless($tracking_id){
      die $self->error_message("\n\nTracking ID not defined\n\n");
    }

    #Get the gene ID from the transcript ID
    unless (defined($ensembl_map->{$tracking_id})){
      die $self->error_message("\n\nCould not map tracking id: $tracking_id to an ensembl gene via ensembl transcript ID!\n\n");
    }
    my $ensg_id = $ensembl_map->{$tracking_id}->{ensg_id};
    my $ensg_name = $ensembl_map->{$tracking_id}->{ensg_name};
    my $gene_biotype = $ensembl_map->{$tracking_id}->{gene_biotype};
    my $transcript_biotype = $ensembl_map->{$tracking_id}->{transcript_biotype};

    #print "\n$rc\t$ensg_id\t$ensg_name\t$gene_biotype\t$transcript_biotype";

    my $locus = $line[$columns{'locus'}{position}];
    my $length = $line[$columns{'length'}{position}];
    my $coverage = $line[$columns{'coverage'}{position}];
    my $FPKM = $line[$columns{'FPKM'}{position}];
    my $FPKM_conf_lo = $line[$columns{'FPKM_conf_lo'}{position}];
    my $FPKM_conf_hi = $line[$columns{'FPKM_conf_hi'}{position}];
    my $FPKM_status;
    if ($columns{'FPKM_status'}){
      $FPKM_status = $line[$columns{'FPKM_status'}{position}];
    }elsif($columns{'status'}){
      $FPKM_status = $line[$columns{'status'}{position}];
    }else{
      die $self->error_message("\n\nRequired column not found: 'FPKM_status' or 'status'");
    }
 
    #Fix gene name and create a new column for this name
    my $fixed_gene_name = $self->fixGeneName('-gene'=>$ensg_name, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);

    $trans{$tracking_id}{record_count} = $rc;

    #Get coords from locus
    my $chr;
    my $chr_start;
    my $chr_end;
    if ($locus =~ /(\w+)\:(\d+)\-(\d+)/){
      $chr = $1;
      $chr_start = $2;
      $chr_end = $3;
    }else{
      die $self->error_message("\n\nlocus format not understood: $locus\n\n");
    }
  
    #Merge down to genes, combining the coverage and FPKM values (cumulatively), coordinates (outer coords), and calculating a new length
    if ($genes{$ensg_id}){
      if ($chr_start <  $genes{$ensg_id}{chr_start}){$genes{$ensg_id}{chr_start} = $chr_start;}
      if ($chr_end >  $genes{$ensg_id}{chr_end}){$genes{$ensg_id}{chr_end} = $chr_end;}
      $genes{$ensg_id}{coverage} += $coverage;
      $genes{$ensg_id}{FPKM} += $FPKM;
      $genes{$ensg_id}{FPKM_conf_lo} += $FPKM_conf_lo;
      $genes{$ensg_id}{FPKM_conf_hi} += $FPKM_conf_hi;
      $genes{$ensg_id}{FPKM_status} = "na";
      $genes{$ensg_id}{transcript_count}++;
    }else{
      $genes{$ensg_id}{ensg_name} = $ensg_name;
      $genes{$ensg_id}{mapped_gene_name} = $fixed_gene_name;
      $genes{$ensg_id}{original_gene_id} = $original_gene_id;
      $genes{$ensg_id}{chr} = $chr;
      $genes{$ensg_id}{chr_start} = $chr_start;
      $genes{$ensg_id}{chr_end} = $chr_end;
      $genes{$ensg_id}{gene_biotype} = $gene_biotype;
      $genes{$ensg_id}{coverage} = $coverage;
      $genes{$ensg_id}{FPKM} = $FPKM;
      $genes{$ensg_id}{FPKM_conf_lo} = $FPKM_conf_lo;
      $genes{$ensg_id}{FPKM_conf_hi} = $FPKM_conf_hi;
      $genes{$ensg_id}{FPKM_status} = "na";
      $genes{$ensg_id}{transcript_count} = 1;
    }
  }
  close(TRANS);

  my $tc = keys %trans;
  unless ($tc == $rc){
    die $self->error_message("\n\nFound $tc distinct transcript entries but $rc data lines - not good...\n\n");
  }

  #If an FPKM status file was defined, use it to add FPKM status value to each gene where possible
  if ($status_file){
    my %columns;
    $header = 1;
    open (STATUS, "$status_file") || die "\n\nCould not open gene file: $status_file\n\n";
    while(<STATUS>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        $header = 0;
        my $p = 0;
        foreach my $head (@line){
          $columns{$head}{position} = $p;
          $p++;
        }
        next();
      }
      my $tracking_id = $line[$columns{'tracking_id'}{position}];
      my $fpkm_status = $line[$columns{'FPKM_status'}{position}];     

      #Update the FPKM status in the %genes hash
      if ($genes{$tracking_id}){
        if (($genes{$tracking_id}{FPKM_status} eq "na" || $genes{$tracking_id}{FPKM_status} eq "OK") && $fpkm_status){
          $genes{$tracking_id}{FPKM_status} = $fpkm_status;
        }
      }
    }
    close (STATUS);
  }

  #Print an outfile sorted on the key
  if ($outfile){
    open (OUT, ">$outfile") || die "\n\nCould not open gene file: $infile\n\n";
    print OUT "tracking_id\tmapped_gene_name\tensg_name\tbiotype\tlocus\tlength\tcoverage\tFPKM\tFPKM_conf_lo\tFPKM_conf_hi\tFPKM_status\n";
    foreach my $gene_id (sort {$genes{$a}{ensg_name} cmp $genes{$b}{ensg_name}} keys %genes){
      my $locus = "$genes{$gene_id}{chr}:$genes{$gene_id}{chr_start}-$genes{$gene_id}{chr_end}";
      my $length = "-";
      print OUT "$gene_id\t$genes{$gene_id}{mapped_gene_name}\t$genes{$gene_id}{ensg_name}\t$genes{$gene_id}{gene_biotype}\t$locus\t$length\t$genes{$gene_id}{coverage}\t$genes{$gene_id}{FPKM}\t$genes{$gene_id}{FPKM_conf_lo}\t$genes{$gene_id}{FPKM_conf_hi}\t$genes{$gene_id}{FPKM_status}\n";
    }
    close(OUT);
  }
  return(\%genes);
}


#################################################################################################################
#Parse a clinseq formated cufflinks expression file, gather basic QC stats and dump them to a Stats.tsv file    #
#################################################################################################################
sub calculateCufflinksStats{
  my $self = shift;
  my %args = @_;
  my $infile = $args{'-infile'};
  my $outfile = $args{'-outfile'};

  #Assume the input file has a column 'FPKM' and base all stats on that.  Write to a Stats.tsv
  #Use the ClinSeq stats standard:
  #Question Answer  Data_Type Analysis_Type Statistic_Type  Extra_Description

  unless (-e $infile){
    die $self->error_message("\n\nInput file to &calculateCufflinksStats() could not be found\n\n");
  }

  #Import FPKM data
  my %fpkms;
  my $grand_fpkm = 0;
  my $detected_feature_count = 0;
  my $header = 1;
  my $l = 0;
  my %columns;
  open (FPKM, "$infile") || die "\n\nCould not open input file: $infile\n\n";
  while(<FPKM>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header){
      my $p = 0;
      foreach my $col (@line){
        $columns{$col}{position} = $p;
        $p++;
      }
      $header = 0;
      next();
    }
    $l++;
    my $fpkm = $line[$columns{'FPKM'}{position}];
    $fpkms{$l}{fpkm} = $fpkm;
    $grand_fpkm += $fpkm;
    if ($fpkm > 0){
      $detected_feature_count++;
    }
  }
  close (FPKM);

  #Summarize:
  #A.) Total Gene/Transcript count
  #B.) Transcripts/Genes observed at > FPKM of: 0, 0.001, 0.01, 0.1, 1, 5, 10, 50, 100, 500
  #C.) Percent of cumulative FPKM consumed by top N% of Transcripts/Genes: 0.01, 0.1, 1, 2, 5, 10, 20, 30, 40, 50
  my $total_feature_count = keys %fpkms;

  my @fpkm_cuts = qw ( 0 0.001 0.01 0.1 1 5 10 50 100 500 );
  my @topp_cuts = qw ( 0.01 0.1 1 2 5 10 20 30 40 50 );

  my %expressed;
  my %consumed;

  foreach my $fc (@fpkm_cuts){
    $expressed{$fc}{count} = 0;
  }
  
  foreach my $tc (@topp_cuts){
    $consumed{$tc}{detected_feature_count} = $detected_feature_count;
    $consumed{$tc}{total_feature_count} = $total_feature_count;
    $consumed{$tc}{percent_consumed} = 0;
    $consumed{$tc}{bin_size} = sprintf("%.0f", ($detected_feature_count * ($tc/100)));
  }

  my $c = 0;
  my $cumulative_fpkm = 0;
  foreach my $l (sort {$fpkms{$b}{fpkm} <=> $fpkms{$a}{fpkm}} keys %fpkms){
    my $fpkm = $fpkms{$l}{fpkm};
    $cumulative_fpkm += $fpkm;
    foreach my $fc (@fpkm_cuts){
      if ($fpkm > $fc){
        $expressed{$fc}{count}++;
      }
    }
    $c++;
    foreach my $tc (@topp_cuts){
      my $bin_size = $consumed{$tc}{bin_size};
      if ($c <= $bin_size){
        $consumed{$tc}{cumulative_fpkm} = $cumulative_fpkm;
        $consumed{$tc}{grand_fpkm} = $grand_fpkm;
        $consumed{$tc}{percent_consumed} = sprintf("%.3f", (($cumulative_fpkm / $grand_fpkm)*100));
      }
    }
    if ($fpkm == 0){
      last();
    }
  }

  #Print out the summary stats
  #Question Answer  Data_Type Analysis_Type Statistic_Type  Extra_Description

  open (STATS, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
  print STATS "Question\tAnswer\tData_Type\tAnalysis_Type\tStatistic_Type\tExtra_Description\n";
  print STATS "Total Cufflinks features count\t$total_feature_count\tRNA-seq\tCufflinks Expression\tCount\tTotal number of genes/transcripts with possible detection by Cufflinks\n";

  #Genes observed at >= 1X	14641	RNA-seq	Alignments	Count	Number of genes (measured by their junctions) observed by >= 1 reads
  #Percent of junctions reads consumed by top 0.01% of detected genes	0.868224775060336	RNA-seq	Alignments	Percent	Percent of junctions reads consumed by top 0.01% of detected genes (>=1 junction read for the gene)
  foreach my $fc (@fpkm_cuts){
    my $count = $expressed{$fc}{count};
    my $percent = sprintf("%.3f", (($count / $total_feature_count)*100));
    print STATS "Number of features observed at FPKM >= $fc\t$count\tRNA-seq\tCufflinks Expression\tCount\tNumber of genes/transcripts with Cufflinks FPKM > $fc\n";
    print STATS "Percent of features observed at FPKM >= $fc\t$percent\tRNA-seq\tCufflinks Expression\tPercent\tPercent of all possible genes/transcripts that were detected by Cufflinks with FPKM > $fc\n";
  }
  foreach my $tc (@topp_cuts){
    print STATS "Percent of FPKM consumed by top $tc% of detected genes\t$consumed{$tc}{percent_consumed}\tRNA-seq\tCufflinks Expression\tPercent\tPercent of cumulative Cufflinks FPKM consumed by top $tc% of detected (FPKM > 0) features\n";
  }
  close(STATS);

  return();
}


1;

