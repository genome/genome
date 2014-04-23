package Genome::Model::ClinSeq::Command::Converge::DocmReport;
use strict;
use warnings;
use Genome;
use Data::Dumper;
use Genome::Info::IUB;
use Spreadsheet::WriteExcel;

class Genome::Model::ClinSeq::Command::Converge::DocmReport {
    is => 'Genome::Model::ClinSeq::Command::Converge::Base',
    has_input => [
        outdir => {
               is => 'FilesystemPath',
               doc => 'Directory where output files will be written',
        },
        docm_variants_file => {
               is => 'FilesystemPath',
               doc => 'Tab delimited variant file with a header. First five columns: chromosome_name,start,stop,reference,variant. Also a column named primary that indicates whether the site is permuted (1) or not (0).',
        },
    ],
    has_optional_input => [
        min_var_support => {
              is => 'Number',
              default => 5,
              doc => 'Minimum number of variant supporting read to consider a variant present',
        },
        min_coverage => {
              is => 'Number',
              default => '10',
              doc => 'Minimum number of aligned read to consider a position covered',
        },
        max_normal_vaf => {
              is => 'Number',
              default => 10,
              doc => 'Variants with a normal VAF greater than this (in any normal sample) will be filtered out'
        },
        min_tumor_vaf => {
              is => 'Number',
              default => 1.0,
              doc => 'Variants with a tumor VAF less than this (in any tumor sample) will be filtered out',
        },
        tmp_space => {
              is => 'Boolean',
              default => 1,
              doc => 'Perform all file I/O in /tmp and copy results to outdir at the end',
        },
        test => {
              is => 'Number',
              doc => 'Only import this many variants (for testing purposes)',
        },
        chromosome => {
              is => 'Text',
              doc => 'Limit analysis to variants on this chromosome only',
        },
    ],
};


sub help_synopsis {
  return <<EOS

genome model clin-seq converge docm-report --builds='id in ["4b7539bb10cc4b9c97577cf11f4c79a2","cdca0edf526c4fe193d3054627a5871b"]' --outdir=/tmp/docm_report/ --docm-variants-file=variants.tsv

genome model clin-seq converge docm-report --builds='model.model_groups.id=9d0fcdca2b5d4f4385b83d2f75addac4,is_last_complete=1' --outdir=/tmp/docm_report/ --docm-variants-file=variants.tsv

genome model clin-seq converge docm-report --builds='model_groups.id=9d0fcdca2b5d4f4385b83d2f75addac4,is_last_complete=1' --outdir=/tmp/docm_report/ --docm-variants-file=variants.tsv

genome model clin-seq converge docm-report --builds='model.id in ["279f50e35d2b479ea3c32486eafd4ad4","7143119a93984056ae3f32c88c9ac2a1"],is_last_complete=1' --outdir=/tmp/docm_report/ --docm-variants-file=variants.tsv

EOS
}

sub help_detail {
  return <<EOS

Create a summary spreadsheet of read counts at canonical mutation sites (SNVs and Indels) for a set of clin-seq builds 

EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);

  unless (-e $self->outdir && -d $self->outdir) {
    push @errors, UR::Object::Tag->create(
                                          type => 'error',
                                          properties => ['outdir'],
                                          desc => "Outdir: " . $self->outdir . " not found or not a directory",
                                        );
  }

  unless (-e $self->docm_variants_file) {
    push @errors, UR::Object::Tag->create(
                                          type => 'error',
                                          properties => ['docm_variants_file'],
                                          desc => "DOCM variants file: " . $self->docm_variants_file . " not found",
                                        );
  }

  return @errors;
}

sub execute {
  my $self = shift;
  my @builds = $self->builds;

  #Add trailing '/' to outdir if needed
  unless ($self->outdir =~ /\/$/){
    my $outdir = $self->outdir . "/";
    $self->outdir($outdir);
  }

  #If the user specified, perform all file I/O in /tmp and copy result over at the end
  my $original_outdir = $self->outdir;
  if ($self->tmp_space){
    my $outdir = Genome::Sys->create_temp_directory;
    $outdir .= "/" unless ($outdir =~ /\/$/);
    $self->outdir($outdir);
  }

  #Get human readable names hash, keyed on build id
  my $subject_labels = $self->resolve_clinseq_subject_labels;

  #Get the case name for this set of builds, if more than one are found, warn the user
  my $case_name = $self->get_case_name;
  $self->status_message("Producing report for individual: $case_name");

  #Get somatic builds associated with the clin-seq builds
  my $somatic_builds = $self->resolve_somatic_builds;

  #Get reference alignment builds associated with the somatic builds
  my $align_builds = $self->get_ref_align_builds('-somatic_builds'=>$somatic_builds);

  #Import the variants, save a local copy for reference, if specified by the user, filter by chr and number
  my $result = $self->gather_variants;
  my $variant_file = $result->{'variant_file'};
  my $variants = $result->{'variants'};
  my $header = $result->{'header'};


  #Get bam-readcounts for all positions for all BAM files
  my $grand_anno_count_file = $self->add_read_counts('-align_builds'=>$align_builds, '-anno_file'=>$variant_file);

  #Parse the variants and readcounts so that they can be summarized and manipulated
  #For each position determine the following by summarizing across data sets (e.g. across samples)
  # - What is the maximum tumor variant read count?
  # - What is the maximum tumor VAF?
  # - What is the maximum normal VAF?
  # - What is the minumum coverage?
  #For each sample, determine the following and produce a table of stats at the data set level (e.g. one row per sample)
  # - How many DOCM variants are supported at min_var_support or greater?
  # - How many DOCM positions are covered at min_coverage or greater?, how many are not?, percent covered?
  # - What is the mean and median coverage across all positions?
  $header = $self->parse_read_counts('-align_builds'=>$align_builds, '-grand_anno_count_file'=>$grand_anno_count_file, '-variants'=>$variants);

  #Apply filters for display of a simplified report
  $self->apply_variant_filters('-variants'=>$variants);

  #Create TSV and Excel Spreadsheet outputs at the variant level
  #- All variants, unfiltered
  #- Primary variants only, unfiltered
  #- All variants, filtered
  #- Primary variants only, filtered
  
  #Create TSV and Excel Spreadsheet outputs as the data sets (e.g. sample) level
  my $result_files = $self->print_final_files('-variants'=>$variants, '-case_name'=>$case_name, '-header'=>$header);



  #print Dumper $align_builds;


  #If the user specified, perform all file I/O in /tmp and copy result over at the end
  if ($self->tmp_space){
    my $cp_cmd = "cp -fr " . $self->outdir . "* $original_outdir";
    Genome::Sys->shellcmd(cmd => $cp_cmd);
  }

  print "\n\n";

  return 1;
};


sub gather_variants{
  my $self = shift;

  my $variant_file = $self->outdir . "docm_variants_used.tsv";
  open(VAR1, $self->docm_variants_file) || die "\n\nCould not open source DOCM variants file: $variant_file\n\n";
  open(VAR2, ">$variant_file") || die "\n\nCould not open new DOCM variants used file: $variant_file\n\n";
  my $h = 1;
  my $c = 0;
  while(<VAR1>){
    if ($h == 1){
      print VAR2 $_;
      $h = 0;
      next;
    }
    my @line = split("\t", $_);
    if ($self->chromosome){
      unless ($self->chromosome eq $line[0]){
        next;
      }
    }
    $c++;
    if ($self->test){
      if ($c > $self->test){
        last;
      }
    }
    print VAR2 $_;
  }
  close(VAR1);
  close(VAR2);

  #Parse all variants into a single hash (keyed on $chr_$start_$end_$ref_$var)
  my $header;
  my %variants;

  open (VAR, $variant_file) || die $self->error_message("could not open file: $variant_file");
  while(<VAR>){
    chomp($_);
    if ($_ =~ /^chromosome\_name/){
      $header = $_;
      next;
    }
    my @line = split("\t", $_);
    my ($chr, $start, $stop, $ref, $var) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
    my $v = $chr . "_$start" . "_$stop" . "_$ref" . "_$var";
    $variants{$v}{anno_line} = $_;
    $variants{$v}{filtered} = 0;
  }
  close(VAR);

  my %result;
  $result{'variant_file'} = $variant_file;
  $result{'variants'} = \%variants;
  $result{'header'} = $header;

  return(\%result);
}


sub parse_read_counts{
  my $self = shift;
  my %args = @_;
  my $align_builds = $args{'-align_builds'};
  my $grand_anno_count_file = $args{'-grand_anno_count_file'};
  my $variants = $args{'-variants'};

  my $header;
  my %columns;
  my $l = 0;
  open (VAR, $grand_anno_count_file) || die $self->error_message("could not open var anno count file: $grand_anno_count_file");
  while(<VAR>){
    chomp($_);
    my @line = split("\t", $_);
    if ($l == 0){
      my $c = 0;
      foreach my $column (@line){
        $columns{$column}{c} = $c;
        $c++;
      }
      #Check for neccessary columns
      foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys  %{$align_builds}){
        my $prefix = $align_builds->{$name}->{prefix};
        my $ref_count_colname = $prefix . "_ref_count";
        my $var_count_colname = $prefix . "_var_count";
        my $vaf_colname = $prefix . "_VAF";
        die $self->error_message("could not find expected columns ($ref_count_colname $var_count_colname $vaf_colname primary) in file $grand_anno_count_file") unless($columns{$ref_count_colname} && $columns{$var_count_colname} && $columns{$vaf_colname} && $columns{'primary'});
      }
      $l++;
      $header = $_;
      next;
    }
    $l++;

    my ($chr, $start, $stop, $ref, $var) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
    my $v = $chr . "_$start" . "_$stop" . "_$ref" . "_$var";
    die $self->error_message("parsed a variant that is not defined in the variant hash") unless $variants->{$v};

    $variants->{$v}->{anno_line} = $_;
    $variants->{$v}->{primary} = $line[$columns{'primary'}{c}];

    #Get max normal VAF (across all 'normal' samples) and min coverage (across all samples)
    my $max_normal_vaf_observed = 0;
    my $max_tumor_vaf_observed = 0;
    my $max_tumor_var_count_observed = 0;
    my $min_coverage_observed = 10000000000000000;
    my $na_found = 0;
    my @covs;
    foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys  %{$align_builds}){
      my $prefix = $align_builds->{$name}->{prefix};
      my $sample_common_name = $align_builds->{$name}->{sample_common_name};
      my $ref_count_colname = $prefix . "_ref_count";
      my $var_count_colname = $prefix . "_var_count";
      my $vaf_colname = $prefix . "_VAF";

      my $vaf = "NA";
      $vaf = $line[$columns{$vaf_colname}{c}] if (defined($line[$columns{$vaf_colname}{c}]));

      if ($vaf eq "NA"){
        $na_found = 1;
        push(@covs, "NA");
        next;
      }
      if ($sample_common_name =~ /normal/){
        my $normal_vaf = $line[$columns{$vaf_colname}{c}];
        $max_normal_vaf_observed = $normal_vaf if ($normal_vaf > $max_normal_vaf_observed);
      }else{
        my $tumor_vaf = $line[$columns{$vaf_colname}{c}];
        $max_tumor_vaf_observed = $tumor_vaf if ($tumor_vaf > $max_tumor_vaf_observed);
        my $tumor_var_count = $line[$columns{$var_count_colname}{c}];
        $max_tumor_var_count_observed = $tumor_var_count if ($tumor_var_count > $max_tumor_var_count_observed);
      }
      my $coverage = $line[$columns{$ref_count_colname}{c}] + $line[$columns{$var_count_colname}{c}];
      push(@covs, $coverage);
      $min_coverage_observed = $coverage if ($coverage < $min_coverage_observed);
    }

    if ($na_found){
      $max_normal_vaf_observed = "NA";
      $max_tumor_vaf_observed = "NA";
      $max_tumor_var_count_observed = "NA";
      $min_coverage_observed = "NA";
    }

    $variants->{$v}->{max_normal_vaf_observed} = $max_normal_vaf_observed;
    $variants->{$v}->{max_tumor_vaf_observed} = $max_tumor_vaf_observed;
    $variants->{$v}->{max_tumor_var_count_observed} = $max_tumor_var_count_observed;
    $variants->{$v}->{min_coverage_observed} = $min_coverage_observed;
    $variants->{$v}->{coverages} = \@covs;
  }
  close(VAR);

  #Get a summary of values observed for each data set (e.g. sample)
  #- How many DOCM variants are supported at min_var_support or greater?
  #- How many DOCM positions are covered at min_coverage or greater?, how many are not?, percent covered?
  #- What is the mean and median coverage across all positions?
  my $summary_file = $self->outdir . "summary.tsv";
  open (SUMMARY, ">$summary_file") || die "\n\nCould not open new summary file: $summary_file\n\n";
  print SUMMARY "sample\tposition_count\tpositions_covered_10x\tpositions_covered_50x\tpositions_covered_100x\tpositions_covered_500x\tpositions_covered_1000x\tsupported_variants\tcoverage_sum\tcoverage_mean\tcoverage_median\n";
  
  my %summary;
  foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys  %{$align_builds}){
    my $prefix = $align_builds->{$name}->{prefix};
    my $libs = $align_builds->{$name}->{libs};

    my $supported_variants = 0;
    my $position_count = 0;
    my $positions_covered_10x = 0;
    my $positions_covered_50x = 0;
    my $positions_covered_100x = 0;
    my $positions_covered_500x = 0;
    my $positions_covered_1000x = 0;
    my $coverage_sum = 0;
    my @grand_cov;

    my $ref_colname = $prefix . "_ref_count";
    my $var_colname = $prefix . "_var_count";
    my $vaf_colname = $prefix . "_VAF";

    open (VAR, $grand_anno_count_file) || die $self->error_message("could not open var anno count file: $grand_anno_count_file");
    my $h = 1;
    while(<VAR>){
      chomp($_);
      if ($h){
        $h = 0;
        next;
      }
      my @line = split("\t", $_);
      my $ref = $line[$columns{$ref_colname}{c}];
      my $var = $line[$columns{$var_colname}{c}];
      my $vaf = $line[$columns{$vaf_colname}{c}];
      $position_count++;
      if($vaf eq "NA"){
        my $coverage = 0;
        push(@grand_cov, $coverage);
      }else{
        my $coverage = $ref + $var;
        $supported_variants++ if ($var >= $self->min_var_support);
        $positions_covered_10x++ if ($coverage >= 10);
        $positions_covered_50x++ if ($coverage >= 50);
        $positions_covered_100x++ if ($coverage >= 100);
        $positions_covered_500x++ if ($coverage >= 500);
        $positions_covered_1000x++ if ($coverage >= 1000);
        push(@grand_cov, $coverage);
        $coverage_sum += $coverage;
      }
    }
    close(VAR);
  
    my $coverage_mean = $coverage_sum/$position_count;
    my @grand_cov_sort = sort {$a <=> $b} @grand_cov;
    my $center = sprintf("%.0f", ($position_count/2));
    my $coverage_median = $grand_cov_sort[$center];
    print SUMMARY "$prefix\t$position_count\t$positions_covered_10x\t$positions_covered_50x\t$positions_covered_100x\t$positions_covered_500x\t$positions_covered_1000x\t$supported_variants\t$coverage_sum\t$coverage_mean\t$coverage_median\n";
  }
  close SUMMARY;

  return $header;
}


sub apply_variant_filters{
  my $self = shift;
  my %args = @_;
  my $variants = $args{'-variants'};
  my $min_var_support = $self->min_var_support;
  my $min_coverage = $self->min_coverage;
  my $max_normal_vaf = $self->max_normal_vaf;
  my $min_tumor_vaf = $self->min_tumor_vaf;

  foreach my $v (keys %{$variants}){
    
    #Make sure at least one sample had a variant readcount greater than $self->min_var_support
    my $max_tumor_var_count_observed = $variants->{$v}->{max_tumor_var_count_observed};
    if ($max_tumor_var_count_observed =~ /\d+/){
      $variants->{$v}->{filtered} = 1 if ($max_tumor_var_count_observed < $min_var_support);
    }elsif($max_tumor_var_count_observed eq "NA"){
      $variants->{$v}->{filtered} = 1;
    }

    #Make sure the coverage observed across all samples was greater than $self->min_coverage
    my $min_coverage_observed = $variants->{$v}->{min_coverage_observed};
    if ($min_coverage_observed =~ /\d+/){
      $variants->{$v}->{filtered} = 1 if ($min_coverage_observed < $min_coverage);
    }elsif($min_coverage_observed eq "NA"){
      $variants->{$v}->{filtered} = 1;
    }

    #Make sure the highest normal VAF observed was less than $self->max_normal_vaf;
    my $max_normal_vaf_observed = $variants->{$v}->{max_normal_vaf_observed};
    if ($max_normal_vaf_observed =~ /\d+/){
      $variants->{$v}->{filtered} = 1 if ($max_normal_vaf_observed > $max_normal_vaf);
    }elsif($max_normal_vaf_observed eq "NA"){
      $variants->{$v}->{filtered} = 1;
    }

    #Make sure the highest tumor VAF observed was greater than $self->min_tumor_vaf;
    my $max_tumor_vaf_observed = $variants->{$v}->{max_tumor_vaf_observed};
    if ($max_tumor_vaf_observed =~ /\d+/){
      $variants->{$v}->{filtered} = 1 if ($max_tumor_vaf_observed < $min_tumor_vaf);
    }elsif($max_tumor_vaf_observed eq "NA"){
      $variants->{$v}->{filtered} = 1;
    }

  }
  return;
}

sub print_final_files{
  my $self = shift;
  my %args = @_;
  my $variants = $args{'-variants'};
  my $case_name = $args{'-case_name'};
  my $header = $args{'-header'};

  my $tsv_all_unfiltered = $self->outdir . "$case_name" . "_all_unfiltered.tsv"; #OUT1
  my $tsv_primary_unfiltered = $self->outdir . "$case_name" . "_primary_unfiltered.tsv"; #OUT2
  my $tsv_all_filtered = $self->outdir . "$case_name" . "_all_filtered.tsv"; #OUT3
  my $tsv_primary_filtered = $self->outdir . "$case_name" . "_primary_filtered.tsv"; #OUT4

  open (OUT1, ">$tsv_all_unfiltered") ||  die "\n\ncould not open new output file: $tsv_all_unfiltered\n\n";
  open (OUT2, ">$tsv_primary_unfiltered") ||  die "\n\ncould not open new output file: $tsv_primary_unfiltered\n\n";
  open (OUT3, ">$tsv_all_filtered") ||  die "\n\ncould not open new output file: $tsv_all_filtered\n\n";
  open (OUT4, ">$tsv_primary_filtered") ||  die "\n\ncould not open new output file: $tsv_primary_filtered\n\n";

  my $full_header = $header . "\tmax_normal_vaf_observed\tmax_tumor_vaf_observed\tmax_tumor_var_count_observed\tmin_coverage_observed\tfiltered\n";

  print OUT1 "$full_header";
  print OUT2 "$full_header";
  print OUT3 "$full_header";
  print OUT4 "$full_header";

  foreach my $v (sort keys %{$variants}){
    my $filtered = $variants->{$v}->{filtered};
    my $primary = $variants->{$v}->{primary};
    my $anno_line = $variants->{$v}->{anno_line};
    my $max_normal_vaf_observed = $variants->{$v}->{max_normal_vaf_observed};
    my $max_tumor_vaf_observed = $variants->{$v}->{max_tumor_vaf_observed};
    my $max_tumor_var_count_observed = $variants->{$v}->{max_tumor_var_count_observed};
    my $min_coverage_observed = $variants->{$v}->{min_coverage_observed};

    my $full_line = "$anno_line\t$max_normal_vaf_observed\t$max_tumor_vaf_observed\t$max_tumor_var_count_observed\t$min_coverage_observed\t$filtered\n";

    #Create TSV and Excel Spreadsheet outputs at the variant level
    #- All variants, unfiltered
    print OUT1 "$full_line";  

    #- Primary variants only, unfiltered
    if ($primary == 1){
      print OUT2 "$full_line";
    }

    #- All variants, filtered
    unless ($filtered == 1){
      print OUT3 "$full_line";
    }

    #- Primary variants only, filtered
    unless ($filtered ==1){
      if ($primary == 1){
        print OUT4 "$full_line";
      }
    }
  }
  close(OUT1);
  close(OUT2);
  close(OUT3);
  close(OUT4);


  #Create an Excel version of all TSV files present
  my @files;
  opendir(DIR, $self->outdir) or die $!;
  while (my $file = readdir(DIR)) {
    next unless ($file =~ /tsv/);
    push(@files, $self->outdir . $file);
  }
  close(DIR);

  foreach my $tsv_file (@files){
    my $xls_file = $tsv_file;
    $xls_file =~ s/tsv/xls/;
    my $workbook  = Spreadsheet::WriteExcel->new("$xls_file");
    my $worksheet = $workbook->add_worksheet();
    open (IN, $tsv_file) || die $self->error_message("Could not open in file: $tsv_file");
    my $row=0;
    while(<IN>){
      chomp($_);
      if ($row == 0){
        $_ =~ s/\_/ /g;
      }
      my @F = split("\t",$_);
      for(my $i=0;$i<@F;$i++){
        $worksheet->write($row, $i, $F[$i]);
      }
      $row++;
    }
    close(IN);
    $workbook->close();
  }

  return;
}


1;
