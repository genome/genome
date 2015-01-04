#Written by Malachi Griffith
package Genome::Model::ClinSeq::Command::Converge::Stats;
use strict;
use warnings;
use Genome;
use Getopt::Long;
use File::Basename;

class Genome::Model::ClinSeq::Command::Converge::Stats {
  is => 'Genome::Model::ClinSeq::Command::Converge::Base',
  has_input => [  
    stats => {
      is => 'Text',
      doc => 'Stats file to summarize',
      valid_values => [qw(all wgs_snv_summary exome_snv_summary
              wgs_exome_snv_summary rnaseq_tumor_cufflinks_isoforms
              rnaseq_tumor_cufflinks_isforms_merged snv_indel_report_stats
              rnaseq_tumor_cufflinks_genes)],
      default => 'all',
    },
    plot => {
      is => 'Boolean',
      doc => 'When set, converge and plot stats. If false, stats are just ' .
        'converged and not plotted.',
      default => 0,
    },
  ],
  doc => 'converge Stats from mutiple clinseq builds'
};

sub help_detail {
  return <<EOS
  For a group of ClinSeq models, converge all Stats.tsv files
  (RNA-seq library quality metrics, SNV concordance, etc.)
  All stats should be pre-calculated.
  Produce an output matrix in which each row is metric and
  each column is a library
  Input:
  A list of Clinseq builds, models, or a Clinseq model-group.
EOS
}

sub help_synopsis {
  return <<EOS
  Example usage:
  genome model clin-seq converge stats  --builds='model_groups.id=786367aa2edc41e1b4a5d33787a8c003,is_last_complete=1' --outdir=/tmp/
  genome model clin-seq converge stats  --builds='model_groups.id=786367aa2edc41e1b4a5d33787a8c003,is_last_complete=1' --outdir=/tmp/ --plot

  Specify *one* of the following as input (each model/build should be a ClinSeq model)
  --build_ids                Comma separated list of specific build IDs
  --model_ids                Comma separated list of specific model IDs
  --model_group_id           A singe genome model group ID

  Test Clinseq model groups:
  32264                      ClinSeq - TechD RNA-seq library comparison - SPIA trimmed data
  32181                      ClinSeq - TechD RNA-seq library comparison
  30176                      ClinSeq - RNA-seq (polyA) vs cDNA Capture (NimbleGen Exome v2) Evaluation
  ";
EOS
}

sub resolve_which_stats_tsv {
  my $self = shift;
  my $b = shift;
  my @stats_files;
  my $stats_file_type  = $self->stats;
  if($stats_file_type =~ m/all/i) {
    push @stats_files, $b->input_summary_stats_file;
    push @stats_files, $b->wgs_snv_summary_stats_file;
    push @stats_files, $b->exome_snv_summary_stats_file;
    push @stats_files, $b->wgs_exome_snv_summary_stats_file;
    push @stats_files, $b->rnaseq_tumor_cufflinks_isoforms_stats_file;
    push @stats_files, $b->rnaseq_tumor_cufflinks_isoforms_merged_stats_file;
    push @stats_files, $b->rnaseq_tumor_cufflinks_genes_stats_file;
    push @stats_files, $b->rnaseq_tumor_tophat_junctions_absolute_stats_file;
    push @stats_files, $b->wgs_cnv_summary_stats_file;
    push @stats_files, $b->sv_stats_file;
    push @stats_files, $b->variant_sc_wgs_stats_file;
    push @stats_files, $b->variant_sc_exome_stats_file;
    $self->add_snv_indel_report_stats_file($b, \@stats_files);
  }
  if($stats_file_type =~ m/wgs_snv_summary/i) {
    push @stats_files, $b->wgs_snv_summary_stats_file;
  }
  if($stats_file_type =~ m/exome_snv_summary/i) {
    push @stats_files, $b->exome_snv_summary_stats_file;
  }
  if($stats_file_type =~ m/wgs_exome_snv_summary/i) {
    push @stats_files, $b->wgs_exome_snv_summary_stats_file;
  }
  if($stats_file_type =~ m/rnaseq_tumor_cufflinks_isoform/i) {
    push @stats_files, $b->rnaseq_tumor_cufflinks_isoforms_stats_file;
  }
  if($stats_file_type =~ m/rnaseq_tumor_cufflinks_isoforms_merged/i) {
    push @stats_files, $b->rnaseq_tumor_cufflinks_isoforms_merged_stats_file;
  }
  if($stats_file_type =~ m/rnaseq_tumor_cufflinks_genes/i) {
    push @stats_files, $b->rnaseq_tumor_cufflinks_genes_stats_file;
  }
  if($stats_file_type =~ m/rnaseq_tumor_tophat_junctions_absolute/i) {
    push @stats_files, $b->rnaseq_tumor_tophat_junctions_absolute_stats_file;
  }
  if($stats_file_type =~ m/wgs_cnv_summary/i) {
    push @stats_files, $b->wgs_cnv_summary_stats_file;
  }
  if($stats_file_type =~ m/input_summary/i) {
    push @stats_files, $b->input_summary_stats_file;
  }
  if($stats_file_type =~ m/sv_summary/i) {
    push @stats_files, $b->sv_stats_file;
  }
  if($stats_file_type =~ m/variant_sc_summary/i) {
    push @stats_files, $b->variant_sc_exome_stats_file;
    push @stats_files, $b->variant_sc_wgs_stats_file;
  }
  if($stats_file_type =~ m/snv_indel_report_stats/i) {
    $self->add_snv_indel_report_stats_file($b, \@stats_files);
  }
  return @stats_files;
}

sub find_stats_tsv {
  my $self = shift;
  my $mb  = shift;
  my $files = shift;
  foreach my $c (keys %$mb){
    my $b = $mb->{$c}->{build};
    my $m = $mb->{$c}->{model};
    my $model_name = $m->name;
    my $data_directory = $b->data_directory;
    my $subject_name = $m->exome_model->tumor_model->subject->common_name;
    $subject_name =~ s/[\s,-]/_/g;
    my $subject_common_name = $b->subject->common_name;
    my $build_id = $b->id;

    #If the subject name is not defined, die
    unless ($subject_name){
      die $self->error_message("\n\nCould not determine subject name for build: $build_id\n\n");
    }

    my $final_name = "Unknown";
    if ($subject_name){$final_name = $subject_name;}
    if ($subject_common_name){$final_name = $subject_common_name;}
    $self->status_message("\n\t$final_name\t$build_id\t$data_directory");

    my %file_info;
    my @stats_files = $self->resolve_which_stats_tsv($b);
    foreach my $stats_file (@stats_files) {
      if($stats_file) {
        if($stats_file =~ /^$data_directory\/.*?\/(.*)/){
          my $distinct_name = $1;
          $file_info{$stats_file}{distinct_name} = $distinct_name;
        } else{
          die $self->error_message("\n\nCould not determine distinct name " . 
              "from path\n$stats_file\n\n");
        }
        $files->{$build_id}{stats_files} = \%file_info;
        $files->{$build_id}{final_name} = $final_name;
        $files->{$build_id}{subject_name} = $subject_name;
        $files->{$build_id}{model_name} = $model_name;
      }
      else {
        $self->warning_message("\n\tWarning.  Could not find any Stats.tsv files for build: $build_id");
      }
    }
  }
}

sub get_column_names {
  my $self = shift;
  my $files = shift;
  my $file_count = shift;
  my $id_choice;
  my %finalnames;
  my %finalnames_subjectnames;
  my %modelnames;
  my %finalnames_subjectnames_builds;

  foreach my $bid (keys %$files){
    my $final_name = $files->{$bid}{final_name};
    my $subject_name = $files->{$bid}{subject_name};
    my $model_name = $files->{$bid}{model_name};
    my $id1 = $final_name;
    my $id2 = "$final_name"."_"."$subject_name";
    my $id3 = $final_name . "_" . $model_name;
    my $id4 = "$final_name"."_"."$subject_name"."_"."$bid";

    $finalnames{$id1}=1;
    $finalnames_subjectnames{$id2}=1;
    $modelnames{$id3}=1;
    $finalnames_subjectnames_builds{$id4}=1;

    $files->{$bid}{id1} = $id1;
    $files->{$bid}{id2} = $id2;
    $files->{$bid}{id3} = $id3;
    $files->{$bid}{id4} = $id4;
  }

  if ($file_count == keys %finalnames){
    $id_choice = "id1";
  }elsif($file_count == keys %finalnames_subjectnames){
    $id_choice = "id2";
  }elsif($file_count == keys %modelnames){
    $id_choice = "id3";
  }elsif($file_count == keys %finalnames_subjectnames_builds){
    $id_choice = "id4";
  }else{
    die $self->error_message("\n\nCould not generate a suitable set of distinct labels for the data files to be joined...\n\n");
  }
  return $id_choice;
}

sub aggregate_stats {
  my $self = shift;
  my %args = @_;
  my $models_builds = $args{'-models_builds'};
  my $verbose = $args{'-verbose'};
  my %files;
  my $fc = 0;
  if ($verbose) { 
    $self->status_message("\n\nGet all Stats files within these builds that match 'Stats.tsv'");
  }
  my %mb = %{$models_builds->{cases}};
  $self->find_stats_tsv(\%mb, \%files);
  my $file_count =  keys %files;
  unless ($file_count > 1){
    die $self->error_message("\n\nFound $file_count files to join (need at least 2 ...)\n\nLooked for files of name: Stats.tsv\n\n");
  }
  my $id_choice = $self->get_column_names(\%files, $file_count); 
  foreach my $fc (keys %files){
    $files{$fc}{column_name} = $files{$fc}{"$id_choice"};
  }
  return(\%files);
}

sub parse_metrics{
  my $self = shift;
  my %args = @_;
  my %files = %{$args{'-files'}};

  my %grand_metrics;
  my %metric_list;
  my $o = 0;
  foreach my $build_id (sort keys %files){
    my %local_metrics;
    my %stats_files = %{$files{$build_id}{stats_files}};
    my $column_name = $files{$build_id}{column_name};
    foreach my $file (sort keys %stats_files){
      my $distinct_name = $stats_files{$file}{distinct_name};

      open (STATS, "$file") || die "\n\nCould not open stats file: $file\n\n";
      my %columns;
      my $header = 1;
      while(<STATS>){
        chomp($_);
        my @line = split("\t", $_);

        #Perform sanity check on header.
        #make sure it conforms to the ClinSeq stats.tsv standard
        if ($header){
          my $p = 0;
          foreach my $col (@line){
            $columns{$col}{position} = $p;
            $p++;
          }
          $header = 0;
          unless ($columns{'Question'} && $columns{'Answer'} &&
            $columns{'Data_Type'} && $columns{'Analysis_Type'} &&
            $columns{'Statistic_Type'} && $columns{'Extra_Description'}){
              die $self->error_message("\n\nRequired column missing from file:".
                "\nRequired columns: Question, Answer, Data_Type," .
                " Analysis_Type, Statistic_Type, Extra_Description\nFile: " .
                "$file\nHeader: @line\n\n");
          }
          next();
        }

        #Parse the metrics data and store in a hash keyed on: build id
        #(one per column in the final output) &&
        #a concatenated string unique to the metric
        my $question = "\"" . $line[$columns{'Question'}{position}] . "\"";
        my $answer = $line[$columns{'Answer'}{position}];
        my $data_type = "\"" . $line[$columns{'Data_Type'}{position}] . "\"";
        my $analysis_type = "\"" . $line[$columns{'Analysis_Type'}{position}]. "\"";
        my $statistic_type = "\"" . $line[$columns{'Statistic_Type'}{position}]. "\"";
        my $extra_description = "\"" . $line[$columns{'Extra_Description'}{position}] . "\"";
        $extra_description =~ s/pre-treatment met/tumor/;
        $extra_description =~ s/recurrence met/tumor/;
        $extra_description =~ s/relapse2/tumor/;
        my $metric_key = $question.$data_type.$analysis_type.$statistic_type.$extra_description.$distinct_name;

        #Clean up the key a bit
        $metric_key =~ s/\'//g;
        $metric_key =~ s/\%//g;
        $metric_key =~ s/\=//g;

        #A metric must not be duplicated within a metrics file!
        if (defined($local_metrics{$metric_key})){
          die $self->error_message("\n\nMetric on this line appears to be a duplicate:\n@line\n\n$file\n\n");
        }else{
          $local_metrics{$metric_key} = 1;

          #Store the descriptive metric info once for all builds in a distinct list
          unless (defined($metric_list{$metric_key})){
            $o++;
            $metric_list{$metric_key}{order} = $o;
            $metric_list{$metric_key}{question} = $question;
            $metric_list{$metric_key}{data_type} = $data_type;
            $metric_list{$metric_key}{analysis_type} = $analysis_type;
            $metric_list{$metric_key}{statistic_type} = $statistic_type;
            $metric_list{$metric_key}{extra_description} = $extra_description;
            $metric_list{$metric_key}{distinct_name} = $distinct_name;
          }
          #Store the actual metric values for each individual build
          $grand_metrics{$build_id}{$metric_key}{answer} = $answer;

        }
      }
      close(STATS);
    }
  }

  my %results;
  $results{metrics} = \%grand_metrics;
  $results{metric_list} = \%metric_list;

  return(\%results);
}

sub get_bq_mq {
  my $self = shift;
  my $bq = $self->bq || die $self->error_message("please enter base quality as a parameter");
  my $mq = $self->mq || die $self->error_message("please enter mapping quality as a parameter");
  return ($bq, $mq);
}

sub add_snv_indel_report_stats_file{
  my $self = shift;
  my $build = shift;
  my $stats_files = shift;
  my $snv_indel_report_stats_file = undef;
  my ($bq, $mq) = $self->get_bq_mq;
  if($bq and $mq) {
    $snv_indel_report_stats_file = $build->snv_indel_report_stats_file($bq, $mq);
  }
  if($snv_indel_report_stats_file) {
    push @$stats_files, $snv_indel_report_stats_file;
  } else {
    $self->warning_message("Unable to add snv_indel_report_stats_file");
  }
}

sub write_output {
  my $self = shift;
  my $metric_list = shift;
  my $metrics = shift;
  my $files = shift;
  my $column_names_s = shift;
  my $outdir = shift;
  my $outfile = $outdir . "/stats.converged.tsv";

  my $header_line = "Question\tData_Type\tAnalysis_Type\tStatistic_Type\tExtra_Description\tFile_Source\t$column_names_s";
  #No go through each build and print out the summary
  open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
  print OUT "$header_line\n";
  foreach my $metric (sort {$metric_list->{$a}->{order} <=> $metric_list->{$b}->{order}} keys %{$metric_list}){
  my $metric_info = $metric_list->{$metric}->{question} . "\t" .
    $metric_list->{$metric}->{data_type} . "\t" . $metric_list->{$metric}->{analysis_type} .
    "\t" . $metric_list->{$metric}->{statistic_type} . "\t" .
    $metric_list->{$metric}->{extra_description} . "\t" . $metric_list->{$metric}->{distinct_name};
  my @build_metrics;
  foreach my $build_id (sort {$files->{$a}->{column_name} cmp $files->{$b}->{column_name}} keys %$files){
      #Watch for undefined metrics
      my $answer = "NA";
      if (defined($metrics->{$build_id}->{$metric})){
        $answer = $metrics->{$build_id}->{$metric}->{answer};
      }
      push(@build_metrics, $answer);
    }
    my $build_metrics_s = join("\t", @build_metrics);
    print OUT "$metric_info\t$build_metrics_s\n";
  }
  close(OUT);
  $self->status_message("Wrote output to $outfile");
}

sub plot_stats {
  my $self = shift;
  my $outdir = shift;
  my $outfile = $outdir . "/stats.converged.tsv";
  my $plot_file = $outdir . "/stats.converged.pdf";
  my $plot_script = __FILE__.".R";
  my $plot_cmd = $plot_script . " " . $outfile . " " . $plot_file;
  Genome::Sys->shellcmd(cmd => $plot_cmd);
  return;
}

sub execute {
  my $self = shift;
  my $outdir = $self->outdir;
  my @builds = $self->builds;
  my @build_ids;
  my $verbose = 1;

  foreach my $build (@builds) {
    push @build_ids, $build->id;
  }

  my $models_builds = $self->getModelsBuilds('-builds'=>\@build_ids,
    '-verbose'=>$verbose);
  my %files = %{$self->aggregate_stats('-models_builds'=>$models_builds,
    '-verbose'=>$verbose)};

  #Build a hash of possible stats values.  Key it on:
  #   Question + Data_Type + Analysis_Type + Statistic_Type + Source_File
  #Make sure these keys are unique within a patient (column_name)!!
  my $result = $self->parse_metrics('-files'=>\%files);
  my $metrics = $result->{'metrics'};
  my $metric_list = $result->{'metric_list'};

  #Build the header line
  my @column_names;
  foreach my $build_id (sort {$files{$a}->{column_name} cmp
        $files{$b}->{column_name}} keys %files){
    my $column_name = $files{$build_id}{column_name};
    push(@column_names, $column_name);
  }
  my $column_names_s = join("\t", @column_names);
  $self->write_output($metric_list, $metrics, \%files, $column_names_s,
      $outdir);
  if($self->plot) {
    $self->plot_stats($outdir);
  }
  return 1;
}

1;

