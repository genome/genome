package Genome::Model::ClinSeq::Command::Converge::MutationSpectrum;
use strict;
use warnings;
use Genome;
use Getopt::Long;
use File::Basename;

class Genome::Model::ClinSeq::Command::Converge::MutationSpectrum {
  is => 'Genome::Model::ClinSeq::Command::Converge::Base',
  has_input => [  
    outfile => {
      is => 'FilesystemPath',
      doc => 'File to write converged-mutation-spectrum results',
    },
  ],
  doc => 'converge Stats from mutiple clinseq builds'
};

sub help_detail {
  return <<EOS
  For a group of ClinSeq models, converge mutation-spectrum results.
  All stats should be pre-calculated. Produce an output matrix in which each row is metric and each column is a sample.
  Input:
  A list of Clinseq builds, models [OR] a Clinseq model-group.
EOS
}

sub help_synopsis {
  return <<EOS
  Example usage: 
  genome model clin-seq converge mutation-spectrum --builds='model_groups.id=6fa120dc0afb400596a1e3d6ecf6167d,is_last_complete=1' --outfile=metris.tsv --outdir=/tmp/
EOS
}

sub resolve_which_stats_tsv {
  my $self = shift;
  my $b = shift;
  my $bq = $self->bq;
  my $mq = $self->mq;
  my @stats_files;
  push @stats_files, $b->mutation_spectrum_exome_summary_file($bq, $mq);
  push @stats_files, $b->mutation_spectrum_wgs_summary_file($bq, $mq);
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
        $self->warning_message("\n\tWarning. Could not find any mutation-spectrum files for build: $build_id");
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
  if ($file_count <= keys %finalnames){
    $id_choice = "id1";
  }elsif($file_count <= keys %finalnames_subjectnames){
    $id_choice = "id2";
  }elsif($file_count <= keys %modelnames){
    $id_choice = "id3";
  }elsif($file_count <= keys %finalnames_subjectnames_builds){
    $id_choice = "id4";
  }else{
    die $self->error_message("\n\nCould not generate a suitable set of distinct labels for the data files to be joined...\n\n");
  }
  return $id_choice;
}

#######################################################################################################################
#Get the desired files from each model                                                                                #
#######################################################################################################################
sub aggregate_stats {
  my $self = shift;
  my %args = @_;
  my $models_builds = $args{'-models_builds'};
  my %files;
  my $fc = 0;
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


#######################################################################################################################
#Build a hash of metrics across the input builds and their stats files                                                #
#######################################################################################################################
sub parse_metrics{
  my $self = shift;
  my %args = @_;
  my %files = %{$args{'-files'}};
  my $results = $args{'-results'};
  foreach my $build_id (sort keys %files){
    my %stats_files = %{$files{$build_id}{stats_files}};
    foreach my $file (sort keys %stats_files){
      my $distinct_name = $stats_files{$file}{distinct_name};
      open (my $STATS, "$file") || die "\n\nCould not open stats file: $file\n\n";
      while(<$STATS>) {
        my @line = split("\t", $_);
        my $mutation_type = $line[0];
        my $mutation_count = $line[1];
        my $mutation_fraction = $line[2];
        $results->{$build_id}{$file}{$mutation_type}{count} = $mutation_count;
        $results->{$build_id}{$file}{$mutation_type}{fraction} = $mutation_fraction;
        my $data_type = $self->_get_data_type($file);
        $results->{$build_id}{$file}{data_type} = $data_type;
      }
      close($STATS);
    }
  }
}

sub write_output {
  my $self = shift;
  my $results = shift;
  my $files = shift;
  my $column_names_s = shift;
  my $outfile = shift;
  my $header_line = "Mutation_Type\tStatistic_Type\tData_Type\t$column_names_s";
  #No go through each build and print out the summary
  open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
  print OUT "$header_line\n";
  my ($a_g_count, $a_g_fraction);
  my ($c_t_count, $c_t_fraction);
  my ($c_g_count, $c_g_fraction);
  my ($c_a_count, $c_a_fraction);
  my ($a_c_count, $a_c_fraction);
  my ($a_t_count, $a_t_fraction);
  my ($transition_count, $transition_fraction);
  my ($transversion_count, $transversion_fraction);
  my $data_type;
  foreach my $build_id (sort {$files->{$a}->{column_name} cmp $files->{$b}->{column_name}} keys %$files){
    my %stats_files = %{$files->{$build_id}{stats_files}};
    foreach my $file (sort keys %stats_files){
      $data_type = $results->{$build_id}{$file}{data_type};
      push @{$a_g_count->{$data_type}}, $results->{$build_id}{$file}{"A->G"}{count};
      push @{$a_g_fraction->{$data_type}}, $results->{$build_id}{$file}{"A->G"}{fraction};
      push @{$c_t_count->{$data_type}}, $results->{$build_id}{$file}{"C->T"}{count};
      push @{$c_t_fraction->{$data_type}}, $results->{$build_id}{$file}{"C->T"}{fraction};
      push @{$a_c_count->{$data_type}}, $results->{$build_id}{$file}{"A->C"}{count};
      push @{$a_c_fraction->{$data_type}}, $results->{$build_id}{$file}{"A->C"}{fraction};
      push @{$c_g_count->{$data_type}}, $results->{$build_id}{$file}{"C->G"}{count};
      push @{$c_g_fraction->{$data_type}}, $results->{$build_id}{$file}{"C->G"}{fraction};
      push @{$c_a_count->{$data_type}}, $results->{$build_id}{$file}{"C->A"}{count};
      push @{$c_a_fraction->{$data_type}}, $results->{$build_id}{$file}{"C->A"}{fraction};
      push @{$a_t_count->{$data_type}}, $results->{$build_id}{$file}{"A->T"}{count};
      push @{$a_t_fraction->{$data_type}}, $results->{$build_id}{$file}{"A->T"}{fraction};
      push @{$transition_count->{$data_type}}, $results->{$build_id}{$file}{"Transitions"}{count};
      push @{$transition_fraction->{$data_type}}, $results->{$build_id}{$file}{"Transitions"}{fraction};
      push @{$transversion_count->{$data_type}}, $results->{$build_id}{$file}{"Transversion"}{count};
      push @{$transversion_fraction->{$data_type}}, $results->{$build_id}{$file}{"Transversion"}{fraction};
    }
  }
  for $data_type (keys %$a_g_count) {
    print OUT "A->G\tcount\t$data_type\t", join("\t", @{$a_g_count->{$data_type}}, "\n");
    print OUT "C->T\tcount\t$data_type\t", join("\t", @{$c_t_count->{$data_type}}, "\n");
    print OUT "A->C\tcount\t$data_type\t", join("\t", @{$a_c_count->{$data_type}}, "\n");
    print OUT "C->G\tcount\t$data_type\t", join("\t", @{$c_g_count->{$data_type}}, "\n");
    print OUT "C->A\tcount\t$data_type\t", join("\t", @{$c_a_count->{$data_type}}, "\n");
    print OUT "A->T\tcount\t$data_type\t", join("\t", @{$a_t_count->{$data_type}}, "\n");
    print OUT "Transition\tcount\t$data_type\t", join("\t", @{$transition_count->{$data_type}}, "\n");
    print OUT "Transversion\tcount\t$data_type\t", join("\t", @{$transversion_count->{$data_type}}, "\n");
    print OUT "A->G\tfraction\t$data_type\t", join("\t", @{$a_g_fraction->{$data_type}}, "\n");
    print OUT "C->T\tfraction\t$data_type\t", join("\t", @{$c_t_fraction->{$data_type}}, "\n");
    print OUT "A->C\tfraction\t$data_type\t", join("\t", @{$a_c_fraction->{$data_type}}, "\n");
    print OUT "C->G\tfraction\t$data_type\t", join("\t", @{$c_g_fraction->{$data_type}}, "\n");
    print OUT "C->A\tfraction\t$data_type\t", join("\t", @{$c_a_fraction->{$data_type}}, "\n");
    print OUT "A->T\tfraction\t$data_type\t", join("\t", @{$a_t_fraction->{$data_type}}, "\n");
    print OUT "Transition\tfraction\t$data_type\t", join("\t", @{$transition_fraction->{$data_type}}, "\n");
    print OUT "Transversion\tfraction\t$data_type\t", join("\t", @{$transversion_fraction->{$data_type}}, "\n");
  }
  close(OUT);
  $self->status_message("Wrote output to $outfile");
}

sub _get_data_type {
  my $self = shift;
  my $stats_file = shift;
  my $data_type = "WGS/Exome";
  if($stats_file =~ /mutation-spectrum.*\/exome\/.*summarize_mutation_spectrum/) {
    $data_type = "exome";
  } elsif($stats_file =~ /mutation-spectrum.*\/wgs\/.*summarize_mutation_spectrum/) {
    $data_type = "wgs";
  }
  return $data_type;
}

sub execute {
  my $self = shift;
  my $outfile = $self->outdir . "/" . basename($self->outfile);
  my @builds = $self->builds;
  my @build_ids;

  foreach my $build (@builds) {
    push @build_ids, $build->id;
  }
  my $models_builds = $self->getModelsBuilds('-builds'=>\@build_ids);
  my %files = %{$self->aggregate_stats('-models_builds'=>$models_builds)};
  my %results;
  $self->parse_metrics('-files'=>\%files, '-results' => \%results);

  #Build the header line
  my @column_names;
  foreach my $build_id (sort {$files{$a}->{column_name} cmp $files{$b}->{column_name}} keys %files){
    my $column_name = $files{$build_id}{column_name};
    push(@column_names, $column_name);
  }

  my $column_names_s = join("\t", @column_names);
  $self->write_output(\%results, \%files, $column_names_s, $outfile);
  return 1;
}

1;

