package Genome::Model::ClinSeq::Command::Converge::CancerRelevanceScore;
use strict;
use warnings;
use Genome;
use Getopt::Long;
use File::Basename;

class Genome::Model::ClinSeq::Command::Converge::CancerRelevanceScore {
    is        => 'Genome::Model::ClinSeq::Command::Converge::Base',
    has_input => [
        outfile => {
            is  => 'FilesystemPath',
            doc => 'File to write converged cancer-relevance-score results',
        },
        gene_category_sum => {
            is          => 'Integer',
            is_optional => 1,
            default     => 5,
            doc         => 'Cutoff for gene_category_sum value',
        },
    ],
    doc => 'converge Stats from mutiple clinseq builds'
};

sub help_detail {
    return <<EOS
  For a group of ClinSeq models, converge cancer-relevance-scores results.
  All stats should be pre-calculated. Produce an output matrix in which
  each row is metric and each column is a sample.
  Input:
    A list of Clinseq builds, models [OR] a Clinseq model-group.
EOS
}

sub help_synopsis {
    return <<EOS
  Example usage: 
  genome model clin-seq converge cancer-relevance-score \\
  --builds='model_groups.id=786367aa2edc41e1b4a5d33787a8c003,is_last_complete=1' \\
  --outfile=metris.tsv --outdir=/tmp/ --gene-category-sum=0
EOS
}

sub resolve_which_stats_tsv {
    my $self = shift;
    my $b    = shift;
    my @stats_files;
    push @stats_files, $b->wgs_snv_tier1_annotated_compact_catanno_file;
    push @stats_files, $b->wgs_indel_tier1_annotated_compact_catanno_file;
    push @stats_files, $b->exome_indel_tier1_annotated_compact_catanno_file;
    push @stats_files, $b->exome_snv_tier1_annotated_compact_catanno_file;
    return @stats_files;
}

sub find_stats_tsv {
    my $self  = shift;
    my $mb    = shift;
    my $files = shift;
    foreach my $c (keys %$mb) {
        my $b              = $mb->{$c}->{build};
        my $m              = $mb->{$c}->{model};
        my $model_name     = $m->name;
        my $data_directory = $b->data_directory;
        my $subject_name   = $m->exome_model->tumor_model->subject->common_name;
        $subject_name =~ s/[\s,-]/_/g;
        my $subject_common_name = $b->subject->common_name;
        my $build_id            = $b->id;

        #If the subject name is not defined, die
        unless ($subject_name) {
            die $self->error_message("\n\nCould not determine subject name for build: $build_id\n\n");
        }
        my $final_name = "Unknown";
        if ($subject_name)        {$final_name = $subject_name;}
        if ($subject_common_name) {$final_name = $subject_common_name;}
        $self->status_message("\n\t$final_name\t$build_id\t$data_directory");
        my %file_info;
        my @stats_files = $self->resolve_which_stats_tsv($b);
        foreach my $stats_file (@stats_files) {

            if ($stats_file) {
                if ($stats_file =~ /^$data_directory\/.*?\/(.*)/) {
                    my $distinct_name = $1;
                    $file_info{$stats_file}{distinct_name} = $distinct_name;
                }
                else {
                    die $self->error_message("\n\nCould not determine distinct name " . "from path\n$stats_file\n\n");
                }
                $files->{$build_id}{stats_files}  = \%file_info;
                $files->{$build_id}{final_name}   = $final_name;
                $files->{$build_id}{subject_name} = $subject_name;
                $files->{$build_id}{model_name}   = $model_name;
            }
            else {
                $self->warning_message("\n\tWarning. Could not find any mutation-spectrum files for build: $build_id");
            }
        }
    }
}

sub get_column_names {
    my $self       = shift;
    my $files      = shift;
    my $file_count = shift;
    my $id_choice;
    my %finalnames;
    my %finalnames_subjectnames;
    my %modelnames;
    my %finalnames_subjectnames_builds;

    foreach my $bid (keys %$files) {
        my $final_name   = $files->{$bid}{final_name};
        my $subject_name = $files->{$bid}{subject_name};
        my $model_name   = $files->{$bid}{model_name};
        my $id1          = $final_name;
        my $id2          = "$final_name" . "_" . "$subject_name";
        my $id3          = $final_name . "_" . $model_name;
        my $id4          = "$final_name" . "_" . "$subject_name" . "_" . "$bid";
        $finalnames{$id1}                     = 1;
        $finalnames_subjectnames{$id2}        = 1;
        $modelnames{$id3}                     = 1;
        $finalnames_subjectnames_builds{$id4} = 1;
        $files->{$bid}{id1}                   = $id1;
        $files->{$bid}{id2}                   = $id2;
        $files->{$bid}{id3}                   = $id3;
        $files->{$bid}{id4}                   = $id4;
    }
    if ($file_count <= keys %finalnames) {
        $id_choice = "id1";
    }
    elsif ($file_count <= keys %finalnames_subjectnames) {
        $id_choice = "id2";
    }
    elsif ($file_count <= keys %modelnames) {
        $id_choice = "id3";
    }
    elsif ($file_count <= keys %finalnames_subjectnames_builds) {
        $id_choice = "id4";
    }
    else {
        die $self->error_message(
            "\n\nCould not generate a suitable set of distinct labels for the data files to be joined...\n\n");
    }
    return $id_choice;
}

#######################################################################################################################
#Get the desired files from each model                                                                                #
#######################################################################################################################
sub aggregate_stats {
    my $self          = shift;
    my %args          = @_;
    my $models_builds = $args{'-models_builds'};
    my %files;
    my $fc = 0;
    my %mb = %{$models_builds->{cases}};
    $self->find_stats_tsv(\%mb, \%files);
    my $file_count = keys %files;

    unless ($file_count > 1) {
        die $self->error_message("\n\nFound $file_count files to join (need at"
                . " least 2 ...)\n\nLooked for files of name: Stats.tsv\n\n");
    }
    my $id_choice = $self->get_column_names(\%files, $file_count);
    foreach my $fc (keys %files) {
        $files{$fc}{column_name} = $files{$fc}{"$id_choice"};
    }
    return (\%files);
}

##########################################################################
#Build a hash of metrics across the input builds and their stats files   #
##########################################################################
sub parse_metrics {
    my $self                     = shift;
    my %args                     = @_;
    my %files                    = %{$args{'-files'}};
    my $results                  = $args{'-results'};
    my $gene_category_sum_cutoff = $self->gene_category_sum;
    foreach my $build_id (sort keys %files) {
        my %stats_files = %{$files{$build_id}{stats_files}};
        foreach my $file (sort keys %stats_files) {
            my $header = 1;
            my @genes;
            my $mutation_count = 0;
            open(my $STATS, "$file")
                || die $self->error_message("\n\nCould not open stats file: $file");
            while (<$STATS>) {
                chomp;
                my @line = split("\t", $_);
                my $ncol = scalar @line;
                if ($header) {
                    if ($line[$ncol - 1] ne "gene_category_sum" or $line[1] ne "gene_name") {
                        die $self->error_message("Last column of $file is not gene_category_sum or "
                                . " second column is not gene_name. ("
                                . $line[$ncol - 1] . ","
                                . $line[1] . ")"
                                . " Please verify file is correct !");
                    }
                    $header = 0;
                    next;
                }
                if ($line[$ncol - 1] >= $gene_category_sum_cutoff) {
                    $mutation_count += 1;
                    push @genes, $line[2] . "(" . $line[4] . ")";
                }
            }
            $results->{$build_id}{$file}{mutation_count} = $mutation_count;
            $results->{$build_id}{$file}{genes} = join(":", @genes);
            my ($data_type, $mutation_type) = $self->_get_data_mutation_type($file);
            $results->{$build_id}{$file}{data_type}     = $data_type;
            $results->{$build_id}{$file}{mutation_type} = $mutation_type;
            close($STATS);
        }
    }
}

sub write_output {
    my $self                     = shift;
    my $results                  = shift;
    my $files                    = shift;
    my $column_names_s           = shift;
    my $outfile                  = shift;
    my $gene_category_sum_cutoff = $self->gene_category_sum;
    my $header_line              = "Mutation_Type\tStatistic_Type\tData_Type\t$column_names_s";
    #No go through each build and print out the summary
    open(OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
    print OUT "#gene_category_sum_cutoff = $gene_category_sum_cutoff\n";
    print OUT "$header_line\n";
    my ($data_type, $mutation_type);
    my ($counts, $genes, $data_types, $mutation_types);
    foreach my $build_id (sort {$files->{$a}->{column_name} cmp $files->{$b}->{column_name}} keys %$files) {
        my %stats_files = %{$files->{$build_id}{stats_files}};
        foreach my $file (sort keys %stats_files) {
            $data_type     = $results->{$build_id}{$file}{data_type};
            $mutation_type = $results->{$build_id}{$file}{mutation_type};
            my $mutation_count = $results->{$build_id}{$file}{mutation_count};
            my $genes1         = $results->{$build_id}{$file}{genes};
            push @{$counts->{$mutation_type . $data_type}}, $mutation_count;
            push @{$genes->{$mutation_type . $data_type}},  $genes1;
            $mutation_types->{$mutation_type . $data_type} = $mutation_type;
            $data_types->{$mutation_type . $data_type}     = $data_type;
        }
    }
    foreach my $variant (keys %$data_types) {
        my $data_type     = $data_types->{$variant};
        my $mutation_type = $mutation_types->{$variant};
        print OUT "$mutation_type\tcount\t$data_type\t" . join("\t", @{$counts->{$variant}}) . "\n";
        print OUT "$mutation_type\tgenes\t$data_type\t" . join("\t", @{$genes->{$variant}}) . "\n";
    }
    close(OUT);
    $self->status_message("Wrote output to $outfile");
}

sub _get_data_mutation_type {
    my $self          = shift;
    my $stats_file    = shift;
    my $data_type     = "WGS/Exome";
    my $mutation_type = "snp/indel";
    if ($stats_file =~ /indel\/wgs/) {
        $data_type     = "wgs";
        $mutation_type = "indel";
    }
    elsif ($stats_file =~ /indel\/exome/) {
        $data_type     = "exome";
        $mutation_type = "indel";
    }
    elsif ($stats_file =~ /snv\/wgs/) {
        $data_type     = "wgs";
        $mutation_type = "snp";
    }
    elsif ($stats_file =~ /snv\/exome/) {
        $data_type     = "exome";
        $mutation_type = "snp";
    }
    return ($data_type, $mutation_type);
}

sub execute {
    my $self    = shift;
    my $outfile = $self->outdir . "/" . basename($self->outfile);
    my @builds  = $self->builds;
    my @build_ids;

    foreach my $build (@builds) {
        push @build_ids, $build->id;
    }
    my $models_builds = $self->getModelsBuilds('-builds' => \@build_ids);
    my %files = %{$self->aggregate_stats('-models_builds' => $models_builds)};
    my %results;
    $self->parse_metrics(
        '-files'   => \%files,
        '-results' => \%results
    );

    #Build the header line
    my @column_names;
    foreach my $build_id (sort {$files{$a}->{column_name} cmp $files{$b}->{column_name}} keys %files) {
        my $column_name = $files{$build_id}{column_name};
        push(@column_names, $column_name);
    }

    my $column_names_s = join("\t", @column_names);
    $self->write_output(\%results, \%files, $column_names_s, $outfile);
    return 1;
}

1;

