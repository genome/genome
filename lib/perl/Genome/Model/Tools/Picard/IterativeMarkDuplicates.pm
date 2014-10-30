package Genome::Model::Tools::Picard::IterativeMarkDuplicates;

use strict;
use warnings;

use Genome;

my $DEFAULT_ASSUME_SORTED = 1;
my $DEFAULT_MAX_RECORDS_IN_RAM = 500000;
my $DEFAULT_STEP_SIZE = 2_000_000;

class Genome::Model::Tools::Picard::IterativeMarkDuplicates {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file => {
            is  => 'String',
            doc => 'The SAM/BAM files to merge.  File type is determined by suffix.',
        },
        output_directory => {
            is  => 'String',
            doc => 'A directory to write intermediate files to when run in distributed mode. In serial, a default file will be written here if output_file is not defined.',
            is_optional => 1,
        },
        output_file => {
            is  => 'String',
            doc => 'The output tsv file of duplication metrics.',
            is_optional => 1,
            is_output => 1,
        },
        # TODO: Resolve this with flagstat, picard alignment summary, or equivalent process
        initial_alignments => {
            is => 'Integer',
            doc => 'The initial number of alignments in the BAM file',
            is_optional => 1,
        },
        distributed => {
            is => 'Boolean',
            doc => 'Run as a distributed workflow.',
            default_value => 0,
            is_optional => 1,
        },
        step_size => {
            is => 'Integer',
            doc => 'The step size of alignments for each iteration',
            default_value => $DEFAULT_STEP_SIZE,
            is_optional => 1,
        },
        alignment_iterations => {
            doc => 'An array reference of alignment iterations to generate.  The default behaviour is to use step_size to determine iterations.',
            is_optional => 1,
        },
        random_seed => {
            is => 'Integer',
            doc => 'Random seed to use if reproducibilty is desired.',
            default_value => 1,
            is_optional => 1,
        },
        assume_sorted => {
            is  => 'Integer',
            valid_values => [1, 0],
            doc => 'Assume the input data is sorted.  default_value='. $DEFAULT_ASSUME_SORTED,
            default_value => $DEFAULT_ASSUME_SORTED,
            is_optional => 1,
        },
        max_sequences_for_disk_read_ends_map => {
            is => 'Integer',
            doc => 'The maximum number of sequences allowed in SAM file.  If this value is exceeded, the program will not spill to disk (used to avoid situation where there are not enough file handles',
            is_optional => 1,
        },
        max_records_in_ram => {
            doc => 'When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.',
            default_value => $DEFAULT_MAX_RECORDS_IN_RAM,
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Tool to mark or remove duplicate reads from a SAM/BAM file.';
}

sub help_detail {
    return <<EOS
    Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules.
    All records are then written to the output file with the duplicate records flagged.
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#MarkDuplicates
EOS
}



sub execute {
    my $self = shift;

    # This may not be ideal when secondary alignments are common
    my $flagstat_file = $self->input_file .'.flagstat';
    if (-e $flagstat_file) {
        my $flagstat_hash_ref = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_file);
        $self->initial_alignments($flagstat_hash_ref->{reads_mapped});
    } else {
        unless(defined($self->initial_alignments)) {
            die('Must define initial_alignments since flagstat file is not found : '. $flagstat_file);
        }
    }
    
    # TODO: Identify the libraries in the BAM file.  Use library names to look for existing mrkdup metrics file for input BAM
    # TODO: We may want to add a limitation to only run on single-library BAM files
    my $cmd = 'samtools view -H '. $self->input_file .' | grep \'@RG\'';
    my $text = `$cmd`;
    my @libraries;
    if ($text eq '') {
        @libraries = ('Unknown Library');
    }
    my @lines = split("\n", $text);
    for my $line (@lines) {
        chomp($line);
        unless ($line =~ /LB:(\S+)/) { die('Failed to match library for line : '. $line); }
        push @libraries, $1;
    }
    if (scalar(@libraries) > 1) {
        die('Please add support for multi-library BAM files!');
    }
    
    unless ($self->output_file) {
        unless ($self->output_directory) {
            die('Either output_file or output_directory must be defined!');
        }
        unless (-e $self->output_directory) {
            unless (Genome::Sys->create_directory($self->output_directory)) {
                die('Failed to create output_directory : '. $self->output_directory);
            }
        }
        $self->output_file($self->output_directory .'/IterativeMarkDuplicates_metrics.tsv');
    }
    if (-e $self->output_file) {
        die('Unable to proceed with existing file : '. $self->output_file);
    }

    my %metrics;
    unless ($self->distributed) {

        my @probabilities;
        # Loop to determine an array of alignment values with the last value being all alignments
        my @alignments;
        if ($self->alignment_iterations) {
            @alignments = @{$self->alignment_iterations};
        } else {
            for (my $i = $self->step_size; $i < $self->initial_alignments; $i += $self->step_size) {
                push @alignments, $i;
            }
            push @alignments, $self->initial_alignments;
        }

        # Loop to determine the probabilities required to randomly revert alignments with the desired step size
        # Since this will run in serial, we must consider the approximate number of alignments in each previous iteration
        my $j = 0;
        my $previous_alignments = $self->initial_alignments;
        for (my $k = scalar(@alignments); $k > 0; $k--) {
            my $goal_alignments = $alignments[$k-1];
            if (exists($alignments[$k+1])) {
                $previous_alignments = $alignments[$k];
            }
            my $ratio = $goal_alignments / $previous_alignments;
            push @probabilities, $ratio;
        }

        # Now that we have all the probabilites required, lets randomly revert if necessary and finally mark duplicates
        my $previous_bam;
        for my $probability (@probabilities) {
            if ($probability < 1) {
                # RandomRevertSam
                my $tmp_revert_bam = Genome::Sys->create_temp_file_path('RandomRevertSam-'. $probability .'.bam');
                unless (Genome::Model::Tools::Picard::RandomRevertSam->execute(
                    use_version => $self->use_version,
                    input_file => $previous_bam,
                    output_file => $tmp_revert_bam,
                    probability => $probability,
                    random_seed => $self->random_seed,
                    remove_alignment_information => 1,
                    remove_duplicate_information => 1,
                    maximum_memory => $self->maximum_memory,
                    maximum_permgen_memory => $self->maximum_permgen_memory,
                    #sort_order => 'coordinate',
                )->result ) {
                    die('Failed to randomly revert BAM : '. $previous_bam);
                }
                # Remove the previous temporary MarkDuplicate BAM file to save /tmp disk space
                # However, check to be sure the original is never removed (just in case...)
                if ($previous_bam ne $self->input_file) {
                    unlink($previous_bam) || die('Failed to remove BAM file : '. $previous_bam);
                }

                # TODO: Remove the sorting if RandomRevertSam can output a sorted BAM file by default
                # SortSam
                my $tmp_sorted_revert_bam = Genome::Sys->create_temp_file_path('RandomRevertSam-'. $probability .'-sorted.bam');
                unless (Genome::Model::Tools::Picard::SortSam->execute(
                    input_file => $tmp_revert_bam,
                    output_file => $tmp_sorted_revert_bam,
                    sort_order => 'coordinate',
                    maximum_memory => $self->maximum_memory,
                    maximum_permgen_memory => $self->maximum_permgen_memory,
                    max_records_in_ram => $self->max_records_in_ram,
                )->result) {
                    die('Failed to sort BAM : '. $tmp_sorted_revert_bam);
                }
                unlink($tmp_revert_bam) || die ('Failed to remove reverted BAM file : '. $tmp_revert_bam);

                # MarkDuplicates
                my $tmp_mrkdup_metrics = Genome::Sys->create_temp_file_path('MarkDuplicates-'. $probability .'-metrics.txt');
                my $tmp_mrkdup_bam = Genome::Sys->create_temp_file_path('MarkDuplicates-'. $probability .'.bam');
                unless (Genome::Model::Tools::Picard::MarkDuplicates->execute(
                    use_version => $self->use_version,
                    input_file => $tmp_sorted_revert_bam,
                    metrics_file => $tmp_mrkdup_metrics,
                    output_file => $tmp_mrkdup_bam,
                    maximum_memory => $self->maximum_memory,
                    maximum_permgen_memory => $self->maximum_permgen_memory,
                    max_sequences_for_disk_read_ends_map => $self->max_sequences_for_disk_read_ends_map,
                    max_records_in_ram => $self->max_records_in_ram,
                    assume_sorted => $self->assume_sorted,
                )->result ) {
                    die('Failed to mark duplicates of BAM : '. $tmp_revert_bam);
                }
                unlink($tmp_sorted_revert_bam) || die ('Failed to remove sorted and reverted BAM file : '. $tmp_sorted_revert_bam);

                # Parse the MarkDuplicates metrics
                my $metrics_hash_ref = Genome::Model::Tools::Picard::MarkDuplicates->parse_file_into_metrics_hashref($tmp_mrkdup_metrics);
                $metrics{$probability} = $metrics_hash_ref;
                $previous_bam = $tmp_mrkdup_bam;
            } else {
                # TODO: shortcut and get existing duplication metrics based on each library name
                my $metrics_hash_ref;
                for my $library (@libraries) {
                    if ($library eq 'Unknown Library') { next; }
                    my ($basename,$dirname,$suffix)  = File::Basename::fileparse($self->input_file,qw/\.bam/);
                    my $glob = $dirname .'/'. $library .'*.metrics';
                    my @metrics_files = glob($glob);
                    if (@metrics_files) {
                        if (scalar(@metrics_files) != 1) {
                            die('Found multiple MarkDuplicates metrics files for library glob : '. $glob);
                        }
                        $metrics_hash_ref = Genome::Model::Tools::Picard::MarkDuplicates->parse_file_into_metrics_hashref($metrics_files[0]);
                    }
                }

                # Well, no existing file was found so we must run MarkDuplicates
                if ($metrics_hash_ref) {
                    $previous_bam = $self->input_file;
                } else {
                    # TODO: Is reseting the duplicate marks necessary?  We have to rerun MarkDuplicates anyway
                    # RevertSam
                    # my $tmp_sorted_revert_bam = Genome::Sys->create_temp_file_path('RevertSam-'. $probability .'-sorted.bam');
                    # unless (Genome::Model::Tools::Picard::RevertSam->execute(
                    #    use_version => $self->use_version,
                    #    input_file => $self->input_file,
                    #    output_file => $tmp_sorted_revert_bam,
                    #    remove_duplicate_information => 1,
                    #    remove_alignment_information => 0,
                    #    restore_original_qualities => 0,
                    #    maximum_memory => $self->maximum_memory,
                    #    maximum_permgen_memory => $self->maximum_permgen_memory,
                    #    sort_order => 'coordinate',
                    # )->result ) {
                    #     die('Failed to revert duplicate marks in BAM : '. $self->input_file);
                    # }

                    # MarkDuplicates
                    my $tmp_mrkdup_metrics = Genome::Sys->create_temp_file_path('MarkDuplicates-'. $probability .'-metrics.txt');
                    my $tmp_mrkdup_bam = Genome::Sys->create_temp_file_path('MarkDuplicates-'. $probability .'.bam');
                    unless (Genome::Model::Tools::Picard::MarkDuplicates->execute(
                        use_version => $self->use_version,
                        input_file => $self->input_file,
                        metrics_file => $tmp_mrkdup_metrics,
                        output_file => $tmp_mrkdup_bam,
                        maximum_memory => $self->maximum_memory,
                        maximum_permgen_memory => $self->maximum_permgen_memory,
                        max_sequences_for_disk_read_ends_map => $self->max_sequences_for_disk_read_ends_map,
                        max_records_in_ram => $self->max_records_in_ram,
                        assume_sorted => $self->assume_sorted,
                    )->result ) {
                        die('Failed to mark duplicates of BAM : '. $self->input_file);
                    }
                    $metrics_hash_ref = Genome::Model::Tools::Picard::MarkDuplicates->parse_file_into_metrics_hashref($tmp_mrkdup_metrics);
                    $previous_bam = $tmp_mrkdup_bam;
                }
                $metrics{$probability} = $metrics_hash_ref;
            }
        }
        # Remove the previous temporary MarkDuplicate BAM file to save /tmp disk space
        # However, check to be sure the original is never removed (just in case...)
        if ($previous_bam ne $self->input_file) {
            unlink($previous_bam) || die('Failed to remove BAM file : '. $previous_bam);
        }
    } else {
        # TODO: Write a distributed workflow version
    }

    my @probabilities = sort {$b <=> $a} keys %metrics;
    #my @libraries = keys %{ $metrics{$probabilities[0]} };
    my @headers= keys %{ $metrics{$probabilities[0]}->{'LIBRARY-'.$libraries[0]} };
    
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->output_file,
        separator => "\t",
        headers => \@headers,
    );
    unless ($writer) {
        die('Failed to open file for writing : '. $self->output_file);
    }

    for my $probability (@probabilities) {
        for my $library (@libraries) {
            my $data = $metrics{$probability}->{'LIBRARY-'. $library};
            $writer->write_one($data);
        }
    }
    $writer->output->close;
    return 1;
}



1;
