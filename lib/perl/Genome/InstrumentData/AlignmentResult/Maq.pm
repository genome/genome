package Genome::InstrumentData::AlignmentResult::Maq;

use strict;
use warnings;

use File::Basename;
use IO::File;
use Genome;

class Genome::InstrumentData::AlignmentResult::Maq {
    is  => ['Genome::InstrumentData::AlignmentResult'],
    has_constant => [
        aligner_name => { value => 'maq', is_param => 1 },
    ],
    has => [
        is_not_run_as_paired_end => {
            type => 'Boolean',
            calculate_from => ['filter_name', 'force_fragment'],
            calculate => q{return $force_fragment
                || ($filter_name && (($filter_name eq 'forward-only')
                || ($filter_name eq 'reverse-only'))); 
            },
        },
    ]
};

sub required_arch_os { 
    'x86_64' 
}

sub required_rusage { 
    "-R 'span[hosts=1] rusage[tmp=50000:mem=12000]' -M 1610612736";
}

sub extra_metrics {
    'contaminated_read_count'
}

#####ALIGNER OUTPUT#####
#the fully quallified file path for aligner output
sub aligner_output_file_path {
    my $self = shift;
    my $lane = $self->instrument_data->subset_name;
    #return $self->temp_scratch_directory . "/alignments_lane_${lane}.map.aligner_output";
    return $self->temp_staging_directory . '/aligner.log';
}


#####UNALIGNED READS LIST#####
#the fully quallified file path for unaligned reads
sub unaligned_reads_list_path {
    my $self        = shift;
    my $subset_name = $self->instrument_data->subset_name;
    return $self->temp_scratch_directory . "/s_${subset_name}_sequence.unaligned";
}

# return list of generated map files for alignments that have completed and been synced to network disk
sub alignment_file_paths {
    my $self = shift;
    my $dir = $self->output_dir;

    return glob($dir . "/*.map");
}

sub _run_aligner {
    die "Maq was deleted - it is obsolete.";
}

sub aligner_params_for_sam_header {
    return undef;
}

sub fillmd_for_sam {
    return 1;
}

sub prepare_reference_sequence_index {
    die "Maq was deleted - it is obsolete.";
}

1;
