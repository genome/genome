package Genome::Model::Tools::Picard::MarkDuplicates;

use strict;
use warnings;

use Genome;

my $DEFAULT_ASSUME_SORTED = 1;
my $DEFAULT_REMOVE_DUPLICATES = 0;

class Genome::Model::Tools::Picard::MarkDuplicates {
    is  => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        input_file => {
            is  => 'String',
            doc => 'The SAM/BAM files to merge.  File type is determined by suffix.',
            picard_param_name => 'INPUT',
        },
        output_file => {
            is  => 'String',
            doc => 'The resulting merged SAM/BAM file.  File type is determined by suffix.',
            picard_param_name => 'OUTPUT',
        },
        metrics_file => {
            is  => 'String',
            doc => 'File to write duplication metrics to.',
            picard_param_name => 'METRICS_FILE',
        },
        assume_sorted => {
            is  => 'Boolean',
            doc => 'Assume the input data is sorted.',
            default_value => $DEFAULT_ASSUME_SORTED,
            is_optional => 1,
            picard_param_name => 'ASSUME_SORTED',
        },
        remove_duplicates => {
            is => 'Boolean',
            doc => 'Merge the seqeunce dictionaries.',
            default_value => $DEFAULT_REMOVE_DUPLICATES,
            is_optional => 1,
            picard_param_name => 'REMOVE_DUPLICATES',
        },
        max_sequences_for_disk_read_ends_map => {
            is => 'Integer',
            doc => 'The maximum number of sequences allowed in SAM file.  If this value is exceeded, the program will not spill to disk (used to avoid situation where there are not enough file handles',
            is_optional => 1,
            picard_param_name => 'MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP ',
        },
        include_comment => {
            is => 'Text',
            doc => 'comment to include as a @CO in the BAM header',
            is_optional => 1,
        },
        read_name_regex => {
            is => 'Text',
            doc => "Set to 'null' to turn off optical duplicate detection",
            is_optional => 1,
            picard_param_name => 'READ_NAME_REGEX',
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
    http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
EOS
}
sub _jar_name {
    return 'MarkDuplicates.jar';
}

sub _java_class {
    my $self = shift;

    if ($self->version_at_least('1.122')) {
        return qw(picard sam markduplicates MarkDuplicates);
    }
    else {
        return qw(picard sam MarkDuplicates);
    }
}

sub _cmdline_args {
    my $self = shift;

    my @args = $self->SUPER::_cmdline_args;
    # MAX_FILE_HANDLES supported in v1.34+
    if ($self->version_at_least('1.34')) {
        # allow picard to use 95% of available file handles for caching reads
        push @args, sprintf("MAX_FILE_HANDLES=%s", $self->calculate_max_file_handles);
    }
    # COMMENT supported in v1.77
    if ($self->include_comment && $self->version_at_least('1.77')) {
        push @args, sprintf("COMMENT=%s", $self->include_comment);
    }

    return @args;
}

sub calculate_max_file_handles {
    return int(0.95 * `sh -ec "ulimit -n"`);
}

1;
