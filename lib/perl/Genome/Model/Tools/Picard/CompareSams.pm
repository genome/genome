package Genome::Model::Tools::Picard::CompareSams;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::CompareSams {
    is  => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        # This picard param uses positional args, so no picard_param_name
        input_file_1 => {
            is  => 'String',
            doc => 'The first SAM file path to compare.',
        },
        input_file_2 => {
            is  => 'String',
            doc => 'The second SAM file path to compare.',
        },
        output_file => {
            is => 'String',
            doc => 'The path to the output file.',
        },
    ],
};

sub help_brief {
    'Tool to compare a BAM or SAM file using Picard';
}

sub help_detail {
    return <<EOS
    Tool to compare a BAM or SAM file using Picard.  For Picard documentation of this command see:
    http://broadinstitute.github.io/picard/command-line-overview.html#CompareSAMs
EOS
}

sub _jar_name {
    return 'CompareSAMs.jar';
}

sub _java_class {
    my $self = shift;
    return qw(samtools apps CompareSAMs) if $self->use_version eq '1.17';
    return qw(picard sam CompareSAMs);
}

sub _validate_params {
    my $self = shift;
    $self->enforce_minimum_version('1.17');
}

sub _redirects {
    my $self = shift;
    return sprintf('> %s', $self->output_file);
}

sub _cmdline_args {
    my $self = shift;
    return (
        $self->input_file_1,
        $self->input_file_2,
        $self->SUPER::_cmdline_args
        );
}

sub _shellcmd_extra_params {
    my $self = shift;
    return (
        input_files => [$self->input_file_1, $self->input_file_2],
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
        );
}

1;
