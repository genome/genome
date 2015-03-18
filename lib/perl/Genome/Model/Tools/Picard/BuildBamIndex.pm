package Genome::Model::Tools::Picard::BuildBamIndex;

use strict;
use warnings;

use Genome;
use File::Spec qw();

class Genome::Model::Tools::Picard::BuildBamIndex {
    is  => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        input_file   => {
            is  => 'String',
            doc => 'A BAM file to process.',
            picard_param_name => 'INPUT',
        },
        output_file => {
            is => 'String',
            doc => 'The BAM index file. Defaults to x.bai if INPUT is x.bam, otherwise INPUT.bai.  If INPUT is a URL and OUTPUT is unspecified, defaults to a file in the current directory.',
            is_optional => 1,
            picard_param_name => 'OUTPUT',
        },
    ],
};

sub help_brief {
    'Generates a BAM index (.bai) file.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://broadinstitute.github.io/picard/command-line-overview.html#BuildBamIndex
EOS
}

sub _jar_name {
    return 'BuildBamIndex.jar';
}

sub _java_class {
    return qw(picard sam BuildBamIndex);
}

sub _validate_params {
    my $self = shift;
    $self->enforce_minimum_version("1.23");
}

sub _shellcmd_extra_params {
    my $self = shift;
    return (
        input_files  => [$self->input_file],
        skip_if_output_is_present => 0,
        );
}

1;
