package Genome::Model::Tools::Picard::BamIndexStats;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::BamIndexStats {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file   => {
            is  => 'String',
            doc => 'A BAM file to process.',
        },
        output_file => {
            is => 'String',
            doc => 'The output file to write stats.',
        },
    ],
};

sub help_brief {
    'Generates BAM index statistics. Input BAM file must have a corresponding index file. ';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#BamIndexStats
EOS
}

sub execute {
    my $self = shift;

    my $jar_path = $self->picard_path .'/BamIndexStats.jar';
    unless (-e $jar_path) {
        if ($self->use_version < 1.29) {
            die('Please define version 1.29 or greater.');
        } else {
            die('Missing jar path: '. $jar_path);
        }
    }
    my $cmd = $jar_path .' net.sf.picard.sam.BamIndexStats INPUT='. $self->input_file .' > '. $self->output_file;

    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => [$self->input_file],
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
