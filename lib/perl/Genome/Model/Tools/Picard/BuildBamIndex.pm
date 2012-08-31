package Genome::Model::Tools::Picard::BuildBamIndex;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::BuildBamIndex {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file   => {
            is  => 'String',
            doc => 'A BAM file to process.',
        },
        output_file => {
            is => 'String',
            doc => 'The BAM index file. Defaults to x.bai if INPUT is x.bam, otherwise INPUT.bai.  If INPUT is a URL and OUTPUT is unspecified, defaults to a file in the current directory.',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Generates a BAM index (.bai) file.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#BuildBamIndex
EOS
}

sub execute {
    my $self = shift;

    my $jar_path = $self->picard_path .'/BuildBamIndex.jar';
    unless (-e $jar_path) {
        if ($self->use_version < 1.23) {
            die('Please define version 1.23 or greater.');
        } else {
            die('Missing jar path: '. $jar_path);
        }
    }
    my $cmd = $jar_path .' net.sf.picard.sam.BuildBamIndex INPUT='. $self->input_file;
    if ($self->output_file) {
        $cmd .= ' OUTPUT='. $self->output_file;
    }

    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => [$self->input_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
