package Genome::Model::Tools::Picard::ReorderSam;

use strict;
use warnings;

use Genome;


class Genome::Model::Tools::Picard::ReorderSam {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file => {
            is  => 'String',
            doc => 'The SAM/BAM file to sort.  File type is determined by suffix.',
        },
        output_file => {
            is  => 'String',
            doc => 'The resulting sorted SAM/BAM file.  File type is determined by suffix.',
        },
        reference_file => {
            is => 'String',
            doc => 'The reference sequence to base order on',
        },
    ],
};

sub help_brief {
    'Tool to reorder a BAM or SAM file using Picard';
}

sub help_detail {
    return <<EOS
    Tool to sort a BAM or SAM file using Picard.  For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#ReorderSam
EOS
}

sub execute {
    my $self = shift;

    my $jar_path = $self->picard_path .'/ReorderSam.jar';
    unless (-e $jar_path) {
        die('Failed to find '. $jar_path .'!  This command may not be available in version '. $self->use_version);
    }
    my $input_file = $self->input_file;
    my $output_file = $self->output_file;
    my $reference_file = $self->reference_file;
    my $sort_cmd = $jar_path .' net.sf.picard.sam.ReorderSam O='. $output_file .' I='. $input_file .' REFERENCE='. $reference_file;
    $self->run_java_vm(
        cmd => $sort_cmd,
        input_files => [$input_file,$reference_file],
        output_files => [$output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
