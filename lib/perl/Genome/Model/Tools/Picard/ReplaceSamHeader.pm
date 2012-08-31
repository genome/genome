package Genome::Model::Tools::Picard::ReplaceSamHeader;

use strict;
use warnings;

use Genome;
use IO::File;
use File::Basename;

class Genome::Model::Tools::Picard::ReplaceSamHeader {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file => {
            is  => 'String',
            doc => 'SAM file from which SAMRecords will be read.  File type is determined by suffix.',
        },
        output_file => {
            is  => 'String',
            doc => 'SAMFileHeader from HEADER file will be written to this file, followed by SAMRecords from INPUT file.  File type is determined by suffix.',
        },
        header_file => {
            is => 'String',
            doc => 'SAM/BAM file from which SAMFileHeader will be read.',
        },
    ],
};

sub help_brief {
    'Tool to replace the header of a BAM or SAM file using Picard';
}

sub help_detail {
    return <<EOS
    Tool to replace the header of a BAM or SAM file using Picard.  For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#SortSam
EOS
}

sub execute {
    my $self = shift;

    my $jar_path = $self->picard_path .'/ReplaceSamHeader.jar';
    unless (-e $jar_path) {
        die('Failed to find '. $jar_path .'!  This command may not be available in version '. $self->use_version);
    }
    my $input_file = $self->input_file;
    my $header_file = $self->header_file;
    my $output_file = $self->output_file;
    my $sort_cmd = $jar_path .' net.sf.picard.sam.ReplaceSamHeader OUTPUT='. $output_file .' INPUT='. $input_file .' HEADER='. $header_file;
    $self->run_java_vm(
        cmd => $sort_cmd,
        input_files => [$input_file,$header_file],
        output_files => [$output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
