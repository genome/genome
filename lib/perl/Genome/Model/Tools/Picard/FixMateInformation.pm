package Genome::Model::Tools::Picard::FixMateInformation;

use strict;
use warnings;

use Genome;
use File::Basename;


class Genome::Model::Tools::Picard::FixMateInformation {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file  => {
            is  => 'String',
            doc => 'The SAM/BAM files to run on.  File type is determined by suffix.',
        },
        output_file => {
            is  => 'String',
            doc => 'The output file',
        },
        sort_order  => {
            is  => 'String',
            doc => 'Optional sort order if the out file should be sorted differently than the INPUT file',
            valid_values => [qw(unsorted queryname coordinate)],
            is_optional  => 1,
        },
    ],
};

sub help_brief {
    'Tool to ensure that all mate-pair information is in sync between each read and its mate pair.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#CollectGcBiasMetrics
EOS
}

sub execute {
    my $self = shift;

    unless ($self->input_file and -s $self->input_file) {
        $self->error_message("Input file is invalid");
        return;
        
    }

    my $cmd = $self->picard_path .'/FixMateInformation.jar net.sf.picard.sam.FixMateInformation';
    $cmd   .= ' OUTPUT='. $self->output_file  .' INPUT='. $self->input_file;
        
    if ($self->sort_order) {
        $cmd .= ' SORT_ORDER=' . $self->sort_order;
    }
    
    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => [$self->input_file],
        #output_files => [$self->output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
