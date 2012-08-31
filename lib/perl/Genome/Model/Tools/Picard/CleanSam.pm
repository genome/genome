package Genome::Model::Tools::Picard::CleanSam;

use strict;
use warnings;

use Genome;
use File::Basename;


class Genome::Model::Tools::Picard::CleanSam {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file  => {
            is  => 'String',
            doc => 'The SAM files to run on',
        },
        output_file => {
            is  => 'String',
            doc => 'The output clean sam file',
        },
    ],
};

sub help_brief {
    'Tool to Read SAM and perform various fix-ups. Currently, the only fix-up it to soft-clip an alignment that hangs off the end of its reference sequence.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#CleanSam
EOS
}

sub execute {
    my $self = shift;

    unless ($self->input_file and -s $self->input_file) {
        $self->error_message("Input file is invalid");
        return;
        
    }

    my $cmd = $self->picard_path .'/CleanSam.jar net.sf.picard.sam.CleanSam';
    $cmd   .= ' OUTPUT='. $self->output_file  .' INPUT='. $self->input_file;
        
    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => [$self->input_file],
        #output_files => [$self->output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
