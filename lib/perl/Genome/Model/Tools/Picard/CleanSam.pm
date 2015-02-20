package Genome::Model::Tools::Picard::CleanSam;

use strict;
use warnings;

use Genome;
use File::Basename;
use Carp qw(confess);


class Genome::Model::Tools::Picard::CleanSam {
    is  => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        input_file  => {
            is  => 'String',
            doc => 'The SAM files to run on',
            picard_param_name => 'INPUT',
        },
        output_file => {
            is  => 'String',
            doc => 'The output clean sam file',
            picard_param_name => 'OUTPUT',
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

sub _jar_name {
    return 'CleanSam.jar';
}

sub _java_class_name {
    return 'net.sf.picard.sam.CleanSam';
}

sub _validate_params {
    my $self = shift;

    unless ($self->input_file and -s $self->input_file) {
        confess $self->error_message("Input file is invalid");
        return;
    }
}

sub _shellcmd_extra_params {
    my $self = shift;
    return (
        input_files  => [$self->input_file],
        #output_files => [$self->output_file],
        skip_if_output_is_present => 0,
        );
}

1;
