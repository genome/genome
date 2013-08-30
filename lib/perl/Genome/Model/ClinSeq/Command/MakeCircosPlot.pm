package Genome::Model::ClinSeq::Command::MakeCircosPlot;

#Written by Ben Ainscough and Scott Smith

use strict;
use warnings;
use Genome; 
use Data::Dumper;

class Genome::Model::ClinSeq::Command::MakeCircosPlot{
    is => 'Command::V2',
    has => [
        build                   => { is => 'Genome::Model::Build::ClinSeq',
                                     doc => 'Clinseq build' },

        output_directory        => { is => 'FilesystemPath',
                                     doc => 'Directory where output will be written', },

    ],
    doc => 'This script attempts to get read counts, frequencies and gene expression values for a series of genome positions',
};


sub sub_command_category { 'pipeline' }

#sub help_detail {
#  return <<EOS
#EOS
#}

#sub help_synopsis {
#  return <<EOS
#EOS
#}

#sub help_usage {
#  return <<EOS
#EOS
#}

sub execute {
    my $self = shift;
    $self->status_message("Hello world.  My build is " . $self->build . " and my output dir is " . $self->output_directory);
    die "Under development";
    return 1;
}

1;

