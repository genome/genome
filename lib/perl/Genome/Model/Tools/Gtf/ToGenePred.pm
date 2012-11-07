package Genome::Model::Tools::Gtf::ToGenePred;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Gtf::ToGenePred {
    is => ['Genome::Model::Tools::Gtf::Base'],
    has => [
        output_file => {
            is => 'Text',
            doc => 'The output genePred format file.',
        },
        use_version => {
            is => 'Text',
            doc => 'The version of the chimerascan software.',
            is_optional => 1,
            default_value => '0.4.5',
        },
    ],
};

sub execute {
    my $self = shift;

    # Is there a better way to find the python script?
    my $cmd = 'gtf_to_genepred.py'. $self->use_version;
    $cmd .= ' '. $self->input_gtf_file .' '. $self->output_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_gtf_file],
        output_files => [$self->output_file],
    );
    return 1;
};
