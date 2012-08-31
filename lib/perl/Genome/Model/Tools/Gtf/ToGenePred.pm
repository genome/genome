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
        extended => {
            is => 'Boolean',
            doc => 'Output the extended format which includes strandedness.',
            default_value => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    my $cmd = '/gscmnt/sata132/techd/solexa/jwalker/bin/gtfToGenePred';
    if ($self->extended) {
        $cmd .= ' -genePredExt';
    }
    $cmd .= ' '. $self->input_gtf_file .' '. $self->output_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_gtf_file],
        output_files => [$self->output_file],
    );
    return 1;
};
