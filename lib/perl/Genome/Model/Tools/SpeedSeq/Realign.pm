package Genome::Model::Tools::SpeedSeq::Realign;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SpeedSeq::Realign {
    is => 'Genome::Model::Tools::SpeedSeq::AlignBase',
    has_input => [
        bams => {
            is => 'Text',
            doc => 'BAM file(s) (must contain read group tags)',
            is_many => 1,
            example_values => ['in.bam'],
            tool_bare_arg_position => '2',
        },
    ],
    has_param => [
        rename_reads => {
            is => 'Boolean',
            doc => 'rename reads for smaller file size',
            is_optional => 1,
            tool_param_name => 'n',
        },
    ],
};

sub _tool_subcommand_name {
    return 'realign';
}

#sub execute {
#    my $self = shift;

    # OUTPUTS
    #unless ($self->output_prefix) {
    #    $self->output_prefix('in.realign');
    #}
#}

1;
