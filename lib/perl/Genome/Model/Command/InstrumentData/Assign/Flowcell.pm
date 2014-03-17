package Genome::Model::Command::InstrumentData::Assign::Flowcell;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::InstrumentData::Assign::Flowcell {
    is => 'Genome::Model::Command::InstrumentData::Assign::Base',
    has_input => [
        flow_cell_id => {
            is => 'Text',
            doc => 'solexa/illumina flow cell id of instrument data to assign',
            shell_args_position => 1,
        },
    ],
};

sub help_brief {
    return "Assign instrument data to a model by flowcell ID";
}

sub help_detail {
    return <<'EOHELP'
This command assigns solexa instrument data to the model based on a flowcell.
EOHELP
;
}

sub _resolve_instrument_data {
    my $self = shift;

    return Genome::InstrumentData::Solexa->get(flow_cell_id => $self->flow_cell_id);
}

1;
