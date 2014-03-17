package Genome::Model::Command::InstrumentData::Assign::AllCompatible;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::InstrumentData::Assign::AllCompatible {
    is => 'Genome::Model::Command::InstrumentData::Assign::Base',
    has_optional_input => [
        include_imported => {
            is => 'Boolean',
            default_value => 0,
            doc => 'include matching imported data',
        },
        maximum_allowed_error => {
            is => 'Number',
            example_values => ['3.0'],
            doc => "If set, the maximum allowed gerald error rate to assign to a model",
        },
    ],
};

sub help_brief {
    return "Assign all compatible instrument data to a model";
}

sub help_detail {
    return <<'EOHELP'
This command asks the model what data is compatible with it and tries to assign all of it.

This should be the same list of instrument data as produced by the command:
`genome model instrument-data list --compatible`.
EOHELP
;
}

sub _resolve_instrument_data {
    my $self = shift;

    my $model = $self->model;
    unless($model->can('compatible_instrument_data')) {
        die $self->error_message('Models of type %s do not currently support finding compatible data.');
    }

    return $model->compatible_instrument_data;
}


sub _assign_instrument_data_if_acceptable {
    my $self = shift;
    my $instrument_data = shift;
    my $filter = shift;

    my @issues = @_;

    if(not $self->include_imported and $instrument_data->isa("Genome::InstrumentData::Imported")) {
        push @issues, 'Data is imported and --include-imported not set.';
    }

    if($self->maximum_allowed_error and not $self->_instrument_data_within_allowed_error($instrument_data, $filter)) {
        push @issues, 'Error rate exceeds maximum allowed.';
    }

    return $self->SUPER::_assign_instrument_data_if_acceptable($instrument_data, $filter, @issues);
}

sub _instrument_data_within_allowed_error {
    my $self = shift;
    my $instrument_data = shift;

    my ($fwd_error, $rev_error);

    if ($instrument_data->fwd_filt_error_rate_avg) {
        $fwd_error = $instrument_data->fwd_filt_error_rate_avg;
    }
    else {
        $fwd_error = $instrument_data->get_default_alignment_metrics('read_1_pct_mismatch');
    }

    if ($instrument_data->filt_error_rate_avg) {
        $rev_error = $instrument_data->filt_error_rate_avg; #this is intentionally not rev_filt_error_rate_avg
    }
    else {
        $rev_error = $instrument_data->get_default_alignment_metrics('read_2_pct_mismatch');
    }

    if($fwd_error > $self->maximum_allowed_error or $rev_error > $self->maximum_allowed_error) {
        return;
    }

    return 1;
}

1;
