# FIXME ebelter
#  Long: Remove or update to use inputs as appropriate.
#
package Genome::Model::Command::InstrumentData::List;

use strict;
use warnings;

use Genome;
use Regexp::Common;

class Genome::Model::Command::InstrumentData::List {
    is => 'UR::Object::Command::List',
    has => [
        model => { is => 'Genome::Model', id_by => 'model_id' },
        model_id => { is => 'Integer', doc => 'ID of the genome model' },
        subject_class_name => {
            is_constant => 1,
            value => 'Genome::InstrumentData',
        },
    ],
    has_optional => [
        assigned => {
            is => 'Boolean',
            default => 0,
            doc => 'List instrument data that has been assigned to the model (default)',
        },
        unassigned => {
            is => 'Boolean',
            default => 0,
            doc => 'List instrument data that is compatable to the model, but not yet assigned',
        },
        compatible => {
            is => 'Boolean',
            default => 0,
            doc => 'List instrument data that has been is compatable to the model',
        },
        show => {
            is => 'Text',
            default => 'id,full_name,full_path,library_name,sample_name,sequencing_platform,subclass_name,target_region_set_name',
            doc => 'Columns of the instrument data to show in the list'
        },
    ],
};

#########################################################

sub help_brief {
    return "List a model's instrument data";
}

sub help_detail {
    return '';
}

#########################################################

sub assignment_types {
    return (qw/ assigned unassigned compatible /);
}

sub _resolve_boolexpr {
    my $self = shift;

    # Verify model
    unless ( $self->model_id ) {
        $self->error_message("No model id given");
        return;
    }

    unless ( $self->model_id =~ /^$RE{num}{int}$/ ) {
        $self->error_message( sprintf('Model id given (%s) is not an integer', $self->model_id) );
        return;
    }

    unless ( $self->model ) {
        $self->error_message( sprintf('Can\'t get model for id (%s) ', $self->model_id) );
        return;
    }

    # Determine what type of inst data to get
    my @types = grep { $self->$_ } assignment_types();
    unless ( @types ) {
        $self->assigned(1);
        @types = (qw/ assigned /);
    }
    elsif ( @types > 1 ) {
        $self->error_message(
            sprintf(
                'Multiple instrument data assignment types requested (%s).  Please select only one: %s',
                join(', ', @types),
                join(', ', assignment_types()),
            )
        );
        return;
    }

    my $method = sprintf('%s_instrument_data', $types[0]);

    return Genome::InstrumentData->define_boolexpr(
        id => [ map { $_->id } $self->model->$method ] ,
    );
}

1;

