package Genome::Model::Command::Define::Base;

use strict;
use warnings;

use Genome;

use Regexp::Common;

class Genome::Model::Command::Define::Base {
    is => 'Genome::Command::Base',
    is_abstract => 1,
    has => [
        name => {
            is => 'Text',
            doc => 'Model name',
        },
        subject => {
            is => 'Genome::Subject',
            doc => 'Model subject. Taxon, individual, sample, etc.',
        },
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            doc => 'Processing profileby id or name.',
        },
        projects => {
            is => 'Genome::Project',
            is_many => 1,
            is_optional => 1,
        },
        params => {
            is_many => 1,
            is_optional => 1,
            shell_args_position => 1,
            doc => 'Params for the model. Use the format: property=value. The value can be an id or filter string. If the property can have many values, use a filter string, or multiple key=value pairs.',
        },
        _model => { is_optional => 1, },
        _model_class_meta => { calculate_from => [qw/ _model_class /], calculate => q| return $_model_class->__meta__; |, },
    ],
};

sub inputs_for_model_class {
    return sort { $a->property_name cmp $b->property_name } 
    grep {
        $_->property_name ne 'instrument_data' and $_->property_name ne 'input_values'
    } grep { 
        $_->via and $_->via eq 'inputs'
    } $_[0]->__meta__->property_meta_for_name('_model_class')->default_value->__meta__->property_metas;
}

sub help_synopsis {
    return;
}

sub help_detail {
    my $class = shift;
    my $model_class = $class->__meta__->property_meta_for_name('_model_class')->default_value;
    my $help = <<HELP;
This will define a model, setting inputs, instrument data as given. List params as space separated bare args. Use the format: property=value. The value can be an id or filter string.  If the property can have many values, use a filter string, or multiple key=value pairs.

Params for all model types
 * auto_assign_inst_data
 * auto_build_alignments
 * instrument_data
HELP
    $help .= "\nSpecific params for $model_class\n";
    for my $input ( $class->inputs_for_model_class ) {
        $help .= ' * '.$input->property_name.' ('.($input->is_many ? 'many' : 'singular').")\n";
    }
    return $help;
}

sub execute {
    my $self = shift;

    $self->status_message('Define model...');

    my $params = $self->_resolve_params;
    return if not $params;
    #print Data::Dumper::Dumper($params);
    my $model = $self->_model_class->create(%$params);
    if ( not $model ) {
        $self->error_message('Failed to create '.$self->_model_class.' model');
        return;
    }
    $self->_model($model);

    $self->status_message('Successfully defined model: '.$model->__display_name__);

    return 1;
}

sub _resolve_params {
    my $self = shift;

    my %params = (
        name => $self->name,
        processing_profile => $self->processing_profile,
        subject_id => $self->subject->id,
        subject_class_name => $self->subject->class,
        projects => [ $self->projects ],
    );

    my $model_meta = $self->_model_class_meta;
    for my $param ( $self->params ) {
        my ($key, $value) = split('=', $param, 2);
        my $property = $model_meta->property_meta_for_name($key);
        if ( not $property ) {
            $self->error_message("Failed to find model property: $key");
            return;
        }

        if ( my ($unallowed) = grep { $property->$_ } (qw/ is_calculated is_constant is_transient /) ){
            $self->error_message("Property ($key) cannot be given on the command line because it is '$unallowed'");
            return;
        }

        my @values = $value;
        my $data_type = $property->data_type;
        if ( not grep { $data_type =~ /^$_$/i } (qw/ boolean integer number string text ur::value /) ) { # hacky...if u kno a better way...
            my $filter = ( $value =~ /^$RE{num}{int}$/ ) ? 'id='.$value : $value;
            my $data_type = $property->data_type;
            my $bx = eval { UR::BoolExpr->resolve_for_string($data_type, $filter); };
            if ( not $bx ) {
                $self->error_message("Failed to create expression for $key ($data_type) from '$value'");
                return;
            }
            @values = $data_type->get($bx);
            if ( not @values ) {
                $self->error_message("Failed to get $key ($data_type) for $value");
                return 1;
            }
        }

        if ( $property->is_many ) {
            push @{$params{$key}}, @values;
        }
        elsif ( @values > 1 or exists $params{$key} ) {
            $self->error_message(
                "Singular property ($key) cannot have more than one value (".join(', ', grep { defined } (@values, $params{$key})).')'
            );
            return;
        }
        else {
            $params{$key} = $values[0];
        }
    }

    return \%params;
}

1;

