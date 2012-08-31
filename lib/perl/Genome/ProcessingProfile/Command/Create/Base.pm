package Genome::ProcessingProfile::Command::Create::Base;

use strict;
use warnings;

use Genome;
use Carp 'confess';
use Regexp::Common;

class Genome::ProcessingProfile::Command::Create::Base {
    is => 'Command::V2',
    has => [
        name => {
            is => 'Text',
            len => 255,
            doc => 'Human readable name.',
        },
        based_on => {
            is => 'Genome::ProcessingProfile',
            doc => "Another profile which is used to specify default values for this new one. To qualify a the based on profile must have params, and at least one must be different. Use --param-name='UNDEF' to indicate that a param that is defined for the based on profile should not be for the new profile.",
            is_optional => 1,
        },
        supersedes => {
            is => 'Text',
            doc => 'The processing profile name that this replaces',
            is_optional => 1,
        },
        describe => {
            is => 'Boolean',
            doc => 'Display the output of `genome processing-profile describe` for the processing profile that is created',
            default => 1,
            is_optional => 1,
        },
        created_processing_profile => {
            is => "Genome::ProcessingProfile",
            is_transient => 1,
            is_optional => 1,
            doc => "the newly created processing profile"
        }
    ],
};

sub _resolve_model_subclass_name {
    my $self = shift;
    my $result = $self->_target_class_name;
    $result =~ s/Genome::ProcessingProfile/Genome::Model/;
    return $result;
}

sub _get_help {
    my ($self, $help_type) = @_;

    my $model_class_name = $self->_resolve_model_subclass_name;
    my $help_fn_name = sprintf("help_%s_for_create_profile", $help_type);
    if ($model_class_name->can($help_fn_name)) {
        my $help = $model_class_name->$help_fn_name();
        return $help;
    }
    elsif ($self->_target_class_name->can($help_fn_name)) {
        my $help = $self->target_class_name->$help_fn_name();
        return $help;
    }
    return;
}

sub help_synopsis {
    my $self = shift;
    return $self->_get_help('synopsis');
}

sub help_detail {
    my $self = shift;
    return $self->_get_help('detail');
}

sub create {
    my $class = shift;
    my $bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);

    if ($bx->specifies_value_for('based_on')) {
        my $other_profile = $bx->value_for('based_on');
        return if not $other_profile;

        my @params = $other_profile->params;

        if (not @params) {
            $class->error_message("In order for a processing profile to used as a 'based on', it must have params that can be copied and at least one of them changed. The based on processing profile (".$other_profile->id." ".$other_profile->name.") does not have any params, and cannot be used.");
            return;
        }
$DB::single=1;
        for my $param (@params) {
            my $param_name = $param->name;
            if ($class->can($param_name)) {
                if ($bx->specifies_value_for($param_name) && $bx->value_for($param_name) eq 'UNDEF') {
                    $bx = $bx->remove_filter($param_name);
                    $bx = $bx->add_filter($param_name => undef);
                }
                unless ($bx->specifies_value_for($param_name)) {
                    my $property_meta = $class->__meta__->property_meta_for_name($param_name);
                    my $value;
                    if ($property_meta->is_many) {
                        my @values = $other_profile->$param_name;
                        $value = \@values;;
                    }
                    else {
                        $value = $other_profile->$param_name;
                    }
                    $bx = $bx->add_filter($param_name => $value);
                }
            } else {
                $class->warning_message("Skipping parameter '$param_name'; It does not exist on " . $class . " perhaps '$param_name' was deprecated or replaced.");
            }
        }
    }

    return $class->SUPER::create($bx);
}

sub execute {
    my $self = shift;

    my $profile_class = $self->_target_class_name;

    my %target_params = (
        name => $self->name,
        $self->_get_target_class_params,
    );
    if ($self->supersedes) {
        $target_params{'supersedes'} = $self->supersedes;
    }

    my $processing_profile = $profile_class->create(%target_params);

    unless ($processing_profile) {
        $self->error_message("Failed to create processing profile.");
        return;
    }

    if (my @problems = $processing_profile->__errors__) {
        $self->error_message("Error(s) creating processing profile\n\t".  join("\n\t", map { $_->desc } @problems));
        return;
    }

    $self->status_message('Created processing profile:');
    if($self->describe) {
        my $describer = Genome::ProcessingProfile::Command::Describe->create(
            processing_profiles => [ $processing_profile] ,
        );
        $describer->execute;
    }

    $self->created_processing_profile($processing_profile);
    return 1;
}

sub _properties_for_class {
    my ($self, $class) = @_;

    my $class_meta = $class->__meta__;
    unless ( $class_meta ) {
        $self->error_message("Can't get class meta object for class ($class)");
        return;
    }

    my %properties;
    for my $param ( $class->params_for_class ) {
        my $property = $class_meta->property_meta_for_name($param);
        unless ( $property ){
            $self->error_message("Can't get property for processing profile param ($param)");
            return;
        }
        my $property_name = $property->property_name;
        no warnings;
        $properties{ $property_name } = {
            is => exists $property->{data_type} ? $property->{data_type} : 'Text',
            is_optional => 1,
            is_many => exists $property->{is_many} ? $property->{is_many} : 0,
            doc =>  $property->doc . ($property->is_optional ? " (optional)" : ""),
        };
        if (defined $property->default_value) {
            $properties{ $property_name }->{'default_value'} = $property->default_value;
        }
    }

    return %properties;
}

sub _target_class_property_names {
    my $self = shift;
    my %properties = $self->_properties_for_class( $self->_target_class_name ) or return;
    return keys %properties;
}

sub _get_target_class_params {
    my $self = shift;
    my %params;
    my $class_name = $self->_target_class_name;
    my $class_meta = $class_name->__meta__;
    for my $property ($class_name->params_for_class) {
        my $meta = $class_meta->property_meta_for_name($property);
        my $property_name = $meta->property_name;
        my @values = grep { defined && length } $self->$property_name;
        if (@values == 0) {
            next;
        }
        elsif ($meta->is_many or @values > 1) {
            $params{$property_name} = \@values;
        }
        else {
            $params{$property_name} = $values[0];
        }
    }

    return %params;
}

1;
