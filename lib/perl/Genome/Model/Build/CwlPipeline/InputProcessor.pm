package Genome::Model::Build::CwlPipeline::InputProcessor;

use strict;
use warnings;

use feature qw(switch);

use Genome;

class Genome::Model::Build::CwlPipeline::InputProcessor {
    is_abstract => 1,
    is => 'UR::Value::Text',
    attributes_have => [
        is_input => { is => 'Boolean', is_optional => 1 },
        is_simple => { is => 'Boolean', is_optional => 1 },
        input_type => { is => 'Text', is_optional => 1 },
    ],
    subclass_description_preprocessor => 'Genome::Model::Build::_preprocess_subclass_description',
    has => [
        id => {
            is => 'UR::Value::Text',
        },
        build => {
            is => 'Genome::Model::Build::CwlPipeline',
            id_by => 'build_id',
        },
        build_id => {
            via => '__self__',
            to => 'id',
        },
    ],
};

sub simple_inputs {
    my $self = shift;

    my %inputs;
    for my $input ($self->__meta__->properties(is_input => 1, is_simple => 1)) {
        my $name = $input->property_name;
        my $type = $input->input_type;

        for ($type) {
            when ('File') {
                $inputs{$name} = {
                    class => 'File',
                    path => '' . $self->$name,
                };
            }
            when ('Text') {
                $inputs{$name} = '' . $self->$name;
            }
            when ('Number') {
                $inputs{$name} = 0 + $self->$name;
            }
            default {
                $self->fatal_message('Unknown input type %s in simple input %s.', $type, $name);
            }
        }
    }

    return \%inputs;
}


sub _preprocess_subclass_description {
    my ($class, $desc) = @_;

    my @names = keys %{ $desc->{has} };
    for my $prop_name (@names) {
        my $prop_desc = $desc->{has}{$prop_name};

        # skip old things for which the developer has explicitly set-up indirection
        next if $prop_desc->{id_by};
        next if $prop_desc->{via};
        next if $prop_desc->{reverse_as};
        next if $prop_desc->{implied_by};

        if ($prop_desc->{'is_input'}) {
            my $assoc = $prop_name . '_association' . ($prop_desc->{is_many} ? 's' : '');
            next if $desc->{has}{$assoc};

            if ($prop_desc->{'data_type'}) {
                my $prop_class = UR::Object::Property->_convert_data_type_for_source_class_to_final_class(
                    $prop_desc->{'data_type'},
                    $class
                );
            }

            $desc->{has}{$assoc} = {
                property_name => $assoc,
                implied_by => $prop_name,
                is => 'Genome::Model::Build::Input',
                via => 'build',
                to => 'inputs',
                where => [ name => $prop_name ],
                is_mutable => $prop_desc->{is_mutable},
                is_optional => $prop_desc->{is_optional},
                is_many => 1, #$prop_desc->{is_many},
            };

            %$prop_desc = (%$prop_desc,
                via => $assoc,
                to => Genome::Model->_resolve_to_for_prop_desc($prop_desc),
            );
        }
    }

    return $desc;
}

1;
