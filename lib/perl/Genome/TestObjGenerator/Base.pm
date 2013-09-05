package Genome::TestObjGenerator::Base;

use strict;
use warnings;
use Genome;
use Genome::TestObjGenerator::Util;

sub fill_in_missing_params {
    my $self = shift;
    my %provided = @_;
    my $required = $self->get_required_params();
    for my $required_param (@$required) {
        unless (defined $provided{$required_param}) {
            my $method_name = "create_$required_param";
            unless ($self->can($method_name)) {
                die("No method available to create default value for required param $required_param");
            }
            $provided{$required_param} = $self->$method_name();
        }
    }
    return \%provided;
}

sub get_required_params {
    #override this method if your object has methods that are required
}

sub generate_obj {
    die ("Abstract class - must override generate_obj");
}

sub setup_object {
    my $self = shift;
    my %params = @_;

    my $filled_in_params = $self->fill_in_missing_params(%params);
    return $self->generate_obj(%$filled_in_params);
}

1;

