package Genome::TestObjGenerator::Base;

use strict;
use warnings;

use Class::ISA;

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
    my $class = shift;
    my @parents = grep { $_ ne __PACKAGE__ } Class::ISA::self_and_super_path($class);
    my @params = map {
        my $var = join('::', $_, 'required_params');
        no strict 'refs';
        @$var;
    } @parents;
    return \@params;
}

# generate_obj and an `our @required_params` are the API
sub generate_obj {
    die "Abstract class - must override generate_obj";
}

sub setup_object {
    my $class = shift;
    my %params = @_;

    my $filled_in_params = $class->fill_in_missing_params(%params);
    return $class->generate_obj(%$filled_in_params);
}

1;
