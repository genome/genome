package Genome::Model::Tools::SvSim::Worklet;

use Genome;

use strict;
use warnings;

class Genome::Model::Tools::SvSim::Worklet {
    has_transient_optional => [
        dag => {},
        _output_directories => {
            is => "HASH",
        },
    ]
};

sub _add_output_directories {
    my ($self, @dirs) = @_;
    @{$self->_output_directories}{@dirs} = 1;
}

sub _unpack_property_names {
    my $property = shift;

    my ($src_name, $dst_name);
    if (ref($property) eq "ARRAY") {
        ($src_name, $dst_name) = @$property;
    }
    else {
        $src_name = $dst_name = $property;
    }

    return ($src_name, $dst_name);
}


sub generate_workflow {
    my $self = shift;
    die $self->error_message(sprintf(
        "Class %s must override generate_workflow function!",
        $self->class));
}

sub create_model {
    my ($self, %params) = @_;
    my $dag = Genome::WorkflowBuilder::DAG->create(%params);
    $self->dag($dag);
    return $dag;
}

sub add_command {
    my ($self, %params) = @_;
    my $dag = $self->dag;
    die "Add command called before create_model" unless $dag;

    my @req = qw(name command);
    for my $req (@req) {
        die "Required parameter $req missing!" unless exists $params{$req};
    }

    my $op = Genome::WorkflowBuilder::Command->create(%params);
    $dag->add_operation($op);
    return $op;
}

sub connect_inputs {
    my ($self, $operation, @params) = @_;
    my $dag = $self->dag;

    for my $p (@params) {
        my ($src_name, $dst_name) = _unpack_property_names($p);
        $dag->connect_input(
            input_property => $src_name,
            destination => $operation,
            destination_property => $dst_name,
            );
    }
}

sub connect_outputs {
    my ($self, $operation, @params) = @_;
    my $dag = $self->dag;

    for my $p (@params) {
        my ($src_name, $dst_name) = _unpack_property_names($p);
        $dag->connect_output(
            source => $operation,
            source_property => $src_name,
            output_property => $dst_name,
            );
    }
}

sub create_links {
    my ($self, $src_op, $dst_op, @property_list) = @_;

    my $dag = $self->dag;

    for my $property (@property_list) {
        my ($src_name, $dst_name) = _unpack_property_names($property);

        printf "ADDING LINK %s/%s -> %s/%s ($property)\n",
            $src_op->name, $src_name,
            $dst_op->name, $dst_name;

        $dag->create_link(
            source => $src_op,
            source_property => $src_name,
            destination => $dst_op,
            destination_property => $dst_name
            );
    }
}

1;
