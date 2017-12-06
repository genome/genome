package Genome::Model::CwlPipeline;

use strict;
use warnings;

use feature qw(switch);
use Genome;

use Set::Scalar;

class Genome::Model::CwlPipeline {
    is => 'Genome::Model',
    has_param => {
        main_workflow_file => {
            is => 'Text',
        },
        primary_docker_image => {
            is => 'Text',
            doc => 'docker image for the main toil worker jobs',
        },
        short_pipeline_name => {
            is => 'Text',
            is_optional => 1,
            doc => 'short name for pipeline to include in default model names',
        },
    },
};

sub create {
    my $class = shift;

    # If create is being called directly on this class or on an abstract subclass, SUPER::create will
    # figure out the correct concrete subclass (if one exists) and call create on it.
    if ($class->__meta__->is_abstract) {
        return $class->SUPER::create(@_);
    }

    my ($bx, %extra) = $class->define_boolexpr(@_);

    my $inputs = delete $extra{input_data};
    unless ($inputs) {
        return $class->SUPER::create(@_);
    }

    $bx = $bx->remove_filter('input_data');

    my $tx = UR::Context::Transaction->begin(commit_validator => sub { 1 });
    my $guard = Scope::Guard->new(sub { local $@; $tx->rollback });

    my $self = $class->SUPER::create($bx);
    $self->process_input_data(%$inputs);

    $tx->commit() && $guard->dismiss();

    return $self;
}

sub get {
    my $class = shift;

    my ($bx, %extra) = $class->define_boolexpr(@_);

    my $input_data = delete $extra{input_data};
    unless($input_data) {
        return $class->SUPER::get(@_);
    }

    $bx = $bx->remove_filter('input_data');
    my $input_info = $class->determine_input_objects(%$input_data);
    my @input_values = values %$input_info;
    for my $v (@input_values) {
        $bx = $bx->add_filter('inputs.value_id', $v->id);
    }

    my @candidates = $class->SUPER::get($bx);

    #now check that we match the same value ID and name
    for my $name (keys %$input_info) {
        my $value = $input_info->{$name};
        my $baseline = (ref $value eq 'ARRAY'?
            Set::Scalar->new(@$value) :
            Set::Scalar->new($value));

        @candidates = grep {
            my @in = grep { $_->name eq $name } $_->inputs();
            my $actual = Set::Scalar->new(map $_->value, @in);

            $baseline->is_subset($actual);
        } @candidates;
    }

    return $class->context_return(@candidates);
}

sub map_workflow_inputs {
    my $self = shift;
    my $build = shift;

    return (build => $build);
}

sub _resolve_workflow_for_build {
    my $self = shift;
    my $build = shift;

    my $wf = Genome::WorkflowBuilder::DAG->create(
        name => $build->workflow_name,
        log_dir => $build->log_directory,
    );

    my $cmd = Genome::WorkflowBuilder::Command->create(
        name => 'Toil Runner',
        command => 'Genome::Model::CwlPipeline::Command::Run',
    );
    $wf->add_operation($cmd);
    $wf->connect_input(
        input_property => 'build',
        destination => $cmd,
        destination_property => 'build'
    );
    $wf->connect_output(
        source => $cmd,
        source_property => 'result',
        output_property => 'result',
    );

    return $wf;
}


#handle specification of inputs such as might be found in an AnP Config YAML
sub process_input_data {
    my $self = shift;
    my %data = @_;

    for my $name (keys %data) {
        my $value = $self->determine_input_object($name, $data{$name});
        unless ($value) {
            $self->fatal_message('Failed to determine input for name "%s" and value "%s".', $name, $data{$name});
        }

        $self->add_input(
            name => $name,
            value_id => $value->id,
            value_class_name => $value->class,
        );
    }

    return 1;
}

sub determine_input_objects {
    my $class = shift;
    my %data = @_;

    my %inputs;
    for my $name (keys %data) {
        my $value = $class->determine_input_object($name, $data{$name});
        unless($value) {
            $class->fatal_message('Unable to determine input for name "%s" and value "%s".', $name, $data{name});
        }

        $inputs{$name} = $value;
    }

    return \%inputs;
}

sub determine_input_object {
    my $class = shift;
    my $name = shift;
    my $value_identifier = shift;

    if (ref $value_identifier eq 'HASH' and exists $value_identifier->{value_class_name} and exists $value_identifier->{value_id}) {
        my $value_class = $value_identifier->{value_class_name};
        unless ($value_class->can('get')) {
            $class->fatal_message('Unknown class name specification: %s', $value_class);
        }

        my $object = $value_class->get($value_identifier->{value_id});
        unless ($object) {
            $class->fatal_message('No %s found with ID %s.', $value_class, $value_identifier->{value_id});
        }

        return $object;
    }

    #This could use something like Command::Dispatch::Shell->resolve_param_value_from_text
    #to allow for non-ID values
    for ($name) {
        when ('instrument_data')          {
            return Genome::InstrumentData->get($value_identifier);
        }
        when ('reference_build')          {
            return Genome::Model::Build::ReferenceSequence->get($value_identifier);
        }
        when ('annotation_build')         {
            return Genome::Model::Build::ImportedAnnotation->get($value_identifier);
        }
        default                           {
            return UR::Value::Text->get($value_identifier);
        }
    }

    return;
}

sub _additional_parts_for_default_name {
    my $self = shift;

    my @parts;
    my $name = $self->processing_profile->short_pipeline_name;

    push @parts, $name if $name;

    return @parts;
}

1;
