package Genome::InstrumentData::Composite::Workflow;

use strict;
use warnings;

use Genome;

use Workflow;
use Workflow::Simple qw(run_workflow_lsf);

class Genome::InstrumentData::Composite::Workflow {
    is => 'Command',
    has => [
        inputs => {
            is => 'HASH',
            doc => 'Contains { name => [values] } for the terms in the alignment strategy',
        },
        strategy => {
            is => 'Text',
            doc => 'The alignment strategy to run',
        },
        merge_group => {
            is => 'Text',
            default_value => 'sample',
            valid_values => ['sample', 'all'],
            doc => 'How to group the instrument data when merging',
        },
    ],
    has_transient_optional => [
        _workflow => {
            is => 'Workflow::Model',
            doc => 'The underlying workflow to run the alignment/merge',
            is_output => 1,
        },
        _result_ids => {
            is_many => 1,
            doc => 'The alignments created/found as a result of running the workflow',
            is_output => 1,
        },
        log_directory => {
            is => 'Text',
            doc => 'Where to write the workflow logs',
        },
        refiners => {
            is => 'ArrayRef',
            doc => 'the refiner strategies used in the workflow',
        },
        first_refiner => {
            is => 'Text',
            doc => 'the first refiner strategy used in the workflow',
        },
        last_refiner => {
            is => 'Text',
            doc => 'the last refiner strategy used in the workflow',
        },
    ],
};

sub execute {
    my $self = shift;

    my $tree = $self->_process_strategy($self->strategy);
    die unless $tree;

    my ($workflow, $inputs) = $self->_generate_workflow($tree);
    $self->_workflow($workflow);

    if($self->log_directory) {
        $workflow->log_dir($self->log_directory);
    }

    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die 'Errors validating workflow';
    }

    my @results = $self->_run_workflow($workflow, $inputs);
    $self->_result_ids(\@results);

    return 1;
}

sub _process_strategy {
    my $self = shift;
    my $strategy_string = shift;

    $self->debug_message('Analyzing strategy...');

    my $strategy = Genome::InstrumentData::Composite::Strategy->create(
        strategy => $strategy_string,
    );
    my $tree = $strategy->execute;

    unless( $tree ) {
        $self->error_message('Failed to analyze strategy: '.$self->strategy);
        return;
    }

    $self->debug_message('Analyzing strategy...OK');
    return $tree;
}

sub _generate_workflow {
    my $self = shift;
    my $tree = shift;

    $self->debug_message('Generating workflow...');

    unless(exists $tree->{data}) {
        die $self->error_message('No data specified in input');
    }

    return Genome::InstrumentData::Composite::Workflow::Generator->generate($tree, $self->inputs, $self->merge_group);
}

sub _run_workflow {
    my $self = shift;
    my $workflow = shift;
    my $inputs = shift;

    Genome::Sys->disconnect_default_handles;

    $self->debug_message('Running workflow...');
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, @$inputs);

    unless($result){
        $self->error_message(join("\n", map {
            ($_->can('name') ? $_->name .': ' : ''). $_->error
        } @Workflow::Simple::ERROR));

        die $self->error_message("Workflow did not return correctly.");
    }

    my @result_ids;
    for my $key (keys %$result) {
        if($key =~ /result_id/) {
            push @result_ids, $result->{$key};
        }
    }

    $self->debug_message('Produced results: ' . join(', ', @result_ids));
    $self->debug_message('Workflow complete.');

    return @result_ids;
}

1;
