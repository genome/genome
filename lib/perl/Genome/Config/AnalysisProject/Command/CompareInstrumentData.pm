package Genome::Config::AnalysisProject::Command::CompareInstrumentData;

use strict;
use warnings;

use Genome;

use Set::Scalar;
use Data::Dumper;

class Genome::Config::AnalysisProject::Command::CompareInstrumentData {
    is => ['Command::V2'],
    has => [
        first_analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'The first Analysis Project to compare Instrument Data.',
            shell_args_position => 1,
        },
        second_analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'The second Analysis Project to compare Instrument Data.',
            shell_args_position => 2,
        },
    ],
    has_optional => [
        list_intersection => {
            is => 'Boolean',
            doc => 'List the instrument data that is shared between Analysis Projects.',
            default_value => 0,
        },
        list_first_diff => {
            is => 'Boolean',
            doc => 'List the Instrument Data that is unique to the first Analysis Project.',
            default_value => 0,
        },
        list_second_diff => {
            is => 'Boolean',
            doc => 'List the instrument data that is unique to the second Analysis Project.',
            default_value => 0,
        },
        show => {
            is => 'Text',
            doc => 'The Instrument Data properties to show when listing. See `genome instrument-data list --help` for more details.',
            example_values => ['sequencing_platform,sample.name,id,ignored,flow_cell_id,projects.id,analysis_projects.id'],
        },
    ],
    doc => 'Compare Instrument Data between two Analysis Projects and optionally list the Instrument Data.',
};

sub execute {
    my $self = shift;

    my $first_analysis_project = $self->first_analysis_project;
    my $second_analysis_project = $self->second_analysis_project;

    my $first_set = Set::Scalar->new($first_analysis_project->instrument_data);
    my $second_set = Set::Scalar->new($second_analysis_project->instrument_data);

    $self->status_message($first_analysis_project->__display_name__ .' set is '. $first_set->compare($second_set) .' as compared to set '. $second_analysis_project->__display_name__);
    $self->status_message($second_analysis_project->__display_name__ .' set is '. $second_set->compare($first_set) .' as compared to set '. $first_analysis_project->__display_name__);

    $self->status_message($first_analysis_project->__display_name__ .' instrument data: '. $first_set->size);
    $self->status_message($second_analysis_project->__display_name__ .' instrument data: '. $second_set->size);

    $self->status_message($first_set->compare($second_set));

    my $set_x = $first_set->intersection($second_set);
    $self->status_message($first_analysis_project->__display_name__ .' x '. $second_analysis_project->__display_name__ .': '.  $set_x->size);
    if ($self->list_intersection) {
        $self->list_set($set_x);
    }

    my $first_set_diff = $first_set->difference($second_set);
    my $second_set_diff = $second_set->difference($first_set);

    $self->status_message($first_analysis_project->__display_name__ .' - '. $second_analysis_project->__display_name__ .': '.  $first_set_diff->size);
    if ($self->list_first_diff) {
         $self->list_set($first_set_diff);
    }

    $self->status_message($second_analysis_project->__display_name__ .' - '. $first_analysis_project->__display_name__  .': '. $second_set_diff->size);
    if ($self->list_second_diff) {
        $self->list_set($second_set_diff);
    }

    return 1;
}

sub list_set {
    my $self = shift;
    my $diff_set = shift;

    my $operator = ':';
    my $separator = '/';

    my @ids = map {$_->id} $diff_set->members;
    if (!@ids) {
        $self->fatal_message('No members in set!');
    } elsif (@ids == 1) {
        $operator = '=';
        $separator = '';
    }
    my %lister_params = (
        filter => 'id'. $operator . join($separator,@ids),
    );
    if ($self->show) {
        $lister_params{'show'} = $self->show;
    }
    my $lister = Genome::InstrumentData::Command::List->create(%lister_params);
    unless ($lister) {
        $self->fatal_message('Failed to create Instrument Data list command with params: '. Data::Dumper::Dumper(%lister_params));
    }
    unless ($lister->execute) {
        $self->fatal_message('Failed to execute Instrument Data list command!');
    }
    return 1;
}

1;
