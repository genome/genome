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
            shell_args_position => 1,
        },
        print_first_diff => {
            is => 'Boolean',
            default_value => 0,
        },
        second_analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            shell_args_position => 2,
        },
        print_second_diff => {
            is => 'Boolean',
            default_value => 0,
        },
    ],
    doc => 'compare instdata between two AnPs.',
};

sub execute {
    my $self = shift;

    my $first_analysis_project = $self->first_analysis_project;
    my $second_analysis_project = $self->second_analysis_project;
    
    my $first_set = Set::Scalar->new($first_analysis_project->instrument_data);
    my $second_set = Set::Scalar->new($second_analysis_project->instrument_data);

    # Compare sets as a status message
    $self->status_message($first_analysis_project->__display_name__ .' set is '. $first_set->compare($second_set) .' as compared to set '. $second_analysis_project->__display_name__);
    $self->status_message($second_analysis_project->__display_name__ .' set is '. $second_set->compare($first_set) .' as compared to set '. $first_analysis_project->__display_name__);

    # show the size of each set 
    $self->status_message($first_analysis_project->__display_name__ .' instrument data: '. $first_set->size);
    $self->status_message($second_analysis_project->__display_name__ .' instrument data: '. $second_set->size);

    $self->status_message($first_set->compare($second_set));

    my $first_set_diff = $first_set->difference($second_set);
    my $second_set_diff = $second_set->difference($first_set);

    $self->status_message($first_analysis_project->__display_name__ .' - '. $second_analysis_project->__display_name__ .': '.  $first_set_diff->size);
    if ($self->print_first_diff) {
         $self->print_diff($first_set_diff);
    }

    $self->status_message($second_analysis_project->__display_name__ .' - '. $first_analysis_project->__display_name__  .': '. $second_set_diff->size);
    if ($self->print_second_diff) {
        $self->print_diff($second_set_diff);
    }
    return 1;
}

sub print_diff {
    my $self = shift;
    my $diff_set = shift;

    for my $diff_member ( $diff_set->members) {
        my @parts = Genome::ProjectPart->get(entity_id => $diff_member->id);
        my @projects = Genome::Project->get([map $_->project_id, @parts]);
        print $diff_member->sequencing_platform ."\t". $diff_member->sample->name ."\t". $diff_member->id ."\t". $diff_member->ignored ."\t". $diff_member->flow_cell_id ."\t". join(',',sort(map {$_->id} @projects)) ."\n";
    }
    return 1;
}

1;
