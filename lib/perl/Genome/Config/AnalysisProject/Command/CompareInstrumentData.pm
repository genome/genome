package Genome::Config::AnalysisProject::Command::CompareInstrumentData;

use strict;
use warnings;

use Genome;

use Set::Scalar;

class Genome::Config::AnalysisProject::Command::CompareInstrumentData {
    is => ['Command::V2'],
    has => [
        analysis_project_a => {
            is => 'Genome::Config::AnalysisProject',
            shell_args_position => 1,
        },
        print_a_diff => {
            is => 'Boolean',
            default_value => 0,
        },
        analysis_project_b => {
            is => 'Genome::Config::AnalysisProject',
            shell_args_position => 2,
        },
        print_b_diff => {
            is => 'Boolean',
            default_value => 0,
        },
    ],
    doc => 'compare instdata between two AnPs.',
};

sub execute {
    my $self = shift;

    my $set_a = Set::Scalar->new($self->analysis_project_a->instrument_data);
    my $set_b = Set::Scalar->new($self->analysis_project_b->instrument_data);

    $self->status_message('A instrument data: '. $set_a->size);
    $self->status_message('B instrument data: '. $set_b->size);

    $self->status_message($set_a->compare($set_b));

    my $set_a_diff = $set_a->difference($set_b);
    my $set_b_diff = $set_b->difference($set_a);

    $self->status_message('A - B: '. $set_a_diff->size);
    if ($self->print_a_diff) {
         $self->print_diff($set_b_diff);
    }

    $self->status_message('B - A: '. $set_b_diff->size);
    if ($self->print_b_diff) {
        $self->print_diff($set_b_diff);
    }
    return 1;
}

sub print_diff {
    my $diff_set = shift;

    for my $diff_member ( $diff_set->members) {
        my @parts = Genome::ProjectPart->get(entity_id => $diff_member->id);
        my @projects = Genome::Project->get([map $_->project_id, @parts]);
        print $diff_member->sequencing_platform ."\t". $diff_member->sample->name ."\t". $diff_member->id ."\t". $diff_member->ignored ."\t". $diff_member->flow_cell_id ."\t". join(',',sort(map {$_->id} @projects)) ."\n";
    }
    return 1;
}

1;
