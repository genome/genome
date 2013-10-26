package Genome::Model::Convergence;

use strict;
use warnings;

use Genome;

class Genome::Model::Convergence {
    is  => 'Genome::ModelDeprecated',
    has => [
        group_id => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'group', value_class_name => 'Genome::ModelGroup' ],
            doc => 'The id for the ModelGroup for which this is the Convergence model',

            is_mutable => 1,
        },
        group => {
            is => 'Genome::ModelGroup',
            id_by => 'group_id',
            doc => 'The ModelGroup for which this is the Convergence model',
        },
        members => {
            is => 'Genome::Model',
            via => 'group',
            is_many => 1,
            to => 'models',
            doc => 'Models that are members of this Convergence model.',
        },
    ],
    doc => <<EODOC
This model type attempts to use the data collected across many samples to generalize and summarize
knowledge about the group as a whole.
EODOC
};


# Updates the convergence model immediately prior to a build starting. In this case, makes sure
# that the convergence model's subject is correct.
# TODO: make this a callback for an observer on Genome::Model::Build::Convergence create()
# so we can get rid of this method.  It is only used for Convergence and refalign.
sub check_for_updates {
    my $self = shift;
    my $subject = $self->subject;
    my $group_subject = $self->group->infer_group_subject;

    if ($group_subject->class eq $subject->class and $group_subject->id eq $subject->id) {
        return 1;
    }

    $self->subject($group_subject);
    return 1;
}

# Return true of the lists have the same values, possibly in a different order
sub _are_list_contents_same {
    my($list1, $list2) = @_;

    return if (@$list1 != @$list2);  # different length

    my @list1 = sort @$list1;
    my @list2 = sort @$list2;

    my $len = scalar(@list1);
    for (my $i = 0; $i < $len; $i++) {
        return unless $list1[$i] eq $list2[$i];
    }

    return 1;
}

sub build_needed {
    my $self = shift;

    my @group_members = $self->members;
    my @b = Genome::Model::Build->get(model_id => [map($_->id, @group_members)], '-hint' => ['the_master_event']); #preload in one query

    my @potential_members = grep(defined $_->last_complete_build, @group_members);
    unless(scalar @potential_members) {
        #$self->status_message('Skipping convergence rebuild--no succeeded builds found among members.');
        return;
    }

    #Check to see if our last build already used all the same builds as we're about to
    if($self->last_complete_build) {
        my $last_build = $self->last_complete_build;

        my @last_members = $last_build->members;

        if (_are_list_contents_same(\@last_members, \@potential_members)) {
            #$self->status_message('Skipping convergence rebuild--list of members that would be included is identical to last build.');
            return;

            #Potentially if some of the underlying builds in the $build->all_subbuilds_closure have changed, a rebuild might be desired
            #For now this will require a manual rebuild (`genome model build start`)
        }
    }

    return 1;
}

sub _resolve_workflow_for_build {
    my $self = shift;
    my $build = shift;

    my $operation = Workflow::Operation->create_from_xml(__FILE__ . '.xml');

    my $log_directory = $build->log_directory;
    $operation->log_dir($log_directory);
    $operation->name($build->workflow_name);

    return $operation;
}

sub map_workflow_inputs {
    my $self = shift;
    my $build = shift;

    my @inputs = ();

    my @members = $build->members;

    #Check that the members are ready for convergence
    for my $member (@members) {
         unless($member->status eq 'Succeeded') {
            $self->error_message("Tried to use non-succeeded build! " . $member->id);
            return;
        }
    } 

    my $data_directory = $build->data_directory;
    unless ($data_directory) {
        $self->error_message("Failed to get a data_directory for this build!");
        return;
    }

    #Assign filenames
    my %default_filenames = $self->default_filenames;
    for my $parameter (keys %default_filenames) {
        push @inputs,
            $parameter => ($data_directory . "/" . $default_filenames{$parameter});
    }

    push @inputs,
        build_id => $build->id,
        skip_if_output_present => 1;

    return @inputs;
}

sub default_filenames{
    my $self = shift;

    my %default_filenames = (
    );

    return %default_filenames;
}

1;
