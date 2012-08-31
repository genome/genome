package Genome::Model::ReferenceAlignment::Command::Convert37To36;
use strict;
use warnings;

class Genome::Model::ReferenceAlignment::Command::Convert37To36 {
    is => 'Genome::Command::Base',
    has_input => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            shell_args_position => 1,
            doc => 'use these build37 models',
        },
    ],
    has_optional_input => [
        region_of_interest_set_name => {
            is => 'Text',
            doc => 'New build36 models will use this ROI',
        },
    ],
    doc => 'Copy models based on build37 to new ones based on build36',
};

sub help_brief {
    'Copy models based on build37 to new ones based on build36'
}

sub help_detail {
    help_brief()
}

sub help_synopsis {
    help_brief()
}

sub execute {
    my $self = shift;
    my @build36_model_ids;
    my %group_to_models;
    my %roi37to36 = reverse Genome::Model::Command::Services::AssignQueuedInstrumentData->get_build36_to_37_rois();
    for my $model ($self->models) {
        my $result = Genome::Model::Command::Copy->execute(
            model => $model,
            overrides => [
                "name=" . $model->name . ".36",
                'reference_sequence_build=101947881', #Build 36 itself
                'dbsnp_build=106227442',
                'annotation_reference_build=113115679',#Model/ImportedAnnotation.pm for build36
            ],
        )->_new_model;

        #Genotype microarray will auto-generate on first build
        map{$_->delete}grep{$_->name eq 'genotype_microarray'}$result->inputs;

        #Set ROI
        if($self->region_of_interest_set_name) {
            $result->region_of_interest_set_name($self->region_of_interest_set_name);
        } elsif ($model->region_of_interest_set_name and exists $roi37to36{$model->region_of_interest_set_name}) {
            print $result->id . " updating region of interest set name to " . $roi37to36{$model->region_of_interest_set_name} . "\n";
            $result->region_of_interest_set_name($roi37to36{$model->region_of_interest_set_name});
        }

        #Track groups to more efficiently add models to groups later
        for my $group ($model->model_groups){
            my $group_id = $group->id;
            push @{$group_to_models{$group_id}}, $result;
        }

        $model->build_requested(1);
        push @build36_model_ids, $result->id;
    }

    while(my ($group_id,$models) = each %group_to_models){
        my $group = Genome::ModelGroup->get($group_id);
        $group->assign_models(@$models);
    }
    print join(',',@build36_model_ids) . "\n";
    return 1;
}
