package Genome::Model::ReferenceAlignment::Command::Convert36To37;
use strict;
use warnings;

class Genome::Model::ReferenceAlignment::Command::Convert36To37 {
    is => 'Genome::Command::Base',
    has_input => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            shell_args_position => 1,
            doc => 'use these build36 models',
        },
    ],
    has_optional_input => [
        region_of_interest_set_name => {
            is => 'Text',
            doc => 'New build37 models will use this ROI',
        },
    ],
    doc => 'Copy models based on build36 to new ones based on build37',
};

sub help_brief {
    'Copy models based on build36 to new ones based on build37'
}

sub help_detail {
    help_brief()
}

sub help_synopsis {
    help_brief()
}

sub execute {
    my $self = shift;
    my @build37_model_ids;
    my %group_to_models;
    my %roi36to37 = Genome::Model::Command::Services::AssignQueuedInstrumentData->get_build36_to_37_rois();
    for my $model ($self->models) {
        my $result = Genome::Model::Command::Copy->execute(
            model => $model,
            overrides => [
                "name=" . $model->name . ".37",
                'reference_sequence_build=106942997', #Build 37 itself
                'dbsnp_build=106375969',
                'annotation_reference_build=106409296', #build37 ensembl only
            ],
        )->_new_model;

        #Genotype microarray will auto-generate on first build
        map{$_->delete}grep{$_->name eq 'genotype_microarray'}$result->inputs;

        #Set ROI
        if($self->region_of_interest_set_name) {
            $result->region_of_interest_set_name($self->region_of_interest_set_name);
        } elsif ($model->region_of_interest_set_name and exists $roi36to37{$model->region_of_interest_set_name}) {
            print $result->id . " updating region of interest set name to " . $roi36to37{$model->region_of_interest_set_name} . "\n";
            $result->region_of_interest_set_name($roi36to37{$model->region_of_interest_set_name});
        }

        #Track groups to more efficiently add models to groups later
        for my $group ($model->model_groups){
            my $group_id = $group->id;
            push @{$group_to_models{$group_id}}, $result;
        }

        $model->build_requested(1);
        push @build37_model_ids, $result->id;
    }

    while(my ($group_id,$models) = each %group_to_models){
        my $group = Genome::ModelGroup->get($group_id);
        $group->assign_models(@$models);
    }
    print join(',',@build37_model_ids) . "\n";
    return 1;
}
