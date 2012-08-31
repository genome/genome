package Genome::Model::Command::Admin::CheckRoi;

class Genome::Model::Command::Admin::CheckRoi {
    is => 'Genome::Command::Base',
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            shell_args_position => 1,
        },
        auto_action => {
            is => 'Boolean',
            default => 0,
        },
    ],
};

use strict;
use warnings;
use Genome;

sub execute {
    my $self = shift;

    my %build36_to_37_rois = Genome::Model::Command::Services::AssignQueuedInstrumentData->get_build36_to_37_rois();

    # if($reference_sequence_build and $reference_sequence_build->name eq 'GRCh37-lite-build37') {
    #     $wuspace_roi_list = 'NCBI-human.combined-annotation-58_37c_cds_exon_and_rna_merged_by_gene';
    # } else {
    #     $wuspace_roi_list = 'NCBI-human.combined-annotation-54_36p_v2_CDSome_w_RNA';
    # }
    # 
    # $self->assign_capture_inputs($model, $capture_target, $roi_list)

    my @models = $self->models;

    my $root_build37_ref_seq = Genome::Model::Build::ImportedReferenceSequence->get(name => 'GRCh37-lite-build37') || die;

    for my $model (@models) {
        my $reference_sequence_build = $model->reference_sequence_build || die;

        unless ($reference_sequence_build->is_compatible_with($root_build37_ref_seq)) {
            print $model->id . "\tis not a build37 model\n";
            next;
        }

        #if ($model->name =~ /\.wu-space$/) {
        #    print $model->id . "\tis a wu-space model\n";
        #    next;
        #}

        my $old_roi_set_name = $model->region_of_interest_set_name;
        if (exists $build36_to_37_rois{$old_roi_set_name}) {
            my $new_roi_set_name = $build36_to_37_rois{$old_roi_set_name};
            if ($self->auto_action) {
               my $old_roi_input = Genome::Model::Input->get(model_id => $model->id, name => 'region_of_interest_set_name');
               my $new_roi_input = $model->add_input(
                   name             => "region_of_interest_set_name",
                   value_class_name => "UR::Value",
                   value_id         => $new_roi_set_name
               );
               if ($new_roi_input) {
                   $old_roi_input->delete;
               }
               $model->build_requested(1);
           }
           else {
               print $model->id . "\tneeds region_of_interest_set_name updated\n";
           }
       }
       else {
           print $model->id . "\thas correct region_of_interest_set_name\n";
       }
   }

   return 1;
}
