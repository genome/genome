package Genome::Model::SomaticValidation::Command::UpdateSamples;

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::Model::SomaticValidation::Command::UpdateSamples {
    is => 'Command::V2',
    has => [
        models => { 
            is => 'Genome::Model::SomaticValidation', 
            is_many => 1,
            doc => 'the models to update',
            shell_args_position => 1,

        },
    ],
    doc => 'switch the sample on a model to an equivalent sample based on available instrument data (and eventually pending workorders)',
};

sub help_synopsis {
    return <<EOS;

# one model
genome model somatic-validation update-samples 12345

# all models using a given set of capture probes
genome model somatic-validation update-samples model_groups.name="AMLx150 Validation (RT#78086)"

# all models in a given group
genome model somatic-validation update-samples target_region_set_name="AML_KP - OID36117 capture chip set"

EOS
}

sub sub_command_category { 'analyst tools' }

sub execute {
    my $self = shift;
    my @models = $self->models;
    my @changes;
    for my $model (@models) {
        my $tumor_sample = $model->tumor_sample;
        my $normal_sample = $model->normal_sample;

        my $patient = $tumor_sample->patient();
        unless ($patient == $normal_sample->patient()) {
            die "patients do not match for tumor and normal on model " . $model->__display_name__;
        }

        my @patient_instdata = Genome::InstrumentData::Solexa->get(
            target_region_set_name => $model->target_region_set_name,
            'sample.patient.id' => $patient->id,
        );

        unless (@patient_instdata) {
            $self->debug_message("No instrument data for patient " . $patient->__display_name__ . " on the target set yet.  Cannot update the model until we have logic to check the workorder");
            next;
        }

        my @instdata_tumor;
        my @instdata_normal;
        my @instdata_unknown;
        for my $instdata (@patient_instdata) {
            if ($instdata->sample == $tumor_sample) {
                push @instdata_tumor, $instdata,
            }
            elsif ($instdata->sample == $normal_sample) {
                push @instdata_normal, $instdata;
            }
            else {
                push @instdata_unknown, $instdata;
            }
        }

        if (@instdata_tumor and @instdata_normal) {
            $self->debug_message("Some data found for both samples.  Any other instata may be for other models.");
            next;
        }

        if (@instdata_tumor == 0) {
            my %tumor_equiv_samples = map { $_->id => $_ } Genome::Sample->get(
                source => $patient,
                common_name => $tumor_sample->common_name,
                tissue_label => $tumor_sample->tissue_label,
                tissue_desc => $tumor_sample->tissue_desc,
            );

            if (%tumor_equiv_samples) {
                my @has_instdata = grep { $tumor_equiv_samples{$_->sample->id} } @instdata_unknown;
                if (@has_instdata == 0) {
                    $self->debug_message("unknown instrument data in are not suitable swaps for the tumor");
                }
                else {
                    my %sample_ids = map { $_->sample_id => 1 } @has_instdata;
                    if (keys(%sample_ids) > 1) {
                        $self->debug_message("ambiguous replacement inst data for the tumor");
                    }
                    else {
                        # probably add the instrument data too
                        my ($new_sample_id) = keys(%sample_ids);
                        my $new_sample = Genome::Sample->get($new_sample_id);
                        $model->tumor_sample($new_sample);

                        for my $instdata (@has_instdata){
                            my $assign = Genome::Model::Command::InstrumentData::Assign::Expression->create(
                                instrument_data => [$instdata], 
                                model => $model,
                                );
                            unless($assign->execute){
                                $self->error_message('Failed to execute instrument data assign for model ' . $model->id . ' and instrument data ' . $instdata->id);
                            }
                        }

                        my %change = (old_sample => $tumor_sample, 
                                      new_sample => $new_sample,
                                      instrument_data => \@has_instdata,
                                      model => $model,
                                      );
                        push @changes, \%change;
                    }
                }
            }
            else {
                $self->debug_message("No equivalent tumor data found for model " . $model->__display_name__);
            }
        }

        if (@instdata_normal == 0) {
            my %normal_equiv_samples = map { $_->id => $_ } Genome::Sample->get(
                source => $patient,
                common_name => $normal_sample->common_name,
                tissue_label => $normal_sample->tissue_label,
                tissue_desc => $normal_sample->tissue_desc,
            );

            if (%normal_equiv_samples) {
                my @has_instdata = grep { $normal_equiv_samples{$_->sample->id} } @instdata_unknown;
                if (@has_instdata == 0) {
                    $self->debug_message("unknown instrument data in are not suitable swaps for the normal");
                }
                else {
                    my %sample_ids = map { $_->sample_id => 1 } @has_instdata;
                    if (keys(%sample_ids) > 1) {
                        $self->debug_message("ambiguous replacement instdata for the normal");
                    }
                    else {
                        my ($new_sample_id) = keys(%sample_ids);
                        my $new_sample = Genome::Sample->get($new_sample_id);
                        $model->normal_sample($new_sample);
                        for my $instdata (@has_instdata){
                            my $assign = Genome::Model::Command::InstrumentData::Assign::Expression->create(
                                instrument_data => [$instdata], 
                                model => $model,
                                );
                            unless($assign->execute){
                                $self->error_message('Failed to execute instrument data assign for model ' . $model->id . ' and instrument data ' . $instdata->id);
                            }
                        }

                        my %change = (old_sample => $normal_sample, 
                                      new_sample => $new_sample,
                                      instrument_data => \@has_instdata,
                                      model => $model,
                                      );
                        push @changes, \%change;
                    }
                }
            }
            else {
                $self->debug_message("No equivalent normal data found for model " . $model->__display_name__);
            }
        }
    }
    $self->print_report(@changes);

    return 1;
}

sub print_report {
    my $self = shift;
    my @changes = @_;
    $self->_print_headers();
    for my $change_ref (@changes){
        my %change = %$change_ref;
        my $model = $change{model};
        my $old_sample = $change{old_sample};
        my $new_sample = $change{new_sample};
        my @instrument_data = @{$change{instrument_data}};
        my $output = join("\t", $model->id, $model->name, $old_sample->name, $new_sample->name, join('/', map($_->id, @instrument_data))); 
        print $output, "\n";
    }
    return 1;
}

sub _print_headers {
    my $self = shift;
    print join("\t", 'MODEL_ID', 'MODEL_NAME', 'OLD SAMPLE NAME', 'NEW SAMPLE NAME', 'INSTRUMENT DATA ASSIGNED'), "\n";
    return 1;
}

1;

