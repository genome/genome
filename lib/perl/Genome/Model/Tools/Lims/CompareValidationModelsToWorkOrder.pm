package Genome::Model::Tools::Lims::CompareValidationModelsToWorkOrder;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Lims::CompareValidationModelsToWorkOrder {
    is => 'Command::V2',
    has => [
        models => {
            is => 'Genome::Model::SomaticValidation',
            doc => 'SomaticValidation models whose samples to check against the work order',
            is_many => 1,
            shell_args_position => 2,
        },
        work_order => {
            is => 'Number',
            doc => 'id of a LIMS work order',
            shell_args_position => 1,
        },
        show => {
            is => 'Text',
            is_optional => 1,
            doc => 'Controls which properties of the samples are displayed',
            default_value => 'id,name',
        },
    ],
};

sub execute {
    my $self = shift;

    my $cmd = Genome::Model::Tools::Lims::SamplesForWorkOrder->create(
        work_order => $self->work_order,
        show => undef,
    );
    unless($cmd->execute) {
        die $self->error_message('Failed to get samples for the work order.');
    }

    my @wo_sample_ids = $cmd->_sample_ids;
    my @wo_samples = Genome::Sample->get(\@wo_sample_ids);

    my @model_samples;
    for my $m ($self->models) {
        push @model_samples, $m->tumor_sample if $m->tumor_sample;
        push @model_samples, $m->normal_sample if $m->normal_sample;
    }

    my %sample_matcher = map { $_->id => 0 } @model_samples;
    my @unmatched_wo_samples;
    for my $s (@wo_samples) {
        if(exists $sample_matcher{$s->id}) {
            $sample_matcher{$s->id}++;
        } else {
            push @unmatched_wo_samples, $s;
        }
    }

    my @unmatched_model_sample_ids = grep { not $sample_matcher{$_} } (keys %sample_matcher); #break this out to avoid NULL in-clause warning
    my @unmatched_model_samples;
    @unmatched_model_samples = Genome::Sample->get(\@unmatched_model_sample_ids)
        if @unmatched_model_sample_ids;

    $self->display_missing_sample_list('model samples', @unmatched_model_samples);
    $self->display_missing_sample_list('work order items', @unmatched_wo_samples);

    return 1; #!(scalar(@unmatched_model_samples) + scalar(@unmatched_wo_samples));
}

sub display_missing_sample_list {
    my $self = shift;
    my $type = shift;
    my @missing_samples = @_;

    if(scalar(@missing_samples) == 0) {
        $self->debug_message("All $type were accounted for.");
    } else {
        if($self->show) {
            $self->warning_message("These $type were not accounted for:");
            Genome::Sample::Command::List->execute(
                show => $self->show,
                filter => 'id:' . join('/', map($_->id, @missing_samples)),
            );
        } else {
            $self->warning_message("Some $type were not accounted for.");
        }
    }
}

1;
