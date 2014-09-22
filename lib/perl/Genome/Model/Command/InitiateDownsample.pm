package Genome::Model::Command::InitiateDownsample;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::InitiateDownsample {
    is => 'Genome::Command::Base',
    doc => "Initiate the downsampled import of all instrument data from a model group.",
    has_input => [
        models => {
            shell_args_position => 1,
            is => 'Genome::Model',
            is_many => 1,
            doc => 'the models to downsample',
        },
        upper_bound => {
            is => 'Float',
            doc => 'The upper bound of the property used to determine the downsample ratio.',
            example_values => ['0.50000','500000000'],
        },
        lower_bound => {
            is => 'Float',
            doc => 'The lower bound of the property used to determine the downsample ratio.',
            example_values => ['0.03125','1000000'],
        },
        property => {
            is => 'Text',
            doc => 'The property used to determine downsample ratios used.',
            valid_values => ['clusters','ratio'],
        },
        step_fraction => {
            is => 'Float',
            doc => 'The fraction (decimal value) of the property used in each downsample step.',
            is_optional => 1,
            default_value => 0.5,
        },
    ],
};

sub execute {
    my $self = shift;
    
    for my $model ($self->models) {
        $self->debug_message('Model: '. $model->id);
        my $model_count = $self->_resolve_model_count($model);
        my $ratios = $self->_resolve_model_ratios($model_count);
        my $cmd = 'genome instrument-data import generate-file-for-import --instrument-data='. join(',', map{$_->id} $model->instrument_data) .' --file='. $model->id .'.tsv --downsample-ratios='. join(',',@{$ratios});
        print $cmd ."\n";
    }
    return 1;
}

sub _resolve_model_ratios {
    my $self = shift;
    my $model_count = shift;

    my $start_ratio = ($self->upper_bound / $model_count);

    my @ratios;
    for ( my $ratio = $start_ratio; ($model_count * $ratio) >= $self->lower_bound; $ratio = ($ratio * $self->step_fraction) ) {
        my $ds_count = $model_count * $ratio;
        $self->debug_message($model_count .' * '. $ratio .' = '. $ds_count);
        push @ratios, $ratio;
    }
    return \@ratios;
} 

sub _resolve_model_count {
    my $self = shift;
    my $model = shift;

    my $model_count = 0;
    if ($self->property eq 'clusters') {
        my @instrument_data = $model->instrument_data;
        unless (@instrument_data) { die('No instrument data found for model: '. $model->id); }
        
        for my $instrument_data (@instrument_data) {
            $model_count += $instrument_data->clusters;
        }
    } elsif ($self->property eq 'ratio') {
        $model_count = 1;
    }
    unless ($model_count > $self->upper_bound) {
        die('Not enough '. $self->property .' found in model '. $model->id .'! Upper bound is \''. $self->upper_bound .'\' but only counted \''. $model_count .'\'');
    }
    return $model_count;
}

1;

