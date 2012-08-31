package Genome::Model::MetagenomicComposition16s::Command::CreateFromInstrumentData; 

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicComposition16s::Command::CreateFromInstrumentData {
    is => 'Command::V2',
    has => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_input => 1,
            is_many => 1,
            require_user_verify => 1,
            where => [ sequencing_platform => [qw/ sanger 454 /] ],
            doc => 'Instrument data to use to create models. Resolved via text string on the command line.',
        },
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            is_input => 1,
            doc => 'Processing profile to use when creating models. Resolved via text string on the command line.',
        },
        project => {
            is => 'Genome::Project',
            is_input => 1,
            is_optional => 1,
            doc => 'Project to assign models. Resolved via text string on the command line.',
        },
        namer => {
            is => 'Perl',
            is_input => 1,
            is_optional => 1,
            doc => 'Regex or perl code to set the name of new models. The subject name is $s and the processing profile name is $pp.'
        },
        _models => { is_transient => 1, is_many => 1, is_optional => 1, },
    ],
};

sub help_brief { 
    return 'Create models from 454 instrument data.';
}

sub help_detail {
    return <<HELP;
This command will take instrument data, create models based on the subject of the instrument data and request a build (if needed). It can then add the models to an existing project.
HELP
}

sub execute {
    my $self = shift;

    $self->status_message('Create models from instrument data....');

    my $namer;
    if ( $self->namer ) {
        $self->status_message('Model namer: '.$self->namer);
        $namer = _wrap_perl_expr($self->namer);
        return if not $namer;
    }

    my $processing_profile = $self->processing_profile;
    if ( not $processing_profile ) {
        $self->error_message('No processing profile given!');
        return;
    }
    $self->status_message('Processing profile: '.$processing_profile->name);

    my @instrument_data = $self->instrument_data;
    if ( not @instrument_data ) {
        $self->error_message('No instrument data given!');
        return;
    }

    my %models;
    my @subject_regexp_to_skip = map { qr/$_/ } (qw/ ^nctrl$ ^n\-cntrl$ ^Pooled_Library /);
    my %metrics = (
        assigned => 0,
        skipped => 0,
        already_assigned => 0,
    );
    INST_DATA: for my $instrument_data ( @instrument_data ) {
        my $sample = $instrument_data->sample;
        for my $subject_regexp_to_skip ( @subject_regexp_to_skip ) {
            if ( $sample->name =~ $subject_regexp_to_skip ) {
                $metrics{skipped}++;
                next INST_DATA;
            }
        }
        my %model_params =  (
            subject => $instrument_data->sample,
            processing_profile => $processing_profile,
        );
        if ( $namer ) {
            my  $name = $namer->($instrument_data->sample->name, $processing_profile->name);
            if ( not $name ) {
                $self->error_message('Failed to get mode lname from namer! '.$self->namer);
                return;
            }
            $model_params{name} = $name;
        }
        my $model = Genome::Model->get(%model_params);
        if ( not $model ) {
            $model = Genome::Model->create(
                auto_assign_inst_data => 0,
                auto_build_alignments => 0,
                %model_params,
            );
            if ( not $model ) {
                $self->error_message('Failed to create model for subject ('.$model->subject.') and processing profile ('.$processing_profile->name.')');
                return;
            }
        }
        $models{$model->id} = $model;
        my $existing_input = $model->inputs(name => 'instrument_data', value => $instrument_data);
        if ( $existing_input) {
            $metrics{already_assigned}++;
            next;
        }
        $metrics{assigned}++;
        $model->add_instrument_data($instrument_data);
    }

    for my $model ( values %models ) {
        $self->status_message($model->name);
        my $build_needed = $model->build_needed;
        next if not $build_needed;
        $model->build_requested(1);
        $metrics{build_needed}++;
    }

    $self->status_message('Models created/found: '.keys(%models));
    $self->status_message('Subject name skips: '.$metrics{skipped});
    $self->status_message('Instrument data count: '.@instrument_data);
    $self->status_message('Instrument data assigned: '.$metrics{assigned});
    $self->status_message('Instrument data already assigned: '.$metrics{already_assigned});
    $self->status_message('Models needing builds: '.$metrics{build_needed});

    my $project = $self->project;
    if ( $project ) {
        $self->status_message('Adding models to project: '.$project->name);
        for my $model ( values %models ) {
            my $existing_part = $project->parts(entity => $model);
            next if $existing_part;
            $project->add_part(entity => $model);
        }
    }

    $self->_models([ values %models ]);
    $self->status_message('Create models from instrument data....DONE');

    return 1;
}

sub _wrap_perl_expr {
    my $expr = shift;
    Carp::confess('No model name expression given!') if not defined $expr;
    my $wrapped = qq|
        sub {
            my (\$s,\$pp) = \@_;
            my \$c = sub { $expr };
            return \$c->();
        }
    |;
    my $sub = eval $wrapped;
    if ($@) {
        Carp::confess("Error in code to set the name of models ($expr): $@");
    }
    return $sub;
}

1;

