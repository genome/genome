package Genome::InstrumentData::Command::FindAssociatedBuilds;

use warnings;
use strict;

use above "Genome";

class Genome::InstrumentData::Command::FindAssociatedBuilds {
    is => 'Command::V2',
    has => {
        instrument_data => {
            is => 'Genome::InstrumentData',
            shell_args_position => 1,
            doc => 'The instrument-data of interest.',
        },
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
            is_optional => 1,
            is_transient => 1,
        },
        max_related_iterations => {
            is_constant => 1,
            default_value => 20,
        }
    },
};

sub help_detail {
    return 'Find the builds that use this instrument-data, as well as builds that use those builds as inputs, ect.';
}

sub execute {
    my $self = shift;

    my $instrument_data = $self->instrument_data;
    my @builds = _find_builds($instrument_data, $self->max_related_iterations);

    for my $build (@builds) {
        $self->status_message($build->id);
    }

    $self->builds(\@builds);
    return 1;
}

sub _find_builds {
    my $instrument_data = shift;
    my $max_related_iterations = shift;

    my %builds;
    my @models = $instrument_data->models;
    for my $model (@models) {
        for my $build ($model->builds) {
            $builds{$build->id} = $build;
        }
    }

    my $prev_num_builds = 0;
    my $count = 0;
    while ($prev_num_builds != keys(%builds)) {
        if ($count > $max_related_iterations) {
            die "Exceeded maximum number of related build iterations."
        }
        $prev_num_builds = keys(%builds);
        $count += 1;

        my @related_builds = _find_related_builds(keys(%builds));
        for my $related_build (@related_builds) {
            $builds{$related_build->id} = $related_build;
        }
    }
    return values(%builds);
}

sub _find_related_builds {
    my @build_ids = @_;
    return unless @build_ids;

    my @inputs = Genome::Model::Build::Input->get('value_id in' => \@build_ids);
    my %related_builds;
    for my $input (@inputs) {
        if ($input->value->isa('Genome::Model::Build')) {
            my $related_build = $input->build;
            $related_builds{$related_build->id} = $related_build;
        }

    }
    return values(%related_builds);
}
