package Genome::Model::Build::Command::ListVolumes;

use Genome;

class Genome::Model::Build::Command::ListVolumes {
    is => 'Command::V2',
    has_optional => [
        status => {
            is                  => 'String',
            doc => 'Build status to filter on',
            valid_values => Genome::Model::Build->__meta__->property('status')->valid_values,
        },
        build_ids => {
            is                  => 'String',
            is_many             => 1,
            shell_args_position => 1,
            doc => 'Build ID(s) to use.'
        },
    ],
};

sub help_detail {
    <<HELP
Prints the mount paths for volumes used by the given builds directly, and by
their software results and instrument data.

Builds may be filtered by build_id or by status.
HELP
}

sub execute {
    my $self = shift;

    print join("\n", $self->unique_volume_names),"\n";
    1;
}

sub unique_volume_names {
    my $self = shift;

    my %filters= $self->_build_filters_for_builds;
    unless (%filters) {
        $self->fatal_message('Refusing to run with no filters');
    }

    my $iter = Genome::Model::Build->create_iterator(%filters, -hints => ['disk_allocations', 'results.disk_allocations', 'instrument_data.disk_allocations']);

    my %seen_vols;
    while( my $guard = UR::Context::AutoUnloadPool->create()
             and
           my $build = $iter->next
    ) {
        foreach my $alloc ( $build->disk_allocations ) {
            $seen_vols{$alloc->mount_path} = undef;
        }

        foreach my $sr ( $build->results ) {
            $seen_vols{$_->mount_path} = undef foreach $sr->disk_allocations;
        }

        foreach my $instr_data ( $build->instrument_data ) {
            $seen_vols{$_->mount_path} = undef foreach $instr_data->disk_allocations;
        }
    }
    return keys %seen_vols;
}

sub _build_filters_for_builds {
    my $self = shift;

    my %filters;
    if ($self->status) {
        $filters{status} = $self->status;
    }

    if ($self->build_ids) {
        $filters{build_id} = [ $self->build_ids ];
    }

    return %filters;
}

1;
