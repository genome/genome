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

    my %filters= $self->_build_filters_for_builds;
    unless (%filters) {
        $self->fatal_message('Refusing to run with no filters');
    }

    my $iter = Genome::Model::Build->create_iterator(%filters, -hints => ['disk_allocations', 'results.disk_allocations', 'instrument_data.disk_allocations']);

    $self->_print_volumes_for_builds($iter);
    1;
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

sub _print_volumes_for_builds {
    my($self, $iter) = @_;

    my %seen_vols;
    my $print_new_vols = sub {
        my $path = shift;
        unless (exists $seen_vols{$path}) {
            $seen_vols{$path} = undef;
            print $path,"\n";
        }
    };

    BUILD:
    while(1) {
        my $guard = UR::Context::AutoUnloadPool->create();
        my $build = $iter->next();
        last BUILD unless $build;

        foreach my $alloc ( $build->disk_allocations ) {
            $print_new_vols->($alloc->mount_path);
        }

        foreach my $sr ( $build->results ) {
            $print_new_vols->($_->mount_path) foreach $sr->disk_allocations;
        }

        foreach my $instr_data ( $build->instrument_data ) {
            $print_new_vols->($_->mount_path) foreach $instr_data->disk_allocations;
        }
    }
}

1;
