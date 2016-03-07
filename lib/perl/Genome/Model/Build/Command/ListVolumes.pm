package Genome::Model::Build::Command::ListVolumes;

use Genome;
use strict;
use warnings;

class Genome::Disk::InstrumentDataAllocation {
    is => 'Genome::Disk::Allocation',
    has => [
        owner => { is => 'Genome::InstrumentData', id_by => 'owner_id' },
    ],
};

class Genome::BuildInstrumentData {
    is => 'Genome::InstrumentData',
    has_many => [
        disk_allocations => { is => 'Genome::Disk::InstrumentDataAllocation', reverse_as => 'owner' },
        build_inputs => { is => 'Genome::Model::Build::Input', reverse_as => 'value_inst_data' },
        builds => { is => 'Genome::Model::Build', via => 'build_inputs', to => 'build' },
    ],
};

class Genome::Disk::BuildAllocation {
    is => 'Genome::Disk::Allocation',
    has => [
        build => { is => 'Genome::Model::Build', id_by => 'owner_id' },
        software_result => { is => 'Genome::SoftwareResult', id_by => 'owner_id' },
        instrument_data => { is_many => 1, is => 'Genome::BuildInstrumentData', reverse_as => 'disk_allocations' },
        instrument_data_builds => { is => 'Genome::Model::Build', is_many => 1, via => 'instrument_data', to => 'builds' },
    ],
};

class Genome::Model::Build::Command::ListVolumes {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name => {
            is_constant => 1,
            value => 'Genome::Model::Build',
        },
    ],
    has_input => [
        filter_class_name => {
            is_constant => 1,
            value => 'Genome::Model::Build'
        },
        show_class_name => {
            is_constant => 1,
            value => 'Genome::Disk::BuildAllocation::Set',
        },
        show => {
            is => 'Text',
            is_optional => 1,
            default_value => 'mount_path',
            doc => 'Specify which columns to show, in order.  Prefix with "+" or "^" to append/prepend to the default list.',
        },
    ],
};

sub _resolve_boolexpr {
    my $self = shift;

    no warnings 'redefine';
    # initially filter on Builds
    local *Genome::Model::Build::Command::ListVolumes::subject_class_name = \&filter_class_name;
    # no hints when getting these builds
    local *Genome::Model::Build::Command::ListVolumes::_hint_string = sub { };
    $self->super_can('_resolve_boolexpr')->($self);
}

sub resolve_show_column_names {
    my $self = shift;

    no warnings 'redefine';
    local *Genome::Model::Build::Command::ListVolumes::subject_class_name = \&show_class_name;
    $self->super_can('resolve_show_column_names')->($self);
}

sub create_iterator_for_results_from_boolexpr {
    my($self, $bx) = @_;

    my @fields = $self->resolve_show_column_names;

    my @build_ids = map { $_->id } $self->filter_class_name->get($bx);
    my @mount_paths = Genome::Disk::BuildAllocation->get(
        -or => [
            [ owner_id => \@build_ids ],
            [ 'software_result.builds.id' => \@build_ids ],
            [ 'instrument_data_builds.id' => \@build_ids ],
        ],
        -group_by => \@fields,
    );

    my(%seen, @results);
    foreach my $set (@mount_paths) {
        my $key = join(':', map { $set->$_ } @fields);
        next if $seen{$key}++;
        push @results, $set;
    }

    return UR::Object::Iterator->create_for_list(@results);
}

1;
