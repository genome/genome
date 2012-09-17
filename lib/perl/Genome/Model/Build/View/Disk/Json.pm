package Genome::Model::Build::View::Disk::Json;

use strict;
use warnings;
require UR;

use XML::Simple;
use JSON;

class Genome::Model::Build::View::Disk::Json {
    is => 'UR::Object::View::Default::Json',
    has_constant => [
        toolkit     => { value => 'json' },
    ],
    has_optional => [
        encode_options => { is => 'ARRAY', default_value => ['ascii', 'pretty', 'allow_nonref', 'canonical'], doc => 'Options to enable on the JSON object; see the documentation for the JSON Perl module' },
    ],
};


sub _generate_content {
    my $self = shift;

    my $build = $self->subject();

    if (!$build) {
        Carp::confess('This JSON view couldnt get the subject of the view. class='
                    , $self->subject_class_name
                    . ' id='
                    . $self->subject_id);
    }

    my $json = {};
    my @allocs = $build->all_allocations();
    my $total = {};

    for my $a (@allocs) {
        my $j = {};

        for my $property (qw/id kilobytes_requested absolute_path owner_class_name owner_id disk_group_name archivable is_archived/) {
            $j->{$property} = $a->$property;
        }

        push @{$json->{'allocations'}}, $j;

        $total->{'kilobytes_requested'} += $a->kilobytes_requested;
        $total->{'archivable'} += $a->archivable ? 1 : 0;
        $total->{'archived'} += $a->is_archived ? 1 : 0;
    }

    $json->{'total_allocations'} = scalar(@allocs);
    $json->{'total_archivable'} = $total->{'archivable'};
    $json->{'total_archived'} = $total->{'archived'};

    my $total_kilobytes = int($total->{'kilobytes_requested'});
    my $total_megabytes = int(($total_kilobytes / 1024) + 0.5);
    my $total_gigabytes = int(($total_kilobytes / 1024 / 1024) + 0.5);
    $json->{'total_kilobytes'} = $total_kilobytes;
    $json->{'total_megabytes'} = $total_megabytes;
    $json->{'total_gigabytes'} = $total_gigabytes;
    $json->{'total_diskspace'} = $total_gigabytes > 1 ? $total_gigabytes . 'GB'
                                : $total_megabytes > 1 ? $total_megabytes . 'MB'
                                : $total_kilobytes . 'KB';

    return $self->_json->encode($json);
}

1;
