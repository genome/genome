#!/usr/bin/env lims-perl

use strict;
use warnings;

use Data::Dumper;
use Carp qw/confess/;

use GSCApp;
App->init();

unless( defined($ARGV[0]) and defined $ARGV[1] ) {
    die "usage: $0 <analysis_id> <mount_path>";
}

my $ii = GSC::IndexIllumina->get(analysis_id => $ARGV[0]);
unless ($ii) {
    Carp::confess('Failed to find the LIMS IndexIllumina by analysis_id: '. $ARGV[0]);
}

my $disk_volume = GSC::DiskVolume->get(mount_path => $ARGV[1]);
unless ($disk_volume) {
    Carp::confess('Failed to find LIMS DiskVolume by mount path: '. $ARGV[1]);
}
unless ($disk_volume->can_allocate) {
    Carp::confess('Unable to use unallocatable disk volume: '. $disk_volume->mount_path);
}

my $csf_id = $ii->get_copy_sequence_files_pse->pse_id;
unless ($csf_id) {
    Carp::confess('Failed to find LIMS CopySequenceFiles pse_id for IndexIllumina analysis_id: '. $ii->analysis_id);
}

my $archive = GSC::DiskArchive->get(
    allocation_path => { operator => 'LIKE', value => '%csf_'. $csf_id },
    status => 'active',
);
unless ($archive) {
    Carp::confess('Failed to find LIMS DiskArchive for CopySequenceFiles pse_id: '. $csf_id);
}

print 'disk_unarchive --keep-archive --disk-volume '. $disk_volume->mount_path .' '. $archive->source_path ."\n";

exit;
