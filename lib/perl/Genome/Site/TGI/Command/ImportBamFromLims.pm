package Genome::Site::TGI::Command::ImportBamFromLims;

use strict;
use warnings;

use File::Basename qw();
use File::Spec qw();
use Try::Tiny qw(try catch);

use Genome;

class Genome::Site::TGI::Command::ImportBamFromLims {
    is => 'Command::V2',
    has_input => {
        instrument_data => {
            is => 'Genome::InstrumentData::Solexa',
            doc => 'The data to look up in the LIMS',
            is_many => 1,
        },
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'The Analysis Project for which this instrument data is being imported',
        },
    },
    doc => 'Import BAM data into the GMS out of the LIMS system',
};

sub help_detail {
    return <<EOHELP
This command looks up the current BAM path in the LIMS system and copies it to an allocation.
EOHELP
}

sub execute {
    my $self = shift;

    for my $i ($self->instrument_data) {
        $self->_process_instrument_data($i);
        UR::Context->commit();
    }

    return 1;
}

sub _process_instrument_data {
    my $self = shift;
    my $data = shift;

    my $bam_path = $data->bam_path;
    unless ($bam_path) {
        $self->error_message('Skipping instrument data %s with no bam_path.', $data->__display_name__);
        return;
    }

    if (-e $bam_path) {
        $self->warning_message('Skipping instrument data %s because %s currently exists.', $data->__display_name__, $bam_path);
        return 1;
    }

    my @alloc = $data->disk_allocations;
    if (@alloc) {
        $self->error_message('Skipping instrument data %s because it already has allocated disk: %s', $data->__display_name__, join(" ", map $_->absolute_path, @alloc));
        return;
    }

    my $lims_path = $self->_resolve_lims_bam_path($data);
    unless ($lims_path) {
        $self->error_message('Skipping instrument data %s because no LIMS BAM path could be found.', $data->__display_name__);
        return;
    }

    if ($lims_path =~ m!^/gscarchive!) {
        $self->warning_message('Skipping instrument data %s because it appears to be in the old archive.', $data->__display_name__);
        return 1;
    }

    my $allocation = $self->_create_allocation($data);
    unless ($allocation) {
        $self->error_message('Failed to allocate space for instrument data %s.', $data->__display_name__);
        return;
    }


    try {
        my ($bam_file) = File::Basename::fileparse($data->bam_path);
        Genome::Sys->rsync_directory(
            source_directory => $lims_path,
            target_directory => $allocation->absolute_path,
        );

        Genome::Sys->shellcmd(
            cmd => ['chgrp', '-R', $self->_user_group, $allocation->absolute_path],
        );

        my $new_path = File::Spec->join($allocation->absolute_path, $bam_file);
        $data->bam_path($new_path);

        $allocation->reallocate;

        $self->status_message('Updated instrument data %s to path: %s.', $data->__display_name__, $new_path);
    }
    catch {
        my $error = $_;
        $allocation->deallocate;
        $self->error_message('Failed to unarchive instrument data %s. -- %s', $data->__display_name__, $error);
    }

    return 1;
}

sub _resolve_lims_bam_path {
    my $self = shift;
    my $data = shift;

    my $id = $data->id;

    my $docker_image = `lims-config docker_images.lims_perl_environment`;
    chomp $docker_image;

    my $guard = Genome::Config::set_env('lsb_sub_additional', "docker($docker_image)");
    my $cmd = [qw(db ii analysis_id), $data->id, qw(-mp get_disk_archive->archive_path)];

    my $log_allocation = Genome::Disk::Allocation->get(owner_class_name => $self->class);
    my $log_dir = $log_allocation->absolute_path;
    my $log_file = File::Spec->join($log_dir, $data->id);

    #not allowed to `docker run`, so `bsub` this query
    #can't nest interactive jobs, so write the output to a file and then read it in
    Genome::Sys::LSF::bsub::bsub(
        cmd => $cmd,
        queue => Genome::Config::get('lsf_queue_build_worker'),
        wait_for_completion => 1,
        log_file => $log_file,
    );

    my @data = Genome::Sys->read_file($log_file);
    unlink $log_file;

    my $path;
    while (!$path and @data) {
        my $next = shift @data;
        $path = $next if $next =~ m!^/gscmnt/!;
    }

    chomp $path if $path;
    return $path if -e $path;

    return;
}

sub _create_allocation {
    my $self = shift;
    my $data = shift;

    my %params = (
        disk_group_name => $self->_disk_group->disk_group_name,
        allocation_path => File::Spec->join('instrument_data',$data->id),
        kilobytes_requested => $data->calculate_alignment_estimated_kb_usage,
        owner_class_name => $data->class,
        owner_id => $data->id,
    );

    my $create_cmd = Genome::Disk::Command::Allocation::Create->create(%params);
    unless ($create_cmd->execute) {
        $self->error_message('Could not create allocation for instrument data: %s', $data->__display_name__);
        return;
    }

    return Genome::Disk::Allocation->get(allocation_path => $params{allocation_path});
}

sub _disk_group {
    my $self = shift;

    unless (exists $self->{_disk_group}) {
        $self->{_disk_group} = $self->_resolve_disk_group;
    }

    return $self->{_disk_group};
}

sub _resolve_disk_group {
    my $self = shift;

    my $anp = $self->analysis_project;
    my $guard = $anp->set_env;

    my $dg = Genome::Config::get('disk_group_alignments');

    return Genome::Disk::Group->get(disk_group_name => $dg);
}

sub _user_group {
    my $self = shift;

    unless (exists $self->{_user_group}) {
        $self->{_user_group} = $self->_resolve_user_group;
    }

    return $self->{_user_group};
}

sub _resolve_user_group {
    my $self = shift;

    my $anp = $self->analysis_project;
    my $guard = $anp->set_env;

    my $group = Genome::Config::get('sys_group');

    return $group;
}

1;
