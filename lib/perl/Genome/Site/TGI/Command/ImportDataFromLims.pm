package Genome::Site::TGI::Command::ImportDataFromLims;

use strict;
use warnings;

use File::Basename qw();
use File::Spec qw();
use Sub::Override;
use Try::Tiny qw(try catch);

use Genome;

class Genome::Site::TGI::Command::ImportDataFromLims {
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
    has_optional_param => {
        remote_host => {
            is => 'Text',
            doc => 'ssh-enabled remote host for transfer--required for remote volumes',
        },
        remote_user => {
            is => 'Text',
            doc => 'remote username for ssh--required for remote volumes',
        },
        remote_port => {
            is => 'Number',
            doc => 'remote port to use for SSH--required for remote volumes',
            default_value => 22,
        },
        authorization_key => {
            is => 'Text',
            doc => 'path to the SSH key to use for rsync transfer--required for remote volumes',
        },
    },
    doc => 'Import instrument data files into the GMS out of the LIMS system',
};

sub help_detail {
    return <<EOHELP
This command looks up the current data path in the LIMS system and copies it to an allocation.
EOHELP
}

sub execute {
    my $self = shift;

    my $remote_params_defined =  grep { defined $self->$_ } (qw(remote_host remote_user remote_port authorization_key));
    if ($remote_params_defined != 0 and $remote_params_defined != 4) {
        $self->fatal_message('All four remote parameters must be defined if any are, but found only %s defined.', $remote_params_defined);
    }

    for my $i ($self->instrument_data) {
        $self->_process_instrument_data($i);
        UR::Context->commit();
    }

    return 1;
}

sub _process_instrument_data {
    my $self = shift;
    my $data = shift;

    my @alloc = $data->disk_allocations;
    if (@alloc) {
        $self->error_message('Skipping instrument data %s because it already has allocated disk: %s', $data->__display_name__, join(" ", map $_->absolute_path, @alloc));
        return;
    }

    my $lims_path = $self->_resolve_lims_path($data);
    unless ($lims_path) {
        $self->error_message('Skipping instrument data %s because no LIMS path could be found.', $data->__display_name__);
        return;
    }

    $self->debug_message('Found LIMS path: %s', $lims_path);

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
        my ($bam_file, $lims_source_dir);
        if ($lims_path =~ /\.bam$/) {
            ($bam_file, $lims_source_dir) = File::Basename::fileparse($lims_path);
        } elsif (-d $lims_path) {
            $lims_source_dir = $lims_path;
        } else {
            die $self->error_message('Unknown LIMS filetype: %s', $lims_path);
        }

        if ($allocation->volume->is_remote_volume) {
            my $source_path = $lims_source_dir . '/';
            my $destination_path = $allocation->absolute_path;
            my $user = $self->remote_user;
            my $host = $self->remote_host;

            my @ssh_opts = ('ssh', '-p', $self->remote_port, '-i', $self->authorization_key, '-o', 'StrictHostKeyChecking=no');

            my $remote_destination = sprintf('%s@%s:%s', $user, $host, $destination_path);

            #first ensure directory exists as target for rsync
            Genome::Sys->shellcmd(
                cmd => [@ssh_opts, sprintf('%s@%s', $user, $host), 'mkdir', '-p', $destination_path],
            );
            Genome::Sys->shellcmd(
                cmd => ['rsync', '-av', '-e', join(' ', @ssh_opts), $source_path, $remote_destination],
            );
        } else {
            Genome::Sys->rsync_directory(
                source_directory => $lims_source_dir,
                target_directory => $allocation->absolute_path,
                chmod => 'Dug=rx,Fug=r',
                chown => ':' . $self->_user_group,
            );
        }

        if ($bam_file) {
            my $new_path = File::Spec->join($allocation->absolute_path, $bam_file);
            $data->bam_path($new_path);
            $self->status_message('Updated instrument data %s to path: %s.', $data->__display_name__, $new_path);
        } else {
            $self->status_message('Data imported for %s to path: %s.', $data->__display_name__, $allocation->absolute_path);
        }

        $allocation->reallocate;
    }
    catch {
        my $error = $_;
        $allocation->deallocate;
        $self->error_message('Failed to unarchive instrument data %s. -- %s', $data->__display_name__, $error);
    };

    return 1;
}

sub _resolve_lims_path {
    my $self = shift;
    my $data = shift;

    my $id = $data->id;

    my $docker_image = `lims-config docker_images.lims_perl_environment`;
    chomp $docker_image;

    my $guard = Genome::Config::set_env('lsb_sub_additional', "docker($docker_image)");
    my $cmd = [qw(db ii analysis_id), $data->id, qw(-mp absolute_path)];

    local $ENV{LSF_DOCKER_PRESERVE_ENVIRONMENT} = 'false';
    local $ENV{LSB_DOCKER_MOUNT_GSC} = 'false';
    local $ENV{LSF_DOCKER_VOLUMES} = undef; #lims-env breaks if /gsc is present.

    my $log_allocation = Genome::Disk::Allocation->get(owner_class_name => $self->class);
    my $log_dir = $log_allocation->absolute_path;
    my $log_file = File::Spec->join($log_dir, $data->id);

    #not allowed to `docker run`, so `bsub` this query
    #can't nest interactive jobs, so write the output to a file and then read it in
    Genome::Sys->bsub_and_wait(
        cmd => $cmd,
        queue => Genome::Config::get('lsf_queue_build_worker'),
        user_group => Genome::Config::get('lsf_user_group'),
        log_file => $log_file,
    );

    my @data = Genome::Sys->read_file($log_file);
    unlink $log_file;

    my $path;
    while (!$path and @data) {
        my $next = shift @data;
        $path = $next if ($next =~ m!^/gscmnt/! and $next !~ m!^/storage./!);
    }

    chomp $path if $path;
    return $path if -e $path;

    return;
}

sub _create_allocation {
    my $self = shift;
    my $data = shift;

    my %params = (
        disk_group_name => $self->_disk_group,
        allocation_path => File::Spec->join('instrument_data',$data->id),
        kilobytes_requested => $data->calculate_alignment_estimated_kb_usage,
        owner_class_name => $data->class,
        owner_id => $data->id,
    );

    if ($self->remote_host) {
        my @volumes = Genome::Disk::Detail::Allocation::Creator->get_candidate_volumes(disk_group_name => $params{disk_group_name});
        if (my @bad = grep !$_->is_remote_volume, @volumes) {
            $self->fatal_message('Remote options specified but disk group contains non-remote volumes: %s', join(',', map $_->mount_path, @bad));
        }

        my $creator_override = Sub::Override->new('Genome::Disk::Detail::Allocation::Creator::create_directory_or_delete_allocation', sub { return 1; }); #we're handling directory creation remotely
        my $volume_override = Sub::Override->new('Genome::Disk::Volume::is_mounted', sub { return 1; }); #mounted remotely

        Genome::Disk::Allocation->_create(
            id => Genome::Disk::Allocation->__meta__->autogenerate_new_object_id,
            %params,
        );
    } else {
        my $create_cmd = Genome::Disk::Command::Allocation::Create->create(%params);
        unless ($create_cmd->execute) {
            $self->error_message('Could not create allocation for instrument data: %s', $data->__display_name__);
            return;
        }
    }

    return Genome::Disk::Allocation->get(allocation_path => $params{allocation_path});
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

    return $dg;
}

1;
