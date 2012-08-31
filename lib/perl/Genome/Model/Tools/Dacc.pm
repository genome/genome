package Genome::Model::Tools::Dacc;

use strict;
use warnings;

use Genome;

require Carp;
require File::Basename;
require File::Copy;
use Data::Dumper 'Dumper';
require File::Temp;
require IPC::Run;

class Genome::Model::Tools::Dacc {
    is  => 'Command',
    #is_abstract => 1,
    has => [
        dacc_directory => {
            is => 'Text',
            is_input => 1,
            shell_args_position => 1,
            doc => 'The directory on the DACC to use.',
        },
        user => { value => 'jmartin', is_constant => 1, },
        site => { value => 'aspera.hmpdacc.org', is_constant => 1, },
        user_and_site => {
            calculate_from => [qw/ user site /],
            calculate => q| return $user.'@'.$site; |,
        },
        dacc_remote_directory => {
            calculate_from => [qw/ user_and_site dacc_directory /],
            calculate => q| return $user_and_site.':'.$dacc_directory; |,
        },
        certificate => { value => '/gsc/scripts/share/certs/dacc/dacc.ppk', is_constant => 1, },
        ssh_key => { value => '/gsc/scripts/share/certs/dacc/dacc.sshkey', is_constant => 1, },
        base_command => {
            calculate_from => [qw/ certificate /],
            calculate => q| return 'ascp -Q -l100M -i '.$certificate; |,
        },
    ],
};

sub help_brief {
    return 'Manipulate files on the DACC';
}

sub help_detail {
    return help_brief();
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $dacc_directory = $self->dacc_directory;
    if ( not $dacc_directory ) {
        $self->error_message('No DACC directory given.');
        return;
    }

    if ( $dacc_directory !~ m#^/# ) {
        $dacc_directory = '/'.$dacc_directory;
    }

    if ( $dacc_directory !~ m/\/$/ ) {
        $dacc_directory .= '/';
    }

    $self->dacc_directory($dacc_directory);

    return $self;
}

sub temp_dir {
    my $self = shift;

    return $self->{_temp_dir} if $self->{_temp_dir};

    $self->{_temp_dir} = File::Temp::tempdir(CLEANUP => 1);

    return $self->{_temp_dir};
}

sub temp_ssh_key {
    my $self = shift;

    return $self->{_temp_ssh_key} if $self->{_temp_ssh_key};

    my $ssh_key = $self->ssh_key;
    my $temp_dir = $self->temp_dir;
    my $temp_ssh_key = $temp_dir.'/dacc.sshkey';
    $self->status_message("Temp ssh key: $temp_ssh_key");

    my $copy_ok = File::Copy::copy($ssh_key, $temp_ssh_key);
    if ( not $copy_ok or not -e $temp_ssh_key ) {
        Carp::confess('Failed to copy the ssh key to temp location.');
    }
    chmod 0400, $temp_ssh_key;

    $self->{_temp_ssh_key} = $temp_ssh_key;

    return $self->{_temp_ssh_key};
}

sub ls_dacc_directory {
    my $self = shift;

    my ($in, $out);
    my $user_and_site = $self->user_and_site;
    my $temp_ssh_key = $self->temp_ssh_key;
    my $harness = IPC::Run::harness([ 'ssh', '-i', $temp_ssh_key, $user_and_site ], \$in, \$out);
    $harness->pump until $out;

    $out = '';
    my $directory = $self->dacc_directory;
    $in = "ls -l $directory\n";
    $harness->pump until $out; 
    $harness->finish;
    $out =~ s/aspsh>\s+$//;

    return $out;
}

sub available_files_and_sizes {
    my $self = shift;

    my $out = $self->ls_dacc_directory;
    return if not $out;

    my %files_and_sizes;
    for my $line ( split("\n", $out) ) {
        my @tokens = split(/\s+/, $line);
        next if not $tokens[8];
        $files_and_sizes{ $tokens[8] } = $tokens[4];
    }

    return %files_and_sizes;
}

sub is_host_a_blade {
    my $self = shift;

    my $hostname = `hostname`;
    if ( not defined $hostname ) {
        $self->error_message('Cannot get hostname');
        return;
    }
    $self->status_message('Host is: '.$hostname);

    return $hostname =~ /blade/ ? 1 : 0;
}

sub is_running_in_lsf {
    my $self = shift;

    if ( $ENV{LSB_JOBID} ) {
        $self->status_message('LSF Job Id: '.$ENV{LSB_JOBID});
        return 1;
    }
    else {
        return;
    }
}

sub validate_running_in_lsf_and_on_a_blade {
    my $self = shift;

    if ( not $self->is_running_in_lsf ) {
        $self->error_message('Must run in LSF.');
        return;
    }

    if ( not $self->is_host_a_blade ) {
        $self->error_message('Must run on a blade.');
        return;
    }

    return 1;
}

sub _launch_to_lsf {
    my ($self, @params) = @_;

    Carp::confess('No params to launch to LSF') if not @params;

    $self->status_message("Launch to LSF");

    my $logging = '-u '. Genome::Config->user_email;
    my @rusage = $self->rusage;
    my $cmd = sprintf(
        'bsub -q long %s -R \'rusage[%s]\' gmt dacc %s %s %s', 
        $logging,
        join(',', @rusage),
        $self->command_name_brief,
        $self->dacc_directory,
        join(' ', @params),
    );

    my $rv = eval { Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        $self->error_message("Failed to launch to LSF: $@");
        return;
    }

    $self->status_message("Launch to LSF...OK");

    return 1;
}

1;

