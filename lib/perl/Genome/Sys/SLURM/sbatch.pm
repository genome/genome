package Genome::Sys::SLURM::sbatch;

use strict;
use warnings;

use Genome::Sys;
use Exporter qw(import);
use Params::Validate qw(:types);

our @EXPORT = qw(sbatch);
our @EXPORT_OK = qw(sbatch);

sub run {
    my $executable = shift;
    my @args = args_builder(@_);

    my %args = @_;

    if (ref($executable) ne 'ARRAY') {
        $executable = [$executable];
    }

    my $docker_image = Genome::Config::get('lsb_sub_additional') || $ENV{LSB_SUB_ADDITIONAL} || $ENV{PYXIS_CONTAINER_IMAGE};
    if( (my $image_name) = $docker_image =~ /^docker\d?\((.*?)\)/) {
        $docker_image = $image_name;
    }
    local $ENV{LSB_SUB_ADDITIONAL} = $docker_image;
    local $ENV{PYXIS_CONTAINER_IMAGE} = $docker_image;

    if (!defined $ENV{LSF_DOCKER_PRESERVE_ENVIRONMENT} or $ENV{LSF_DOCKER_PRESERVE_ENVIRONMENT} eq 'true') {
        unshift @args,
            '--container-env=PATH'; #most env is getting preserved automatically but this one needs special care
    }

    my $docker_volumes = $ENV{LSF_DOCKER_VOLUMES} || $ENV{PYXIS_CONTAINER_MOUNTS} || $ENV{APPTAINER_BIND} || '';
    $docker_volumes =~ s/ /,/g;
    if (my $config_docker_volumes = Genome::Config::get('docker_volumes')) {
        $config_docker_volumes =~ s/ /,/g;
        if ($docker_volumes) {
            $docker_volumes .= ',' . $config_docker_volumes;
        } else {
            $docker_volumes = $config_docker_volumes;
        }
    }
    local $ENV{LSF_DOCKER_VOLUMES} = $docker_volumes if $docker_volumes;
    local $ENV{PYXIS_CONTAINER_MOUNTS} = $docker_volumes if $docker_volumes;

    if (exists $args{interactive} and $args{interactive}) {
        if($executable->[0] eq 'sbatch') {
            $executable->[0] = 'srun';
        }
        Genome::Sys->shellcmd(
            cmd => [@$executable, @args],
        );

        return -1; #not a valid job ID.
    } else {
        my @output = Genome::Sys->capture(@$executable, @args);

        my $job_id = ($output[-1] =~ /^Submitted batch job (\d+)/)[0];
        unless ($job_id) {
            die "Could not get job id from sbatch output!";
        }

        return $job_id;
    }
}

sub sbatch {
    return run('sbatch', @_);
}

sub args_builder {
    my %args = _args(@_);
    my $command = delete $args{cmd};
    my $resource_string = delete $args{resource_string};

    #remove LSF-specific flags with no equivalent
    for my $skipped_key ('job_group', 'never_rerunnable') {
        delete $args{$skipped_key};
    }

    my @args;
    for my $flag (_flags()) {
        if ($args{$flag}) {
            push @args, _option_mapper($flag, \%args);
        }
        delete $args{$flag};
    }
    if (keys %args) {
        push @args, map { _option_mapper($_, \%args), $args{$_} } sort keys %args;
    }

    if (ref($command) ne 'ARRAY') {
        $command = [$command];
    }

    unless ($args{interactive}) {
        my @propagate_conf;
        if (defined $ENV{SLURM_CONF}) {
            push @propagate_conf, 'SLURM_CONF="' . $ENV{SLURM_CONF} . '"';
        }
        $command = ['--wrap', join(' ', @propagate_conf, @$command)];
    }

    if ($resource_string) {
        my $parsed = _parse_slurm_params($resource_string);
        if (exists $parsed->{options}{partition}) {
            push @args, '--partition', $parsed->{options}{partition};
        }
        if (exists $parsed->{options}{cpus}) {
            push @args, '--cpus-per-task', $parsed->{options}{cpus};
        }
        if (exists $parsed->{options}{mem}) {
            push @args, '--mem', $parsed->{options}{mem};
        }
    }

    return (@args, @$command);
}

sub _args {
    return Params::Validate::validate(@_, _args_spec());
}

sub _option_mapper {
    my $arg = shift;
    my $args_ref = shift;
    my $spec = _args_spec();
    if (not exists $spec->{$arg}{option_flag}) {
        die qq(Could not find option flag for '$arg' in args spec);
    }
    my $flag = $spec->{$arg}{option_flag};
    if (ref($flag) eq 'CODE') {
        return $flag->($arg, $args_ref);
    }
    return $flag;
}

sub _flags {
    my $spec = _args_spec();
    return grep { $spec->{$_}{type} == BOOLEAN } sort keys %{$spec};
}

sub _args_spec {
    return {
        rerunnable => {
            optional => 1,
            type => BOOLEAN,
            option_flag => '--requeue',
        },
        never_rerunnable => {
            optional => 1,
            type => BOOLEAN,
            option_flag => sub { return () },
        },
        hold_job => {
            optional => 1,
            option_flag => '--hold',
            type => BOOLEAN,
        },
        send_job_report => {
            optional => 1,
            option_flag => '--mail-type=END',
            type => BOOLEAN,
        },
        queue => {
            optional => 1,
            option_flag => '--partition',
            type => SCALAR,
            callbacks => {
                'valid SLURM partition' => \&_valid_slurm_partition,
            },
        },
        project => {
            optional => 1,
            option_flag => '--comment',
            type => SCALAR,
        },
        email => {
            optional => 1,
            option_flag => '--mail-user',
            type => SCALAR,
        },
        err_file => {
            optional => 1,
            option_flag => '--error',
            type => SCALAR,
        },
        log_file => {
            optional => 1,
            option_flag => '--output',
            type => SCALAR,
        },
        job_group => {
            optional => 1,
            option_flag => sub { return () },
            type => SCALAR,
        },
        user_group => {
            optional => 1,
            option_flag => '--account',
            type => SCALAR,
        },
        job_name => {
            optional => 1,
            option_flag => '--job-name',
            type => SCALAR,
        },
        depend_on => {
            optional => 1,
            option_flag => sub {
                my ($arg, $args_ref) = @_;
                my $lsf_dep = $args_ref->{depend_on};
                my $slurm_dep = _translate_dependency($lsf_dep);
                return ('--dependency', $slurm_dep);
            },
            type => SCALAR,
        },
        wait_for_completion => {
            optional => 1,
            option_flag => '--wait',
            type => BOOLEAN,
        },
        interactive => {
            optional => 1,
            option_flag => '--pty',
            type => BOOLEAN,
        },
        post_exec_cmd => {
            optional => 1,
            option_flag => '--epilog',
            type => SCALAR,
        },
        cmd => {
            type => SCALAR | ARRAYREF,
        },
        resource_string => {
            optional => 1,
            type => SCALAR,
        },
    };
}

sub _translate_dependency {
    my $lsf_dep = shift;

    my $slurm_dep = $lsf_dep;

    $slurm_dep =~ s/ended\((\d+)\)/after:$1/g;
    $slurm_dep =~ s/done\((\d+)\)/afterok:$1/g;
    $slurm_dep =~ s/started\((\d+)\)/afterany:$1/g;
    #there are other types of dependency that we're ignoring at the moment.

    if($slurm_dep =~ m/\s*&&\s*/ or $slurm_dep =~ m/\s*\|\|\s*/) {
        die 'complex LSF dependency translation not supported';
    }

    return $slurm_dep;
}

sub _parse_slurm_params {
    my $resource_string = shift;

    my %result = (options => {});

    #We're going to ignore the '-R' in LSF resource strings here. The short options in the below
    #match the parts we care about.
    if ($resource_string =~ /-q\s+(\S+)/) {
        $result{options}{partition} = $1;
    }
    if ($resource_string =~ /-n\s+(\d+)/) {
        $result{options}{cpus} = $1;
    }
    if ($resource_string =~ /-M\s+(\S+)/) {
        my $val = $1;
        if ($val =~ /\d$/) {
            $val .= 'K'; #memlimit in LSF defaulted to KB
        }
        $result{options}{mem} = $val;
    }

    return \%result;
}

sub _valid_slurm_partition {
    my $requested_partition = shift;

    system(qq(scontrol show partition $requested_partition 1> /dev/null 2>&1));
    return $? == 0;
}

1;
