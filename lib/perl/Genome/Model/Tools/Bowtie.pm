package Genome::Model::Tools::Bowtie;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use File::Temp;
use IPC::Cmd;

my $BOWTIE_DEFAULT = '0.12.5';

class Genome::Model::Tools::Bowtie {
    is => ['Command'],
    has_optional => [
        use_version => {
                    is    => 'string',
                    doc   => 'version of Bowtie application to use',
                    default_value => $BOWTIE_DEFAULT
                },
        _tmp_dir => {
                    is => 'string',
                    doc => 'a temporary directory for storing files',
                },
    ],
    doc => 'tools to work with the Bowtie aliger'
};

sub help_detail {                           # This is what the user will see with --help <---
    return <<EOS

EOS
}

my %BOWTIE_VERSIONS = (
    '0.12.7' => $ENV{GENOME_SW} . '/bowtie/bowtie-0.12.7',
    '0.12.5' => $ENV{GENOME_SW} . '/bowtie/bowtie-0.12.5',
    '0.12.1' => $ENV{GENOME_SW} . '/bowtie/bowtie-0.12.1',
    '0.10.0.2' => $ENV{GENOME_SW} . '/bowtie/bowtie-0.10.0.2',
    '0.9.9.2' => $ENV{GENOME_SW} . '/bowtie/bowtie-0.9.9.2',
    '0.9.8' => $ENV{GENOME_SW} . '/bowtie/bowtie-0.9.8',
    '0.9.4' => $ENV{GENOME_SW} . '/bowtie/bowtie-0.9.4',
);

sub path_for_bowtie_version {
    my ($class, $version, $subcommand) = @_;
    $version ||= $BOWTIE_DEFAULT;
    my $path = $BOWTIE_VERSIONS{$version};
    if (defined($path)) {
        if (Genome::Config->arch_os =~ /64/) {
            $path .= '-64/bowtie';
        }
        elsif (Genome::Config->arch_os =~ /i686/) {
            # There is a x86 version at /gsc/pkg... but memory likely will be insufficient.
            # In any case, the Genome::InstrumentData::AlignmentResult::Bowtie will enforce X86_64
            $path .= '/bowtie';
        }
        if($subcommand && $subcommand ne ''){
            $path .= "-$subcommand";
        }
        return $path;
    } else {
        my $command = 'bowtie';
        if($version =~ /^2/){
            $command .= "2";
        }
        if($subcommand && $subcommand ne ''){
            $command .= "-$subcommand";
        }
        $path = IPC::Cmd::can_run($command . $version);
        return $path if $path;
    }
    die 'No path found for bowtie version: '.$version;
}

#some programs (e.g. tophat) use `which` to decide where to call bowtie
#and so we have to prepend to PATH to make sure they call the right version
sub path_variable_for_bowtie_version {
    my ($class, $version) = @_;

    my $path = $class->path_for_bowtie_version($version);
    unless($path) { die 'No path found for bowtie version: '.$version; }

    if($version =~ /^2/) {
        $path =~ s!^/usr/bin/bowtie2!/usr/lib/bowtie!;
        $path .= '/bin';
    } else {
        $path =~ s'/bowtie$'';
    }

    return $path;
}

sub default_bowtie_version {
    die "default bowtie version: $BOWTIE_DEFAULT is not valid" unless $BOWTIE_VERSIONS{$BOWTIE_DEFAULT};
    return $BOWTIE_DEFAULT;
}

sub default_version { return default_bowtie_version; }

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless (Genome::Config->arch_os =~ /64/) {
        $self->error_message('Most Bowtie tools must be run from 64-bit architecture');
        return;
    }
    unless ($self->temp_directory) {
        my $base_temp_directory = Genome::Sys->base_temp_directory;
        my $temp_dir = File::Temp::tempdir($base_temp_directory .'/Bowtie-XXXX', CLEANUP => 1);
        Genome::Sys->create_directory($temp_dir);
        $self->_tmp_dir($temp_dir);
    }
    return $self;
}


1;

