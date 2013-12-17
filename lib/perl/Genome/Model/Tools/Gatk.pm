package Genome::Model::Tools::Gatk;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use File::Temp;
use POSIX qw(floor);

my $DEFAULT_VERSION = '5336';
my $GATK_BASE = 'GenomeAnalysisTK';
my $GATK_COMMAND = "$GATK_BASE.jar";

class Genome::Model::Tools::Gatk {
    is => ['Command'],
    has_input => [
        version => {
            is    => 'string',
            doc   => 'version of Gatk application to use',
            default => $DEFAULT_VERSION,
        },
        max_memory => {
            # Accessor is overridden so that it can be limited based on available memory or LSF limit.
            is => 'Text',
            doc => 'The maximum memory (GB) to use when running Java VM. Limited to environmental constraints.',
            default => '4',
        },
    ],
    has_optional => [
        tmp_dir => {
            is => 'Text',
            doc => 'Temporary directory for Java.',
        },
    ],
    has_param => [
        lsf_resource => {
            default => '-M 16777216 rusage[mem=16384] select[type==LINUX64 & mem > 16384] span[hosts=1]',
        },
    ],
};

sub help_brief {
    "tools to work with Gatk output"
}

sub help_detail {                           # This is what the user will see with --help <---
    return <<EOS

EOS
}

sub gatk_versions {
    my %GATK_VERSIONS = (
        'v1' => $ENV{GENOME_SW} . '/gatk/GenomeAnalysisTK-1.0.5336/' . $GATK_COMMAND, # This is temporary... "v1" is what is in the first set of somatic-variation processing profiles
        '2986' => $ENV{GENOME_SW} . '/gatk/GenomeAnalysisTK-1.0.2986/' . $GATK_COMMAND,
        '3362' => $ENV{GENOME_SW} . '/gatk/GenomeAnalysisTK-1.0.3362/' . $GATK_COMMAND,
        '3362P' => $ENV{GENOME_SW} . '/gatk/GenomeAnalysisTK-1.0.3362P/' . $GATK_COMMAND,
        '3423' => $ENV{GENOME_SW} . '/gatk/GenomeAnalysisTK-1.0.3423/' . $GATK_COMMAND,
        '3471' => $ENV{GENOME_SW} . '/gatk/GenomeAnalysisTK-1.0.3471/' . $GATK_COMMAND,
        '4168' => $ENV{GENOME_SW} . '/gatk/GenomeAnalysisTK-1.0.4168/' . $GATK_COMMAND,
        '5336' => $ENV{GENOME_SW} . '/gatk/GenomeAnalysisTK-1.0.5336/' . $GATK_COMMAND,
        '5777' => $ENV{GENOME_SW} . '/gatk/GenomeAnalysisTK-1.0.5777/' . $GATK_COMMAND,
        '2.4' => $ENV{GENOME_SW} . '/gatk/GenomeAnalysisTK-2.4/' . $GATK_COMMAND,
        #'2.4' => Genome::Sys->jar_path($GATK_BASE, "2.4"),
    );
    return %GATK_VERSIONS;
}

our @legacy_versions = qw(v1 2986 3362 3362P 3423 3471 4168 5336 5777);

sub is_legacy_version {
    my $self = shift;
    my $version = shift;
    return grep {$_ eq $version} @legacy_versions;
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    unless (Genome::Config->arch_os =~ /64/) {
        $self->error_message('We recommend running GATK from 64-bit architecture');
        return;
    }

    if ( not $self->tmp_dir ) { 
        my $tempdir = Genome::Sys->create_temp_directory;
        $self->tmp_dir($tempdir);
    }

    return $self;
}

sub gatk_path {
    my $self = $_[0];
    return $self->path_for_gatk_version($self->version);
}

sub available_gatk_versions {
    my $self = shift;
    my %versions = $self->gatk_versions;
    return keys %versions;
}

sub path_for_gatk_version {
    my $class = shift;
    my $version = shift;
    my %versions = $class->gatk_versions;
    if (defined $versions{$version}) {
        return $versions{$version};
    }
    die('No path for gatk version '. $version);
}

sub default_gatk_version {
    my $class = shift;
    my %versions = $class->gatk_versions;
    die "default gatk version: $DEFAULT_VERSION is not valid" unless $versions{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}

sub has_version {
    my $self = shift;
    my $version = shift;
    unless(defined($version)){
        $version = $self->version;
    }
    my %versions = $self->gatk_versions;
    if(exists($versions{$version})){
        return 1;
    }
    return 0;
}

sub max_memory {
    my $self = shift;
    my $max_memory = $self->__max_memory(@_);
    my $max_memory_kb = $max_memory * 1_048_576;
    my $mem_limit_kb = Genome::Sys->mem_limit_kb;
    if ($mem_limit_kb) {
        my $safe_mem_limit_kb = int(0.8 * $mem_limit_kb);
        if ($max_memory_kb > $safe_mem_limit_kb) {
            my $safe_mem_limit_gb = floor($safe_mem_limit_kb / 1_048_576);
            if ($safe_mem_limit_gb == 0) {
                die "Does not work on systems with less than 1GB of memory.\n";
            }
            $max_memory = $safe_mem_limit_gb;
            $self->__max_memory($max_memory);
            warn "Overriding max_memory due to environmental limitations.";
        }
    }
    return $max_memory;
}

sub base_java_command {
    my $self = shift;

    my $gatk_path = $self->gatk_path;
    my $java_cmd = "java";
    if (defined $self->max_memory) {
        $java_cmd .= sprintf(" -Xmx%dg", $self->max_memory);
    }
    $java_cmd .= " -Djava.io.tmpdir=" . $self->tmp_dir;

    $java_cmd .= " -jar $gatk_path -et NO_ET";

    return $java_cmd;
}

1;

