package Genome::Model::Tools::Gatk;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use File::Temp;

my $DEFAULT_VERSION = '5336';
my $GATK_COMMAND = 'GenomeAnalysisTK.jar';

class Genome::Model::Tools::Gatk {
    is => ['Command'],
    has_optional => [
        version => {
            is    => 'string',
            doc   => 'version of Gatk application to use',
            default => $DEFAULT_VERSION,
        },
        _tmp_dir => {
            is => 'string',
            doc => 'a temporary directory for storing files',
        },
        max_memory => {
            is => 'Text',
            doc => 'Parameter to provide to the Java -Xmx argument for maximum memory. Should be something like "3000m" or "16g".',
            is_optional => 1,
        },
    ]
};

sub help_brief {
    "tools to work with Gatk output"
}

sub help_detail {                           # This is what the user will see with --help <---
    return <<EOS

EOS
}

my %GATK_VERSIONS = (
    'v1' => '/gsc/scripts/pkg/bio/gatk/GenomeAnalysisTK-1.0.5336/' . $GATK_COMMAND, # This is temporary... "v1" is what is in the first set of somatic-variation processing profiles
    '2986' => '/gsc/scripts/pkg/bio/gatk/GenomeAnalysisTK-1.0.2986/' . $GATK_COMMAND,
    '3362' => '/gsc/scripts/pkg/bio/gatk/GenomeAnalysisTK-1.0.3362/' . $GATK_COMMAND,
    '3362P' => '/gsc/scripts/pkg/bio/gatk/GenomeAnalysisTK-1.0.3362P/' . $GATK_COMMAND,
    '3423' => '/gsc/scripts/pkg/bio/gatk/GenomeAnalysisTK-1.0.3423/' . $GATK_COMMAND,
    '3471' => '/gsc/scripts/pkg/bio/gatk/GenomeAnalysisTK-1.0.3471/' . $GATK_COMMAND,
    '4168' => '/gsc/scripts/pkg/bio/gatk/GenomeAnalysisTK-1.0.4168/' . $GATK_COMMAND,
    '5336' => '/gsc/scripts/pkg/bio/gatk/GenomeAnalysisTK-1.0.5336/' . $GATK_COMMAND,
    '5777' => $ENV{GENOME_SW} . '/gatk/GenomeAnalysisTK-1.0.5777/' . $GATK_COMMAND,
);

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    unless (Genome::Config->arch_os =~ /64/) {
        $self->error_message('We recommend running GATK from 64-bit architecture');
        return;
    }
    my $tempdir = Genome::Sys->create_temp_directory;
    $self->_tmp_dir($tempdir);

    return $self;
}

sub gatk_path {
    my $self = $_[0];
    return $self->path_for_gatk_version($self->version);
}

sub available_gatk_versions {
    my $self = shift;
    return keys %GATK_VERSIONS;
}

sub path_for_gatk_version {
    my $class = shift;
    my $version = shift;

    if (defined $GATK_VERSIONS{$version}) {
        return $GATK_VERSIONS{$version};
    }
    die('No path for gatk version '. $version);
}

sub default_gatk_version {
    die "default gatk version: $DEFAULT_VERSION is not valid" unless $GATK_VERSIONS{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}

sub has_version {
    my $self = shift;
    my $version = shift;
    unless(defined($version)){
        $version = $self->version;
    }
    if(exists($GATK_VERSIONS{$version})){
        return 1;
    }
    return 0;
}

sub base_java_command {
    my $self = shift;

    my $gatk_path = $self->gatk_path;
    my $java_cmd = "java";
    if (defined $self->max_memory) {
        $java_cmd .= " -Xmx".$self->max_memory;
    }
    if (defined $self->_tmp_dir) {
        $java_cmd .= " -Djava.io.tmpdir=" . $self->_tmp_dir;
    }

    $java_cmd .= " -jar $gatk_path -et NO_ET";

    return $java_cmd;
}

1;

