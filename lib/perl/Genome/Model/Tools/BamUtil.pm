package Genome::Model::Tools::BamUtil;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use File::Temp;
use POSIX qw(floor);

my $DEFAULT_VERSION = '1.0.11';
my $BAMUTIL_COMMAND = "bam";

class Genome::Model::Tools::BamUtil {
    is => ['Command'],
    has_input => [
        version => {
            is    => 'string',
            doc   => 'version of BamUtil application to use',
            default => $DEFAULT_VERSION,
        },
    ],
    # has_param => [
        # lsf_resource => {
            # default => '-M 16777216 rusage[mem=16384] select[type==LINUX64 & mem > 16384] span[hosts=1]',
        # },
    # ],
};

sub help_brief {
    "tools to work with BamUtil"
}

sub help_detail {                           # This is what the user will see with --help <---
    return <<EOS

EOS
}

sub bamutil_versions {
    my %BAMUTIL_VERSIONS = (
        '1.0.11' => File::Spec->join('', 'gscuser', 'ssiebert', 'bamUtil', $BAMUTIL_COMMAND)
        # '2986' => $ENV{GENOME_SW} . '/gatk/GenomeAnalysisTK-1.0.2986/' . $GATK_COMMAND,
    );
    return %BAMUTIL_VERSIONS;
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    return $self;
}

sub bamutil_path {
    my $self = $_[0];
    return $self->path_for_bamutil_version($self->version);
}

sub available_bamutil_versions {
    my $self = shift;
    my %versions = $self->bamutil_versions;
    return keys %versions;
}

sub path_for_bamutil_version {
    my $class = shift;
    my $version = shift;
    my %versions = $class->bamutil_versions;
    if (defined $versions{$version}) {
        return $versions{$version};
    }
    die('No path for bamutil version '. $version);
}

sub default_bamutil_version {
    my $class = shift;
    my %versions = $class->bamtutil_versions;
    die "default BamUtil version: $DEFAULT_VERSION is not valid" unless $versions{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}

sub has_version {
    my $self = shift;
    my $version = shift;
    unless(defined($version)){
        $version = $self->version;
    }
    my %versions = $self->bamutil_versions;
    if(exists($versions{$version})){
        return 1;
    }
    return 0;
}

sub base_command {
    my $self = shift;

    my $bamutil_path = $self->bamutil_path;

    return $bamutil_path;
}

1;

