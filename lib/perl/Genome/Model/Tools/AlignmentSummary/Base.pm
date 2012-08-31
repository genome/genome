package Genome::Model::Tools::AlignmentSummary::Base;

use strict;
use warnings;

use Genome;
use YAML;

my $DEFAULT_VERSION = '1.2.6';

class Genome::Model::Tools::AlignmentSummary::Base {
    is  => 'Command::V2',
    is_abstract => 1,
    has_input => [
        use_version => {
            is => 'Version',
            doc => 'The version of the C++ alignment_summary command to use.',
            default_value => $DEFAULT_VERSION,
            is_optional => 1,
        },
    ],
};

my %VERSIONS = (
    '1.2.1'  => '/gsc/pkg/alignment_summary_cpp/alignment_summary_cpp-1.2.1/Debug/alignment_summary_cpp',
    '1.2.6' => '/gsc/scripts/opt/genome_legacy_code/bin/alignment-summary-v1.2.6',
);

sub available_versions {
    my $self = shift;
    my @versions = keys %VERSIONS;
    return @versions;
}

sub path_for_version {
    my $class = shift;
    my $version = shift;

    if (defined $VERSIONS{$version}) {
        return $VERSIONS{$version};
    }
    die('No path for version '. $version);
}

sub default_version {
    die "default version: $DEFAULT_VERSION is not valid" unless $VERSIONS{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}

sub run_in_bash_env {
    my $self = shift;
    my %params = @_;

    my $as_cmd = delete($params{cmd});
    
    my $bash_cmd = "bash -c 'LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:/gsc/scripts/opt/genome_legacy_code/lib:/gsc/pkg/boost/boost_1_42_0/lib $as_cmd'";
    $params{cmd} = $bash_cmd;
    return Genome::Sys->shellcmd(%params);
}

sub parse_metrics_file_into_hashref {
    my ($class,$metrics_file) = @_;
    
    my $res = YAML::LoadFile($metrics_file);
    unless (ref($res) eq 'HASH') {
        die('Failed to parse YAML hash from alignment_summary output: '. $metrics_file);
    }
    return $res;
}



1;
