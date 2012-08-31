package Genome::Model::Tools::Varscan;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use File::Temp;

my $DEFAULT_VERSION = '2.2.6';

class Genome::Model::Tools::Varscan {
    is => ['Command'],
    has_optional_input => [
        version => {
            is    => 'String',
            doc   => 'version of Varscan application to use',
            default_value => 'latest',
        },
        no_headers => {
            is => 'Boolean',
            doc => 'Stop varscan from putting headers on its output files',
            default_value => '0',
        },
    ],
};

sub help_brief {
    "tools to work with Varscan output"
}

sub help_detail {                           # This is what the user will see with --help <---
    return <<EOS

EOS
}

my %VARSCAN_VERSIONS = (
    '2.2.9' => '/gsc/scripts/lib/java/VarScan/VarScan.v2.2.9.jar',
    '2.2.6' => '/gsc/scripts/lib/java/VarScan/VarScan.v2.2.6.jar',
    '2.2.4' => '/gsc/scripts/lib/java/VarScan/VarScan.v2.2.4.jar',
);

sub java_command_line {
    my $self = shift;
    my $parameter_string = shift;

    my $path = $self->path_for_version($self->version);
    my $headers = $self->no_headers ? "--no-headers 1" : "";
    my $command_line = 'bash -c "java -jar ' . $path . ' ' . $parameter_string . ' '.$headers.' "';

    return $command_line;
}

sub path_for_version {
    my $class = shift;
    my $version = shift || $DEFAULT_VERSION;

    if($version eq 'latest') {
        return $class->path_for_latest_version;
    }

    unless(exists $VARSCAN_VERSIONS{$version}) {
        $class->error_message('No path found for Varscan Version ' . $version);
        die $class->error_message;
    }

    return $VARSCAN_VERSIONS{$version};
}

sub pileup_command_for_reference_and_bam {
    my $self = shift;
    my $reference = shift;
    my $bam = shift;
    my $mapqual = shift;

    $mapqual = 10 if(!$mapqual);

    my $command;
    # TODO this should be made a little cleaner, but it works for now.
    if ($self->version eq "2.2.4") {
        my $samtools_path = Genome::Model::Tools::Sam->path_for_samtools_version("r963"); # The last version of samtools that supports pileup
        $command = "$samtools_path view -b -u -q $mapqual $bam | $samtools_path pileup -f $reference -";
    } else {
        my $samtools_path = Genome::Model::Tools::Sam->path_for_samtools_version("r963"); # The latest version of samtools installed should go here
        $command = "$samtools_path mpileup -f $reference -q $mapqual $bam";
    }

    return $command;
}

sub path_for_latest_version {
    my $class = shift;
    my $link = '/gsc/scripts/lib/java/VarScan/VarScan.jar';

    unless(-e $link and -l $link) {
        $class->error_message('Link to latest version not found or not a link!');
    }

    return $link;
}

sub default_version {
    my $class = shift;

    unless(exists $VARSCAN_VERSIONS{$DEFAULT_VERSION}) {
        $class->error_message('Default Varscan version (' . $DEFAULT_VERSION . ') is invalid.');
        die $class->error_message;
    }

    return $DEFAULT_VERSION;
}

sub available_varscan_versions {
    return keys(%VARSCAN_VERSIONS);
}

1;

