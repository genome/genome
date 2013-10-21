package Genome::Model::Tools::Varscan;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use File::Temp;

my $DEFAULT_VERSION = '2.3.6';

class Genome::Model::Tools::Varscan {
    is => ['Command'],
    has => [
        samtools_version => {
            is    => 'String',
            doc => 'version of samtools to use when doing pileup or mpileup. Must be r963 or earlier if using samtools 2.2.4, which only supports pileup.',
            default_value => "r963",
        },
        samtools_use_baq => {
            is => 'Boolean',
            doc => 'When doing pileup/mpileup, should we enable baq (-B) option',
            is_input => 1,
            default_value => 1,
        },
        samtools_params => {
            is    => 'String',
            doc => 'Additional parameters to pass to samtools when doing pileup/mpileup',
            is_optional => 1,
        },
    ],
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
    '2.3.6' => $ENV{GENOME_SW_LEGACY_JAVA} . '/VarScan/VarScan.v2.3.6.jar',
    '2.3.5' => $ENV{GENOME_SW_LEGACY_JAVA} . '/VarScan/VarScan.v2.3.5.jar',
    '2.3.2' => $ENV{GENOME_SW_LEGACY_JAVA} . '/VarScan/VarScan.v2.3.2.jar',
    '2.3.1' => $ENV{GENOME_SW_LEGACY_JAVA} . '/VarScan/VarScan.v2.3.1.jar',
    '2.2.9' => $ENV{GENOME_SW_LEGACY_JAVA} . '/VarScan/VarScan.v2.2.9.jar',
    '2.2.6' => $ENV{GENOME_SW_LEGACY_JAVA} . '/VarScan/VarScan.v2.2.6.jar',
    '2.2.4' => $ENV{GENOME_SW_LEGACY_JAVA} . '/VarScan/VarScan.v2.2.4.jar',
);

sub java_command_line {
    my $self = shift;
    my $parameter_string = shift;

    my $path = $self->path_for_version($self->version);
    my $headers = $self->no_headers ? "--no-headers 1" : "";
    my $command_line = 'java -jar ' . $path . ' ' . $parameter_string . ' '.$headers;

    return $command_line;
}

sub command_line {
    my $self = shift;
    my $parameter_string = shift;
    return sprintf('bash -c "%s "', $self->java_command_line($parameter_string));
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

# Bams are now passed in as an arrayref so we can formulate a mpileup command with multiple bams
sub pileup_command_for_reference_and_bam {
    my $self = shift;
    my $reference = shift;
    my $bams = shift;
    my $mapqual = shift;

    $mapqual = 10 if(!$mapqual);

    my $command;
    my $samtools_version = $self->samtools_version;
    my $samtools_params = $self->samtools_params || "";
    my $samtools_use_baq = $self->samtools_use_baq;
    unless ($samtools_use_baq) {
        $samtools_params = join(" ", ($samtools_params, "-B")); #turn off baq
    }
    my $samtools_path = Genome::Model::Tools::Sam->path_for_samtools_version($samtools_version);

    # Use pileup for legacy varscan v2.2.4, because it could not handle mpileup output
    if ($self->version eq "2.2.4") {
        if (scalar(@$bams) > 1) {
            die $self->error_message("Multiple bams not allowed in samtools pileup");
        }
        my $bam = shift @$bams;
        $command = "$samtools_path view -b -u -q $mapqual $bam | $samtools_path pileup $samtools_params -f $reference -";
    } else {
        my $bam_string = join(" ", @$bams);
        $command = "$samtools_path mpileup -f $reference -q $mapqual $samtools_params $bam_string";
    }

    return $command;
}

sub path_for_latest_version {
    my $class = shift;
    my $link = $ENV{GENOME_SW_LEGACY_JAVA} . '/VarScan/VarScan.jar';

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

sub samtools_path {
    my $self = shift;
    return Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);
}

1;

