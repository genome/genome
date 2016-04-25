package Genome::Model::Tools::Varscan;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use File::Temp;
use POSIX (qw(floor));

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
            default_value => $DEFAULT_VERSION,
        },
        no_headers => {
            is => 'Boolean',
            doc => 'Stop varscan from putting headers on its output files',
            default_value => '0',
        },
        java_interpreter => {
            is => 'Text',
            doc => 'The java interpreter to use',
            #LSF:  Not able to use the Genome::Sys->java_executable_path('1.7')
            #      since it depends on java version below.
            #
            #  java version "1.6.0_18"
            #  Java(TM) SE Runtime Environment (build 1.6.0_18-b07)
            #  Java HotSpot(TM) 64-Bit Server VM (build 16.0-b13, mixed mode)
            #
            default_value => 'java',
        },
        max_memory => {
            # Accessor is overridden so that it can be limited based on available memory or LSF limit.
            is => 'Text',
            doc => 'The maximum memory (GB) to use when running Java VM. Limited to environmental constraints.',
            default => '4',
        },
        tmp_dir => {
            is => 'Text',
            doc => 'Temporary directory for Java.',
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

my $sw_legacy_java = Genome::Config::get('sw_legacy_java');
my %VARSCAN_VERSIONS = (
    '2.3.6' => $sw_legacy_java . '/VarScan/VarScan.v2.3.6.jar',
    '2.3.5' => $sw_legacy_java . '/VarScan/VarScan.v2.3.5.jar',
    '2.3.2' => $sw_legacy_java . '/VarScan/VarScan.v2.3.2.jar',
    '2.3.1' => $sw_legacy_java . '/VarScan/VarScan.v2.3.1.jar',
    '2.2.9' => $sw_legacy_java . '/VarScan/VarScan.v2.2.9.jar',
    '2.2.6' => $sw_legacy_java . '/VarScan/VarScan.v2.2.6.jar',
    '2.2.4' => $sw_legacy_java . '/VarScan/VarScan.v2.2.4.jar',
);

sub java_command_line {
    my $self             = shift;
    my $parameter_string = shift;

    my @java_cmd = ( $self->java_interpreter );

    my $path = $self->path_for_version( $self->version );
    if ( defined $self->max_memory ) {
        push @java_cmd, sprintf( "-Xmx%dg", $self->max_memory );
    }
    if ( not $self->tmp_dir ) {
        my $tempdir = Genome::Sys->create_temp_directory;
        $self->tmp_dir($tempdir);
    }
    push @java_cmd, "-Djava.io.tmpdir=" . $self->tmp_dir, '-jar' => $path;
    push @java_cmd, $parameter_string;
    push @java_cmd, ( "--no-headers" => 1 ) if ( $self->no_headers );

    return join( ' ', @java_cmd );
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
    my $link = Genome::Config::get('sw_legacy_java') . '/VarScan/VarScan.jar';

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

1;

