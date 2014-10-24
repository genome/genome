package Genome::Model::Tools::Sam::Readcount;

use strict;
use warnings;

use Genome;
use File::Touch qw(touch);

my $DEFAULT_VER = '0.6';

class Genome::Model::Tools::Sam::Readcount {
    is  => 'Command',
    has_input => [
        bam_file => {
            is => 'String',
            doc => "The bam file from which to obtain readcounts",
        },
        reference_fasta => {
            is => 'String',
            doc => "The reference fasta to be used. This corresponds to the -f parameter.",
        },
        output_file => {
            is => 'String',
            doc => "The output file containing readcounts",
            is_output => 1,
        },
        region_list => {
            is => 'String',
            doc => "list of regions to report readcounts within. This should be in a tab delimited format with 'chromosome start stop'. This is the -l parameter.",
        },
    ],
    has_optional_input => [
        use_version => {
            is  => 'Version',
            doc => "bam-readcount version to be used.",
            default_value => $DEFAULT_VER,
        },
        minimum_mapping_quality => {
            is => 'Integer',
            default => 0,
            doc => "filter reads with mapping quality less than this. This is the -q parameter.",
        },
        minimum_base_quality => {
            is => 'Integer',
            default => 0,
            doc => "don't include reads where the base quality is less than this. This is the -b parameter. This is only available in versions 0.3 and later.",
        },
        comma_separated_format => {
            is  => 'Boolean',
            default => 0,
            doc => "report the mapping qualities as a comma separated list. This is the -d parameter. This is only available in versions 0.3 to 0.4.",
        },
        max_count => {
            is  => 'Integer',
            doc => "max depth to avoid excessive memory. This is the -d parameter in version 0.5.",
        },
        per_library => {
            is  => 'Bool',
            default => 0,
            doc => "report results per library. This is the -p parameter in version 0.5.",
        },
        insertion_centric => {
            is  => 'Bool',
            default => 0,
            doc => "do not include reads containing insertions after the current position in per-base counts. This is the -i parameter in version 0.5.",
        },
    ],
    doc => "Tool to get readcount information from a bam",
};

sub help_synopsis {
    "gmt readcount bam --minimum-base-quality 15 --use-version 0.4 --reference-fasta /path/to/fasta.fa --region-list /path/to/regions --bam-file /path/to/file.bam";
}

sub help_detail {
    "used to get readcount information from a bam";
}

sub default_version {
    return $DEFAULT_VER;
}

# The -w option was introduced in 0.5
sub version_has_warning_suppression {
    my ($self, $version) = @_;
    return ($version >= "0.5");
}

sub readcount_path {
    my $self = shift;
    my $version = $self->use_version;

    my $path = "/usr/bin/bam-readcount$version";
    if (! -x $path) {
        die $self->error_message("Failed to find executable bam-readcount version $version at $path!");
    }
    return $path;
}

sub command {
    my $self = shift;

    my $version = $self->use_version;

    my $command = sprintf "%s %s -f %s -l %s", $self->readcount_path, $self->bam_file, $self->reference_fasta, $self->region_list;

    if ($self->minimum_mapping_quality) {
        $command .= " -q " . $self->minimum_mapping_quality;
    }
    if ($self->minimum_base_quality) {
        $command .= " -b " . $self->minimum_base_quality;
    }
    if ($self->comma_separated_format && ($version eq "0.3" or $version eq "0.4")) {
        $command .= " -d";
    }
    if(defined $self->max_count) {
        $command .= " -d " . $self->max_count;
    }
    if($self->per_library) {
        $command .= " -p";
    }
    if($self->insertion_centric) {
        $command .= " -i";
    }
    if($self->version_has_warning_suppression($version)) {
        $command .= " -w 1";    # suppress error messages to a single report
    } else {
        $command .= ' 2> /dev/null';
    }

    return $command;
}

sub validate {
    my $self = shift;

    Genome::Sys->validate_file_for_reading($self->bam_file);
    Genome::Sys->validate_file_for_reading($self->reference_fasta);

    my $output_file = $self->output_file;
    if ($output_file =~ /[\(\)]/) {
        $output_file =~ s{\(}{\\(}g;
        $output_file =~ s{\)}{\\)}g;
        $self->output_file($output_file);
    }
    Genome::Sys->validate_file_for_writing($output_file);

    return 1;
}

sub execute {
    my $self = shift;

    $self->validate;

    unless (-s $self->region_list) {
        $self->warning_message("Region list provided (%s) does not exist or does not have size. Skipping run and touching output (%s).", $self->region_list, $self->output_file);
        touch($self->output_file);
        return 1;
    }

    Genome::Sys->shellcmd(
        cmd => $self->command . " > " . $self->output_file,
        input_files => [$self->bam_file, $self->reference_fasta, $self->region_list],
        output_files => [$self->output_file],
        allow_zero_size_output_files => 1,
    );
    $self->debug_message('Done running BAM Readcounts.');

    return 1;
}

1;

