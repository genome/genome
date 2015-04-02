package Genome::Model::Tools::Picard::CollectRnaSeqMetrics;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::CollectRnaSeqMetrics {
    is  => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        input_file => {
            is => 'String',
            doc => 'The SAM/BAM files to run on. File type is determined by suffix.',
            picard_param_name => 'INPUT',
        },
        output_file => {
            is => 'String',
            doc => 'The output summary file',
            picard_param_name => 'OUTPUT',
        },
        chart_output => {
            is => 'String',
            doc => 'The output PDF chart.',
            is_optional => 1,
            picard_param_name => 'CHART_OUTPUT',
        },
        refseq_file => {
            is => 'String',
            doc => 'The reference sequence file',
            picard_param_name => 'REFERENCE_SEQUENCE',
        },
        ribosomal_intervals_file => {
            is => 'Text',
            doc => 'Location of rRNA sequences in genome, in interval_list '
                 . 'format. If not specified no bases will be identified as '
                 . 'being ribosomal. Interval files can be created with gmt '
                 . 'picard bed-to-interval-list. More details here: '
                 . 'http://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/util/IntervalList.html',
            picard_param_name => 'RIBOSOMAL_INTERVALS',
        },
        ref_flat_file => {
            is => 'Text',
            doc => 'Gene annotations in refFlat form. Format described here: '
                 . 'http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat',
            picard_param_name => 'REF_FLAT',
        },
        assume_sorted => {
            is => 'Boolean',
            doc => 'coordinate sorted beforehand or not, default 1',
            default_value => 1,
            is_optional => 1,
            picard_param_name => 'ASSUME_SORTED',
        },
        stop_after => {
            is => 'Integer',
            doc => 'Stop after processing N reads, mainly for debugging.',
            is_optional => 1,
            picard_param_name => 'STOP_AFTER',
        },
        strand_specificity => {
            is => 'String',
            doc => 'For strand-specific library prep. For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND if the reads are expected to be on the transcription strand.',
            valid_values => ['NONE', 'FIRST_READ_TRANSCRIPTION_STRAND', 'SECOND_READ_TRANSCRIPTION_STRAND'],
            default_value => 'NONE',
            is_optional => 1,
            picard_param_name => 'STRAND_SPECIFICITY',
        },
    ],
};

sub help_brief {
    'Program to collect metrics about the alignment of RNA to various functional classes of loci in the genome: coding, intronic, UTR, intergenic, ribosomal.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics
EOS
}

sub _jar_name {
    return 'CollectRnaSeqMetrics.jar';
}

sub _java_class {
    return qw(picard analysis CollectRnaSeqMetrics);
}

sub _validate_params {
    my $self = shift;

    if ($self->chart_output and $self->version_older_than('1.52')) {
        die $self->error_message('Picard CollectRnaSeqMetrics version ',
            $self->use_version, ' does not support chart_output');
    }
}

sub _shellcmd_extra_params {
    my $self = shift;

    my @inputs = grep defined, map $self->$_, map "${_}_file", qw(
        input ribosomal_intervals ref_flat refseq
    );

    return (
        input_files => \@inputs,
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
    );
}

sub parse_file_into_metrics_hashref {
    my ($class,$metrics_file) = @_;
    my $as_fh = Genome::Sys->open_file_for_reading($metrics_file);
    my @headers;
    my %data;
    while (my $line = $as_fh->getline) {
        chomp($line);
        if ($line =~ /^## METRICS CLASS/) {
            my $next_line = $as_fh->getline;
            chomp($next_line);
            @headers = split("\t",$next_line);
            next;
        }
        if (@headers) {
            if ($line =~ /^\s+$/) {
                last;
            } else {
                my @values = split("\t",$line);
                for (my $i = 0; $i < scalar(@values); $i++) {
                    my $header = $headers[$i];
                    my $value = $values[$i];
                    $data{$header} = $value;
                }
            }
        }
    }
    return \%data;
}

1;
