package Genome::Model::Tools::Picard::CollectRnaSeqMetrics;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::CollectRnaSeqMetrics {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file => {
            is  => 'String',
            doc => 'The SAM/BAM files to run on.  File type is determined by suffix.',
        },
        output_file => {
            is  => 'String',
            doc => 'The output summary file',
        },
        chart_output => {
            is_optional => 1,
            is => 'String',
            doc => 'The output PDF chart.',
        },
        refseq_file => {
            is  => 'String',
            doc => 'The reference sequence file',
        },
        ribosomal_intervals_file => {
            is  => 'Text',
            doc => 'Location of rRNA sequences in genome, in interval_list format. If not specified no bases will be identified as being ribosomal. Format described here: http://picard.sourceforge.net/javadoc/net/sf/picard/util/IntervalList.html',
        },
        ref_flat_file => {
            is  => 'Text',
            doc => 'Gene annotations in refFlat form. Format described here: http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat',
        },
        assume_sorted => {
            is  => 'Boolean',
            doc => 'coordinate sorted beforehand or not, default 1',
            default_value => 1,
            is_optional   => 1,
        },
        stop_after => {
            is  => 'Integer',
            doc => 'Stop after processing N reads, mainly for debugging.',
            is_optional   => 1,
        },
        strand_specificity => {
            is  => 'String',
            doc => 'For strand-specific library prep. For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND if the reads are expected to be on the transcription strand.',
            valid_values =>['NONE','FIRST_READ_TRANSCRIPTION_STRAND','SECOND_READ_TRANSCRIPTION_STRAND'],
            default_value => 'NONE',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Program to collect metrics about the alignment of RNA to various functional classes of loci in the genome: coding, intronic, UTR, intergenic, ribosomal.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#CollectRnaSeqMetrics
EOS
}

sub execute {
    my $self = shift;

    my @input_files;
    for my $type qw(input ribosomal_intervals ref_flat) {
        my $property_name = $type.'_file';
        unless ($self->$property_name and -s $self->$property_name) {
            $self->error_message("$property_name is invalid");
            return;
        }
        push @input_files, $self->$property_name;
    }

    my $cmd = $self->picard_path .'/CollectRnaSeqMetrics.jar net.sf.picard.analysis.CollectRnaSeqMetrics';
    $cmd   .= ' OUTPUT='. $self->output_file  .' INPUT='. $self->input_file  .' RIBOSOMAL_INTERVALS='. $self->ribosomal_intervals_file .' REF_FLAT='. $self->ref_flat_file;
    if ($self->refseq_file) {
        $cmd .= ' REFERENCE_SEQUENCE='. $self->refseq_file;
        push @input_files, $self->refseq_file;
    }
    if ($self->assume_sorted) {
        $cmd .= ' ASSUME_SORTED=true';
    } else {
        $cmd .= ' ASSUME_SORTED=false';
    }
    if ($self->stop_after) {
        $cmd .= ' STOP_AFTER='. $self->stop_after;
    }
    if ($self->strand_specificity) {
        $cmd .= ' STRAND_SPECIFICITY='. $self->strand_specificity;
    }
    if ($self->chart_output) {
        if ($self->use_version >= '1.52') {
            $cmd .= ' CHART_OUTPUT='. $self->chart_output;
        } else {
            $self->warning_message('Picard version '. $self->use_version .' does not support CHART_OUTPUT!');
        }
    }
    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => \@input_files,
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
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
