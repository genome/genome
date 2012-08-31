package Genome::Model::Tools::Fastq::QualityTrim;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Fastq::QualityTrim {
    is => ['Command'],
    has => [
        input_fastq_file => {
            is => 'Text',
            doc => 'The input fastq file',
            is_input => 1,
        },
        output_fastq_file => {
            is => 'Text',
            doc => 'The trimmed output fastq file',
            is_input => 1,
        },
        summary_file => {
            is => 'Text',
            doc => 'A summary file containing histogram of trimmed read lengths',
            is_input => 1,
        },
        phred_quality_value => {
            doc => 'the sanger quality value to start read trim',
            default_value => 2,
            is_optional => 1,
            is_param => 1,
        },
        input_quality_format => {
            doc => 'the quality format of the input file',
            default_value => 'sanger',
            valid_values => ['sanger','illumina'],
        }
    ],
};

sub help_brief {
    'A tool to remove all bases and quality values after the defined quality_value param',
}

sub help_synopsis {
    'remove all bases from illumina fastq at quality_value and to the end of the read',
}

sub execute {
    my $self = shift;

    my $in_fh = Genome::Sys->open_file_for_reading($self->input_fastq_file);
    my $out_fh = Genome::Sys->open_file_for_writing($self->output_fastq_file);
    my $qual_chr;
    if ($self->input_quality_format eq 'illumina') {
        $qual_chr = chr(64 + $self->phred_quality_value);
    } elsif ($self->input_quality_format eq 'sanger') {
        $qual_chr = chr(33 + $self->phred_quality_value);
    } else {
        die('Failed to convert '. $self->input_quality_format .' for phred_quality_value '. $self->phred_quality_value);
    }
    my %histogram;
    while (my $desc_line = $in_fh->getline) {
        my $seq_line = $in_fh->getline;
        chomp($seq_line);
        my $opt_line = $in_fh->getline;
        my $qual_line = $in_fh->getline;
        chomp($qual_line);
        my $qual_len;
        if ($qual_line =~ /$qual_chr/) {
            $qual_line = $`;
            $qual_len = length($qual_line);
            if ($qual_len) {
                $seq_line = substr($seq_line,0,$qual_len);
            } else {
                $seq_line = 'N';
                $qual_line = $qual_chr;
            }
        } else {
            $qual_len = length($qual_line);
        }
        $histogram{$qual_len}++;
        print $out_fh $desc_line . $seq_line ."\n". $opt_line . $qual_line ."\n";
    }
    $in_fh->close;
    $out_fh->close;
    my $summary_fh = Genome::Sys->open_file_for_writing($self->summary_file);
    for my $length (sort {$a <=> $b} keys %histogram) {
        print $summary_fh $length ."\t". $histogram{$length} ."\n";
    }
    $summary_fh->close;
    return 1;
}
