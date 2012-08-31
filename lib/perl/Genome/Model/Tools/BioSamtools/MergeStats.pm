package Genome::Model::Tools::BioSamtools::MergeStats;

use strict;
use warnings;

use Genome;
use Statistics::Descriptive;

class Genome::Model::Tools::BioSamtools::MergeStats {
    is => ['Genome::Model::Tools::BioSamtools'],
    has_input => [
        stats_file => {
            is => 'Text',
            doc => 'A STATS file output from refcov',
        },
        output_file => {
            is => 'Text',
            doc => 'The output file to write summary stats',
        },
        merge_by => {
            is => 'Text',
            is_optional => 1,
            default_value => 'gene',
            valid_values => ['exome','gene','transcript','exon'],
        },
        delimiter => {
            is => 'Text',
            doc => 'The character that delimits GENE, TRANSCRIPT, TYPE, and ORDINAL',
            default_value => ':',
            valid_values => [':','_','.'],
        },
        minimum_breadth_filters => {
            is => 'Text',
            is_optional => 1,
            doc => 'A comma-delimited list of minimum breadths as percentages.',
            default_value => '100,90,80',
        },
    ],
};

sub help_detail {
    'This command takes the STATS file from refcov and generates summary statistics based on all regions/targets.  This command will run on 32 or 64 bit architechture'
}

sub execute {
    my $self = shift;
    my $stats_fh = Genome::Sys->open_file_for_reading($self->stats_file);
    unless ($stats_fh) {
        die('Failed to read file stats file: '. $self->stats_file);
    }
    my @minimum_breadths = split(',',$self->minimum_breadth_filters);
    my %coverage;
    while (my $line = $stats_fh->getline) {
        chomp($line);
        my @entry = split("\t",$line);
        my $name = $entry[0];
        my ($gene,$transcript,$type,$ordinal) = split($self->delimiter,$name);
        my $key;
        if ($self->merge_by eq 'exome') {
            $key = 'exome';
        } elsif ($self->merge_by eq 'gene') {
            $key = $gene;
        } elsif ($self->merge_by eq 'transcript') {
            $key = $transcript;
        } elsif ($self->merge_by eq 'exon') {
            $key = $gene . $self->delimiter . $transcript . $self->delimiter . $type;
            if (defined $ordinal) {
                $key .= $self->delimiter . $ordinal;
            }
        }
        my $size = $entry[2];
        my $covered = $entry[3];
        my $breadth = $entry[1];
        $coverage{$key}{size} += $size;
        $coverage{$key}{covered} += $covered;
        $coverage{$key}{intervals}++;
        for my $min_breadth (@minimum_breadths) {
            if ($breadth >= $min_breadth) {
                $coverage{$key}{$min_breadth}{covered_intervals}++;
            }
        }
    }
    $stats_fh->close;
    my @headers = qw/NAME LENGTH COVERED MEAN_BREADTH MINIMUM_BREADTH INTERVALS COVERED_INTERVALS MISSING_INTERVALS/;
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    print $output_fh join("\t",@headers) ."\n";
    for my $min_breadth (@minimum_breadths) {
        for my $key (keys %coverage) {
            my $length = $coverage{$key}{size};
            my $covered = $coverage{$key}{covered} ;
            my $mean_breadth = (($covered / $length) * 100);
            my $intervals = $coverage{$key}{intervals};
            my $covered_intervals = 0;
            if ($coverage{$key}{$min_breadth}{covered_intervals}) {
                $covered_intervals = $coverage{$key}{$min_breadth}{covered_intervals};
            }
            my $missing_intervals = $intervals - $covered_intervals;
            print $output_fh $key ."\t". $length ."\t". $covered ."\t". $mean_breadth ."\t". $min_breadth ."\t". $intervals ."\t". $covered_intervals ."\t".  $missing_intervals ."\n";
        }
    }
    return 1;
}


1;
