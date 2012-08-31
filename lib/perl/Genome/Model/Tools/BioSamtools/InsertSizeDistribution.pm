package Genome::Model::Tools::BioSamtools::InsertSizeDistribution;

use strict;
use warnings;

use Genome;
use Statistics::Descriptive;

class Genome::Model::Tools::BioSamtools::InsertSizeDistribution {
    is => 'Genome::Model::Tools::BioSamtools',
    has_input => [
        bam_file => {
            is => 'Text',
            doc => 'A path to a BAM format file of aligned capture reads',
        },
        output_file => {
            is => 'Text',
            doc => 'The output file to write histogram',
        },
    ],
};

sub execute {
    my $self = shift;
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    my $refcov_bam  = Genome::Model::Tools::RefCov::Bam->create(bam_file => $self->bam_file );
    unless ($refcov_bam) {
        die('Failed to load bam file '. $self->bam_file);
    }
    my $bam  = $refcov_bam->bio_db_bam;
    my $header = $bam->header();
    my $header_text = $header->text;
    my @lines = split/\n/, $header_text;
    my @read_group_lines = grep { $_ =~ /^\@RG/ } @lines;
    my %read_groups;
    my %library_metrics;
    my $default_library = 'library';
    my $default_rg = 'read_group';
    unless (@read_group_lines) {
        $read_groups{$default_rg} = $default_library;
        $library_metrics{$default_library} = {};
    }
    for my $read_group_line (@read_group_lines) {
        unless ($read_group_line =~ /ID:(\d+)/) {
            die('Failed to parse read group line for ID:  '. $read_group_line);
        }
        my $id = $1;
        unless ($read_group_line =~ /LB:(\S+)/) {
            die('Failed to parse read group line for LB:  '. $read_group_line);
        }
        my $lib = $1;
        unless (defined($read_groups{$id})) {
            $read_groups{$id} = $lib;
            $library_metrics{$lib} = {};
        } else {
            die('Found same read group twice with ID:  '. $id);
        }
    }
    
    for my $lib (keys %library_metrics) {
        $library_metrics{$lib}{all_stats} = Statistics::Descriptive::Sparse->new();
        $library_metrics{$lib}{histogram} = {};
    }
    
    my %mate_pairs;
    while (my $align = $bam->read1()) {
        my $rg = $align->aux_get('RG');
        unless ($rg) {
            $rg = $default_rg;
        }
        my $lib = $read_groups{$rg};
        my $flag = $align->flag;
        # paired in sequencing and proper pair
        if (($flag & 1) && ($flag & 2)) {
            my $qname = $align->qname;
            my $isize = $align->isize;
            # If the second read-end is found skip
            # TODO: This probably won't work for all aligners
            # MAQ just happens to truncate the read-end from read names
            if ($mate_pairs{$qname}) {
                delete($mate_pairs{$qname});
                next;
            }
            # TODO: Try to calculate inner distance for Tophat
            if (defined($isize)) {
                if ( $isize < 1 ) {
                    print $output_fh "#\t". $qname ."\t". $align->pos ."\t". $align->calend ."\t". $align->mate_start ."\t". $align->mate_end ."\t". $isize ."\n";
                }
                $library_metrics{$lib}{all_stats}->add_data($isize);
                $library_metrics{$lib}{histogram}->{$isize}++;
            }
            $mate_pairs{$qname} = 1;
        }
    }

    for my $lib (sort keys %library_metrics) {
        my $stats = $library_metrics{$lib}{all_stats};
        print $output_fh "#\t". $lib ."\t". $stats->mean ."\t". $stats->standard_deviation ."\n";
        for my $bin ( sort { $a <=> $b } keys %{$library_metrics{$lib}{histogram}}) {
            print $output_fh $bin ."\t". $library_metrics{$lib}{histogram}{$bin} ."\n";
        }
    }
    $output_fh->close;
    return 1;
}

1;
