package Genome::Model::Tools::Nimblegen::DesignFromLowexoncov;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;

class Genome::Model::Tools::Nimblegen::DesignFromLowexoncov {
    is => 'Command',
    has => [
    low_exon_cov_file => {
        type => 'String',
        is_optional => 1,
        doc => "An file with low exon coverage. Assumes STDIN if not specified",
    },
    output_file => {
        type => 'String',
        is_optional => 1,
        doc => "Output file. A list of probe regions for capture validation. Assumes STDOUT if not specified",
    },
    span => {
        type => 'Integer',
        is_optional => 1,
        default => 200,
        doc => "The number of bases to span the region upstream and downstream of a variant locus",
    },
    exclude_non_canonical_sites => {
        type => 'Bool',
        is_optional => 1,
        default => 1,
        doc => "Whether or not to remove sites on the mitochondria or non-chromosomal contigs",
    },
    include_y => {
        type => 'Bool',
        is_optional => 1,
        default => 1,
        doc => "Whether or not to include sites on the Y chromosome in the output",
    },
    reference_index => {
        type => 'String',
        is_optional => 0,
        doc => "Samtools index of the reference sequence (To check chromosomal bounds of regions)",
        default => "/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa.fai",
    },
    
    ]
};

sub help_brief {

    "The low exon coverage input file should have the format 'gene   Strand    Chromosome    start    end    Avg_Coverage";

}


sub execute {
    my $self=shift;

    my $reference_index = $self->reference_index;
    my $span = $self->span;

    my $fh = IO::File->new($reference_index,"r");
    unless($fh) {
        $self->error_message("Unable to open the reference sequence index: $reference_index");
        return;
    }

    #read in the index to get the chromosome lengths
    my %chromosome_lengths;
    while(my $line = $fh->getline) {
        chomp $line;
        my ($chr, $length) = split /\t/, $line;
        $chromosome_lengths{$chr} = $length;
    }
    $fh->close;

    my $output_fh;
    if(defined $self->output_file) {
        $output_fh = IO::File->new($self->output_file,"w");
        unless($output_fh) {
            $self->error_message("Unable to open file " . $self->output_file . " for writing.");
            return;
        }
    }
    else {
        $output_fh = IO::File->new_from_fd(fileno(STDOUT),"w");
        unless($output_fh) {
            $self->error_message("Unable to open STDOUT for writing.");
            return;
        }
    }

    my $input_fh;
    if(defined $self->low_exon_cov_file) {
        $input_fh = IO::File->new( $self->low_exon_cov_file,"r");
        unless($input_fh) {
            $self->error_message("Unable to open file ". $self->low_exon_cov_file . " for reading.");
            return;
        }
    }
    else {
        $input_fh = IO::File->new_from_fd(fileno(STDIN),"r");
        unless($input_fh) {
            $self->error_message("Unable to open STDIN for reading.");
            return;
        }
    }

    my %canonical_chrs = map{ $_ => 1 } (1..22, qw( X Y x y));
    while(my $line = $input_fh->getline) {
        next if( $line =~ /^(#|chromosome_name|Chr\t|Gene)/ ); #Skip headers
        next if( !$self->include_y && $line =~ /(^y\t|^chry\t)/i );
        chomp $line;
	my ($gene,$strand,$chr,$start,$stop,$avg_cov) = split /\t/, $line;
        $chr =~ s/chr//i;
        next if( $self->exclude_non_canonical_sites && ( !defined $canonical_chrs{$chr} ));

        if( $start < 1 || $stop < 1 || $start > $chromosome_lengths{$chr} || $stop > $chromosome_lengths{$chr}) {
            $self->error_message("Coordinates out of bounds. Skipping line: $line");
            next;
        }

        #If the span goes out of bounds of the chromosome, then clip it
        my $new_start = (( $start - $span < 1 ) ? 0 : $start - $span );
        my $new_stop = (( $stop + $span > $chromosome_lengths{$chr} ) ? $chromosome_lengths{$chr} : $stop + $span );
        printf $output_fh "chr%s\t%d\t%d\t%d\t%s\n", $chr, $new_start, $new_stop, ($new_stop - $new_start), $line;
    }
    $input_fh->close;
    $output_fh->close;

    return 1;

}

sub help_detail {
    return <<EOS

Takes a low_exon_coverage input file with tab-delimited column 
'gene'   'Strand'    'Chromosome'    'start'    'end'    'Avg_Coverage'
and produces a list of regions to target for validation via 
Nimblegen Solid Phase Capture Array.

EOS

}

1;
