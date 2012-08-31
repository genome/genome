package Genome::Model::Tools::Nimblegen::DesignFromAnnotation;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;

class Genome::Model::Tools::Nimblegen::DesignFromAnnotation {
    is => 'Command',
    has => [
    annotation_file => {
        type => 'String', is_optional => 1,
        doc => "An annotation file of sites to validate. Assumes STDIN if not specified",
    },
    output_file => {
        type => 'String', is_optional => 1,
        doc => "Output file. A list of probe regions for capture validation. Assumes STDOUT if not specified",
    },
    span => {
        type => 'Integer', is_optional => 1, default => 200,
        doc => "The number of bases to span the region upstream and downstream of a variant locus",
    },
    include_mitochondrial_sites => {
        type => 'Boolean', is_optional => 1, default => 0,
        doc => "Whether or not to remove sites on the mitochondria or non-chromosomal contigs",
    },
    include_unplaced_contig_sites => {
        type => 'Boolean', is_optional => 1, default => 1,
        doc => "Whether or not to remove sites on the unplaced contigs of the chromosome",
    },
    include_y_chrom_sites => {
        type => 'Boolean', is_optional => 1, default => 1,
        doc => "Whether or not to include sites on the Y chromosome in the output (if cases are all female)",
    },
    reference_index => {
        type => 'String', is_optional => 0,
        default => "/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa.fai",
        doc => "Samtools index of the reference sequence (To check chromosomal bounds of regions)",
    },
    ]
};

sub execute {
    my $self=shift;
    my $span = $self->span;
    my $reference_index = $self->reference_index;
    my $include_mitochondrial_sites = $self->include_mitochondrial_sites;
    my $include_unplaced_contig_sites = $self->include_unplaced_contig_sites;
    my $include_y_chrom_sites = $self->include_y_chrom_sites;

    # Hash to help convert hg19 unplaced contig names to the UCSC equivalent
    my %ucsc_unplaced_contigs = (( 'GL000207.1', 'chr18_gl000207_random' ), ( 'GL000226.1', 'chrUn_gl000226' ), ( 'GL000229.1', 'chrUn_gl000229' ), ( 'GL000231.1', 'chrUn_gl000231' ), ( 'GL000210.1', 'chr21_gl000210_random' ), ( 'GL000239.1', 'chrUn_gl000239' ), ( 'GL000235.1', 'chrUn_gl000235' ), ( 'GL000201.1', 'chr9_gl000201_random' ), ( 'GL000247.1', 'chrUn_gl000247' ), ( 'GL000245.1', 'chrUn_gl000245' ), ( 'GL000197.1', 'chr8_gl000197_random' ), ( 'GL000203.1', 'chr17_gl000203_random' ), ( 'GL000246.1', 'chrUn_gl000246' ), ( 'GL000249.1', 'chrUn_gl000249' ), ( 'GL000196.1', 'chr8_gl000196_random' ), ( 'GL000248.1', 'chrUn_gl000248' ), ( 'GL000244.1', 'chrUn_gl000244' ), ( 'GL000238.1', 'chrUn_gl000238' ), ( 'GL000202.1', 'chr11_gl000202_random' ), ( 'GL000234.1', 'chrUn_gl000234' ), ( 'GL000232.1', 'chrUn_gl000232' ), ( 'GL000206.1', 'chr17_gl000206_random' ), ( 'GL000240.1', 'chrUn_gl000240' ), ( 'GL000236.1', 'chrUn_gl000236' ), ( 'GL000241.1', 'chrUn_gl000241' ), ( 'GL000243.1', 'chrUn_gl000243' ), ( 'GL000242.1', 'chrUn_gl000242' ), ( 'GL000230.1', 'chrUn_gl000230' ), ( 'GL000237.1', 'chrUn_gl000237' ), ( 'GL000233.1', 'chrUn_gl000233' ), ( 'GL000204.1', 'chr17_gl000204_random' ), ( 'GL000198.1', 'chr9_gl000198_random' ), ( 'GL000208.1', 'chr19_gl000208_random' ), ( 'GL000191.1', 'chr1_gl000191_random' ), ( 'GL000227.1', 'chrUn_gl000227' ), ( 'GL000228.1', 'chrUn_gl000228' ), ( 'GL000214.1', 'chrUn_gl000214' ), ( 'GL000221.1', 'chrUn_gl000221' ), ( 'GL000209.1', 'chr19_gl000209_random' ), ( 'GL000218.1', 'chrUn_gl000218' ), ( 'GL000220.1', 'chrUn_gl000220' ), ( 'GL000213.1', 'chrUn_gl000213' ), ( 'GL000211.1', 'chrUn_gl000211' ), ( 'GL000199.1', 'chr9_gl000199_random' ), ( 'GL000217.1', 'chrUn_gl000217' ), ( 'GL000216.1', 'chrUn_gl000216' ), ( 'GL000215.1', 'chrUn_gl000215' ), ( 'GL000205.1', 'chr17_gl000205_random' ), ( 'GL000219.1', 'chrUn_gl000219' ), ( 'GL000224.1', 'chrUn_gl000224' ), ( 'GL000223.1', 'chrUn_gl000223' ), ( 'GL000195.1', 'chr7_gl000195_random' ), ( 'GL000212.1', 'chrUn_gl000212' ), ( 'GL000222.1', 'chrUn_gl000222' ), ( 'GL000200.1', 'chr9_gl000200_random' ), ( 'GL000193.1', 'chr4_gl000193_random' ), ( 'GL000194.1', 'chr4_gl000194_random' ), ( 'GL000225.1', 'chrUn_gl000225' ), ( 'GL000192.1', 'chr1_gl000192_random' ));

    # Depending on the bool flags the user sets, sites on some chromosomes will be excluded
    my %valid_chrs = map{ chomp; $_ => 1 } `cut -f 1 $reference_index`; # Used to check input files for valid ref names
    my %include_chrs = ( $include_unplaced_contig_sites ? %valid_chrs : ( map{ $_ => 1 } ( 1..22, qw( X Y MT ))));
    delete $include_chrs{MT} unless( $include_mitochondrial_sites );
    delete $include_chrs{Y} unless( $include_y_chrom_sites );

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
    if(defined $self->annotation_file) {
        $input_fh = IO::File->new($self->annotation_file,"r");
        unless($input_fh) {
            $self->error_message("Unable to open file ". $self->annotation_file . " for reading.");
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

    while(my $line = $input_fh->getline) {
        next if( $line =~ /^(#|chromosome_name|Chr\t)/ ); #Skip headers
        chomp $line;
        my ($chr,$start,$stop,) = split /\t/, $line;
        $chr =~ s/chr//i;
        next unless( defined $include_chrs{$chr} );

        if( $start < 1 || $stop < 1 || $start > $chromosome_lengths{$chr} || $stop > $chromosome_lengths{$chr}) {
            $self->error_message("Coordinates out of bounds. Skipping line: $line");
            next;
        }

        #If the span goes out of bounds of the chromosome, then clip it
        #Nimblegen design files are one based, so don't clip to 0 and don't calculate the start as a 0-based coordinate
        my $new_start = (( $start - $span < 1 ) ? 1 : $start - $span);
        my $new_stop = (( $stop + $span > $chromosome_lengths{$chr} ) ? $chromosome_lengths{$chr} : $stop + $span );
        $chr = $ucsc_unplaced_contigs{$chr} if( defined $ucsc_unplaced_contigs{$chr} );
        $chr = "chr$chr" unless( $chr =~ m/^chr/ );
        printf $output_fh "%s\t%d\t%d\t%d\t%s\n", $chr, $new_start, $new_stop, ($new_stop - $new_start + 1), $line;
    }
    $input_fh->close;
    $output_fh->close;

    return 1;

}

sub help_brief {
    "Takes an annotation file and produces a Design file for capture validation.";
}

sub help_detail {
    return <<EOS
Takes a file in annotation format and produces a list of regions to target for validation via
Nimblegen Solid Phase Capture Array.
EOS
}

1;
