package Genome::Model::Tools::Bed::Convert::Snv::CosmicToBed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::Convert::Snv::CosmicToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Snv'],
    has => [
        reference_fasta => {
            is => 'String',
            doc => 'Fasta file of the reference chromosomes',
            is_optional => 0,
            default_value => '/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa',
        },
    ],
};

sub help_brief {
    "Tools to convert Cosmic format to BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert snv cosmic-to-bed --source snps_all_sequences --output snps_all_sequences.bed
EOS
}

sub help_detail {                           
    return <<EOS
    This is a small tool to take SNV calls in cosmic format and convert them to a common BED format.
EOS
}

sub process_source {
    my $self = shift;

    my $input_fh = $self->_input_fh;

    my $reference_build = Genome::Model::Build->get($self->reference_build_id);


    while(my $line = <$input_fh>) {
        chomp $line;

        next if ($line =~ /Gene/);
        my @fields = split (/\t/, $line);
        my $cDNA = $fields[2];
        my $reference;
        my $variant;

        #extract positions
        my ($chr,$start,$end) = split /[:-]/,$fields[6];

        unless(defined $chr && defined $start && defined $end) {
            print STDERR "No build37 coordinates found for: $line\n";
            next;
        }

        if ($chr eq "23") {
            $chr = "X";
        }

        if ($cDNA =~ /c.\d+([+-]\d+)*(_\d+([+-]\d+)*)*ins[ACTG]+/) {
            ($variant) = $cDNA =~ /\d+ins([ACTG]+)/;
            $reference = 0;
            $end = $start;
        }
        elsif ($cDNA =~ /c.\d+([+-]\d+)*(_\d+([+-]\d+)*)*del[ACTG]+/) {
            ($reference) = $cDNA =~ /\d+del([ACTG]+)/;
            $variant = 0;
            #check strand
            my $cmd = "samtools faidx ".$self->reference_fasta." $chr:$start-$end";
            my $out = `$cmd`;
            chomp $out;
            my @lines = split(/\n/, $out);
            if ($reference ne $lines[1]) {
                $reference = reverse($reference);
                $reference =~ tr/ACGT/TGCA/;
                if ($reference ne $lines[1]) {
                    print STDERR "Not equal del $reference ".$lines[1]." $line\n";
                    next;
                }
            }
            $start--;
        }
        elsif ($cDNA =~ /\d+\D>\D/) {
            ($reference,$variant) = $cDNA =~ /(\D)>(\D)/;
            if ($end - $start + 1 != length($variant)) {
                next;
            }
            #check strand
            my $cmd = "samtools faidx ".$self->reference_fasta." $chr:$start-$end";
            my $out = `$cmd`;
            chomp $out;
            my @lines = split(/\n/, $out);
            if ($reference ne $lines[1]) {
                $reference = reverse($reference);
                $variant = reverse($variant);
                $reference =~ tr/ACGT/TGCA/;
                $variant =~ tr/ACGT/TGCA/;
                if ($reference ne $lines[1]) {
                    print STDERR "Not equal $reference ".$lines[1]." $line\n";
                    next;
                }
            }
            $start--;
        }
        else {
            print STDERR "Unable to extract alleles for $line\n";
            next;
        }

        #position => 1-based position of the SNV
        #BED uses 0-based position of and after the event
        $self->write_bed_line($chr, $start, $end, $reference, $variant);
    }

    return 1;
}

1;
