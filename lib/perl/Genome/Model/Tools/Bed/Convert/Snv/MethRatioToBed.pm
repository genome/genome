package Genome::Model::Tools::Bed::Convert::Snv::MethRatioToBed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::Convert::Snv::MethRatioToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Snv'],
};

sub help_brief {
    "Tools to convert meth-ratio SNV format to BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert snv meth-ratio-to-bed --source snps_all_sequences --output snps_all_sequences.bed
EOS
}

sub help_detail {
    return <<EOS
    This is a small tool to take SNV calls in meth-ratio format and convert them to a common BED format (using the first five columns).
EOS
}

sub process_source {
    my $self = shift;

    my $input_fh = $self->_input_fh;

    while(my $line = <$input_fh>) {
        chomp($line);
        my @fields = split("\t", $line);

        #skip header
        next if $line =~/context/;

        #ratio is fixed-width column, so convert to integer by removing decimal and leading zero(es)
        my $score = $fields[4];
        $score =~ s/\.//;
        $score =~ s/^0+(\d)/$1/;

        $self->write_bed_line(
            $fields[0],   #chr
            $fields[1]-1, #st
            $fields[1],   #sp
            "C",          #ref (will always be C->T)
            "T",          #var
            $score,       #ratio as score
            $fields[2],   #strand
            $fields[5],   #depth
            $fields[3],   #context
            $fields[7],   #CI_low
            $fields[8],   #CI_high
            )
    }

    return 1;
}

1;
