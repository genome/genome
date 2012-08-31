package Genome::Model::Tools::Bed::Convert::Snv::VarscanSomaticToBed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::Convert::Snv::VarscanSomaticToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Snv'],
};

sub help_brief {
    "Tools to convert varscan-somatic SNV format to BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert snv varscan-somatic-to-bed --source snps_all_sequences --output snps_all_sequences.bed
EOS
}

sub help_detail {                           
    return <<EOS
    This is a small tool to take SNV calls in var-scan format and convert them to a common BED format (using the first four columns + quality).
EOS
}

sub process_source {
    my $self = shift;
    
    my $input_fh = $self->_input_fh;
    
    while(my $line = $input_fh->getline) {
        my ($chromosome, $position, $reference, undef,$depth1, $depth2, undef, undef, undef,$qual , undef,$consensus, @extra) = split("\t", $line);
        no warnings qw(numeric);
        next unless $position eq int($position); #Skip header line(s)
        use warnings qw(numeric);
        
        #position => 1-based position of the SNV
        #BED uses 0-based position of and after the event
        my $depth = $depth1+$depth2;
        $self->write_bed_line($chromosome, $position - 1, $position, $reference, $consensus, $qual, $depth);
    }
    
    return 1;
}

1;
