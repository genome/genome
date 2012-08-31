package Genome::Model::Tools::Bed::Convert::Indel::MaqToBed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::Convert::Indel::MaqToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Indel'],
};

sub help_brief {
    "Tools to convert MAQ indel format to BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert indel maq-to-bed --source indelpe.out --output indels_all_sequences.bed
EOS
}

sub help_detail {                           
    return <<EOS
    This is a small tool to take indel calls in MAQ format and convert them to a common BED format (using the first four columns).
EOS
}

sub process_source {
    my $self = shift;
    
    my $input_fh = $self->_input_fh;
    
    while(my $line = <$input_fh>) {
        my ($chromosome, $position, $quality, $depth, $length, $bases, @extra) = split(/[\t:]/, $line);
        
        unless($quality =~ m/[*+-.]/) {
            $self->error_message('The file does not appear to be the output from `maq indelpe`. (Encountered unexpected quality value: ' . $quality . ') This converter does not support `maq indelsoa` output.');
            return;
        }
        
        my ($minus) = $length =~ s/-//;
        
        my ($reference, $variant, $start, $stop);
        
        $start = $position - 1; #Convert to 0-based coordinate
        
        if($minus) {
            $reference = $bases;
            $variant = '*';
            $stop = $start + $length;
        } else {
            $start -= 1; #Want to output the base before the insertion
            $reference = '*';
            $variant = $bases;
            $stop = $start; #Two positions are included-- but an insertion has no "length" so stop and start are the same
        }
        
        $self->write_bed_line($chromosome, $start, $stop, $reference, $variant, $quality, $depth);
    }
    
    return 1;
}

1;
