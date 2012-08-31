package Genome::Model::Tools::Bed::Convert::Indel::SniperToBed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::Convert::Indel::SniperToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Indel'],
};

sub help_brief {
    "Tools to convert sniper indel format to BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert indel sniper-to-bed --source indels_all_sequences --output indels_all_sequences.bed
EOS
}

sub help_detail {                           
    return <<EOS
    This is a small tool to take indel calls in sniper format and convert them to a common BED format (using the first four columns).
EOS
}

sub process_source {
    my $self = shift;
    
    my $input_fh = $self->_input_fh;
    
    while(my $line = <$input_fh>) {
        my ($chromosome, $position, $star, $score, $indel_call_1, $indel_call_2, $indel_length_1, $indel_length_2, @extra) = split("\t", $line);
        # tumor read depth is the 13th field in the file
        my $depth = $extra[4];

        $self->_process_indel($chromosome, $position, $indel_call_1, $indel_length_1, $score, $depth)
            or return;
        $self->_process_indel($chromosome, $position, $indel_call_2, $indel_length_2, $score, $depth)
            or return;
    }
    
    return 1;
}

sub _process_indel {
    my $self = shift;
    my ($chromosome, $position, $indel, $length, $score, $depth) = @_;
    
    return 1 if $indel eq '*'; #Indicates only one indel call...and this isn't it!
    
    #position => 1-based position of the start of the indel
    #BED uses 0-based position of and after the event

    my ($reference, $variant, $start, $stop);
    
    $start = $position; #Convert to 0-based coordinate
    
    if($length > 0) {
        $reference = '*';
        $variant = $indel;
        $stop = $start; #Two positions are included-- but an insertion has no "length" so stop and start are the same
    } elsif($length < 0) {
        $reference = $indel;
        $variant = '*';
        $stop = $start + abs($length);
    } else {
        $self->error_message('Zero length indel encountered.');
        return;
    }

    $self->write_bed_line($chromosome, $start, $stop, $reference, $variant, $score, $depth);
    
    return 1;
}

1;
