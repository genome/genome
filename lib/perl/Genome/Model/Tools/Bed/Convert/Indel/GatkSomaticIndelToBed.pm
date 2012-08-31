package Genome::Model::Tools::Bed::Convert::Indel::GatkSomaticIndelToBed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::Convert::Indel::GatkSomaticIndelToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Indel'],
};


sub process_source {
    my $self = shift;
    
    my $input_fh = $self->_input_fh;
    
    while(my $line = $input_fh->getline) {
        chomp $line;
        my ($chr,$start,$stop, $refvar) = split("\t", $line);
        my ($ref, $var, $type);
        ($refvar) = split ":", $refvar;
        if ($refvar =~ m/\-/) {
            $type = '-';
        }
        else {
            $type = '+';
        }
        $refvar =~ s/[\+,\-]//;
        if ($type eq '+') {
            $ref = '0';
            $var = $refvar;
        }
        else {
            $var = '0';
            $ref = $refvar;
        }
        $self->write_bed_line($chr, $start, $stop, $ref, $var);
    }
    $input_fh->close;
    return 1;
}

1;
