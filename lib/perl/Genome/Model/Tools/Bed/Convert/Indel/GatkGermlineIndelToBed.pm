package Genome::Model::Tools::Bed::Convert::Indel::GatkGermlineIndelToBed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::Convert::Indel::GatkGermlineIndelToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Indel'],
};

sub process_source {
    my $self = shift;

    my $input_fh = $self->_input_fh;

    while(my $line = $input_fh->getline) {
        chomp $line;
        my ($chr,$start,$stop, $vardepth) = split("\t", $line);
        my ($refvar, $ref, $var, $depth_string, $type);
        ($refvar, $depth_string) = split ":", $vardepth;

        # The column containing depth will look like this: OBS_COUNTS[C/A/T]:8/8/9 ... the final number is the total depth
        my ($read_support, $depth) = split "/", $depth_string;

        # Score is not available from this mode of running gatk, which kind of sucks.
        my $score = $read_support;
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
        $self->write_bed_line($chr, $start, $stop, $ref, $var, $score, $depth);
    }
    $input_fh->close;
    return 1;
}

1;
