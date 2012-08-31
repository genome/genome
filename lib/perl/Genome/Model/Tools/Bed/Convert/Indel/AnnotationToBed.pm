package Genome::Model::Tools::Bed::Convert::Indel::AnnotationToBed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::Convert::Indel::AnnotationToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Indel'],
};

sub help_brief {
    "Tools to convert Annotation indel format to BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert indel annotation-to-bed --source annotated_indels --output annotated_indels.bed
EOS
}

sub help_detail {                           
    return <<EOS
    This is a small tool to take indel calls in annotation format and convert them to a common BED format (using the first five columns).
EOS
}

sub process_source {
    my $self = shift;

    my $input_fh = $self->_input_fh;

    while(my $line = <$input_fh>) {
        chomp $line;
        my ($chromosome, $start, $stop, $reference, $consensus, $type, @extra) = split("\t", $line);
        # convert from 1-based to 0 based
        if(defined($type)){
            if ($type eq "INS") {
                --$stop;
            } elsif ($type eq "DEL") {
                --$start;
            } else {
                next;
            }
        } else {
            if ($reference =~ /0|\-|\*/){ ##INS
                --$stop;
            } elsif ($consensus =~ /0|\-|\*/){ ##DEL
                --$start;
            } else {
                $self->error_message("not recognized as an indel: $reference $consensus");
                next
            }
        }
        #my $depth = 0;
        #my $qual = 0;
        #$self->write_bed_line($chromosome, $start, $stop, $reference, $consensus, $depth, $qual);
        $self->write_bed_line($chromosome, $start, $stop, $reference, $consensus);
    }

    return 1;
}

1;
