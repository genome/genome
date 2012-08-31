package Genome::Model::Tools::Bed::Convert::Snv::AnnotationToBed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::Convert::Snv::AnnotationToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Snv'],
};

sub help_brief {
    "Tools to convert Annotation SNV format to BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert snv annotation-to-bed --source snps_all_sequences --output snps_all_sequences.bed
EOS
}

sub help_detail {
    return <<EOS
    This is a small tool to take SNV calls in annotation format and convert them to a common BED format (using the first five columns).
EOS
}

sub process_source {
    my $self = shift;

    my $input_fh = $self->_input_fh;

    while(my $line = <$input_fh>) {
        chomp $line;
        my ($chromosome, $start, $stop, $reference, $consensus, @extra) = split("\t", $line);
        unless($start =~ /^\d+$/) {
            next; #skip header lines
        }
        unless($start eq $stop) {
            $self->error_message('Start and stop positions do not match in line: ' . $line);
            return;
        }

        #position => 1-based position of the SNV
        #BED uses 0-based position of and after the event
        $self->write_bed_line($chromosome, $start - 1, $stop, $reference, $consensus);
    }

    return 1;
}

1;
