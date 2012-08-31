package Genome::Model::Tools::Fastq::Sol2phred;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Fastq::Sol2phred {
    is => 'Genome::Model::Tools::Fastq',
    has => [
        phred_fastq_file => {
            is => 'Text',
            is_optional => 1,
            shell_args_position => 1,
            doc => 'The output fastq file for phred quality sequences',
        },
    ],
    doc => 'convert an illumina/solexa fastq into a sanger/phred fastq'
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless $self;

    unless ($self->phred_fastq_file) {
        $self->phred_fastq_file($self->fastq_file .'.phred');
    }

    return $self;
}

# prints solexa ascii characters for 0-40
#perl -e 'for (0 .. 40) { print  chr($_ + 64); }'
# prints phred ascii characters for 0-40
#perl -e 'for (0 .. 40) { print  chr($_ + 33); }'

#For the new Solexa-Phred qualities with an offset of 64, the equation
#simplifies to
#  $fastq = chr(ord($solq) - 64 + 33);
#or just
#  $fastq = chr(ord($solq) - 31);

sub execute {
    my $self = shift;
    my $reader = Genome::Sys->open_file_for_reading($self->fastq_file);
    binmode $reader, ":utf8";
    my $writer = Genome::Sys->open_file_for_writing($self->phred_fastq_file);
    binmode $writer, ":utf8";
    while (my $line = $reader->getline) {
        chomp($line);
        if ($line =~ /^\+/) {
            # print the quality read name line
            print $writer $line ."\n";

            # parse the solexa quality data line and interpolate to phred quality values
            my $qual_line = $reader->getline;
            chomp($qual_line);
            my @sol_quals = unpack("C*", $qual_line);
            print $writer (join "", map {chr($_-31)} @sol_quals), "\n";
        } else {
            # print the sequence read name and sequence data lines
            print $writer $line ."\n";
        }
    }
    $writer->close;
    $reader->close;
    return 1;
};
