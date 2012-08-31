package Genome::Model::Tools::Fastq::Sol2sanger;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Fastq::Sol2sanger {
    is => 'Genome::Model::Tools::Fastq',
    has => [
        sanger_fastq_file => {
            is => 'Text',
            is_optional => 1,
            shell_args_position => 1,
            doc => 'The output fastq file for sanger quality sequences',
        },
    ],
    doc => 'convert an old solexa fastq into a sanger fastq (using the method of `maq sol2sanger`'
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless $self;

    unless ($self->sanger_fastq_file) {
        $self->sanger_fastq_file($self->fastq_file .'.phred');
    }

    return $self;
}

## Conversion table copied from MAQ's fq_all2std.pl

# Solexa->Sanger quality conversion table
my @conv_table;
for (-64..64) {
  $conv_table[$_+64] = chr(int(33 + 10*log(1+10**($_/10.0))/log(10)+.499));
}

sub execute {
    my $self = shift;
    my $reader = Genome::Sys->open_file_for_reading($self->fastq_file);
    binmode $reader, ":utf8";
    my $writer = Genome::Sys->open_file_for_writing($self->sanger_fastq_file);
    binmode $writer, ":utf8";
    while (my $line = $reader->getline) {
        chomp($line);
        if ($line =~ /^\+/) {
            # print the quality read name line
            print $writer $line ."\n";

            # parse the solexa quality data line and interpolate to sanger quality values
            my $qual_line = $reader->getline;
            chomp($qual_line);
            my @sol_quals = unpack("C*", $qual_line);
            print $writer (join "", map { $conv_table[$_] } @sol_quals), "\n";
        } else {
            # print the sequence read name and sequence data lines
            print $writer $line ."\n";
        }
    }
    $writer->close;
    $reader->close;
    return 1;
};

1;
