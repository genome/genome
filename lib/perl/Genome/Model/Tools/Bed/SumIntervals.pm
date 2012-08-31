package Genome::Model::Tools::Bed::SumIntervals;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::SumIntervals {
    is => ['Genome::Model::Tools::Bed'],
    has_input => [
        bed_file => {
            is => 'Text',
            doc => 'The BED format file to sum interval sizes of.',
        },
        one_based_start => {
            is => 'Boolean',
            doc => 'A true BED file uses zero-based start positions.  However, typically one-based start positions are used.',
            default_value => 1,
        }
    ],
};

sub execute {
    my $self = shift;

    my $bed_fh = Genome::Sys->open_file_for_reading($self->bed_file);
    my $offset = 0;
    unless ($self->one_based_start) {
        $offset = 1;
    }
    my $sum;
    while (my $line = $bed_fh->getline) {
        chomp($line);
        if ($line =~ /^#/) { next; }
        my @entry = split("\t", $line);
        $sum += ($entry[2] - $entry[1] - $offset) ;
    }
    print $sum ."\n";
}
