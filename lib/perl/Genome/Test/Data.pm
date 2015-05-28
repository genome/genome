package Genome::Test::Data;

use strict;
use warnings;
use Genome;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    get_test_file
);

sub get_test_file {
    my $data_set = shift;
    my $file_name = shift;
    my $file = File::Spec->join(__FILE__ . '.d', $data_set, $file_name);
    unless (-e $file) {
        die "Test file $file doen't exist";
    }
    return $file;
}

1;
