package Genome::Model::Tools::TestHelpers::General;

use strict;
use warnings;

use Params::Validate qw(:types);
use above 'Genome';
use Genome::Utility::Test 'compare_ok';
use Test::More;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    get_test_dir
    get_test_data
);

sub get_test_data {
    my ($class, $data_path, $version) = Params::Validate::validate_pos(@_,
        1, 1, 1);

    my $test_dir = get_test_dir($class, $version, 1);
    my $test_data = File::Spec->join($test_dir, $data_path);

    if (-e $test_data) {
        note "Found test data: ($test_data)";
    } else {
        die "Couldn't find test data: ($test_data)";
    }

    return $test_data;
}

sub get_test_dir {
    my ($class, $version, $silent) = Params::Validate::validate_pos(@_,
        1, 1, {default => 0});

    my $test_dir = _test_dir($class, "v$version");
    note("Test data located at ($test_dir)\n") unless $silent;

    return $test_dir;
}

sub _test_dir {
    my ($package, $test_version) = @_;

    (my $dirname = $package) =~ s/::/-/g;
    my $test_dir = File::Spec->join(Genome::Config::get('test_inputs'),
        $dirname, $test_version);

    return $test_dir;
}


