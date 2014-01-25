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
    get_blessed_file
    ensure_file
    compare_to_blessed_file
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
    my $test_dir = File::Spec->join($ENV{GENOME_TEST_INPUTS},
        $dirname, $test_version);

    return $test_dir;
}

# Find the blessed version of an output:
# i.e. /tmp/foo/bar/baz.out -> /test_data_directory/foo/bar/baz.out
# with test_dir => /test_data_directory    and   output_dir => /tmp
sub get_blessed_file {
    my %params = _validate_for_blessed_lookup(@_);

    my $test_dir = $params{test_dir};
    my $output_path = $params{output_path};
    my $output_dir = $params{output_dir};

    (my $blessed_path = $output_path) =~ s/$output_dir/$test_dir/;

    ensure_file($blessed_path);
    return $blessed_path;
}

sub _validate_for_blessed_lookup {
    my %params = Params::Validate::validate(@_,{
        output_path => {
            type => SCALAR,
        },
        test_dir => {
            type => SCALAR,
        },
        output_dir => {
            type => SCALAR,
            default => '::INFERRED::',
        },
    });
    if ($params{output_dir} eq '::INFERRED::') {
        $params{output_dir} = basename($params{output_path});
    }
    return %params;
}

sub ensure_file {
    my $path = shift;

    unless (-f $path) {
        die "Couldn't find file at ($path)";
    }
}

sub compare_to_blessed_file {
    my %params = _validate_for_blessed_lookup(@_);

    my $blessed_path = get_blessed_file(@_);
    compare_ok($blessed_path, $params{output_path},
        sprintf("blessed-file (%s) matches output-file (%s)",
                $blessed_path, $params{output_path}),
    );
}
