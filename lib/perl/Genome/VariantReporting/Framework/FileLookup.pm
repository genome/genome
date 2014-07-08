package Genome::VariantReporting::Framework::FileLookup;

use strict;
use warnings FATAL => 'all';
use Genome;
use Cwd qw(abs_path);

use Exporter 'import';

our @EXPORT_OK = qw(
    is_file
    calculate_lookup
);

sub is_file {
    my $maybe_file = shift;

    return if !defined($maybe_file);

    if (ref $maybe_file) {
        return;
    }

    {
        no warnings 'newline';
        return -f $maybe_file;
    }
}

sub calculate_lookup {
    my $file = shift;

    my $fullpath = abs_path($file);

    if (my $allocation = Genome::Disk::Allocation->get_allocation_for_path($fullpath)) {
        my $allocation_absolute_path = $allocation->absolute_path;
        my $allocation_id = $allocation->id;
        (my $modified_path = $fullpath) =~ s/$allocation_absolute_path/$allocation_id/;

        return $modified_path;
    } else {
        return Genome::Sys->md5sum($fullpath);
    }
}

1;
