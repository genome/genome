package Genome::Utility::Copy;

use strict;
use warnings FATAL => 'all';

use Exporter qw(import);
our @EXPORT_OK = qw(copy);

sub copyable_properties {
    my $source_class = shift;
    if ($source_class->can('copyable_properties')) {
        return $source_class->copyable_properties();
    }

    my $meta = $source_class->__meta__;
    return grep { !$_->is_delegated && !$_->is_id } $meta->properties;
}

sub copy {
    my $source = shift;

    my @copyable_properties = copyable_properties($source->class);
    my %params = map {
        my $name = $_->property_name;
        my $value = $source->$name;
        defined $value ? ($name => $value) : ();
    } @copyable_properties;

    return $source->class->create(%params);
}
