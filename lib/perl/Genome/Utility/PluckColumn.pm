package Genome::Utility::PluckColumn;

use strict;
use warnings;

use List::Util qw(first);

sub pluck_column_from_class {
    my $class_name = shift;

    my $column_name = _get_column_name_from_args($class_name, @_);
    my $table_name = _get_table_name_from_class($class_name);
    my $dbh = _get_dbh_from_class($class_name);

    my $query = _build_query($column_name, $table_name, $dbh);

    return $dbh->selectcol_arrayref($query);
}

sub _get_column_name_from_args {
    my $class_name = shift;
    my %opts = @_;

    my $column_name = delete $opts{column_name};
    my $property_name = delete $opts{property_name};

    die("You must specifiy either a column name or a property name!") unless $column_name || $property_name;
    my $query_column = $column_name ? $column_name : _get_column_name_from_property_name_and_class_name($property_name, $class_name);
    die("You specified a property ($property_name) that doesn't correspond to a database column!") unless $query_column;

    return $query_column;
}

sub _get_column_name_from_property_name_and_class_name {
    my $property_name = shift;
    my $class_name = shift;

    die("You must specify a property name and a class name!") unless $property_name && $class_name;

    my $column_name = $class_name->__meta__->property($property_name)->column_name;

    unless ($column_name) {
        my $item = first { $_->property($property_name) && $_->property($property_name)->column_name }
            $class_name->__meta__->ancestry_class_metas;

        $column_name =  $item->property($property_name)->column_name if $item;
    }

    die("There isn't a column for $property_name") unless $column_name;
    return $column_name;
}

sub _get_table_name_from_class {
    my $class_name = shift;

    die("You must specify a class name!") unless $class_name;
    my $table_name = $class_name->__meta__->table_name;

    unless($table_name) {
        my $item = first { $_->table_name }
            $class_name->__meta__->ancestry_class_metas;
        $table_name = $item->table_name;
    }

    die("There isn't a table name for $class_name") unless $table_name;
    return $table_name;
}

sub _get_dbh_from_class {
    my $class_name = shift;

    die("You must specify a class name!") unless $class_name;
    my $dbh = $class_name->__meta__->data_source->get_default_handle;
    die("Failed to get a database handle!") unless $dbh;

    return $dbh;
}

sub _build_query {
    my $column = shift;
    my $table_name = shift;
    my $dbh = shift;

    die('You must specify a column, table name and database handle!')
        unless $column && $table_name && $dbh;

    return $dbh->prepare(sprintf('SELECT %s FROM %s', $column, $table_name));
}

1;
