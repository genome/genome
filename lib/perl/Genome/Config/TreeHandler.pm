package Genome::Config::TreeHandler;
use warnings;
use strict;

use Genome;
use File::Basename;

class Genome::Config::TreeHandler {
    is => 'UR::Object',
    has => [
        base_path => {
            is => 'Text'
        },
        directory_content_hash => {
            is => 'HASH',
            is_optional => 1,
            is_transient => 1,
        },
    ],
    doc => 'handles a directory tree of YML configuration files',
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    eval {
        die("Given base_path " . $self->base_path . " doesn't exist!") unless (-d $self->base_path);
        $self->directory_content_hash(_build_hash_for_directory($self->base_path));
    };
    if (my $error = $@) {
        $self->delete();
        die($error);
    }
    return $self;
}

sub get_leaf {
    my ($self, @args) = @_;
    my $path = _traverse_hash($self->directory_content_hash, @args);
    die('You have entered a query for non-existent configuration!') unless $path;
    return $path
}

sub parameters_exist {
    my ($self, @args) = @_;
    return _traverse_hash($self->directory_content_hash, @args);
}

sub _traverse_hash {
    my ($hash_ref, @args) = @_;

    foreach my $i (0..$#args) {
        my $item = $hash_ref->{$args[$i]};
        if (ref($item) eq 'HASH') {
            splice(@args, $i, 1);
            return _traverse_hash($item, @args);
        } elsif ($item) {
            return $item;
        }
    }
    return;
}

sub _build_hash_for_directory {
    my $directory = shift;
    die('No directory given!') unless $directory;

    my $hash = {};

    my @item_list = (glob("$directory/*"));
    foreach my $item (@item_list) {
        if (-f $item) {
            die("More than one item found in a directory with a config file! $item") if @item_list > 1;
            return $item;
        } else {
            die("$item isn't a directory or a config file!") unless (-d $item);
            my $current_basename = basename($item);
            $hash->{$current_basename} = _build_hash_for_directory($item)
        }
    }
    return $hash;
}

1;
