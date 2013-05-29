package Genome::Object::Command::List;

use strict;
use warnings;

use Genome;
use List::Util qw(reduce);
use Term::ReadKey qw(GetTerminalSize);
use Text::ParseWords qw(quotewords);

# This is mostly here so that we aren't making UR dependent on Term::ReadKey,
# which is an XS module. It currently just extends it's parent class by
# wrapping some text lines better. To do it more completely will require
# further refactoring of the parent class (actually all the way up to
# Command::V1 and Command::V2).

class Genome::Object::Command::List {
    is => 'UR::Object::Command::List',
};

# this is a copy of what Tom has put into UR::Object::Command::List
# remove it after that goes stable (redundant)
sub _properties_for_class_to_document {
    my $self = shift;
    my $target_class_name = shift;

    my $target_class_meta = $target_class_name->__meta__;
    my @id_by = $target_class_meta->id_properties;

    my @props = $target_class_meta->properties;

    no warnings;
    # These final maps are to get around a bug in perl 5.8 sort
    # involving method calls inside the sort sub that may
    # do sorts of their own
    return 
        map { $_->[1] }
        sort { $a->[1]->position_in_module_header <=> $b->[1]->position_in_module_header or $a->[0] cmp $b->[0] }
        map { [ $_->property_name, $_ ] }
        grep {
            substr($_->property_name, 0, 1) ne '_'
            and not $_->implied_by
            and not $_->is_transient
            and not $_->is_deprecated
        }
        @props;
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless $self;
    my $show = $self->show;
    if (defined($show) and $show =~ /\bdelete\b/) {
        $self->delete;
        die "please don't use this to delete things\n";
    }
    return $self;
}

sub _format_property_doc_data {
    my ($class, @data) = @_;

    my @lines = $class->SUPER::_format_property_doc_data(@data);

    my @names = map { $_->[0] } grep { ref $_ } @data;
    my $longest_name = reduce { length($a) > length($b) ? $a : $b } @names;
    my $indent = length($longest_name) + 5;

    my ($term_width) = GetTerminalSize();
    @lines = map { _reformat_doc_line_for_term_width($term_width, $indent, $_) } @lines;

    return @lines;
}

# Split a string at given width without breaking words.
sub _split_width {
    my ($width, $line) = @_;

    my @w = quotewords(' ', 1, $line);

    my ($left, $right);
    while (@w) {
        my $w = shift @w;
        $w = '' unless (defined $w);
        if ( ! defined $left) {
            $left = $w
        } elsif ((length($left) + length($w)) >= $width) {
            no warnings 'uninitialized';  # some remaining values of @w might be undef
            $right = join(' ', $w, @w);
            last;
        } else {
            $left = join(' ', $left, $w);
        }
    }

    return ($left, $right);
}

sub _reformat_doc_line_for_term_width {
    my $term_width = shift;
    my $indent = shift;
    my $doc_line = shift;

    return $doc_line unless ($term_width && $indent);

    my $doc_width = $term_width - $indent - 2;

    my ($left, $right) = _split_width($term_width - 2, $doc_line);
    my @lines = ($left);
    while (defined $right) {
        ($left, $right) = _split_width($doc_width, $right);
        push @lines, (' ' x $indent) . $left;
    }

    return @lines;
}

1;
