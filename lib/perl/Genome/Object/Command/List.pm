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
        if ( ! defined $left) {
            $left = $w
        } elsif ((length($left) + length($w)) >= $width) {
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

    my $doc_width = $term_width - $indent - 1;

    my ($left, $right) = _split_width($term_width - 1, $doc_line);
    my @lines = ($left);
    while (defined $right) {
        ($left, $right) = _split_width($doc_width, $right);
        push @lines, (' ' x $indent) . $left;
    }

    return @lines;
}

1;
