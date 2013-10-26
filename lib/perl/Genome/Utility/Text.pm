package Genome::Utility::Text;

use strict;
use warnings;

use POSIX "floor";
use List::Util 'max';
require Carp;
use Params::Validate qw(:types);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(camel_case_to_string
                    capitalize_words
                    hash_to_string
                    justify
                    param_string_to_hash
                    sanitize_string_for_filesystem
                    find_diff_pos
                    side_by_side
                    string_to_camel_case
                    strip_color
                    tree_to_condensed_string
                    tree_to_string
                    width
                    );

#< Camel Case >#
sub string_to_camel_case {
    my $string = shift;
    unless ( $string ) {
        Carp::cluck('No string to convert to camel case');
        return;
    }
    my $split_chars = $_[0] ? $_[0] : qr([\s_]+);
    return join('', map { ucfirst } split($split_chars, $string));
}

sub camel_case_to_string {
    my $camel_case = shift;
    unless ( $camel_case ) {
        Carp::cluck('No camel case to get words');
        return;
    }
    my @words = split( /(?=(?<![A-Z])[A-Z])|(?=(?<!\d)\d)/, $camel_case); #split on the first capital or the start of a number
    my $join = ( @_ ) ? $_[0] : ' ';
    return join($join, map { lc } @words);
}

#< Params as String and Hash >#
sub param_string_to_hash {
    my ($param_string, $value_split) = @_;

    unless ( $param_string ) {
        Carp::cluck('No param string to convert to hash');
        return;
    }

    unless ($param_string =~ m#^-#) {
        Carp::cluck('Param string must start with a dash (-)');
        return;
    }

    my %params;
    my @params = split(/\s?(\-{1,2}\D[\w\-]*)[\s=]?/, $param_string);
    shift @params;
    for ( my $i = 0; $i < @params; $i += 2 ) {
        my $key = $params[$i];
        $key =~ s/^\-{1,2}//;
        Carp::cluck("Malformed param string ($param_string).  Found empty dash (-).") if $key eq '';
        my $value = $params[$i + 1];
        #$params{$key} = ( $value ne '' ? $value : 1 );
        $value =~ s/\s*$// if defined $value;
        if ( not defined $value or $value eq '' ) {
            $params{$key} = 1;
        }
        elsif ( defined $value_split ) {
            $params{$key} = [ split($value_split, $value) ];
        }
        else {
            $params{$key} = $value;
        }
    }

    return %params;
}

sub hash_to_string {
    my ($hash, $delimiter) = Params::Validate::validate_pos(
        @_, { type => HASHREF, },
            {
                type => SCALAR,
                default => q('),
                regex => qr/^('|"|q|qq)$/,
            }
    );

    my @params;
    for my $key (sort keys %$hash) {
        next unless defined $hash->{$key};
        push @params, "$key => " . wrap_as_string($hash->{$key}, $delimiter);
    }
    return join(",", @params);
}

sub wrap_as_string {
    my ($value, $delimiter) = Params::Validate::validate_pos(
        @_, 1,
            {
                type => SCALAR,
                default => q('),
                regex => qr/^('|"|q|qq)$/,
            }
    );

    if ($delimiter =~ /q/) {
        return sprintf("%s(%s)", $delimiter, $value);
    } else {
        return sprintf("%s%s%s", $delimiter, $value, $delimiter);
    }
}

#< Sanitize for File System >#
sub sanitize_string_for_filesystem {
    my $string = shift;
    unless ( $string ) {
        Carp::cluck('No string to sanitize for filesystem.');
        return;
    }
    my $OK_CHARS = '-a-zA-Z0-9_.';
    $string =~ s/[^$OK_CHARS]/_/go;
    return $string;
}

#< Capitalize >#
sub capitalize_words {
    my $string = shift;

    unless ( $string ) {
        Carp::confess('No string to capitalize words.');
    }

    my $seps = join('', ' ', @_); # allow other separators
    my $regexp = qr/[$seps]+/;

    return join(' ', map { ucfirst } split($regexp, $string));
}

sub find_diff_pos {
    my ($s1, $s2) = @_;

    die "You must supply two strings to 'find_diff_pos" unless defined($s1) and defined($s2);

    my $mask = $s1 ^ $s2;

    my @diff_pos;
    while ($mask =~ /[^\0]/g) {
        push @diff_pos, $-[0];
    }
    return @diff_pos;
}

# put text side_by_side
# args:
#     strings => ARRAY of strings
# kwargs:
#     separator => String
#     justification => String or ARRAY of Strings
#           must be one of 'left', 'right' or 'center'
#     fill => String or ARRAY of Strings
#     max_width => Integer, lines longer than this will be tucked.
#           if = -1 no tucking will occur.
#     stack => Boolean, if true side_by_side will stack instead of tuck.
# Tucking:
#     If putting the strings side-by-side would result in strings longer
#     than <max_width>, then 1 or more strings will be tucked under
#     their leftward neighbor instead of placed to the right of them.
sub side_by_side {
    my $strings = shift;
    my @strings = @{$strings};
    my %params = @_;

    my $separator = delete $params{separator} || " ";
    my $justification = delete $params{justification} || "left";
    my $fill = delete $params{fill} || " ";
    my $max_width = delete $params{max_width} || -1;
    my $stack = delete $params{stack} || 0;
    if(keys %params) {
        Carp::croak("Unsupported optional argument(s) supplied to " .
            "side_by_side:" . join(", ", keys %params));
    }

    my @justifications;
    if(ref($justification)) {
        @justifications = @{$justification};
        if(scalar(@justifications) != scalar(@strings)) {
            Carp::croak("You must supply as many justifications as strings.");
        }
    } else {
        @justifications = map {$justification || 'left'} (0..$#strings);
    }

    my @fills;
    if(ref($fill)) {
        @fills = @{$fill};
        if(scalar(@fills) != scalar(@strings)) {
            Carp::croak("You must supply as many fills as strings.");
        }
    } else {
        @fills = map {$fill || ' '} (0..$#strings);
    }

    my @line_arrays;
    my @widths;
    for my $i (0..$#strings) {
        my $string = $strings[$i];

        my $width = width($string);
        my @lines = split("\n", $string);
        push(@line_arrays, \@lines);

        my $fill = $fills[$i];
        $width += width($fill)*2 + 1 if $fill ne ' ';

        push(@widths, $width);
    }

    my @combined_lines;
    my @num_lines = map {scalar(@{$_})} @line_arrays;
    my $max_i = (max(@num_lines) || 1) - 1;
    for my $i (0..$max_i) {
        my $combined_line = '';
        for my $j (0..$#line_arrays) {
            my $array_ref = $line_arrays[$j];
            my $this_part = $array_ref->[$i];
            $this_part = '' unless defined($this_part);

            my $width = $widths[$j];
            my $justification = $justifications[$j];
            my $fill = $fills[$j];

            my $tucking = '';
            my $new_length = $width + width($combined_line) +
                             width($separator);
            if($new_length > $max_width and $max_width != -1) {
                if($stack) {
                    return join("\n\n", @strings);
                }
                $tucking = "\n" . '  'x($new_length/$max_width) . '`-> ';
            }

            my $next_final_length = $width + width($combined_line) +
                                    ($widths[$j+1] || 0) +
                                    2*width($separator);
            if($next_final_length > $max_width) {
                $fill = undef; #don't fill if will tuck next line.
            }

            $this_part = justify($this_part, $justification, $width, $fill);
            $combined_line .= $tucking . $this_part;
            $combined_line .= $separator unless $j == $#line_arrays;
        }
        push(@combined_lines, $combined_line);
    }
    return join("\n", @combined_lines);
}

sub strip_color {
    my ($string) = @_;
    $string =~ s/\e\[[\d;]*[a-zA-Z]//g;
    return $string;
}

# return the VISIBLE width of the of a string.
sub width {
    my ($string) = @_;
    $string = strip_color($string);
    my @lines = split(/\n/, $string);

    my $width = length($string);
    if(scalar(@lines)) {
        $width = max(map {length($_ . '')} @lines);
    }
    return $width;
}

sub justify {
    my ($string, $kind, $field_width, $fill, $spacer) = @_;
    $fill = $fill || " ";
    $spacer = " " unless defined($spacer);

    my $num_spaces_needed = $field_width - width($string);
    my $fill_string = "$fill"x$num_spaces_needed;

    if($num_spaces_needed > 0) {
        if(width($spacer) > $num_spaces_needed) {
            $spacer = substr($spacer, 0, $num_spaces_needed);
        }
        $num_spaces_needed -= width($spacer);
        my $result = $string . $spacer;
        if($num_spaces_needed == 0) {
            return $result;
        }
    } else {
        return $string;
    }

    my $format;
    if($kind eq 'left') {
        $format = "%s" . $spacer .
                substr($fill_string, -1 * $num_spaces_needed);
    } elsif($kind eq 'right') {
        $format = substr($fill_string, 0, $num_spaces_needed ) .
                $spacer . "%s";
    } elsif($kind eq 'center') {
        my $left_spaces = floor($num_spaces_needed/2);
        my $right_spaces = $num_spaces_needed - $left_spaces;
        $format = substr($fill_string, 0, $left_spaces) .  "%s" .
                  substr($fill_string, -1 * $right_spaces)
    } else {
        Carp::croak("kind argument must be one of 'left',".
                    " 'right', or 'center', not '$kind'");
    }
    return sprintf($format, $string);
}

sub _next_foundation {
    my ($foundation, $prefix) = @_;

    return '' unless defined($prefix);
    if($prefix =~ m/\|/) {
        return $foundation . "| ";
    } else {
        return $foundation . "  ";
    }
}

sub tree_to_condensed_string {
    return tree_to_string(shift, shift || '', shift || '', 1);
}

sub tree_to_string {
    my ($tree, $foundation, $prefix, $condensed) = @_;
    $foundation = $foundation || '';
    $prefix = $prefix || '';

    my $result = '';
    my $ref = ref($tree);
    if($ref eq 'HASH' or $ref eq 'ARRAY') {
        my $root;
        my @items;
        if($ref eq 'HASH') {
            $root = "%";
            @items = sort(keys %{$tree});
        } else {
            $root = "@";
            @items = @{$tree};
        }
        my $next_foundation = $foundation;
        unless($condensed) {
            $result .= tree_to_string($root, $foundation, $prefix, $condensed);
            $next_foundation = _next_foundation($foundation, $prefix);
        }
        my $count = 0;
        for my $item (@items) {
            my $next_prefix = $count == $#items ? '`-' : '|-';
            $count += 1;
            $next_prefix = '  ' if $condensed;

            $result .= tree_to_string($item, $next_foundation,
                    $next_prefix, $condensed);

            if($ref eq 'HASH') {
                my $child_foundation = _next_foundation($next_foundation,
                        $next_prefix);
                # can only have a single 'value' for each key.
                $result .= tree_to_string($tree->{$item}, $child_foundation,
                        '`-', $condensed);
            }
        }
    } else {
        $result .= sprintf("%s%s%s\n", $foundation, $prefix, $tree);
    }
    return $result;
}

sub table_to_delim_string {
    my ($delimiter, $array) = Params::Validate::validate_pos(
        @_, { type => SCALAR },
            { type => ARRAYREF },
    );

    return join("\n",
        map { join($delimiter, @$_) } @$array
    );
}

sub table_to_tab_string {
    table_to_delim_string("\t", @_);
}

1;
