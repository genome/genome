use strict;
use warnings;

package Sort::strverscmp;

use Exporter 'import';
our @EXPORT = qw(strverscmp);
our @EXPORT_OK = qw(strverssort);

sub strverscmp {
    my ($l, $r) = @_;

    my ($ls, $ln, $rs, $rn);
    ($ls, $ln, $l) = decompose_version($l);
    ($rs, $rn, $r) = decompose_version($r);

    my $scmp = ($ls cmp $rs);
    return $scmp if ($scmp != 0);

    # TODO: How can this be refactored?
    if ($ln eq '' && $rn eq '') {
        return 0;
    } elsif ($ln eq '' || $rn eq '') {
        return ($ln eq '' ? -1 : 1);
    }

    my $ncmp = ncmp($ln, $rn);
    return $ncmp if ($ncmp != 0);

    if (length($l) || length($r)) {
        return strverscmp($l, $r);
    } else {
        return 0;
    }
}

sub strverssort {
    return sort { strverscmp($a, $b) } @_;
}

sub decompose_version {
    my ($string, $number, $remainder) = shift =~ /^(\D*)(\d*)(.*)$/;
    return ($string, $number, $remainder);
}

sub ncmp {
    my ($l, $r) = @_;

    if (!is_fractional($l) && !is_fractional($r)) {
        return $l <=> $r;
    } elsif (is_fractional($l) && is_fractional($r)) {
        return fcmp($l, $r);
    } else {
        return (is_fractional($l) ? -1 : 1);
    }
}

sub is_fractional {
    my $n = shift;
    return (index($n, '0') == 0);
}

sub fcmp {
    my ($l, $r) = @_;

    my ($lz, $ln, $rz, $rn);
    ($lz, $ln) = decompose_fractional($l);
    ($rz, $rn) = decompose_fractional($r);

    if (length($lz) == length($rz)) {
        return $ln <=> $rn;
    } else {
        return (length($lz) > length($rz) ? -1 : 1);
    }
}

sub decompose_fractional {
    my ($zeroes, $number) = shift =~ /^(0*)(\d+)$/;
    return ($zeroes, $number);
}

1;

__END__

=head1 NAME

Sort::strverscmp -- Compare strings while treating digits characters numerically.

=head1 SYNOPSIS

  my @list = qw(a A beta9 alpha9 alpha10 alpha010 1.0.5 1.05);
  my @them = strverssort(@list);
  print join(' ', @them), "\n";

Prints:

  1.05 1.0.5 A a alpha010 alpha9 alpha10 beta9

=head1 DESCRIPTION

Pure Perl implementation of GNU strverscmp.

=head1 COPYRIGHT AND DISCLAIMER

Copyright 2013, The Genome Institute at Washington University C<nnutter@cpan.org>, all rights
reserved.  This program is free software; you can redistribute it
and/or modify it under the same terms as Perl itself.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.

=head1 AUTHOR

Nathaniel Nutter C<nnutter@cpan.org>
