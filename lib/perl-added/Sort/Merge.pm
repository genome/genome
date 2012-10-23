#!/usr/bin/perl
# Copyright 2003, James Mastros, <james@mastros.biz>
# Released under the same terms as perl itself.
# See bottom for full copyright info.
package Sort::Merge;

use warnings;
use strict;
use Data::Dumper;

our $VERSION = 0.01;

# "Big" data structure is sources.
# This is an array(ref) of arrayrefs.
# Each inner array has three elements.
# The first is a coderef, which should return the next two elements when called.
# The second is the sort key for the next data point.
# The third is the data point itself.

# Given a sources array, makes sure that each has it's next datapoint loaded,
# and returns the new sources list.
sub prime {
  my @sources=@{shift @_};
  my $indexes=shift;
  $indexes ||= [0..$#sources];
  my @indexes=@$indexes;
  
#  print "In prime: ", Dumper(\@sources);
  foreach my $source (@sources[@$indexes]) {
	if (not defined $source->[1]) {
#	  print "Priming.\n";
	  my ($key, $data) = $source->[0]->();
	  if (!defined $key) {
#		print "Source dead.\n";
		# We can't just delete it here, so mark it undef, and filter them out later.
		$source=undef;
	  }
#	  print "Key: $key, Data: $data\n";
	  $source->[1] = $key;
	  $source->[2] = $data;
	}
  }
  
#  print "Before purge: ", Dumper(\@sources);
  @sources = grep {defined $_->[0]} @sources;
#  print "After purge: ", Dumper(\@sources);

  return @sources;
}

# Given a sources arrayref, and a coderef to call with each data point, sort.
sub sort_inner {
  my $sources=shift;
  my $output=shift;

  my @sources=@$sources;
  my $lastsource;
  while (@sources) {
	@sources=prime(\@sources, $lastsource);
	last if !@sources;

#	print Dumper \@sources;
	
	my($firstsource,$earliesttime)=(0,~0);
	
	foreach (0..$#sources) {
	  if ($earliesttime > $sources[$_][1]) {
		$earliesttime=$sources[$_][1];
		$firstsource=$_;
	  }
	}
	$output->($sources[$firstsource]);
	$sources[$firstsource][1] = undef;
	$sources[$firstsource][2] = undef;
	$lastsource=[$firstsource];
  }
}

# Given an arrayref of input coderefs, and a single output coderef, sort.
sub sort_coderefs {
  my $inputs=shift;
  my $output=shift;
  my @sources;

  @sources = map {[$_]} @$inputs;
  return sort_inner(\@sources, $output);
}

=head1 NAME

Sort::Merge - general merge sort

=head1 SYNOPSIS

   use Sort::Merge;
   use IO::File;
   my @output;
   Sort::Merge::sort_coderefs([sub{our $i=0; return ($i++ x 2)},
                               sub{our $j=0; return ($j++*5 x 2)}],
                              sub{print shift->[1});

=head1 DESCRIPTION

Merge sorting is a technique for merging together multiple input
streams, each of which is sorted, into one larger output stream, still
sorted.  This module implements that algorithm.

Note the large caveat above: the inputs must already be sorted.  If
the inputs aren't sorted, you'll need to sort them yourself first, in
which case a merge sort may not be the right algorithm for the job.

Sort::Merge was designed to give you, the programmer, control over
where your data comes from, and where it goes -- it's simply supposed
to sort it, and not fetch it from files, parse out sort keys, etc,
etc.  There's at least one merge-sorting module already on CPAN,
File::MergeSort, but I wanted one that I could use to sort together
multiple files that were of different formats, and not line-oriented
formats either.

The module I got out doesn't even require that the inputs are files at
all, or that the output gets printed.  It also streams the output, so
you needn't wait until the sort is complete before printing the first
piece of output.  (The fact that this is possible at all is one of
more useful properties of merge sorts vs. other sorts.)

Most (OK, at the moment all) of the available interfaces require you to
write your input and output handlers in terms of coderefs.  If this is
too much work for you, I encourage you to not use this module, or,
alternatively, to propose another interface to me.  I'm actively
looking for more and better interfaces.

=head2 Sources

There are two types of coderefs Sort::Merge makes you write, sources
and the output.

Sources are sources of input data.  They are called with no
parameters.  They are expected to return a two-element list if they
have more data.

The first element should be a sort key, which will be compared
lexographicly (using cmp).  If you wish to think in terms of numerical
sort keys, simply run it through chr, pack, or sprintf to get a
representation that will sort properly.  (Passing it to chr requires
that the number is a nonnegative integer, and may warn for some
values, but is quite fast.)

The second element may be any arbitrary scalar.  This is your
datapoint.  It is passed to the output subroutine without being
examined.  It can be a line of text, an arrayref or hashref, or even
an object, or undef for all Sort::Merge cares.

If the source has run out of data, it should return nothing (an empty
list, not C<undef>).  If it wishes to signal an error and abort
processing, it can simply C<die>; the error will be propagated.

=head2 Output

Because I'm lazy, the output callback gets passed the source
structure.  That is, it gets an arrayref of three values.  The first
is the coderef, the second is the key, the third is the data.

To put it another way, it gets an arrayref of the coderef, followed by
what the coderef has returned.

=head2 C<sort_coderefs>

C<sort_coderefs> takes an array-ref of source coderefs, and a single output coderef.

The output coderef is called for each element, in sorted order, and
the function returns empty-list.

=head1 SEE ALSO

L<File::MergeSort>, the source code.

=head1 COPYRIGHT AND DISCLAIMERS

Copyright (c) 2003 James M. Mastros.  All rights reserved.

This library is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.

=head1 AUTHOR

James Mastros E<lt>james@mastros.bizE<gt>, L<theorbtwo on perlmonks|http://perlmonks.org/?node=theorbtwo>

=cut
