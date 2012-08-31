package UNIVERSAL::can;

use strict;
use warnings;

use 5.006;

use vars qw( $VERSION $recursing );
$VERSION = '1.12';

use Scalar::Util 'blessed';
use warnings::register;

my $orig;
use vars '$always_warn';

BEGIN
{
	$orig = \&UNIVERSAL::can;

	no warnings 'redefine';
	*UNIVERSAL::can = \&can;
}

sub import
{
	my $class = shift;
	for my $import (@_)
	{
		$always_warn = 1 if $import eq '-always_warn';
		no strict 'refs';
		*{ caller() . '::can' } = \&can if $import eq 'can';
	}
}

sub can
{
	# can't call this on undef
	return _report_warning() unless defined $_[0];

	# don't get into a loop here
	goto &$orig if $recursing;

	# call an overridden can() if it exists
	local $@;
	my $can = eval { $_[0]->$orig('can') || 0 };

	# but not if it inherited this one
	goto &$orig if $can == \&UNIVERSAL::can;

	# make sure the invocant is useful
	unless ( _is_invocant( $_[0] ) )
	{
		_report_warning();
		goto &$orig;
	}

	# redirect to an overridden can, making sure not to recurse and warning
	local $recursing = 1;
	my $invocant = shift;

	_report_warning();
	return $invocant->can(@_);
}

sub _report_warning
{
	if ( $always_warn || warnings::enabled() )
	{
		my $calling_sub = ( caller(2) )[3] || '';
        #warnings::warn("Called UNIVERSAL::can() as a function, not a method")
		#	if $calling_sub !~ /::can$/;
	}

	return;
}

sub _is_invocant
{
	my $potential = shift;
	return unless length $potential;
	return 1 if blessed($potential);

	my $symtable = \%::;
	my $found    = 1;

	for my $symbol ( split( /::/, $potential ) )
	{
		$symbol .= '::';
		unless ( exists $symtable->{$symbol} )
		{
			$found = 0;
			last;
		}

		$symtable = $symtable->{$symbol};
	}

	return $found;
}

1;
__END__

=head1 NAME

UNIVERSAL::can - Hack around people calling UNIVERSAL::can() as a function

=head1 VERSION

Version 1.01

=head1 SYNOPSIS

To use this module, simply:

  use UNIVERSAL::can;

=head1 DESCRIPTION

The UNIVERSAL class provides a few default methods so that all objects can use
them.  Object orientation allows programmers to override these methods in
subclasses to provide more specific and appropriate behavior.

Some authors call methods in the UNIVERSAL class on potential invocants as
functions, bypassing any possible overriding.  This is wrong and you should not
do it.  Unfortunately, not everyone heeds this warning and their bad code can
break your good code.

This module replaces C<UNIVERSAL::can()> with a method that checks to see if
the first argument is a valid invocant (whether an object -- a blessed referent
-- or the name of a class).  If so, and if the invocant's class has its own
C<can()> method, it calls that as a method.  Otherwise, everything works as you
might expect.

If someone attempts to call C<UNIVERSAL::can()> as a function, this module will
emit a lexical warning (see L<perllexwarn>) to that effect.  You can disable it
with C<no warnings;> or C<no warnings 'UNIVERSAL::isa';>, but don't do that;
fix the code instead.

Some people argue that you must call C<UNIVERSAL::can()> as a function because
you don't know if your proposed invocant is a valid invocant.  That's silly.
Use C<blessed()> from L<Scalar::Util> if you want to check that the potential
invocant is an object or call the method anyway in an C<eval> block and check
for failure.

Just don't break working code.

=head1 EXPORT

This module can I<optionally> export a C<can()> subroutine that works exactly
as described.  It's a convenient shortcut for you.  This actually works in
version 1.11.

Also, if you pass the C<-always_warn> flag on the import line, this module will
warn about all incorrect uses of C<UNIVERSAL::can()>.  This can help you change your code to be correct.

=head2 can()

The C<can()> method takes two arguments, a potential invocant and the name of a
method that that invocant may be able to call.  It attempts to divine whether
the invocant is an object or a valid class name, whether there is an overridden
C<can()> method for it, and then calls that.  Otherwise, it calls
C<UNIVERSAL::can()> directly, as if nothing had happened.

=head1 AUTHOR

chromatic, C<< <chromatic@wgz.org> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-universal-can@rt.cpan.org>,
or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=UNIVERSAL-can>.  This will
contact me, hold onto patches so I don't drop them, and will notify you of
progress on your request as I make changes.

=head1 ACKNOWLEDGEMENTS

Inspired by L<UNIVERSAL::isa> by Yuval Kogman, Autrijus Tang, and myself.

Adam Kennedy has tirelessly made me tired by reporting potential bugs and
suggesting ideas that found actual bugs.

Mark Clements helped to track down an invalid invocant bug.

=head1 COPYRIGHT & LICENSE

Copyright (c) 2005 - 2006 chromatic. All rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut
