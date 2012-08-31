package Finfo::ClassUtils;

use strict;
use warnings;

use Carp;
use Data::Dumper;

require Scalar::Util;

my @exported_subs = (qw/ anon_scalar class ident use_class /);

sub import
{
    my $self = shift;
    my $caller = caller;

    my @subs_to_export = ( @_ )
    ? split(/\s+/, $_[0])
    : @exported_subs;

    for my $sub ( @subs_to_export ) 
    {
        croak "Invalid sub ($sub) to export to $caller" 
        unless grep { $sub eq $_ } @exported_subs;
        
        #next if $caller->can($sub);
        no strict 'refs';
        *{ $caller . '::' . $sub } = \&{ $sub };
    }

    return 1;
}

sub anon_scalar { return \my $scalar; }

sub class
{
    my $self = shift;

    return ref($self) || $self;
}

*ident = \&Scalar::Util::refaddr;

sub use_class {
    my $class = shift;

    croak "No class given" unless $class;

    eval("use $class") if defined $class;
    croak $@ if $@; 

    return 1;
}

1;

=pod

=head1 Name

Finfo::ClassUtils

=head1 Synopsis

Some generic methods exported to your class.  Some methods inspired by Conway's Class::Std.  The importer will warn if overwriting methods already created for a class.

=head1 Usage

I<To use it, use it!>

B<Exports all methods>
 use Finfo::ClassUtils;

B<Exports class method>
 use Finfo::ClassUtils 'class';

B<Exports anon_scalar and class methods>
 use Finfo::ClassUtils 'anon_scalar class';

I<...or require then import!>

=head1 Exported Methods

=head2 anon_scalar

 my $class = shift;
 my $obj = bless anon_scalar(), $class;

=over

=item I<Synopsis>   returns a reference to an anonymous scalar, suitable for blessing as an object

=item I<Params>     none

=item I<Returns>    a reference to an anonymous scalar 

=back

=head2 class

 my $class = $obj->class;

=over

=item I<Synopsis>   Gets the class of the object/class

=item I<Params>     none

=item I<Returns>    class of object 

=back

=head2 ident

 my $ident = $obj->ident;

=over

=item I<Synopsis>   Gets the ref address of the object

=item I<Params>     none

=item I<Returns>    ref address of the object 

=back

=head1 See Also

I<Class::Std::Utils>, I<Scalar::Util>

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

Eddie Belter <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
