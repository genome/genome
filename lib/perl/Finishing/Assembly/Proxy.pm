package Finishing::Assembly::Proxy;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;

my %source :name(source:r) :isa('object');
my %obj_constructor :name(object_constructor:r) :isa(code);

sub get_method
{
    my ($self, $requested_method, @args) = @_;

	if ( $requested_method eq 'source' )
	{
		 return sub
         {
            $self->source->source(@args);
         };
	}

    if ( $self->can($requested_method) )
    {
        return sub{ $self->$requested_method(@args) };
    }

    return $self->_check_methods($requested_method, @args);
}

sub _check_methods : RESTRICTED
{
    my ($self, $requested_method, @args) = @_;

    my @methods = $self->methods;
    return unless @methods;
    foreach my $method ( @methods )
    {
        my $methods_method = 'methods_for_' . $method;
        my %valid_methods = $self->$methods_method;
        if ( grep { $requested_method eq $_ } keys %valid_methods )
        {
            return sub
            {
                $self->$method
                (
                    method => $requested_method, 
                    params => $valid_methods{$requested_method}, 
                    args => \@args
                );
            };
        }
    }

    return;
}

#- This base proxy will look for source methods that can be executed
#- directly.  This can be an arry or an hasref.  If hashref, the value
#- of the valid method key should be the object type to contruct.  If undef,
#- the value will be returned as is.
sub methods : CUMULATIVE
{
    return (qw/ source_method /);
}

sub methods_for_source_method : CUMULATIVE
{
    return;
}

sub source_method : RESTRICTED
{
    my ($self, %p) = @_;

    my $method = $p{method};
    my $val = $self->source->$method( @{ $p{args} });

    return unless defined $val;

    return $val unless $p{params};

    return $self->_construct_objects($p{params}, $val);
}

sub _construct_objects : RESTRICTED
{
    my ($self, $type, $source) = @_;

    my $aryref = ref($source) eq 'ARRAY';
    my @sources = ( $aryref )
    ? @$source
    : $source;

    my @objects;
    foreach my $source ( @sources )
    {
        push @objects, $self->_construct_object($type, $source);
    }

    return ( $aryref ) ? \@objects : $objects[0];
}

sub _construct_object : RESTRICTED
{
    my ($self, $type, $source) = @_;

    $self->fatal_msg
    (
        "No type to construct object",
        { caller_level => 2 }
    ) unless $type;

    $self->fatal_msg
    (
        "No source(s) to construct object ($type)",
        { caller_level => 2 }
    ) unless $source;

    return $self->object_constructor->
    (
        type => $type, 
        source => $source, 
    );
}

1;

=pod

=head1 Name

Finishing::Assembly::Proxy

=head1 Synopsis

=head1 Usage

=head1 Methods

=head1 See Also

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Proxy.pm $
#$Id: Proxy.pm 31442 2008-01-03 23:47:59Z adukes $

