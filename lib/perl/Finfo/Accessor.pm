package Finfo::Accessor;
use base 'Class::Accessor';
use strict;

=head1 NAME

Finfo::Accessor::Fast - Faster, but less expandable, accessors, modified from Class::Accessor::Fast to include an accessors method, which returns an array of accessor names for a class.

=head1 SYNOPSIS

  package Foo;
  use base qw(Finfo::Accessor::Fast);

  # The rest as Class::Accessor except no set() or get().

=head1 DESCRIPTION

This is a somewhat faster, but less expandable, version of
Class::Accessor.  Class::Accessor's generated accessors require two
method calls to accompish their task (one for the accessor, another
for get() or set()).  Class::Accessor::Fast eliminates calling
set()/get() and does the access itself, resulting in a somewhat faster
accessor.

The downside is that you can't easily alter the behavior of your
accessors, nor can your subclasses.  Of course, should you need this
later, you can always swap out Class::Accessor::Fast for
Class::Accessor.

=cut

use Finfo::Logging;

my %accessor_names;

sub new
{
    my ($class, %p) = @_;

    return bless \%p, $class;
}

sub make_accessor {
    my($class, $field) = @_;
   
    push @{$accessor_names{$class}}, $field;

    return sub {
        my $self = shift;
        return $self->{$field} unless @_;
        $self->{$field} = (@_ == 1 ? $_[0] : [@_]);
    };
}


sub make_ro_accessor {
    my($class, $field) = @_;

    push @{$accessor_names{$class}}, $field;

    return sub {
        return $_[0]->{$field} unless @_ > 1;
        my $caller = caller;
        require Carp;
        Carp::croak("'$caller' cannot alter the value of '$field' on ".
                    "objects of class '$class'");
    };
}


sub make_wo_accessor {
    my($class, $field) = @_;

    push @{$accessor_names{$class}}, $field;

    return sub {
        my $self = shift;

        unless (@_) {
            my $caller = caller;
            require Carp;
            Carp::croak("'$caller' cannot access the value of '$field' on ".
                        "objects of class '$class'");
        }
        else {
            return $self->{$field} = (@_ == 1 ? $_[0] : [@_]);
        }
    };
}

sub accessors {
    my $class = shift;
    return @{$accessor_names{$class}};
}


=head1 EFFICIENCY

L<Class::Accessor/EFFICIENCY> for an efficiency comparison.

=head1 CURRENT AUTHOR

Marty Pauley <marty+perl@kasei.com>

=head1 ORIGINAL AUTHOR

Michael G Schwern <schwern@pobox.com>

=head1 SEE ALSO

L<Class::Accessor>

=cut

1;

#$HeadURL$
#$Id$

