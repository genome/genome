package Finfo::Singleton;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;

my %instance;

sub instance 
{
    my $class = shift;
    my @caller = caller;


    if ( ref($class) )
    {
        $class->warn_msg("Attempting to get re-get instance at line $caller[2] in $caller[1]");
        return $class;
    }

    unless ( exists $instance{$class} ) 
    {
        $instance{$class} = Finfo::Std::new($class, @_);
    }
    elsif ( @_ )
    {
        $instance{$class}->fatal_msg
        (
            "Attempting to set params to an existing singleton thru new/instance method at line $caller[2] in $caller[1], use accessor instead."
        );
        return;
    }

    return $instance{$class};
}

sub _enforce_instance : RESTRICTED
{
    my $self = shift;

    unless ( ref $self )
    {
        $self->fatal_msg
        (
            "Need to get instance to continue", 
            { 'caller' => [ caller] }
        );
        return;
    }

    return 1;
}

no warnings 'redefine';

# remove this??
sub new
{
    my $class = shift;
    $class->warn_msg("Use method 'instance' to create/get singleton object.  In the future, this may be a fatl error");
    return $class->instance(@_);
}

1;

=pod

=head1 Name

Finfo::Singleton
 
=head1 Synopsis

A Singleton describes an object class that can have only one instance in any system.  An example of a Singleton might be a print spooler or system registry.  This module implements a Singleton class from which other classes can be derived.  This class provides the management of the instantiation of a single object, attribute handling (see Finfo::Std) and logging (see Finfo::logging).  In deriving a class from Finfo::Singleton, your module will inherit the instance and enforce_instance methods and then can implement whatever specific functionality is required. 

For a description and discussion of the Singleton class, see "Design Patterns", Gamma et al, Addison-Wesley,
1995, ISBN 0-201-63361-2.

=head1 Usage

B<In your class...>

 pakage Good;
 
 use base 'Finfo::Singleton';

 # see Finfo::Std for attribute handling
 my %foo :name(foo:r) :type(string);
 my %bar :name(bar:o) :type(int) :default(20);

 sub say_hello
 {
    my $self = shift;

    $self->_enforce_instance; # make sure we have the instance, not class
    
    $self->info_msg("hello, foo is".$self->foo);

    return 1
 }
 
B<In the code...>

 # Create the good object, sets foo to 'cool', bar to 20 (default)
 my $good = Good->instance(foo => 'cool');
 $good->say_hello;
 make_foo_hot();
 raise_the_bar();

 # add an info msg - prints the msg to the screen and will execute any logging appenders...
 $good->instance->info_msg( sprintf('foo is %s, bar is %d', $good->instance->foo, $good->instance->bar) );
 # or do the same throught the instance method...
 Good->instance->info_msg( sprintf('foo is %s, bar is %d', Good->instance->foo, Good->instance->bar) );
 
 # foo is 'hot', bar is 30

 ## SUBS ##
 sub make_foo_hot
 {
    # get the Good instance, don't send params. it will error out
    # because the instance is already created. use the accessors to
    # change the attributes
    my $good = Good->instance;
    
    return $good->foo('hot'); # foo is 'hot'
 }

 sub raise_the_bar
 {
    my $val = shift;

    # set bar to 30, using the instance
    my $bar = Good->instance->bar;
    $bar += $val;

    return Good->instance->bar($bar);
 }
   
=head1 Public Methods

=head2 instance

 my $obj = Good->instance;

=over

=item I<Synopsis>   Creates/gets the current instance of the class

=item I<Params>     depends on class

=item I<Returns>    current class instance

=back

=head1 Private Methods

=head2 _enforce_instance

 $obj->_enforce_instance;

=over

=item I<Synopsis>   Verifies that the object is indeed an object.

=item I<Params>     none

=item I<Returns>    true on success, fatal on failure

=back

=head1 See Also

=over

=item Finfo::Std

=item Finfo::Logging

=back

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
