package Finfo::Object;

use strict;
use warnings;
no warnings 'reserved';

use Data::Dumper;
use Finfo::ClassUtils 'class';
use Finfo::Logging;
use Finfo::Validate;
use Finfo::OldValidate;

sub new
{
    my $self = bless {}, shift;

    return unless $self->_add_attributes;

    return unless $self->_add_params(@_);

    return unless $self->_init;

    return $self;
}

sub _mk_accessor 
{
    my ($self, $attr) = @_;
    
    my $accessor =  sub 
    {
        my ($self, $value) = @_;

        # Check if private
        if ( $self->attribute_is_private($attr) )
        {
            my $caller = caller;
            $self->error_msg
            (
                sprintf('Cannot set a private attribute (%s) from outside %s', $attr, $self->class)
            ) 
                and return unless UNIVERSAL::isa($self, $caller);
        }

        if ( defined $value ) #set
        {
            my $access = $self->access_for_attribute($attr);
            $self->error_msg("Attribute ($attr) is read only, cannot set")
                and return if $access eq 'ro';

            return unless $self->_validate($attr, $value);

            my $old_value = $self->{$attr};

            $self->{$attr} = $value;

            if ( $self->can('_post_set') )
            {
                $self->_post_set($attr, $value, $old_value);
            }
            
            return $self->{$attr};
        }
        else # get
        {
            return $self->{$attr};
        }
    };

    my $class = $self->class;

    no strict 'refs';
    *{$class."\:\:$attr"}  = $accessor unless defined &{$class."\:\:$attr"};

    return 1;
}

sub undef_attribute
{
    my ($self, $attr) = @_;

    return unless Finfo::Validate->validate
    (
        attr => 'attr to set undef',
        value => $attr,
        msg => 'fatal_msg',
    );

    $self->error_msg("Invalid attr to set undef")
        and return unless $self->can($attr);

    $self->{$attr} = undef;

    return 1;
}

sub _init
{
    return 1;
}

# attrs and params
sub _add_params
{
    my ($self, %p) = @_;

    foreach my $attr ( $self->attributes )
    {
        $self->error_msg("Can't find attribute ($attr)")
            and return unless $self->can($attr);

        my $value = delete $p{$attr};

        # If no value, check if attr is required, and set by a default
        # Note - moved setting default vals to method _add_attributes
        unless ( defined $value )
        {
            if ( $self->attribute_is_required($attr) and not defined $self->$attr )
            {
                $self->error_msg("$attr is required")
                    and return unless defined $value;
            }

            next;
        }

        return unless $self->_validate($attr, $value);

        $self->{$attr} = $value;
    }

    # check if there are left over keys
    $self->error_msg
    (
        sprintf
        ( 
            'Invalid params passed in (%s)',
            join(', ', keys %p),
        )
    )
        and return if %p;

    return 1;
}

# attributes - also sets defaults
sub _add_attributes
{
    my $self = shift;

    my $attrs = $self->_attrs;

    return unless Finfo::Validate->validate
    (
        attr => 'attribute hashref',
        value => $attrs,
        ds => 'hashref',
        empty_ok => 1,
        msg => 'fatal',
    );

    foreach my $attr ( keys %$attrs )
    {
        my ($name) = $attr =~ /(.+):[rop]/;
        $self->error_msg("Invalid setup for attribute ($attr)")
            and return unless $name;

        if ( not defined $self->{$name} and exists $attrs->{$attr}->{default} )
        {
            $self->{$name} = $attrs->{$attr}->{default} 
        }

        $self->_mk_accessor($name);
    }
    
    return 1;
}

sub _attrs
{
    my $self = shift;
    
    return $self->{_attrs} if exists $self->{_attrs};

    my $attrs = {};
    foreach my $attr_method (qw/ _reqs _opts _private_attributes /)
    {
        next unless $self->can($attr_method);

        my $class = $self->class;
        no strict 'refs';
        unless ( ${"$class".'::OLD_OBJECT_WARNED'} )
        {
            $self->info_msg
            (
                "Finfo::Object attribute api has changed, and this class has not been updated"
            );
        }
        ${"$class".'::OLD_OBJECT_WARNED'} = 1;
        
        use strict;

        my ($attr_type) = $attr_method =~ /^_(\w)/;

        my $old_attrs = $self->$attr_method;
        if ( ref $old_attrs eq 'HASH' )
        {
            foreach my $attr ( keys %$old_attrs )
            {
                my $attr_name = "$attr\:$attr_type";
                $attrs->{$attr_name}->{type} = $old_attrs->{$attr}->[0] || 'defined';
                $attrs->{$attr_name}->{default} = $old_attrs->{$attr}->[1] || undef;
                $attrs->{$attr_name}->{options} = $old_attrs->{$attr}->[2] || undef;
            }
        }
        else
        {
            foreach my $attr ( $self->$attr_method )
            {
                my $attr_name = "$attr\:$attr_type";
                $attrs->{$attr_name}->{type} = 'defined';
            }
        }
    }

    $self->{_attrs} = $attrs;

    return $attrs;
}

sub attributes
{
    my $self = shift;

    return ( $self->required_attributes, $self->optional_attributes );
}

sub required_attributes
{
    my $self = shift;

    my $attrs = $self->_attrs;

    my @reqs = grep { s/:r$// } keys %$attrs;
    
    return @reqs;
}

sub optional_attributes
{
    my $self = shift;

    my $attrs = $self->_attrs;

    my @opts = grep {  s/:o$// } keys %$attrs;

    return @opts; 
}

sub attribute_is_required
{
    my ($self, $attr) = @_;

    $self->error_msg("No attribute to test if required")
        and return unless defined $attr;

    my $attrs = $self->_attrs;
    
    return $attrs->{"$attr\:r"};
}

sub attribute_is_private
{
    my ($self, $attr) = @_;

    $self->error_msg("No attribute to test if required")
        and return unless defined $attr;

    my $attrs = $self->_attrs;
    
    return $attrs->{"$attr\:p"};
}

# validate
sub _validate
{
    my ($self, $attr, $value) = @_;

    return Finfo::OldValidate->validate
    (
        attr => $attr,
        value => $value,
        type => $self->validation_type_for_attribute($attr),
        ref_type => $self->validation_ref_type_for_attribute($attr),
        options => $self->valid_options_for_attribute($attr),
        err_cb => $self,
    );
}

# defaults
sub default_for_attr { attributes_attribute(@_) }
sub _default_for_attr { attributes_attribute(@_) }

sub attributes_attribute
{
    my ($self, $attr, $attrs_attr) = @_;
    
    $self->error_msg("No attribute to get default for attribute")
        and return unless defined $attr;
    
    return unless Finfo::Validate->validate
    (
        attr => "attributes ($attr) attribute",
        value => $attrs_attr,
        isa => 'in_list attr_type type access options ref_type clo desc',
        msg => 'fatal'
    );

    my $attrs = $self->_attrs;

    return $attrs->{"$attr\:r"}->{$attrs_attr} if exists $attrs->{"$attr\:r"};
    return $attrs->{"$attr\:o"}->{$attrs_attr} if exists $attrs->{"$attr\:o"};
    return $attrs->{"$attr\:p"}->{$attrs_attr} if exists $attrs->{"$attr\:p"};
    return;
}

sub validation_ref_type_for_attribute
{
    my ($self, $attr) = @_;

    return $self->_default_for_attr($attr, 'ref_type');
}

sub validation_type_for_attribute
{
    my ($self, $attr) = @_;

    return $self->_default_for_attr($attr, 'type') || 'defiend';
}

sub default_value_for_attribute
{
    my ($self, $attr) = @_;

    return $self->_default_for_attr($attr, 'default');
}

sub valid_options_for_attribute
{
    my ($self, $attr) = @_;

    return $self->_default_for_attr($attr, 'options');
}

sub access_for_attribute
{
    my ($self, $attr) = @_;

    return $self->_default_for_attr($attr, 'access') || 'rw';
}

sub description_for_attribute
{
    my ($self, $attr) = @_;

    return $self->_default_for_attr($attr, 'desc');
}

sub command_line_option_for_attribute
{
    my ($self, $attr) = @_;

    return $self->default_for_attr($attr, 'clo');
}

1;

=pod

=head1 Name

Finfo::Object

** abstract base class ** 

=head1 Synopsis

This hash based abstract base class provides a contsructor which will auto build accessors, validate and set the incoming params and execute the _init method.  Uses Finfo::Logging to handle messaging and logging.  If the object fails the construction, an error message will be set on the class.
  
=head1 Methods to Overwrite in Your Class

=head2 _attrs

 sub _attrs
 {
    return
    {
        'attr1:r' => # a required attribute
        {
            type => 'inlist', # validation type, optional - default is defined
            default => 'foo', # default value to set to attribute if none given, optional - no default
            options => [qw/ foo bar /], # aryref of valid options for the validation type, optional - no default
            access => 'ro', # access level, optional - default is rw (read-write)
        },
        'attr2:o' => # an optional attribute
        {
            type => 'inlist', # validation type, optional - default is defined
            default => 'foo', # default value to set to attribute if none given, optional - no default
            options => [qw/ foo bar /], # aryref of valid options for the validation type, optional - no default
            access => 'ro', # access level, optional - default is rw (read-write)
        },
        '_attr3:p' => # private attribute, access only from the class and it's children
        {
            type => 'int', # validation type, optional - default is defined
            default => 'foo', # default value to set to attribute if none given, optional - no default
            options => [qw/ foo bar /],        # aryref of valid options for the validation type, optional - no default
        },
    };
}

=over

=item I<Attribute's Attributes>

=over

=item B<hash key>       REQUIRED! This is the name of the attribute folowed by a ':' and r, o or p (indicating the type)

=over

=item I<r>       Required attribute - must be supplied in the constructor

=item I<o>       Optional attribute - not necessary in constructor, may have a default

=item I<p>       Private attribute - only accessible to the class and it's children

=back

=item B<type>       The validation type to use (optional, default is 'defined')

=item B<options>    The validation options to use (required for some  validation types)

=item B<ref_type>   The validation ref type (optional, no default see Finfo::Validate)

=item B<default>    The default value to set to the attribute (optional, no default)

=item B<access>     The access of the atribute (optional default is 'rw')

=item B<clo>        The command line option of the atribute (optional, no default)

=item B<desc>       The desciption for the attribute of the atribute (optional, no default)

=back

=back

=head2 _init

 sub _init
 {
    my $self = shift;

    my $fh =  IO::File->new('<' . $self->file);
    $self->error_msg( sprintf('Can\'t open file (%s)', $self->file) )
        and return unless $fh;
    
    $self->file_handle($fh)
    
    return 1; # return true for success, false for failure
 }

=item I<Synopsis>   This is method is called after the attributes have been added and the params have been set.  This allows the object to have additional intialization, like opening files, etc.
This method does not need to be overwritten.

=head1 Caveats

=head2 Calling SUPER in _attrs and _init

If inheritting from a class that already uses the _attrs and/orthe _init methods, it will be necessary to call SUPER::_attrs or SUPER::_init resppectfully to ensure that the object has been setup properly.

=head1 See Also

I<Finfo::Validate> I<Finfo::Logging> I<Finfo::CommandLineOptions>

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> <ebelter@watson.wustl.edu>
B<Adam Dukes> <adukes@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
