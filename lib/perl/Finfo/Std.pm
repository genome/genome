package Finfo::Std;
use version; $VERSION = qv('0.0.1');

use strict;
use warnings;
no warnings 'reserved';

use Carp;
use Data::Dumper;
use Finfo::ClassUtils 'anon_scalar class ident';
use Finfo::Logging 'info_msg fatal_msg';
use Finfo::Validate;
use Scalar::Util 'blessed';
use Storable 'dclone';
$Storable::Deparse = 1;
$Storable::Eval = 1;
use overload;

my (%attribute, %attributes_of, %cumulative, %anticumulative, %restricted, %private, %overload);

my @exported_subs = qw(
    new
    undef_attribute
    attributes
    required_attributes
    optional_attributes
    attributes_attribute 
    DESTROY
    AUTOLOAD
    DUMP
    STORABLE_freeze
    STORABLE_thaw
);

my @exported_extension_subs = qw(
    MODIFY_HASH_ATTRIBUTES
    MODIFY_CODE_ATTRIBUTES
);

sub import # eb modified
{
    my $caller = caller;

    if ( $caller eq 'main' )
    {
        __PACKAGE__->fatal_msg("Can't use in main");
        return;
    }

    # import logging methods
    Finfo::Logging->import({ class => $caller });

    no strict 'refs';
    for my $sub ( @exported_subs ) 
    {
        *{ $caller . '::' . $sub } = \&{$sub};
    }

    for my $sub ( @exported_extension_subs ) 
    {
        my $target = $caller . '::' . $sub;
        my $real_sub = *{ $target }{CODE} || sub { return @_[2..$#_] };
        no warnings 'redefine';
        *{ $target } = sub 
        {
            my ($package, $referent, @unhandled) = @_;
            for my $handler ($sub, $real_sub) {
                next if !@unhandled;
                @unhandled = $handler->($package, $referent, @unhandled);
            }
            return @unhandled;
        };
    }
}

sub _find_sub
{
    my ($package, $sub_ref) = @_;
    
    my $caller = caller;
    unless ( $caller eq __PACKAGE__ )
    {
        __PACKAGE__->fatal_msg
        (
            "Method (_find_sub) is private, for use only from " . __PACKAGE__
        );
        return;
    }

    no strict 'refs';
    for my $name (keys %{$package.'::'}) 
    {
        my $candidate = *{$package.'::'.$name}{CODE};
        return $name if $candidate && $candidate == $sub_ref;
    }
    
    $package->fatal_msg(q{Can't make anonymous subroutine cumulative});
}

sub _get_value_from_attr_config
{
    my ($class, $attr_config) = @_;

    my $caller = caller;
    unless ( $caller eq __PACKAGE__ )
    {
        __PACKAGE__->fatal_msg
        (
            "Method (_get_value_from_attr_config) is private, for use only from " . __PACKAGE__
        );
        return;
    }

    $attr_config =~ s/__PACKAGE__/$class/g;

    my $value = eval $attr_config;

    if ( $@ )
    {
        $value = $attr_config;
    }

    #print Data::Dumper::Dumper({in => \@_, val => $value});
    return $value;
}

sub undef_attribute
{
    my ($self, $attr) = @_;

    unless ( $attr )
    {
        $self->fatal_msg("Need attribute to set undef");
        return;
    }
    
    my $id = ident($self);

    unless ( $id )
    {
        __PACKAGE__->fatal_msg("Need object for " . class($self));
        return;
    }

    unless ( UNIVERSAL::can($self, $attr) )
    {
        $self->fatal_msg("Can't find attribute ($attr)");
        return;
    }
    
    my $class = class($self);
    for my $base_class ( _reverse_hierarchy_of($class) ) 
    {
        for my $attr_ref ( @{ $attribute{$base_class} } ) 
        {
            next unless $attr_ref->{name} eq $attr;
            $attr_ref->{ref}{$id} = undef;
            return 1;
        }
    }

    # Shouldn't get here
    __PACKAGE__->fatal_msg("$class can do attribute ($attr), but can't find storage var");
    return;
}

sub attributes
{
    my $self = shift;

    return ($self->required_attributes, $self->optional_attributes);
}

sub required_attributes
{
    my $self = shift;

    my $class = class($self);
    
    my @attrs;
    foreach my $class ( _hierarchy_of($class) )
    {
        push @attrs, map { $_->{name} } grep { $_->{attr_type} eq 'r' } @{ $attribute{ $class } };
    }

    return @attrs;
}

sub optional_attributes
{
    my $self = shift;

    my $class = class($self);

    my @attrs;
    foreach my $class ( _hierarchy_of($class) )
    {
        push @attrs, map { $_->{name} } grep { $_->{attr_type} eq 'o' } @{ $attribute{ $class } };
    }

    return @attrs;
}

sub attributes_attribute
{
    my ($self, $attr, $attrs_attr) = @_;

    Finfo::Validate->validate
    (
        attr => 'attribute',
        value => $attr,
        msg => 'fatal',
    );

    Finfo::Validate->validate
    (
        attr => "attribute's ($attr) attribute",
        value => $attrs_attr,
        isa => 'in_list name attr_type access default isa ds desc clo',
        msg => 'fatal',
    );

    my $attr_ref;
    foreach my $class ( _hierarchy_of( class($self) ) )
    {
        ($attr_ref) = grep { $attr eq $_->{name} } @{ $attribute{ $class } };
        last if $attr_ref;
    }

    unless ( $attr_ref )
    {
        $self->fatal_msg("Can't find attribute ($attr) for " . class($self) . " to get attributes");
        return;
    }

    if ( ref($attr_ref->{$attrs_attr}) eq 'CODE' )
    {
        return $attr_ref->{$attrs_attr}->();
    }

    return $attr_ref->{$attrs_attr};
}

sub MODIFY_HASH_ATTRIBUTES
{
    my ($package, $referent, @attrs) = @_;

    my $caller = [ caller ];

    my %pkg_attrs_attrs = 
    (
        isa => 'defined',
        ds => 'scalar',
        empty_ok => 0,
        default => undef,
        access => 'rw', 
        clo => '', 
        desc => '',
    );

    my $name_config = shift @attrs;
    @pkg_attrs_attrs{qw/ name attr_type /} = $name_config =~ /^name\(([\w\d_]*):([rop])\)$/i;

    # Validate name
    Finfo::Validate->validate
    (
        attr => "attribute name in attribute config ($name_config)",# in $package",
        value => $pkg_attrs_attrs{name},
        obj => $package,
        msg => 'fatal',
        caller_level => 3,
    );

    # Validate attr type
    Finfo::Validate->validate
    (
        attr => "attribute type for attribute ($pkg_attrs_attrs{name})",# in $package",
        value => $pkg_attrs_attrs{attr_type},
        isa => 'in_list r o p',
        obj => $package,
        msg => 'fatal',
        caller_level => 3,
    );

    foreach my $attr_config ( @attrs )
    {
        my ($attr_name, $attr_value) = $attr_config =~ /^([\w\d_]+)\((.+)\)$/i;
        #print Data::Dumper::Dumper([$attr_config, $attr_name, $attr_value]);

        unless ( $attr_name )
        {
            __PACKAGE__->fatal_msg("Error in attribute's ($pkg_attrs_attrs{name}) config ($attr_config) n$package");
        }

        unless ( grep { $attr_name eq $_ } keys %pkg_attrs_attrs )
        {
            __PACKAGE__->fatal_msg
            (
                "Unknown attribute's ($pkg_attrs_attrs{name}) attribute ($attr_name) in $package"
            );
        }
        
        Finfo::Validate->validate
        (
            attr => "attribute's ($pkg_attrs_attrs{name}) attribute value ($attr_value)",# in $package",
            value => $attr_value,
            obj => $package,
            msg => 'fatal',
            caller_level => 3,
        );

        $attr_value =~ s/__PACKAGE__/$package/g;

        my $value = eval $attr_value;

        #print Data::Dumper::Dumper([$attr_name,$attr_value,$value]);
        
        if ( $@ )
        {
            $value = $attr_value;
        }

        $pkg_attrs_attrs{$attr_name} = $value;
    }

    # Validate type
    Finfo::Validate->is_isa
    (
        attr => "$pkg_attrs_attrs{name} in $package",
        value => $pkg_attrs_attrs{isa},
        obj => $package,
        msg => 'fatal',
        caller_level => 3,
    );

    no strict 'refs';
    # Accessor
    *{ $package . '::'. $pkg_attrs_attrs{name} } = sub 
    {
        my ($self, $value) = @_;

        # Check if private attribute
        if ( $pkg_attrs_attrs{attr_type} eq 'p' )
        {
            my $caller = caller;
            unless ( $self->isa($caller) )
            {
                $self->fatal_msg("Can't access private attribute ($pkg_attrs_attrs{name}) from outside $package");
                return;
            }
        }

        if ( defined $value )
        {
            # Check access
            if ( $pkg_attrs_attrs{access} eq 'ro' )
            {
                $self->fatal_msg("Can't set read only attribute ($pkg_attrs_attrs{name})");
                return;
            }

            # Validate
            Finfo::Validate->validate
            ( 
                attr => $pkg_attrs_attrs{name},
                value => $value,
                ds => $pkg_attrs_attrs{ds}, 
                isa => $pkg_attrs_attrs{isa}, 
                empty_ok => $pkg_attrs_attrs{empty_ok}, 
                obj => $package,
                msg => 'fatal',
                caller_level => 1,
            );

            $referent->{ ident($self) } = $value;
        }
        else
        {
            return $referent->{ ident($self) };
        }
    };

    use strict 'refs';

    push @{ $attributes_of{$package} },
    {
        ref      => $referent,
        name     => $pkg_attrs_attrs{name},
    };

    push @{ $attribute{$package} }, 
    {
        ref       => $referent,
        %pkg_attrs_attrs,
    };

    return;
}

sub _DUMP
{
    my (%params) = @_;

    my $caller = caller;
    unless ( $caller eq __PACKAGE__ )
    {
        __PACKAGE__->fatal_msg
        (
            "Method (_DUMP) is private, for use only from " . __PACKAGE__
        );
        return;
    }

    my $object = delete $params{object};
    __PACKAGE__->fatal_msg("Need object to DUMP")
        and return unless $object;
    my $separate = delete $params{separate} || 0;
    my $recursive = delete $params{recursive} || 0;

    $object->fatal_msg("Unknown params passed into DUMP: " . join(', ', keys %params))
        and return if %params;

    my $id = ident($object);
    unless ( $id )
    {
        $object->fatal_msg("Can't call DUMP on a class");
        return;
    }

    my $class = class($object);
    my @package_list = $class;
    my %package_seen = ( $class => 1 ); #ignore diamond/looped base classes :-)

    my %obj_attrs;
    PACKAGE:
    while( my $package = shift @package_list) 
    {
        #make sure we add any base classes to the list of
        #packages to examine for attributes.
        { no strict 'refs';
            for my $base_class ( @{"${package}::ISA"} ) 
            {
                push @package_list, $base_class unless $package_seen{$base_class}++;
            }
        }
        #examine attributes from known packages only
        my $attr_list_ref = $attributes_of{$package} or next PACKAGE;

        #look for any attributes of this object for this package
        ATTR:
        for my $attr_ref ( @{$attr_list_ref} ) 
        {
            my $value = $attr_ref->{ref}{$id};

            #nothing to do if attr not set for this object
            next ATTR unless defined $value;

            if ( $recursive and my $value_class = blessed($attr_ref->{ref}{$id}) )
            {
                # ok, it's an object, test whether it is stored w/in Finfo::Std
                if ( exists $attribute{$value_class} )
                {
                    $value = _DUMP
                    (
                        object => $value,
                        separate => $separate, 
                        recursive => $recursive
                    );
                }
            }

            #save the attr by name into the package or class hash
            ( $separate )
            ? $obj_attrs{$package}{ $attr_ref->{name} } = $value
            : $obj_attrs{$class}{ $attr_ref->{name} } = $value;
        }
    }

    return \%obj_attrs;
}

sub DUMP 
{
    my $self = shift;

    my $obj = _DUMP(object => $self, @_);
    
    my $dump = Data::Dumper::Dumper($obj);
    $dump =~ s/^.{8}//gxms;

    return $dump;
}

my $STD_OVERLOADER
    = q{ package %%s;
         use overload (
            q{%s} => sub { $_[0]->%%s($_[0]->ident()) },
            fallback => 1
         );
       };

my %OVERLOADER_FOR = (
    STRINGIFY => sprintf( $STD_OVERLOADER, q{""}   ),
    NUMERIFY  => sprintf( $STD_OVERLOADER, q{0+}   ),
    BOOLIFY   => sprintf( $STD_OVERLOADER, q{bool} ),
    SCALARIFY => sprintf( $STD_OVERLOADER, q{${}}  ),
    ARRAYIFY  => sprintf( $STD_OVERLOADER, q{@{}}  ),
    HASHIFY   => sprintf( $STD_OVERLOADER, q{%%{}} ),  # %% to survive sprintf
    GLOBIFY   => sprintf( $STD_OVERLOADER, q{*{}}  ),
    CODIFY    => sprintf( $STD_OVERLOADER, q{&{}}  ),
);

sub MODIFY_CODE_ATTRIBUTES 
{
    my ($package, $referent, @attrs) = @_;

    for my $attr (@attrs) 
    {
        if ($attr eq 'CUMULATIVE') 
        {
            push @{$cumulative{$package}}, $referent;
        }
        elsif ($attr =~ m/\a cumulative \s* [(] \s* base \s* first \s* [)] \z/xms) 
        {
            push @{$anticumulative{$package}}, $referent;
        }
        elsif ($attr =~ m/\A RESTRICTED \z/xms) 
        {
            push @{$restricted{$package}}, $referent;
        }
        elsif ($attr =~ m/\A PRIVATE \z/xms) 
        {
            push @{$private{$package}}, $referent;
        }
        elsif (exists $OVERLOADER_FOR{$attr}) 
        {
            push @{$overload{$package}}, [$referent, $attr];
        }
        undef $attr;
    }
    
    return grep {defined} @attrs;
}

my %_hierarchy_of;

sub _hierarchy_of
{
    my ($class) = @_;
    
    my $caller = caller;
    unless ( $caller eq __PACKAGE__ )
    {
        __PACKAGE__->fatal_msg
        (
            "Method (_hierarchy_of) is private, for use only from " . __PACKAGE__
        );
        return;
    }

    return @{$_hierarchy_of{$class}} if exists $_hierarchy_of{$class};

    no strict 'refs';

    my @hierarchy = $class;
    my @parents   = @{$class.'::ISA'};

    while (defined (my $parent = shift @parents)) {
        push @hierarchy, $parent;
        push @parents, @{$parent.'::ISA'};
    }

    for (my $i = @hierarchy - 1; $i >= 0; $i--)
    {
        my $pos = $i;
        for (my $j = 1; $j <= $i; $j++)
        {
            if ( $hierarchy[$j]->isa($hierarchy[$j - 1]) )
            {
                my $temp = $hierarchy[$j - 1];
                $hierarchy[$j - 1] = $hierarchy[$j];
                $hierarchy[$j] = $temp;
            }
        }
    }
                           
    #print Dumper(['hierarchy', $class,\@hierarchy]);

    return @{$_hierarchy_of{$class}} = @hierarchy;
    
    # old sort
    #my %seen;
    #return @{$_hierarchy_of{$class}}
    #= sort { $a->isa($b) ? -1
    #    : $b->isa($a) ? +1
    #    :                0
        #} grep !$seen{$_}++, @hierarchy;
}

my %_reverse_hierarchy_of;
sub _reverse_hierarchy_of 
{
    my ($class) = @_;

    my $caller = caller;
    unless ( $caller eq __PACKAGE__ )
    {
        __PACKAGE__->fatal_msg
        (
            "Method (_reverse_hierarchy_of) is private, for use only from " . __PACKAGE__
        );
        return;
    }

    return @{$_reverse_hierarchy_of{$class}}
    if exists $_reverse_hierarchy_of{$class};

    no strict 'refs';

    my @hierarchy = $class;
    my @parents   = reverse @{$class.'::ISA'};

    while (defined (my $parent = shift @parents)) {
        push @hierarchy, $parent;
        push @parents, reverse @{$parent.'::ISA'};
    }

    for (my $i = @hierarchy - 1; $i >= 0; $i--)
    {
        my $pos = $i;
        for (my $j = 1; $j <= $i; $j++)
        {
            if ( $hierarchy[$j - 1]->isa($hierarchy[$j]) )
            {
                my $temp = $hierarchy[$j - 1];
                $hierarchy[$j - 1] = $hierarchy[$j];
                $hierarchy[$j] = $temp;
            }
        }
    }

    #print Dumper(['reverse_hierarchy', $class, \@hierarchy]);:
    
    return @{$_reverse_hierarchy_of{$class}} = @hierarchy;

    # old sort
    #my %seen;
    #return @{$_reverse_hierarchy_of{$class}} = reverse sort { $a->isa($b) ? -1
    #    : $b->isa($a) ? +1
    #    :                0
    #} grep !$seen{$_}++, @hierarchy;
}

{
    no warnings qw( void );
    CHECK { initialize() }
}

sub initialize # eb modified
{
    # Short-circuit if nothing to do...
    return if keys(%restricted) + keys(%private)
            + keys(%cumulative) + keys(%anticumulative)
            + keys(%overload)
                == 0;

    my (%cumulative_named, %anticumulative_named);

    # Implement restricted methods (only callable within hierarchy)...
    for my $package (keys %restricted) 
    {
        for my $sub_ref (@{$restricted{$package}}) {
            my $name = _find_sub($package, $sub_ref);
            no warnings 'redefine';
            no strict 'refs';
            my $sub_name = $package.'::'.$name;
            my $original = *{$sub_name}{CODE};
            unless ( $original )
            {
                $package->fatal_msg("Restricted method ${package}::$name() declared but not defined");
                return;
            }
            *{$sub_name} = sub 
            {
                my $caller;
                my $level = 0;
                while ($caller = caller($level++)) 
                {
                    last if $caller !~ /^(?: Finfo::Std | attributes )$/xms;
                }
                goto &{$original} if !$caller || $caller->isa($package) || $package->isa($caller);
                                              
                $package->fatal_msg("Can't call restricted method $sub_name() from class $caller");
                return;
            }
        }
    }

    # Implement private methods (only callable from class itself)...
    for my $package (keys %private)
    {
        for my $sub_ref (@{$private{$package}}) {
            my $name = _find_sub($package, $sub_ref);
            no warnings 'redefine';
            no strict 'refs';
            my $sub_name = $package.'::'.$name;
            my $original = *{$sub_name}{CODE};
            unless ( $original )
            {
                $package->fatal_msg("Private method ${package}::$name() declared but not defined");
                return;
            }
            
            *{$sub_name} = sub {
                my $caller = caller;
                goto &{$original} if $caller eq $package;
                $package->fatal_msg
                (
                    "Can't call private method $sub_name() from class $caller",
                    { 'caller' => [ caller ] },
                );
                return;
            }
        }
    }

    for my $package (keys %cumulative) {
        for my $sub_ref (@{$cumulative{$package}}) {
            my $name = _find_sub($package, $sub_ref);
            $cumulative_named{$name}{$package} = $sub_ref;
            no warnings 'redefine';
            no strict 'refs';
            *{$package.'::'.$name} = sub {
                my @args = @_;
                my $class = ref($_[0]) || $_[0];
                my $list_context = wantarray; 
                my (@results, @classes);
                for my $parent (_hierarchy_of($class)) {
                    my $sub_ref = $cumulative_named{$name}{$parent} or next;
                    ${$parent.'::AUTOLOAD'} = our $AUTOLOAD if $name eq 'AUTOLOAD';
                    if (!defined $list_context) {
                        $sub_ref->(@args);
                        next;
                    }
                    push @classes, $parent;
                    if ($list_context) {
                        push @results, $sub_ref->(@args);
                    }
                    else {
                        push @results, scalar $sub_ref->(@args);
                    }
                }
                return if !defined $list_context;
                return @results if $list_context;
                return Finfo::Std::SCR->new({
                    values  => \@results,
                    classes => \@classes,
                });
            };
        }
    }

    for my $package (keys %anticumulative) {
        for my $sub_ref (@{$anticumulative{$package}}) {
            my $name = _find_sub($package, $sub_ref);
            if ($cumulative_named{$name}) {
                for my $other_package (keys %{$cumulative_named{$name}}) 
                {
                    next unless $other_package->isa($package) || $package->isa($other_package);
                    __PACKAGE__->fatal_msg
                    (
                        "Conflicting definitions for cumulative method ($name) ",
                        "(specified as :CUMULATIVE in $other_package ",
                        " but declared :CUMULATIVE(BASE FIRST) in $package)\n"
                    );
                }
            }
            $anticumulative_named{$name}{$package} = $sub_ref;
            no warnings 'redefine';
            no strict 'refs';
            *{$package.'::'.$name} = sub {
                my $class = ref($_[0]) || $_[0];
                my $list_context = wantarray; 
                my (@results, @classes);
                for my $parent (_reverse_hierarchy_of($class)) {
                    my $sub_ref = $anticumulative_named{$name}{$parent} or next;
                    if (!defined $list_context) {
                        &{$sub_ref};
                        next;
                    }
                    push @classes, $parent;
                    if ($list_context) {
                        push @results, &{$sub_ref};
                    }
                    else {
                        push @results, scalar &{$sub_ref};
                    }
                }
                return if !defined $list_context;
                return @results if $list_context;
                return Finfo::Std::SCR->new({
                    values  => \@results,
                    classes => \@classes,
                });
            };
        }
    }

    for my $package (keys %overload) {
        foreach my $operation (@{ $overload{$package} }) {
            my ($referent, $attr) = @$operation;
            local $^W;
            my $method = _find_sub($package, $referent);
            eval sprintf $OVERLOADER_FOR{$attr}, $package, $method;
            die "Internal error: $@" if $@;
        }
    }

    # Remove initialization data to prevent re-initializations...
    %restricted     = ();
    %private        = ();
    %cumulative     = ();
    %anticumulative = ();
    %overload       = ();
}

sub new # eb modified
{
    my ($class, %params) = @_;

    #print Dumper([$class, \%params]);

    Finfo::Std::initialize();   # Ensure run-time (and mod_perl) setup is done

    no strict 'refs';
    unless ( keys %{$class.'::'} )
    {
        __PACKAGE__->fatal_msg("Can't find class $class");
        return;
    }

    my $new_obj = bless anon_scalar(), $class;
    #my $new_obj = bless \my($anon_scalar), $class;
    my $new_obj_id = ident($new_obj);

    for my $base_class (_reverse_hierarchy_of($class)) 
    {
        # BUILD - called before params are added
        {
            no warnings 'once';
            if ( my $build_ref = *{$base_class.'::BUILD'}{CODE} ) 
            {
                $build_ref->($new_obj, $new_obj_id, \%params)
                    or __PACKAGE__->fatal_msg("Error in BUILDing $class");
            }
        }

        # Apply params and defaults to attributes
        ATTR: for my $attr_ref ( @{ $attribute{$base_class} } ) 
        {
            my $value = delete $params{ $attr_ref->{name} };

            # Check if initializing a private attribute
            if ( defined $value and $attr_ref->{attr_type} eq 'p' )
            {
                $new_obj->fatal_msg
                (
                    sprintf
                    (
                        "Can't initialize a private attribute (%s)",
                        $attr_ref->{name},
                    ),
                    { 'caller' => [ caller ] },
                );
            }

            # Get default
            if ( not defined $value and defined $attr_ref->{default} )
            {
                my $default = $attr_ref->{default};
                my $default_ref = ref($default);
                $value = ( $default_ref and grep { $default_ref eq $_ } (qw/ ARRAY HASH /) )
                ? dclone($default) 
                : $default;
            }

            # Check if atribute is required
            if ( not defined $value and $attr_ref->{attr_type} eq 'r' )
            {
                $new_obj->fatal_msg
                (
                    $attr_ref->{name} . " is required",
                    { caller_level => 1 },
                );
            }

            next ATTR unless defined $value;

            # Validate value
            return unless Finfo::Validate->validate
            ( 
                attr => $attr_ref->{name},
                value => $value,
                ds => $attr_ref->{ds}, 
                isa => $attr_ref->{isa}, 
                empty_ok => $attr_ref->{empty_ok}, 
                obj => $new_obj,
                msg => 'fatal',
                'caller' => [ caller ],
            );

            # Set value
            $attr_ref->{ref}{$new_obj_id} = $value;
        }
    }

    # Check if any unrecognized params passed in
    if ( %params )
    {
        $class->fatal_msg
        (
            sprintf
            ( 
                'Invalid params (%s) passed in',
                join(', ', keys %params),
            ),
            { caller_level => 1 },
        );
        return;
    }

    START:
    for my $base_class ( _reverse_hierarchy_of($class) )
    {
        # Apply START() methods...
        {
            no warnings 'once';
            if (my $init_ref = *{$base_class.'::START'}{CODE}) 
            {
                unless ( $init_ref->($new_obj, $new_obj_id) )
                { 
                    $base_class->fatal_msg
                    (
                        "Error initializing $base_class in START method",
                        { 'caller' => [ caller ] },
                    );
                    return;
                }
            }
        }
    }

    return $new_obj;
}

sub uniq (@) {
    my %seen;
    return grep { $seen{$_}++ } @_;
}

sub DESTROY 
{
    my ($self) = @_;

    my $id = ident($self);
    
    push @_, $id;

    DEMOLISH: 
    for my $base_class (_hierarchy_of(ref $_[0])) {
        no strict 'refs';
        if (my $demolish_ref = *{$base_class.'::DEMOLISH'}{CODE}) {
            &{$demolish_ref};
        }

        for my $attr_ref ( @{$attribute{$base_class}} ) {
            delete $attr_ref->{ref}{$id};
        }
    }

    return 1;
}

sub AUTOLOAD 
{
    my ($invocant) = @_;
    
    my $invocant_class = class($invocant);
    my ($package_name, $method_name) = our $AUTOLOAD =~ m/ (.*) :: (.*) /xms;

    my $ident = ident($invocant);
    if (!defined $ident) { $ident = $invocant }

    AUTOMETHOD:
    for my $parent_class ( _hierarchy_of($invocant_class) )
    {
        no strict 'refs';
        if (my $automethod_ref = *{$parent_class.'::AUTOMETHOD'}{CODE}) {
            local $CALLER::_ = $_;
            local $_ = $method_name;
            if ( my $method_impl = $automethod_ref->($invocant, $ident, @_[1..$#_]) ) 
            {
                goto &$method_impl;
            }
        }
    }

    my $type = ref $invocant ? 'object' : 'class';

    $package_name->fatal_msg
    (
        "Can't locate $type method ($method_name) via package ($package_name)",
        { 'caller' => [ caller ] }
    );

    return;
}

# eb added storable functions
sub STORABLE_freeze 
{
    my $caller = caller;
    #__PACKAGE__->fatal_msg("Freeze must be called from Storable", { caller_level => 1 })
    #    and return unless $caller eq 'Storable';
    
    my($self, $cloning) = @_;

    $self->STORABLE_freeze_pre($cloning)
        if $self->can("STORABLE_freeze_pre");
        
    my $id = ident($self);
    
    require Storable;
    
    my $serialized = Storable::freeze( \ (my $anon_scalar) );
    
    my $frozen_attr = _DUMP(object => $self, separate => 1);

    $self->STORABLE_freeze_post($cloning, $frozen_attr) if $self->can("STORABLE_freeze_post");

    #print Dumper($frozen_attr);
    
    return ($serialized, $frozen_attr );
}

sub STORABLE_thaw 
{
    my $caller = caller;
    #__PACKAGE__->fatal_msg("Thaw must be called from Storable", { caller_level = 1})
    #    and return unless $caller eq 'Storable';

    my($self, $cloning, $serialized, $frozen_attr_ref) = @_;

    #we can ignore $serialized, as we know it's an anon_scalar.

    $self->STORABLE_thaw_pre($cloning, $frozen_attr_ref)
        if $self->can("STORABLE_thaw_pre");

    my $id = ident($self);

    PACKAGE:
    while( my ($package, $pkg_attr_ref) = each %$frozen_attr_ref ) 
    {
        $self->fatal_msg("Unknown base class '$package' seen while thawing ")
        unless UNIVERSAL::isa($self, $package);

        my $attr_list_ref = $attributes_of{$package};
        ATTR:
        for my $attr_ref ( @{$attr_list_ref} ) #for known attrs...
        { 
            #nothing to do if frozen attr doesn't exist
            next ATTR unless exists $pkg_attr_ref->{ $attr_ref->{name} };

            #block attempts to meddle with existing objects
            $self->fatal_msg("Trying to modify existing attributes")
            if exists $attr_ref->{ref}{$id};

            #ok, set the attribute
            $attr_ref->{ref}{$id} = delete $pkg_attr_ref->{ $attr_ref->{name} };
        }

        if ( my @extra_keys = keys %$pkg_attr_ref ) 
        {
            #this is probably serious enough to throw an exception.
            #however, TODO: it would be nice if the class could somehow
            #indicate to ignore this problem.
            $self->fatal_msg
            (
                sprintf
                (
                    "Unknown attributes (%s) seen while thawing",
                    join(', ', @extra_keys),
                )
            );
        }
    }

    $self->STORABLE_thaw_post($cloning) if $self->can("STORABLE_thaw_post");
}


###################################################

package Finfo::Std::SCR;

use base qw( Finfo::Std );

use Finfo::ClassUtils 'ident';

my %values_of  :name(values:o);
my %classes_of :name(classes:o);

sub new {
    my ($class, $opt_ref) = @_;
    my $new_obj = bless \do{my $scalar}, $class;
    my $new_obj_id = ident($new_obj);
    $values_of{$new_obj_id}  = $opt_ref->{values};
    $classes_of{$new_obj_id} = $opt_ref->{classes};
    return $new_obj;
}

use overload (
    q{""}  => sub { return join q{}, grep { defined $_ } @{$values_of{ident($_[0])}}; },
    q{0+}  => sub { return scalar @{$values_of{ident($_[0])}};    },
    q{@{}} => sub { return $values_of{ident($_[0])};              },
q{%{}} => sub {
my ($self) = @_;
my %hash;
@hash{@{$classes_of{ident($self)}}} = @{$values_of{ident($self)}};
return \%hash;
    },
    fallback => 1,
);

1;

=pod 

=head1 Name

Finfo::Std

=head1 Synopsis

 Support for creating standard "inside-out" classes.  Please see Damian Conway's Class::Std for more information about inside out classes.
 
=head1 Usage

 package MyClass;
    
 use Finfo::Std;

 use IO::File;

 # Create storage for object attributes...
 my %file 
    :name(file:r) # this is the accessor name separated by a colon for the attribute type (r,o,p)
    :isa(file_w) # what the attribute is (validation type) - optional, default is defined
    :clo('in=s') # command line option - optional, no default
    :desc('Input file to parse'); # description for the clo - optional, no default
 my %fh
    :name(_fh:p)
    :isa('object IO::File');

 # The START method is optional.  It is called for the reverse
 # hierarchy of the class after the incoming parameters are set and
 # validated.  Use theis method to initialize the object
 sub START 
 {
     my ($self, $obj_ID, $arg_ref) = @_;

     # open a file
     my $fh = IO::File->new('<'.$self->file);
     $self->fatal_msg( sprintf('Can\'t open file (%s): %s', $self->file, $!) ) unless $fh;
     $self->_fh($fh);
     
     return 1; # always return true for success!
 }

 
 # The DEMOLISH method is optional.  It is called for the hierarchy of
 # the class when an object is DESTROYED.  Use it to close file handles, DB
 # connections etc.  
 sub DEMOLISH 
 {
     my ($self, $obj_id) = @_;

     # close file handle
     $self->_fh->close;
     
     return 1;
 }

 # The AUTOMETHOD is also optional.  This will be called for the
 # hierarchy of the class when a method can't be found.
 sub AUTOMETHOD 
 {
     my ($self, $obj_ID, @other_args) = @_;

     # Return any public data...
     if ( m/\A get_(.*)/ ) {  # Method name passed in $_
         my $get_what = $1;
         return sub {
             return $public_data{$obj_ID}{$get_what};
         }
     }

     warn "Can't call $method_name on ", ref $self, " object";

     return;   # The call is declined by not returning a sub ref
 }

=head1 Methods created automatically

=head2 new

Every class that loads the Finfo::Std module automatically has a C<new()>
constructor, which returns an inside-out object (i.e. a blessed scalar).

    $obj = MyClass->new();

The constructor can be passed a single argument to initialize the
object. This argument must be a hash reference. 

 $obj = MyClass->new
 (
    name => 'Foo', 
    location => 'bar',
 );

It is almost always an error to implement your own C<new()> in any class
that uses Finfo::Std. You almost certainly want to write a C<BUILD()> or
C<START()> method instead. See below.

=head2 undef_attribute

 $obj->undef_attribute('attribute');

Each time an attribute is set, it expects a value.  To explicitly undef an attribute, use this method.

=head2 required_attributes

 my @required_attributes = $obj->required_attributes;

Returns a list of the attributes that are required on creation of the object.

=head2 optional_attributes

 my @optional_attributes = $obj->optional_attributes;

Returns a list of the attributes that are optional on creation of the object.

=head2 attributes_attribute 
 
 my $attrs_attr = $obj->attributes_attribute('attribute', 'isa');

Gets the attribute of an object's attribute.

=head2 DESTROY

Every class that loads the Finfo::Std module automatically has a C<DESTROY>
destructor, which automatically cleans up any attributes declared.

It is almost always an error to write your own C<DESTROY> in any class that
uses Finfo::Std. You almost certainly want to write your own C<DEMOLISH>
instead. See below.

=head2 AUTOLOAD

Every class that loads the Finfo::Std module automatically has an
C<AUTOLOAD> method, which implements the C<AUTOMETHOD()> mechanism
described below. 

It is almost always an error to write your own C<AUTOLOAD()> in any class that
uses Finfo::Std. You almost certainly want to write your own C<AUTOMETHOD()>
instead.

=head2 DUMP(%params)

This method returns a string that represents the internal state (i.e. the
attribute values) of the object on which it's called.

B<Params:>

=over

=item I<separate>   boolean - separate the attributes into their respective classes

=item I<recursive>  boolean - recursively DUMP Finfo::Std objects in attributes (only works on scalar attributes)

=back

=head2 STORABLE_freeze

Allows for use of Storable.

=head2 STORABLE_thaw

Allows for use of Storable.

=head1 Methods that can be supplied by the developer

The following subroutines can be specified as standard methods of a
Finfo::Std class.

=head2 BUILD

When the C<new> constructor of a Class::Std class is called, it
automatically calls every method named C<BUILD> in I<all> the classes
in the new object's hierarchy. That is, when the constructor is called,
it walks the class's inheritance tree (from base classes downwards) and
calls every C<BUILD> method it finds along the way.

This means that, to initialize any class, you merely need to provide a
C<BUILD> method for that class. You don't have to worry about ensuring
that any ancestral C<BUILD> methods also get called; the constructor
will take care of that.

Each C<BUILD> method is called with three arguments: the invocant object,
the identifier number of that object, and a reference to the hash of 
arguments that was originally passed to the constructor:

    sub BUILD {
        my ($self, $ident, $args_ref) = @_;
        ...
    }

=head2 START

Once the initialization values or defaults have been subsequently applied to
attributes, Finfo::Std arranges for any C<START()> methods
in the class's hierarchy to be called befre the constructor finishes.
That is, after the build and default initialization processes are
complete, the constructor walks down the class's inheritance tree a
second time and calls every C<START()> method it finds along the way.

Each C<START()> method is called with three arguments:
the invocant object, the identifier number of that object, and a
reference to (a customized version of) the hash of arguments that was
originally passed to the constructor.

=head2 DEMOLISH

The C<DESTROY()> method that is automatically provided by Finfo::Std ensures
that all the marked attributes of an object, from all the classes in its 
inheritance hierarchy, are automatically cleaned up.

But, if a class requires other destructor behaviours (e.g. closing
filehandles, decrementing allocation counts, etc.) then you may need to
specify those explicitly.

Whenever an object of a Finfo::Std class is destroyed, the C<DESTROY()>
method supplied by Finfo::Std automatically calls every method named
C<DEMOLISH> in I<all> the classes in the new object's hierarchy. That
is, when the destructor is called, it walks the class's inheritance
tree (from derived classes upwards) and calls every C<DEMOLISH()> method it
finds along the way.

This means that, to clean up any class, you merely need to provide a
C<DEMOLISH> method for that class. You don't have to worry about ensuring
that any ancestral C<DEMOLISH> methods also get called; the destructor
will take care of that.

Each C<DEMOLISH> method is called with two arguments: the invocant object,
and the identifier number of that object. For example:

    sub DEMOLISH {
        my ($self, $ident) = @_;

        $filehandle_of{$ident}->flush();
        $filehandle_of{$ident}->close();
    }

Note that the attributes of the object are cleaned up I<after> the
C<DEMOLISH()> method is complete, so they may still be used within
that method.

=head2 AUTOMETHOD

There is a significant problem with Perl's built-in C<AUTOLOAD> mechanism:
there's no way for a particular C<AUTOLOAD> to say "no".

If two or more classes in a class hierarchy have separate C<AUTOLOAD>
methods, then the one belonging to the left-most-depth-first class in
the inheritance tree will always be invoked in preference to any others.
If it can't handle a particular call, the call will probably fail
catastrophically. This means that derived classes can't always be used
in place of base classes (a feature known as "Liskov substitutability")
because their inherited autoloading behaviour may be pre-empted by some
other unrelated base class on their left in the hierarchy.

Finfo::Std provides a mechanism that solves this problem: the
C<AUTOMETHOD> method. An AUTOMETHOD is expected to return either a
handler subroutine that implements the requested method functionality,
or else an C<undef> to indicate that it doesn't know how to handle the
request. Finfo::Std then coordinates every C<AUTOMETHOD> in an object's
hierarchy, trying each one in turn until one of them produces a
suitable handler.

The advantage of this approach is that the first C<AUTOMETHOD> that's
invoked doesn't have to disenfranchise every other C<AUTOMETHOD> in the
hierarchy. If the first one can't handle a particular method call, it
simply declines it and Finfo::Std tries the next candidate instead.

Using C<AUTOMETHOD> instead of C<AUTOLOAD> makes a class
cleaner, more robust, and less disruptive in class hierarchies.
For example:

    package Phonebook;
    use Finfo::Std;
    {
        my %entries_of : ATTR;

        # Any method call is someone's name:
        # so store their phone number or get it...
        sub AUTOMETHOD {
            my ($self, $ident, $number) = @_;

            my $subname = $_;   # Requested subroutine name is passed via $_

            # Return failure if not a get_<name> or set_<name>
            # (Next AUTOMETHOD() in hierarchy will then be tried instead)...
            my ($mode, $name) = $subname =~ m/\A ([gs]et)_(.*) \z/xms
                or return;

            # If get_<name>, return a handler that just returns the old number...
            return sub { return $entries_of{$ident}->{$name}; }
                if $mode eq 'get';

            # Otherwise, set_<name>, so return a handler that
            # updates the entry and then returns the old number...
            return sub {
                $entries_of{$ident}->{$name} = $number;
                return;
            };
        }
    }

    # and later...

    my $lbb = Phonebook->new();

    $lbb->set_Jenny(867_5309);
    $lbb->set_Glenn(736_5000);

    print $lbb->get_Jenny(), "\n";
    print $lbb->get_Glenn(), "\n";

Note that, unlike C<AUTOLOAD()>, an C<AUTOMETHOD()> is called with both the
invocant and the invocant's unique C<ident> number, followed by the actual
arguments that were passed to the method.

Note too that the name of the method being called is passed as C<$_>
instead of C<$AUTOLOAD>, and does I<not> have the class name prepended
to it, so you don't have to strip that name off the front like almost
everyone almost always does in their C<AUTOLOAD()>. If your C<AUTOMETHOD()>
also needs to access the C<$_> from the caller's scope, that's still
available as C<$CALLER::_>.

=head1 Variable traits that can be ascribed

The following markers can be added to the definition of any hash used as an attribute storage within a Finfo::Std class.  These attribute's attributes can be on the same line or on different lines.

Example:

 my %name
    :name(name:r)
    :isa(string)
    :clo('name=s')
    :desc('Your name');
 my %age
    :name(age:r)
    :isa('int pos')
    :clo('age=i')
    :desc('Your age');
 my %color
    :name(favorite_color:o)
    :isa([ 'inlist', __PACKAGE__->colors ])
    :default(__PACKAGE__->default_color)
    :clo(color=s')
    :desc('Favorite color';
 my %friends
    :name(friends:o)
    :ds(aryref)
    :isa(string)
    :clo('fr=@s{,}')
    :desc('Your friends');
 
 sub colors
 {
    return (qw/ red blue green yellow /);
 }

 sub default_color
 {
    return (colors())[0];
 }

=head2 name

 :name(attribute_name:[rop])

Required.  No need to quote.  The name of the attribute.  It is followed by a colon (:) and an r, o or p.  The letter after the name tells Std how to treat the attribute: 

=over

=item B<r> - required on object creation

=item B<o> - optional on object creation

=item B<p> - private attribute 

=back

=head2 ds (data structure)

 :ds(aryref)

Optional, default is 'scalar.'  Indicats the data structure of the attribute.

=over

=item scalar - a single value

=item aryref - a reference to an array 

=item scalar - a reference to a hash

=back

=head2 isa

 :isa('int pos')
 
Optional, default is 'defined.'  Use quotes when using int.  Also accepts aryrefs for types in_list and object.  What the attribute is.  If there is a ds (data structure) indicated, each value will also be checked.  See Finfo::Validate for a list of valid isas. 

=head2 default

Optional, no default.  Can be a string, array ref, etc.  Sets this value to the attribute if none is given oin object creation.  Can also be used to set private attributes with a default.

=head2 clo

Optional, no default.  Use quotes.  The commandline option for the attribute.  Use quotes around the value.  This allows Finfo::CommandLineOptions use a class and derive commanline options form the attributes, then get them when a command is executed.  Use only when wanting to allow the attribute to be assigned from the commandline.  See GetOpt:Long for formatting.

=head2 desc

Optional, no default.  The description of the attribute.

=head1 Method attributes that can be ascribed

=head2 :RESTRICTED and :PRIVATE

Occasionally, it is useful to be able to create subroutines that can only be
accessed within a class's own hierarchy (that is, by derived classes). And
sometimes it's even more useful to be able to create methods that can only be
called within a class itself.

Typically these types of methods are I<utility> methods: subroutines
that provide some internal service for a class, or a class hierarchy.
Finfo::Std supports the creation of these kinds of methods by providing two
special markers: C<:RESTRICTED()> and C<:PRIVATE()>.

Methods marked C<:RESTRICTED()> are modified at the end of the
compilation phase so that they throw an exception when called from
outside a class's hierarchy. Methods marked C<:PRIVATE()> are modified
so that they throw an exception when called from outside the class in
which they're declared.

For example:

    package DogTag;
    use Finfo::Std;
    {
        my %ID_of   : ATTR;
        my %rank_of : ATTR;

        my $ID_num = 0;

        sub _allocate_next_ID : RESTRICTED {
            my ($self) = @_;
            $ID_of{ident $self} = $ID_num++;
            return;
        }

=head2 :CUMULATIVE

One of the most important advantages of using the C<BUILD()> and C<DEMOLISH()>
mechanisms supplied by Finfo::Std is that those methods don't require
nested calls to their ancestral methods, via the C<SUPER> pseudo-class. The
constructor and destructor provided by Finfo::Std take care of the
necessary redispatching automatically. Each C<BUILD()> method can focus
solely on its own responsibilities; it doesn't have to also help
orchestrate the cumulative constructor effects across the class
hierarchy by remembering to call C<< $self->SUPER::BUILD() >>.

Moreover, calls via C<SUPER> can only ever call the method of exactly one
ancestral class, which is not sufficient under multiple inheritance.

Finfo::Std provides a different way of creating methods whose effects
accumulate through a class hierarchy, in the same way as those of
C<BUILD()> and C<DEMOLISH()> do. Specifically, the module allows you to define
your own "cumulative methods".

An ordinary non-cumulative method hides any method of the same name
inherited from any base class, so when a non-cumulative method is
called, only the most-derived version of it is ever invoked. In
contrast, a cumulative method doesn't hide ancestral methods of the same
name; it assimilates them. When a cumulative method is called, the
most-derived version of it is invoked, then any parental versions, then any
grandparental versions, etc. etc, until every cumulative method of the
same name throughout the entire hierarchy has been called.

For example, you could define a cumulative C<describe()> method to the various
classes in a simple class hierarchy like so:

    package Wax::Floor;
    use Finfo::Std;
    {
        my %name_of    :name(name:r);
        my %patent_of  :name(patent:r);

        sub describe :CUMULATIVE {
            my ($self) = @_;

            print "The floor wax $name_of{ident $self} ",
                  "(patent: $patent_of{ident $self})\n";

            return;
        }
    }

    package Topping::Dessert;
    use Finfo::Std;
    {
        my %name_of     :name(name:r);
        my %flavour_of  :name(flavour:r);

        sub describe :CUMULATIVE {
            my ($self) = @_;

            print "The dessert topping $name_of{ident $self} ",
                  "with that great $flavour_of{ident $self} taste!\n";

            return;
        }
    }

    package Shimmer;
    use base qw( Wax::Floor  Topping::Dessert );
    use Finfo::Std;
    {
        my %name_of    :name(name:r);
        my %patent_of  :name(patent:r);

        sub describe :CUMULATIVE {
            my ($self) = @_;

            print "New $name_of{ident $self} ",
                  "(patent: $patent_of{ident $self})\n",
                  "Combining...\n";

            return;
        }
    }

Because the various C<describe()> methods are marked as being cumulative, a
subsequent call to:

    my $product 
        = Shimmer->new({
              name    => 'Shimmer',
              patent  => 1562516251,
              flavour => 'Vanilla',
          });

    $product->describe();

will work its way up through the classes of Shimmer's inheritance tree
(in the same order as a destructor call would), calling each C<describe()>
method it finds along the way. So the single call to C<describe()> would
invoke the corresponding method in each class, producing:

    New Shimmer (patent: 1562516251)
    Combining...
    The floor wax Shimmer (patent: 1562516251)
    The dessert topping Shimmer with that great Vanilla taste!

Note that the accumulation of C<describe()> methods is hierarchical, and
dynamic in nature. That is, each class only sees those cumulative
methods that are defined in its own package or in one of its ancestors.
So calling the same C<describe()> on a base class object:

    my $wax 
        = Wax::Floor->new({ name=>'Shimmer ', patent=>1562516251 });

    $wax->describe();

only invokes the corresponding cumulative methods from that point on up
the hierarchy, and hence only prints:

    The floor wax Shimmer (patent: 1562516251)

Cumulative methods also accumulate their return values. In a list
context, they return a (flattened) list that accumulates the lists
returned by each individual method invoked.

In a scalar context, a set of cumulative methods returns an object that,
in a string context, concatenates individual scalar returns to produce a
single string. When used as an array reference that same scalar-context-return
object acts like an array of the list context values. When used as a hash
reference, the object acts like a hash whose keys are the classnames from the
object's hierarchy, and whose corresponding values are the return values of
the cumulative method from that class.

For example, if the classes each have a cumulative method that returns
their list of sales features:

    package Wax::Floor;
    use Finfo::Std;
    {
        sub feature_list :CUMULATIVE {
            return ('Long-lasting', 'Non-toxic', 'Polymer-based');
        }
    }

    package Topping::Dessert;
    use Finfo::Std;
    {
        sub feature_list :CUMULATIVE {
            return ('Low-carb', 'Non-dairy', 'Sugar-free');
        }
    }

    package Shimmer;
    use Finfo::Std;
    use base qw( Wax::Floor  Topping::Dessert );
    {
        sub feature_list :CUMULATIVE {
            return ('Multi-purpose', 'Time-saving', 'Easy-to-use');
        }
    }

then calling feature_list() in a list context:

    my @features = Shimmer->feature_list();
    print "Shimmer is the @features alternative!\n";

would produce a concatenated list of features, which could then be
interpolated into a suitable sales-pitch:

    Shimmer is the Multi-purpose Time-saving Easy-to-use
    Long-lasting Non-toxic Polymer-based Low-carb Non-dairy
    Sugar-free alternative!

It's also possible to specify a set of cumulative methods that
start at the base class(es) of the hierarchy and work downwards, the way
BUILD() does. To get that effect, you simply mark each method with
:CUMULATIVE(BASE FIRST), instead of just :CUMULATIVE. For example:

    package Wax::Floor;
    use Finfo::Std;
    {
        sub active_ingredients :CUMULATIVE(BASE FIRST) {
            return "\tparadichlorobenzene, cyanoacrylate, peanuts\n";
        }
    }

    package Topping::Dessert;
    use Finfo::Std;
    {
        sub active_ingredients :CUMULATIVE(BASE FIRST) {
            return "\tsodium hypochlorite, isobutyl ketone, ethylene glycol\n";
        }
    }

    package Shimmer;
    use Finfo::Std;
    use base qw( Wax::Floor  Topping::Dessert );

    {
        sub active_ingredients :CUMULATIVE(BASE FIRST) {
            return "\taromatic hydrocarbons, xylene, methyl mercaptan\n";
        }
    }

So a scalar-context call to active_ingredients():

    my $ingredients = Shimmer->active_ingredients();
    print "May contain trace amounts of:\n$ingredients";

would start in the base classes and work downwards, concatenating base-
class ingredients before those of the derived class, to produce:

    May contain trace amounts of:
        paradichlorobenzene, cyanoacrylate, peanuts
        sodium hypochlorite, isobutyl ketone, ethylene glycol
        aromatic hydrocarbons, xylene, methyl mercaptan

Or, you could treat the return value as a hash:

    print Dumper(\%$ingredients);

and see which ingredients came from where:

    $VAR1 = {
       'Shimmer'
            => 'aromatic hydrocarbons, xylene, methyl mercaptan',

       'Topping::Dessert'
            => 'sodium hypochlorite, isobutyl ketone, ethylene glycol',

        'Wax::Floor'
            => 'Wax: paradichlorobenzene,  hydrogen peroxide, cyanoacrylate',
    };

Note that you can't specify both C<:CUMULATIVE> and C<:CUMULATIVE(BASE
FIRST)> on methods of the same name in the same hierarchy. The resulting
set of methods would have no well-defined invocation order, so
Finfo::Std throws a compile-time exception instead.


=head2 :STRINGIFY

If you define a method and add the C<:STRINGIFY> marker then that method
is used whenever an object of the corresponding class needs to be
coerced to a string. In other words, instead of:

    # Convert object to a string...
    sub as_str {
        ...
    }

    # Convert object to a string automatically in string contexts...
    use overload (
        q{""}    => 'as_str',
        fallback => 1,
    );

you can just write:

    # Convert object to a string (automatically in string contexts)...
    sub as_str : STRINGIFY {
        ...
    }


=head2 :NUMERIFY

If you define a method and add the C<:NUMERIFY> marker then that method
is used whenever an object of the corresponding class needs to be
coerced to a number. In other words, instead of:

    # Convert object to a number...
    sub as_num {
        ...
    }

    # Convert object to a string automatically in string contexts...
    use overload (
        q{0+}    => 'as_num',
        fallback => 1,
    );

you can just write:

    # Convert object to a number (automatically in numeric contexts)...
    sub as_num : NUMERIFY {
        ...
    }


=head2 :BOOLIFY

If you define a method and add the C<:BOOLIFY> marker then that method
is used whenever an object of the corresponding class needs to be
coerced to a boolean value. In other words, instead of:

    # Convert object to a boolean...
    sub as_bool {
        ...
    }

    # Convert object to a boolean automatically in boolean contexts...
    use overload (
        q{bool}    => 'as_bool',
        fallback => 1,
    );

you can just write:

    # Convert object to a boolean (automatically in boolean contexts)...
    sub as_bool : BOOLIFY {
        ...
    }


=head2 :SCALARIFY :ARRAYIFY :HASHIFY :GLOBIFY :CODIFY

If a method is defined with one of these markers, then it is automatically
called whenever an object of that class is treated as a reference of the
corresponding type.

For example, instead of:

    sub as_hash {
        my ($self) = @_;

        return {
            age      => $age_of{ident $self},
            shoesize => $shoe_of{ident $self},
        };
    }

    use overload (
        '%{}'    => 'as_hash',
        fallback => 1,
    );

you can just write:

    sub as_hash : HASHIFY {
        my ($self) = @_;

        return {
            age      => $age_of{ident $self},
            shoesize => $shoe_of{ident $self},
        };
    }

Likewise for methods that allow an object to be treated as a scalar
reference (C<:SCALARIFY>), a array reference (C<:ARRAYIFY>), a
subroutine reference (C<:CODIFY>), or a typeglob reference
(C<:GLOBIFY>).

=head1 DIAGNOSTICS

=over 

=item Can't find class %s

You tried to call the Finfo::Std::new() constructor on a class 
that isn't built using Finfo::Std. Did you forget to write C<use Finfo::Std>
after the package declaration?

=item Argument to %s->new() must be hash reference

The constructors created by Finfo::Std require all initializer values
to be passed in a hash, but you passed something that wasn't a hash.
Put your constructor arguments in a hash.

=item Missing initializer label for %s: %s

You specified that one or more attributes had initializer values (using the
C<init> argument inside the attribute's C<ATTR> marker), but then failed
to pass in the corresponding initialization value. Often this happens because
the initialization value I<was> passed, but the key specifying the
attribute name was misspelled.

=item Can't make anonymous subroutine cumulative

You attempted to use the C<:CUMULATIVE> marker on an anonymous subroutine.
But that marker can only be applied to the named methods of a class. Convert
the anonymous subroutine to a named subroutine, or find some other way to 
make it interoperate with other methods.

=item Conflicting definitions for cumulative method: %s

You defined a C<:CUMULATIVE> and a C<:CUMULATIVE(BASE FIRST)> method of the
same name in two classes within the same hierarchy. Since methods can only be
called going strictly up through the hierarchy or going strictly down 
through the hierarchy, specifying both directions is obviously a mistake.
Either rename one of the methods, or decide whether they should accumulate
upwards or downwards.

=item Missing new value in call to 'set_%s' method

You called an attribute setter method without providing a new value 
for the attribute. Often this happens because you passed an array that
happened to be empty. Make sure you pass an actual value.

=item Can't locate %s method "%s" via package %s

You attempted to call a method on an object but no such method is defined
anywhere in the object's class hierarchy. Did you misspell the method name, or
perhaps misunderstand which class the object belongs to?

=item %s method %s declared but not defined

A method was declared with a C<:RESTRICTED> or C<:PRIVATE>, like so:

    sub foo :RESTRICTED;
    sub bar :PRIVATE;

But the actual subroutine was not defined by the end of the compilation
phase, when the module needed it so it could be rewritten to restrict or
privatize it.


=item Can't call restricted method %s from class %s

The specified method was declared with a C<:RESTRICTED> marker but
subsequently called from outside its class hierarchy. Did you call the
wrong method, or the right method from the wrong place?


=item Can't call private method %s from class %s

The specified method was declared with a C<:PRIVATE> marker but
subsequently called from outside its own class. Did you call the wrong
method, or the right method from the wrong place?


=item Internal error: %s

Your code is okay, but it uncovered a bug in the Finfo::Std module.
L<BUGS AND LIMITATIONS> explains how to report the problem.

=back


=head1 CONFIGURATION AND ENVIRONMENT

Finfo::Std requires no configuration files or environment variables.


=head1 DEPENDENCIES

Finfo::Std depends on the following modules:

=over

=item *

version

=item *

Scalar::Util

=item *

Data::Dumper

=back

=head1 BUGS AND LIMITATIONS

=over

=item *

Does not handle threading (including C<fork()> under Windows).

=item *

Attribute declarations cannot include variables, since these are not
interpolated into the declaration (a limitation in Perl itself).  But may contain
package methods.  For same package use __PACKAGE__ or get them from another source.

=back

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
