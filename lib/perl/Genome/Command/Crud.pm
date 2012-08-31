package Genome::Command::Crud;

use strict;
use warnings;

use Genome;
      
require Carp;
use Data::Dumper 'Dumper';
require Lingua::EN::Inflect;

class Genome::Command::Crud {
    doc => 'Class for dynamically building CRUD commands',
};

sub camel_case_to_string {
    my $string = shift;
    unless ($string) {
        Carp::confess "No camel case string to convert to string"
    }

    # Split on capital letters
    my @words = split( /(?=(?<![A-Z])[A-Z])|(?=(?<!\d)\d)/, $string);
    my $join = ( @_ ) ? $_[0] : ' '; 
    return join($join, map { lc } @words);
}

sub display_name_for_value {
    my ($class, $value) = @_;

    return 'NULL' if not defined $value;

    my $ref = ref $value;
    if ( $ref eq 'ARRAY' and @$value > 10 ) {
        return scalar @$value.' items';
    }

    my @display_names;
    for my $val ( $ref eq 'ARRAY' ? @$value : $value ) {
        if ( not defined $val ) {
            push @display_names, 'NULL';
        }
        elsif ( not Scalar::Util::blessed($val) ) {
            push @display_names, $val;
        }
        elsif ( my $display_name_sub = $val->can('__display_name__') ) {
            push @display_names, $display_name_sub->($val);
        }
        else {
            push @display_names, $val->id;
        }
    }

    return join(" ", @display_names);
}

our %inited;
sub init_sub_commands {
    my ($class, %incoming_config) = @_;

    # Config: target class, namespace
    Carp::confess('No target class given to init_sub_commands') if not $incoming_config{target_class};
    my %config;
    $config{target_class} = delete $incoming_config{target_class};
    $config{namespace} = ( exists $incoming_config{namespace} )
    ? delete $incoming_config{namespace}
    : $config{target_class}.'::Command';

    # Ok if we inited already
    return 1 if $inited{ $config{namespace} };

    # Names for objects
    my $target_name = ( defined $incoming_config{target_name} )
    ? delete $incoming_config{target_name}
    : join(' ', map { camel_case_to_string($_) } split('::', $config{target_class}));
    Lingua::EN::Inflect::classical(persons => 1);
    $target_name =~ s/[-_]/ /g;
    $config{target_name} = $target_name;
    $config{target_name_ub} = $config{target_name};
    $config{target_name_ub} =~ s/ /_/;
    $config{target_name_pl} = ( exists $incoming_config{target_name_pl} )
    ? delete $incoming_config{target_name_pl} 
    : Lingua::EN::Inflect::PL($target_name);
    $config{target_name_pl_ub} = $config{target_name_pl};
    $config{target_name_pl_ub} =~ s/ /_/;

    # Main tree command
    if ( not $class->_build_command_tree(%config) ) {
        Carp::confess('Failed to create main tree class for '.$config{namespace});
    }

    # Get the current sub commands
    my @namespace_sub_command_classes = $config{namespace}->sub_command_classes;
    my @namespace_sub_command_names = @namespace_sub_command_classes;
    @namespace_sub_command_names = map {
        s/$config{namespace}:://; $_ = lc($_); $_;
    } @namespace_sub_command_names;

    # Create the sub commands
    my @command_names = (qw/ create list update delete /);
    my @command_classes;
    for my $command_name ( @command_names ) {
        # config for this sub command
        my %command_config;
        if ( exists $incoming_config{$command_name} ) {
            my $ref = ref $incoming_config{$command_name};
            if ( not $ref or $ref ne 'HASH' ) {
                Carp::confess("Invalid config for $command_name sub command: ".Dumper($incoming_config{$command_name}));
            }
            %command_config = %{delete $incoming_config{$command_name}};
        }

        # skip existing sub commands, except update
        if ( $command_name ne 'update' and grep { $command_name eq $_ } @namespace_sub_command_names ) {
            next if not %command_config;
            #Carp::confess("Subcommand '$sub_command' for namespace '$config{namespace}' already exists, but there is CRUD config for it. Please correct.");
            next;
        }

        # skip if requested not to init
        next if delete $command_config{do_not_init};

        # build the sub class
        my $method = '_build_'.$command_name.'_command';
        my $sub_class = $class->$method(%config, %command_config);
        if ( not $sub_class ) {
            Carp::confess('Cannot dynamically create class for sub command name: '.$command_name);
        }
        push @command_classes, $sub_class;
        no strict;
        *{ $sub_class.'::_display_name_for_value' } = \&display_name_for_value;
    }

    # Note inited
    $inited{ $config{namespace} } = 1;

    # Check for left over config
    Carp::confess('Unknown config for CRUD commands: '.Dumper(\%incoming_config)) if %incoming_config;

    # Overload sub command classes to return these in memory ones, plus the existing ones
    my %sub_command_classes = map { $_ => 1 } ( @command_classes, @namespace_sub_command_classes );
    no strict;
    *{ $config{namespace}.'::sub_command_classes' } = sub{ return keys %sub_command_classes; };
    
    return 1;
}

sub _build_command_tree {
    my ($class, %config) = @_;

    my $meta = eval{ $config{namespace}->__meta__; };
    return 1 if $meta;

    UR::Object::Type->define(
        class_name => $config{namespace},
        is => 'Command::Tree',
        doc => 'work with '.$config{target_name_pl},
    );

    return 1;
}

sub _build_create_command {
    my ($class, %config) = @_;

    my @exclude = $class->_property_name_ref_to_array($config{exclude});

    my $target_meta = $config{target_class}->__meta__;
    my %properties;
    for my $target_property ( $target_meta->property_metas ) {
        my $property_name = $target_property->property_name;

        next if grep { $property_name eq $_ } @exclude;
        next if $target_property->class_name eq 'UR::Object';
        next if $property_name =~ /^_/;
        next if grep { $target_property->$_ } (qw/ is_calculated is_constant is_transient /);
        next if $target_property->is_id and ($property_name eq 'id' or $property_name =~ /_id$/);
        next if grep { not $target_property->$_ } (qw/ is_mutable /);
        next if $target_property->is_many and $target_property->is_delegated and not $target_property->via; # direct relationship

        my %property = (
            property_name => $property_name,
            data_type => $target_property->data_type,
            is_many => $target_property->is_many,
            is_optional => $target_property->is_optional,
            valid_values => $target_property->valid_values,
            default_value => $target_property->default_value,
            doc => $target_property->doc,
        );

        if ( $property_name =~ s/_id(s)?$// ) {
            $property_name .= $1 if $1;
            my $object_meta = $target_meta->property_meta_for_name($property_name);
            if ( $object_meta and  not grep { $object_meta->$_ } (qw/ is_calculated is_constant is_transient id_class_by /) ) {
                $property{property_name} = $property_name;
                $property{data_type} = $object_meta->data_type;
                $property{doc} = $object_meta->doc if $object_meta->doc;
            }
        }

        $properties{$property{property_name}} = \%property;
    }

    if ( $config{add} ) {
        for my $property_name ( keys %{$config{add}} ) {
            $properties{$property_name} = $config{add}->{$property_name};
        }
    }

    if ( not %properties ) {
        Carp::confess('No properties found for target class: '.$config{target_class});
    }

    my $sub_class = $config{namespace}.'::Create';
    UR::Object::Type->define(
        class_name => $sub_class,
        is => 'Genome::Command::Create',
        has => [ %properties ],
        doc => 'create '.$config{target_name_pl},
    );

    no strict;
    *{$sub_class.'::_target_class'} = sub{ return $config{target_class}; };
    *{$sub_class.'::_target_name'} = sub{ return $config{target_name}; };
    *{$sub_class.'::_before'} = $config{before} if $config{before};

    return $sub_class;
}

sub _build_list_command {
    my ($class, %config) = @_;

    my @has =  (
        subject_class_name  => {
            is_constant => 1,
            value => $config{target_class},
        },
    );
    if ( $config{show} ) {
        push @has, show => { default_value => $config{show}, };
    }
    if ( $config{order_by} ) {
        push @has, order_by => { default_value => $config{order_by}, };
    }

    my $list_command_class_name = $config{namespace}.'::List';
    UR::Object::Type->define(
        class_name => $list_command_class_name,
        is => 'UR::Object::Command::List',
        has => \@has,
    );

    no strict;
    *{$list_command_class_name.'::help_brief'} = sub{ return $config{target_name_pl}; };

    return $list_command_class_name;
}
   
sub _build_update_command {
    my ($class, %config) = @_;

    # Config
    # include only these properties
    my @include_only = $class->_property_name_ref_to_array($config{include_only});
    my @exclude = $class->_property_name_ref_to_array($config{exclude});
    if ( @include_only and @exclude ) {
        Carp::confess('Can not include only and exclude update sub commands!');
    }

    # only if null
    my (%only_if_null, $all_only_if_null);
    if ( my $only_if_null = delete $config{only_if_null} ) {
        my $ref = ref $only_if_null;
        if ( $only_if_null eq 1 ) {
            $all_only_if_null = 1;
        }
        elsif ( not $ref ) {
            Carp::confess("Unknown 'only_if_null' config: $only_if_null");
        }
        else {
            %only_if_null = map { $_ => 1 } map { s/_id$//; $_; } ( $ref eq 'ARRAY' ? @$only_if_null : keys %$only_if_null )
        }
    }

    # Update tree
    my $update_class_name = $config{namespace}.'::Update';
    my $update_class = eval{ $update_class_name->class; };
    my (@update_sub_commands, @update_sub_command_names);
    if ( not $update_class ) {
        UR::Object::Type->define(
            class_name => $update_class_name,
            is => 'Genome::Command::UpdateTree',
            doc => 'properties on '.$config{target_name_pl},
        );
    }
    else {
        @update_sub_commands = $update_class_name->sub_command_classes;
        @update_sub_command_names = $update_class_name->sub_command_names;
    }

    # Properties make a command for each
    my $target_meta = eval{ $config{target_class}->__meta__; };
    my %properties_seen;
    PROPERTY: for my $target_property ( $target_meta->property_metas ) {
        my $property_name = $target_property->property_name;
        next if grep { $property_name eq $_ } @update_sub_command_names;
        next if @include_only and not grep { $property_name =~ /^$_(_id)?$/ } @include_only;
        next if @exclude and grep { $property_name =~ /^$_(_id)?$/ } @exclude;

        next if $target_property->class_name eq 'UR::Object';
        next if $property_name =~ /^_/;
        next if grep { $target_property->$_ } (qw/ is_id is_calculated is_constant is_transient /);
        next if grep { not $target_property->$_ } (qw/ is_mutable /);
        next if $target_property->is_many and $target_property->is_delegated and not $target_property->via; # direct relationship

        my %property = (
            name => $target_property->singular_name,
            name_pl => $target_property->plural_name,
            is_many => $target_property->is_many,
            data_type => $target_property->data_type,
            doc => $target_property->doc,
        );
        $property{valid_values} = $target_property->valid_values if defined $target_property->valid_values;

        if ( $property_name =~ s/_id(s)?$// ) {
            $property_name .= $1 if $1;
            my $object_meta = $target_meta->property_meta_for_name($property_name);
            if ( $object_meta ) {
                next if grep { $object_meta->$_ } (qw/ is_calculated is_constant is_transient id_class_by /);
                $property{name} = $object_meta->singular_name;
                $property{name_pl} = $object_meta->plural_name;
                $property{is_optional} = $object_meta->is_optional;
                $property{data_type} = $object_meta->data_type;
                $property{doc} = $object_meta->doc if $object_meta->doc;
            }
        }
        next if $properties_seen{$property_name};
        $properties_seen{$property_name} = 1;

        $config{property} = \%property;
        $config{only_if_null} = ( $all_only_if_null or exists $only_if_null{$property_name} ) ? 1 : 0;
        my $update_sub_command;
        if ( $property{is_many} ) {
            $update_sub_command = $class->_build_add_remove_property_sub_commands(%config);
        }
        else {
            $update_sub_command = $class->_build_update_property_sub_command(%config);
        }
        push @update_sub_commands, $update_sub_command if $update_sub_command;
    }

    no strict;
    *{$update_class_name.'::sub_command_classes'} = sub{ return @update_sub_commands; };

    return $update_class_name;
}

sub _build_update_property_sub_command {
    my ($class, %config) = @_;

    my $property = $config{property};
    my $update_property_class_name = $config{namespace}.'::Update::'.join('', map { ucfirst } split('_', $property->{name}));
    my $update_property_class = eval{ $update_property_class_name->class; };
    return if $update_property_class; # OK

    UR::Object::Type->define(
        class_name => $update_property_class_name,
        is => 'Genome::Command::UpdateProperty',
        has => [ 
            $config{target_name_pl_ub} => {
                is => $config{target_class},
                is_many => 1,
                shell_args_position => 1,
                doc => ucfirst($config{target_name_pl}).' to update, resolved via string.',
            },
            value => {
                is => $property->{data_type},
                valid_values => $property->{valid_values},
                doc => $property->{doc},
            },
        ],
        doc => 'update '.$config{target_name_pl}.' '.$property->{name},
    );

    no strict;
    *{ $update_property_class_name.'::_target_name_pl' } = sub{ return $config{target_name_pl}; };
    *{ $update_property_class_name.'::_target_name_pl_ub' } = sub{ return $config{target_name_pl_ub}; };
    *{ $update_property_class_name.'::_property_name' } = sub{ return $property->{name}; };
    *{ $update_property_class_name.'::_property_doc' } = sub{ return $property->{doc}; } if $property->{doc};
    *{ $update_property_class_name.'::_only_if_null' } = sub{ return $config{only_if_null}; };
    *{ $update_property_class_name.'::_display_name_for_value' } = \&display_name_for_value;

    return $update_property_class_name;
}

sub _build_add_remove_property_sub_commands { 
    my ($class, %config) = @_;

    my $property = $config{property};
    my $update_tree_class_name = $config{namespace}.'::Update::'.join('', map { ucfirst } split('_', $property->{name_pl}));
    UR::Object::Type->define(
        class_name => $update_tree_class_name,
        is => 'Command::Tree',
        doc => 'add/remove '.$property->{name_pl},
    );

    my @update_sub_command_class_names;
    no strict;
    *{$update_tree_class_name.'::_target_name'} = sub{ return $config{target_name}; };
    *{$update_tree_class_name.'::_target_name_ub'} = sub{ return $config{target_name_ub}; };
    *{$update_tree_class_name.'::_property_name'} = sub{ return $property->{name}; };
    *{$update_tree_class_name.'::_property_name_pl'} = sub{ return $property->{name_pl}; };
    *{$update_tree_class_name.'::_display_name_for_value'} = \&display_name_for_value;
    *{$update_tree_class_name.'::sub_command_classes'} = sub{ return @update_sub_command_class_names; };
    use strict;

    for my $function (qw/ add remove /) {
        my $update_sub_command_class_name = $update_tree_class_name.'::'.ucfirst($function);
        push @update_sub_command_class_names, $update_sub_command_class_name;
        UR::Object::Type->define(
            class_name => $update_sub_command_class_name,
            is => 'Genome::Command::AddRemoveProperty',
            has => [ 
                $config{target_name_pl_ub} => {
                    is => $config{target_class},
                    is_many => 1,
                    shell_args_position => 1,
                    doc => ucfirst($config{target_name_pl}).' to update, resolved via string.',
                },
                'values' => => {
                    is => $property->{data_type},
                    is_many => 1,
                    valid_values => $property->{valid_values},
                    doc => $property->{doc},
                },
            ],
            #doc => $function.' '.$property->{name_pl}.' to '.$config{target_name_pl},
            doc => $config{target_name_pl}.' '.$function.' '.$property->{name_pl},
        );
        no strict;
        *{$update_sub_command_class_name.'::_add_or_remove'} = sub{ return $function; };
        *{$update_sub_command_class_name.'::_target_name'} = sub{ return $config{target_name}; };
        *{$update_sub_command_class_name.'::_target_name_pl'} = sub{ return $config{target_name_pl}; };
        *{$update_sub_command_class_name.'::_target_name_pl_ub'} = sub{ return $config{target_name_pl_ub}; };
        *{$update_sub_command_class_name.'::_property_name'} = sub{ return $property->{name}; };
        *{$update_sub_command_class_name.'::_property_name_pl'} = sub{ return $property->{name_pl}; };
        *{$update_sub_command_class_name.'::_display_name_for_value'} = \&display_name_for_value;
    }

    return $update_tree_class_name;
}

sub _build_delete_command {
    my ($class, %config) = @_;

    my $sub_class = $config{namespace}.'::Delete';
    UR::Object::Type->define(
        class_name => $sub_class,
        is => 'Genome::Command::Delete',
        has => [ 
            $config{target_name_pl_ub} => {
                is => $config{target_class},
                is_many => 1,
                shell_args_position => 1,
                require_user_verify => 1, # needed?
                doc => ucfirst($config{target_name_pl}).' to delete, resolved via text string.',
            },
        ],
        doc => 'delete '.$config{target_name_pl},
    );

    no strict;
    *{ $sub_class.'::_target_name_pl' } = sub{ return $config{target_name_pl}; };
    *{ $sub_class.'::_target_name_pl_ub' } = sub{ return $config{target_name_pl_ub}; };

    return $sub_class;
}

sub _property_name_ref_to_array {
    my ($class, $property_names) = @_;

    return if not $property_names;

    my @property_names;
    my $ref = ref $property_names;
    if ( not $ref ) {
        @property_names = $property_names;
    }
    elsif ( $ref eq 'HASH' ) {
        @property_names = keys %$property_names;
    }
    else { # ARRAY
        @property_names = @$property_names;
    }

    return map { s/_id$//; $_; } @property_names; # rm '_id'
}

1;

