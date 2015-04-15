package Genome::Model::Tools::Sx::Validate;

use strict;
use warnings;

use Genome;

use Params::Validate qw(:types);
use List::MoreUtils;

class Genome::Model::Tools::Sx::Validate { 
    has => {
        command => { is => 'Text', },
    },
    has_transient => {
        command_class => { is => 'Text', },
        command_params => { is => 'HASH', },
    },
};

sub is_valid {
    my ($self) = Params::Validate::validate_pos(@_, {type => HASHREF});

    my $command = $self->command;
    if ( not defined $command ) { 
        $self->error_message('Cannot validate command. None given.');
        return;
    }

    my @tokens = split(/\s+/, $command);
    my $gmt = shift @tokens;
    my $sx = $tokens[0];
    if ( not $gmt or $gmt ne 'gmt' or not $sx or $sx ne 'sx' ) {
        $self->error_message("Command ($command) must start w/ 'gmt sx'");
        return;
    }

    my @subclass_parts;
    while ( my $token = shift @tokens ) {
        if ( $token =~ /^\-/ ) {
            unshift @tokens, $token;
            last;
        }
        push @subclass_parts, $token;
    }

    unless ( @subclass_parts ) {
        $self->error_message("Could not get class from command: $command");
        return;
    }

    my $class = 'Genome::Model::Tools::'.
    join(
        '::', 
        map { Genome::Utility::Text::string_to_camel_case($_) }
        map { s/\-/ /g; $_; }
        @subclass_parts
    );

    my $class_meta = eval{ $class->__meta__; };
    if ( not $class_meta ) {
        $self->error_message("Cannot validate class ($class) for command ($command): $@");
        return;
    }

    my %params;
    if ( @tokens ) {
        my $params_string = join(' ', @tokens);
        eval{
            %params = Genome::Utility::Text::param_string_to_hash(
                $params_string
            );
        };
        unless ( %params ) {
            $self->error_message("Can't get params from params string: $params_string");
            return;
        }
    }

    my %converted_params;
    for my $key ( keys %params ) {
        if ( $key =~ /_/ ) { # underscores not allowed
            $self->error_message("Param ($key) for command ($command) params has an underscore. Use dashes (-) instead");
            return;
        }
        my $new_key = $key; 
        $new_key =~ s/\-/_/g; # sub - for _ to create processor
        my $property = $class_meta->property_meta_for_name($new_key);
        unless ( $property ) {
            $self->error_message("Cannot find property ($new_key) in class ($class)");
            return;
        }
        if ( $property->is_many and not ref($params{$key}) ) {
            $converted_params{$new_key} = [ $params{$key} ];
        }
        else {
            $converted_params{$new_key} = $params{$key};
        }
    }
    $self->command_class($class);
    $self->command_params(\%converted_params);

    my $obj = eval{ $class->create(%converted_params); };
    unless ( $obj ) {
        $self->error_message("Can't validate command ($command) using class ($class)".( $@ ? ": $@" : '') );
        return;
    }

    my @errors = $obj->__errors__;
    if ( @errors ) {
        $self->error_message("Found errors with command:\n".join("\n", map { $_->__display_name__ } @errors));
        return;
    }
    $obj->delete;

    return 1;
}

sub number_of_threads_required {
    my $self = shift;
    my $command_params = $self->command_params;
    die $self->error_message('No command params set! Please make sure the command is_valid') if not $command_params;
    my $param_name = List::MoreUtils::firstval { $_ =~ /thread/ } keys %$command_params;
    return 1 if not $param_name;
    return $command_params->{$param_name};
}

1;

