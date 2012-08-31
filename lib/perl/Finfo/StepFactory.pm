package Finfo::StepFactory;

use strict;
use warnings;
use Data::Dumper;
use Finfo::StepConfiguration;
use base qw(Finfo::Object);

our $AUTOLOAD;

sub _attrs
{
    return
    {
        'configuration_file:r' => 
        {
            type => 'input_file',
        },
        'parameters:r' => 
        {
            type => 'non_empty_hashref',
        },
        'configuration_object:p' => 
        {
            type => 'object'
        },
    };
}

sub _init
{
    my ($self) = @_;
    my $conf = Finfo::StepConfiguration->new( configuration_file => $self->configuration_file );
    return 0 unless defined $conf;
    $self->check_parameters; 
    $self->configuration_object($conf);
    return 1;
}

sub check_parameters{
    my ($self) = @_;
    my @accessors = keys %{$self->parameters};
    my @forbidden_words = qw(_attrs _init next_step get_step create_step_obj execute);
    foreach my $accessor (@accessors){
       $self->warning_message
       (
           "Parameter $accessor conflicts with a reserved name of a method in the factory class.  If use of this varibale name is continued, it will only be accessable through \$factory_obj->parameters->{\$accessor}"
       ) if grep { $accessor eq $_ } @forbidden_words;
    }
}


sub AUTOLOAD
{
    my $call = $AUTOLOAD;
    my ($self, $param ) = @_;
    $call =~ s/(.*):://;
    my $caller = "$1";
    return if $call =~ /^[A-Z]/;
    my $params = $self->{parameters};
    if (grep {$call eq $_} keys %$params){
        if (defined $param){
            $params->{$call} = $param;
            return $param;
        }else{
            return $params->{$call} 
        }
    }
    $self->error_msg("method call $call to $caller is invalid!") and die;
}

sub next_step{
    my ($self) = @_;
    return $self->configuration_object->next_step;
}

sub get_step{
    my ($self, $step) = @_;
    my $reqs = $self->configuration_object->requirements_for_step($step);
    my $class;
    my $success = 0;
    $success = 1 unless scalar @$reqs > 0;
    foreach my $reqset (@$reqs){
        my ($class_name, $reqhash) = @$reqset;
        $class = $class_name;
        unless (defined $reqhash){
            $success = 1;
            last;
        }
        foreach my $key (keys %$reqhash){
            my $result = $reqhash->{$key};
            next unless defined $self->parameters->{$key};
            $success = 1 if $result eq 1;
            $success = 1 if $self->parameters->{$key} eq $result;
        }
        last if $success;
    }
    if ($success){
        $self->debug_msg("Creating object of class $class for step $step");
        return $self->create_step_obj($class);
    }else{
        $self->debug_msg("Step $step doesn't meet requirements in config file")
            and return $self->get_step($self->next_step);
    }
}

sub create_step_obj{
    my ($self, $class) = @_;
    require $class.".pm";
    import $class;
    my %hash;
    foreach my $attr($class->required_attributes,$class->optional_attributes){
        if ($attr eq 'factory'){
            $hash{$attr}=$self;
        }elsif (defined $self->parameters->{$attr}){
            $hash{$attr}=$self->parameters->{$attr};
        }
    }
    my $obj = $class->new(%hash);
    return $obj if defined $obj and $obj->isa($class);
    $self->error_msg("Couldn't create object of class $class!") and return undef;
    
}

sub execute
{
    my ($self) = @_;
    while (my $step_name = $self->next_step){
        my $step = $self->get_step($step_name);
        $self->error_message("Couldn't get a step obj for $step_name. Exiting") and die unless defined $step;
        $self->info_msg("Beginning step: $step_name");
        eval{
            $step->execute;
        }
    }
}

=pod

=head1 NAME
StepFactory - This is a module use to manage a set of Finfo::Object 'steps' in a process.

=head1 SYNOPSIS
The StepFactory module provides a simple interface for executing several Finfo::Object modules in sequence.
my $factory = Finfo::StepFactory->new(configuration_file => $conf_file, parameters => \%parameters);
$factory->execute;

The configuration file determines step order and implementation details, the parameters contain the hashref of parameters to be used by all steps.

Each step should be a finfo object and contain an execute method which calls all necessary code for the step.

=head1 DESCRIPTION 

StepFactory was created to facilitate the creation, design, and maintenance of multi-part operations and applications.  By breaking down a long process into step modules, inheriting from Finfo::Object, different steps can be modified, maintained or swapped out for an entirely new implementation without affecting other parts of the process.  

Using a configuration file, one can define a series of steps to be executed in order.  For each step a class is provided that implements the step.  Multiple classes can be defined for a step, as well as requirements for a certain class to be used for the step.  If no classes satisfy their requirements, that step is skipped.

=head1 SEE ALSO

Finfo::Object for details on building a Finfo::Object

Finfo::StepConfiguration for information on the configuration file syntax

=cut

1;

#$HeadURL$
#$Id$


