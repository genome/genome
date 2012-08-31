package Finfo::StepConfiguration;

use strict;
use warnings;
use Data::Dumper;
use IO::File;

use base qw(Finfo::Object);

sub _attrs
{
    return {
        'configuration_file:r' =>
        {
            type => 'input_file'
        },
        'steps:p' =>
        {
            type => 'hashref'
        },
        'requirements:p' =>
        {
            type => 'hashref'
        },
        'current_step:p' =>
        {
            type => 'defined'
        },
    };
}

sub _init
{
    my $self= shift;
    $self->requirements({});
    $self->steps({});
    $self->current_step('start');
    $self->parse_conf_file;
    return 1;
}

sub parse_conf_file{
    my $self = shift;
    my $fh = IO::File->new("< ".$self->configuration_file);
    while ( my $line = $fh->getline ){
        if ($line =~ /^STEPS\{/){       #Begin parsing conf file
            my $prev_step;          #define prev step for building step hash
            while ( my $line = $fh->getline ){
                chomp $line;
                last if $line =~ /\}/;
                my $step = $line;
                $step =~ s/\[//;
                if ( $line =~ /\[/ ){
                    $self->parse_step_conf($step, $fh);
                }
                $self->steps->{$prev_step} = $step if defined $prev_step;
                $prev_step = $step;
            }
        }
    }
}

sub parse_step_conf{
    my ($self, $step, $fh) = @_;
    $self->requirements->{$step}=[];
    while ( my $line = $fh->getline ){
        chomp $line;
        return if $line =~ /\]/;
        my %reqhash;
        my @class_array;
        my ($class, $requirementstring) = split(/\s+/,$line);
        $self->error_msg("No class info for step!") and return undef unless defined $class;
        if (defined $requirementstring){
            my @reqirements = split(/,/,$requirementstring);
            foreach (@reqirements){
                my ($requirement, $value) = split(/=/,$_);
                if (defined $value){
                    $reqhash{$requirement}=$value;
                }else{
                    $reqhash{$requirement}=1;
                }
            }
        }
        push @class_array, $class;
        push @class_array, \%reqhash if defined \%reqhash;
        push @{$self->requirements->{$step}}, \@class_array;
    }
}

sub next_step{
    my ($self, $step) = @_;
    if (defined $step){
        $self->error_msg("step name provided does not have a next step!") 
            and return undef unless grep {$step eq $_} keys %{$self->steps};
        return $self->steps->{$step};
    }else{
        $step = $self->current_step;
        my $next_step = $self->steps->{$step};
        $self->current_step($next_step);
        return $next_step;
    }
}

sub requirements_for_step{
    my ($self, $step) = @_;
    if (defined $step){
        $self->error_msg("step name provided does not have a class/requirement data!") 
            and return undef unless grep {$step eq $_} keys %{$self->requirements};
    }else{
        $step = $self->current_step;
    }
    return $self->requirements->{$step};
}

=pod

=head1 NAME
StepConfiguration - parser for configuration files 

=head1 SYNOPSIS
in configuration file:

STEPS{
initialize[
Example::Initialize
]
execute[
Example::SpecialExecute             special_flag
Example::Execute
]
cleanup[
Example::Cleanup                    cleanup_flag
]
}

=head1 DESCRIPTION 

=head3 File Format:
STEPS{
<step-block1>
<step-block2>
...
}
Determines step execution order

=head3 Format of step_block
<step_name>[
Classname_for_step     <requirements>
Alternate_classname    <requirements>
]

Multiple class choices are evaluated in order, as soon as a requirement is met, that class is assigned for the step.
If no requirements are given for a class, that class is used for the step and the step is executed.

=head3 Format of requirements
<requirement1>,<requirement2>=<value>
<requirement>=<value>
comma separated requirements are or'd to determine if step will be executed, if the requirement is not associated with a value, it will be checked for definition

Requirements must be valid keys in the Finfo::StepFactory parameters. 
multiple class choices are evaluated in order, as soon as a requirement is met, that class is assigned for the step
if no requirements are given for a class, that class is used for the step and the step is executed

=cut

1;

#$HeadURL$
#$Id$


