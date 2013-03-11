package Genome;
use warnings;
use strict;

our $VERSION = '0.080001';
$DB::deep = 10000;

use UR;
use UR::ObjectV001removed;
use File::Temp;
use IO::String;
use File::Basename;
use Carp;
use Carp::Heavy;

# Standard namespace declaration for a UR namespace
UR::Object::Type->define(
    class_name => 'Genome',
    is => ['UR::Namespace'],
    english_name => 'genome',
);

# Local configuration
require Genome::Site;

# Checks that all variables that start with GENOME_ have a corresponding Genome/Env/* module
# and assigns default values to any variables that have one set.
require Genome::Env;

# If the search engine is installed, configure its hooks
eval {
    local $SIG{__WARN__};
    local $SIG{__DIE__};
    require Genome::Search;
};

require Genome::Site::TGI::LegacyConfig;

# This ensures that the search system is updated when certain classes are updated
# The search system is optional so it skips this if usage above fails
if ($INC{"Genome/Search.pm"}) {
    Genome::Search->register_callbacks('Genome::Searchable');
}

# Account for a perl bug in pre-5.10 by applying a runtime patch to Carp::Heavy
if ($] < 5.01) {
    no warnings;
    *Carp::caller_info = sub {
        package
            Carp;
        our $MaxArgNums;
        my $i = shift(@_) + 1;
        package DB;
        my %call_info;
        @call_info{
            qw(pack file line sub has_args wantarray evaltext is_require)
        } = caller($i);

        unless (defined $call_info{pack}) {
            return ();
        }

        my $sub_name = Carp::get_subname(\%call_info);
        if ($call_info{has_args}) {
            # SEE IF WE CAN GET AROUND THE BIZARRE ARRAY COPY ERROR...
            my @args = ();
            if ($MaxArgNums and @args > $MaxArgNums) { # More than we want to show?
                $#args = $MaxArgNums;
                push @args, '...';
            }
            # Push the args onto the subroutine
            $sub_name .= '(' . join (', ', @args) . ')';
        }
        $call_info{sub_name} = $sub_name;
        return wantarray() ? %call_info : \%call_info;
    };
    use warnings;
}

require Genome::Sys;

# this is only on the publicly-released GMS to handle things we must ship for functionality
# but will not be part of individual package releases.
require Genome::Site::Deprecated;

# fixes which have not yet gone into the main-line master branch
require Genome::Patch; 

# DB::single is set to this value in many places, creating a source-embedded break-point
# set it to zero in the debugger to turn off the constant stopping...
$DB::stopper = 1;

# This is tested by Genome/SoftwareResult/Default.t
# And also by Genome/Command/WithSavedResult-buildwrapper.t
# It is only used by the full GMS install.
# They probably needs a better home than the base module.

sub __extend_namespace__ {
    my ($self,$ext) = @_;

    # Auto-generate a Genome::SoftwareResult::Default class for any Command::V2.
    # The class will have the name ${COMMAND_CLASS_NAME}::Result

    my $meta = $self->SUPER::__extend_namespace__($ext);
    if ($meta) {
        return $meta;
    }

    my $new_class = $self . '::' . $ext;
    
    if ($new_class  =~ /^(.*)::Result$/) {
        my $command_class = $1;
        return $self->_generate_result_class($new_class, $command_class);
    }

    if ($new_class  =~ /^(.*)::BuildStepWrapper$/) {
        my $command_class = $1;
        return $self->_generate_build_extender_class($new_class, $command_class);
    }
    
    return;
}

sub _generate_build_extender_class {
    my ($self, $new_class, $command_class) = @_;

    my $new_class_base = 'Command::V2';
    my %has = $self->_wrap_command_class($command_class,$new_class_base);
    
    Sub::Install::install_sub({
        into => $new_class,
        as => 'execute',
        code => sub {
            my $self = shift;

            my $build;
            my $label;

            my $command;
            my $result;
            my $called_as_class_method;
            if (ref($self)) {
                $called_as_class_method = 0;
                $build = $self->wrapper_build;
                $label = $self->wrapper_build_label;
                unless ($build) {
                    die "no build on $self!";
                }
                my %params = map { 
                        if ($has{$_}{is_many}) {
                            $_ => [$self->$_]
                        }
                        else {
                            my $value = $self->$_;
                            $_ => $value
                        }
                    } keys %has;
                $command = $command_class->create(%params);
                $result = $command->execute(@_);
                $result = $command->result();
            }
            else {
                $called_as_class_method = 1;
                my $bx = $command_class->define_boolexpr(@_);
                $build = $bx->value_for("wrapper_build");
                $label = $bx->value_for("wrapper_build_label");
                unless ($build) {
                    die "no add_to_build in $bx!";
                }
                unless ($label) {
                    die "no add_to_label in $bx!";
                }
                $bx = $bx->remove_filter('wrapper_build')->remove_filter('wrapper_build_label');
                $command = $command_class->execute($bx);
                $result = $command->result;
            }

            $result->add_user(label => $label, user => $build);

            if ($called_as_class_method) {
                return $command;
            }
            else {
                return $result;
            }
        }
    });

    class {$new_class} {
        is => $new_class_base,
        has => [ 
            %has,
            wrapper_build => { is => 'Genome::Model::Build' },
            wrapper_build_label => { is => 'Text' },
        ],
        doc => "wrap $command_class in build workflows",
    };
}

sub _generate_result_class {
    my ($self, $result_class, $command_class) = @_;
    unless ($command_class->isa("Command")) {
        # the prefix of the class name is not a command class
        return;
    }
    unless ($command_class->isa("Command::V2")) {
        $self->warning_message("cannot autogenerate results for $command_class since it is not ::V2");
        return;
    }
    if ($command_class->isa("Command::Tree")) {
        # no results for "trees" :)
        return;
    }

    unless (UR::Object::Type->get("Genome::SoftwareResult::Default")) {
        $self->warnings_message("Cannot find a working Genome::SoftwareResult::Default.  Autogenerating $result_class only works with a full GMS install.");
        return;
    }

    my $new_class_base = 'Genome::SoftwareResult::Default';
    my %has = $self->_wrap_command_class($command_class,$new_class_base);

    class {$result_class} {
        is => $new_class_base,
        has => [ %has ],
        doc => "results for $command_class",
    };
}


sub _wrap_command_class {
    my ($self, $command_class, $new_class_base) = @_;

    my $command_meta = $command_class->__meta__;
    my @properties = $command_meta->properties();
    
    my %has;
    for my $property (@properties) {
        my %desc;
        next unless $property->can("is_param") and $property->can("is_input") and $property->can("is_output");

        my $name = $property->property_name;

        next if $new_class_base->can($name);

        if ($property->is_param) {
            $desc{is_param} = 1;
        }
        elsif ($property->is_input) {
            $desc{is_input} = 1;
        }
        #elsif ($property->can("is_metric") and $property->is_metric) {
        #    $desc{is_metric} = 1;
        #}
        #elsif ($property->can("is_output") and $property->is_output) {
        #    $desc{is_output} = 1;
        #}
        else {
            next;
        }

        $has{$name} = \%desc;
        $desc{is} = $property->data_type;
        $desc{doc} = $property->doc;
        $desc{is_many} = $property->is_many;
        $desc{is_optional} = $property->is_optional;
    }

    return %has;
}

1;

=pod

=head1 NAME

Genome - pipelines, tools, and data managment for genomics

=head1 SYNOPSIS

 use Genome;

 # modules in the Genome namespace will now dynamically load

 @i = Genome::InstrumentData::Illumina->get(...);
 $m = Genome::Model::SomaticVariation->create(...);

=head1 DESCRIPTION

This is the base namespace module for the Genome software tree.

That tree has several primary components:

 Genome::Model:         a data modeling pipeline management system for genomics

 Genome::Model::Tools   a tree of >1000 tools and tool wrappers for genomics

 Genome::*              a variety of sample tracking classes with an RDBMS back-end

Only the tools system is currently released.

See B<genome> for a complete inventory of all tool packages, and for command-line access to
those tools.

=head1 AUTHORS

 This software is developed by the analysis and engineering teams at
 The Genome Center at Washington Univiersity in St. Louis, with funding from
 the National Human Genome Research Institute.  Richard K. Wilson, P.I.

 Scott Abbott
 Travis Abbott
 Edward Belter
 Paul Bender
 Anthony Brummett
 Mark Burnett
 Todd C. Carter
 Matthew Callaway
 C.J. Carey
 Lynn Carmichael
 Ken Chen
 Lei Chen
 Eric Clark
 Adam Coffman
 Kevin Crouse
 Indraniel Das
 Nathan Dees
 Eric deMello
 Brian Derickson
 Alice Diec
 David Dooling
 Feiyu Du
 Adam Dukes
 James Eldred
 Xian Fan
 Ian Ferguson
 Malachi Griffith
 Obi Griffith
 Chris Harris
 Amy Hawkins
 Todd Hepler
 Xin Hong
 Shunfang Hou
 Jasreet Hundal
 Erik Hvatum
 Mark Johnson
 Krisha-Latha Kanchi
 Cyriac Kandoth
 Phil Kimmey
 Michael Kiwala
 Daniel Koboldt
 James Koval
 Karthik Kota
 Kim Kyung
 David Larson
 Sai Lek
 Shawn Leonard
 Shin Leong
 Ling Lin
 Justin Lolofie
 Robert Long
 Charles Lu
 Chris Maher
 John Martin
 Josh McMichael
 Rick Meyer
 Thomas Mooney
 David Morton
 William Nash
 Nathan Nutter
 Ben Oberkfell
 John Osborne
 Josh Peck
 Jerome Peirick
 Craig Pohl
 Allison Regier
 Ryan Richt
 Noorus Sahar Abubucker
 Gabriel Sanderson
 William Schierding
 Jon Schindler
 William Schroeder
 Christopher Schuster
 Xiaoqi Shi
 Scott Smith
 Gary Stiehr
 Sasi Suruliraj
 Kenneth Swanson
 Jason Walker
 John Wallis
 Jim Weible
 Mike Wendl
 Todd Wylie

=head1 LICENSE

Copyright (C) 2007-2012 Washington University in St. Louis.

It is released under the Lesser GNU Public License (LGPL) version 3.  See the
associated LICENSE file in this distribution.

=head1 BUGS

For defects with any software in the genome namespace,
contact gmt ~at~ genome.wustl.edu.

=cut

