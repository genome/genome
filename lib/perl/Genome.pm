package Genome;


use warnings;
use strict;

our $VERSION = '0.080100';

use UR;
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

if ($ENV{GENOME_SYS_UMASK}) {
    my $old_umask = umask oct($ENV{GENOME_SYS_UMASK});
    if (!defined($old_umask)) {
        die 'failed to set umask';
    }
}

# If the search engine is installed, configure its hooks
eval {
    local $SIG{__WARN__};
    local $SIG{__DIE__};
    require Genome::Search;
};
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

# DB::single is set to this value in many places, creating a source-embedded break-point
# set it to zero in the debugger to turn off the constant stopping...
$DB::stopper = 1;

# This is tested by Genome/SoftwareResult/Default.t
# And also by Genome/Command/WithSavedResult-buildwrapper.t
# It is only used by the full GMS install.
# They probably needs a better home than the base module.

our %NAMESPACE_EXTENSIONS = (
    Result => {
        base_class => 'Genome::SoftwareResult::Default',
        doc => 'results for CLASS',
    },
    BuildStepWrapper => {
        base_class => 'Genome::Command::BuildStepWrapper',
        doc => 'wrap CLASS as a step in a build workflow',
    },
);

sub __extend_namespace__ {
    my ($self,$ext) = @_;

    my $meta = $self->SUPER::__extend_namespace__($ext);
    if ($meta) {
        return $meta;
    }

    my $new_class = $self . '::' . $ext;
    my ($command_class,$final_ext) = ($new_class =~ /^(.*)::(.*)$/);

    my $cmd_meta = UR::Object::Type->get($command_class);
    return unless $cmd_meta;

    return unless $command_class->isa("Command");

    return unless $command_class->isa("Command::V2");

    my $spec = $NAMESPACE_EXTENSIONS{$final_ext};
    return unless $spec;

    my $base_class = $spec->{base_class};
    my $doc = $spec->{doc};
    $doc =~ s/CLASS/$new_class/;

    my %has = $self->_wrapper_has($command_class,$base_class);

    class {$new_class} {
        is => $base_class,
        has => [ %has ],
        doc => $doc,
    };
}


# TODO: switch to calling the static _wrapper_has method on Command::V2
sub _wrapper_has {
    my ($self, $command_class, $base) = @_;

    my $command_meta = $command_class->__meta__;
    my @properties = $command_meta->properties();

    my %has;
    for my $property (@properties) {
        my %desc;
        next unless $property->can("is_param") and $property->can("is_input") and $property->can("is_output");

        my $name = $property->property_name;

        next if $base->__meta__->properties(property_name => $name);

        if ($property->is_param) {
            $desc{is_param} = 1;
        }
        elsif ($property->is_input) {
            $desc{is_input} = 1;
        }
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

=head1 VERSION

This document describes Genome 0.8.1.0

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

