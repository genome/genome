package Genome;

use warnings;
use strict;

our $VERSION = '0.080100';

require Genome::Config;

use UR;
use File::Temp;
use IO::String;
use File::Basename;
use Carp;
use Carp::Heavy;

require Sys::Hostname;

# Standard namespace declaration for a UR namespace
UR::Object::Type->define(
    class_name => 'Genome',
    is => ['UR::Namespace'],
    english_name => 'genome',
);

# There is no "dev mode" anymore, but continue to produce a warning if no-commit is on
my $dev_mode = UR::DBI->no_commit;
if ($dev_mode) {
    my $h = Sys::Hostname::hostname;
    warn "***** GENOME_DEV_MODE ($h) *****";
}

sub execution_id {
    unless ($ENV{GENOME_EXECUTION_ID}) {
        $ENV{GENOME_EXECUTION_ID} = UR::Object::Type->autogenerate_new_object_id_uuid();
    }
    return $ENV{GENOME_EXECUTION_ID};
}

require Genome::Site;

if (my $umask = Genome::Config::get('sys_umask')) {
    my $old_umask = umask oct($umask);
    if (!defined($old_umask)) {
        die 'failed to set umask';
    }
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

