package Genome::Site::TGI;

use strict;
use warnings;

BEGIN {
    require Genome::Config;
};

# do this first so we get usage metrics even if something crashes below
use Genome::Site::TGI::UsageLog;

use File::Spec qw();
my $plugins_dir;
BEGIN {
    my $file = __FILE__;
    my $tgi_dir = ($file =~ /(.*)\.pm$/)[0];
    $plugins_dir = File::Spec->join($tgi_dir, 'SiteLib');
};
use lib $plugins_dir;

BEGIN {
    if ($ENV{UR_DBI_NO_COMMIT}) {
        Genome::Config::set_env('db_pause', '');
    }

    my $sys_id = Genome::Config::get('sys_id');
    if (!$sys_id || $sys_id ne 'GMS1') {
        die q(ERROR: At TGI we expect the sys_id to be 'GMS1'.  Other sites should not use the Genome::Site::TGI module.)
    }
};

# this conflicts with all sorts of Finishing/Finfo stuff
# ironicall it is used by Pcap stuff
BEGIN { $INC{"UNIVERSAL/can.pm"} = 'no' };
BEGIN { $INC{"UNIVERSAL/isa.pm"} = 'no' };

BEGIN { $INC{ "UR/Time.pm"} = "no" };

# this keeps available parts of the UR pre-0.01 API we still use
use UR::ObjectV001removed;

# we removed UR::Time, but lots of things still depend on it
# this brings back UR::Time as a namespace, but only for legacy things
use Genome::Site::TGI::LegacyTime;

# bring in the regular Genome::Sys, then extend
use Genome::Sys;
use Genome::Site::TGI::Extension::Sys;      # extensions to Genome::Sys

use DBIx::RetryConnect Pg => sub { return {total_delay => 90} };

BEGIN {
    unless ($ENV{UR_DBI_NO_COMMIT}) {
        require Genome::Site::TGI::Extension::Logger;
        Genome::Site::TGI::Extension::Logger->import();
    }
};

use Genome::Sys::Lock;
use Genome::Site::TGI::Extension::Sys::Lock;

# configuration for internal WUGC network software & LIMS
# this module is called by Genome::Config::edu::wustl::gsc right now on all *.gsc.wustl.edu hosts
# print STDERR "using " . __PACKAGE__ . "\n";

# ensure we can get to legacy modules
use Class::Autouse;
Class::Autouse->autouse(qr/Bio.*/);

# Loads site-specific observers
use Genome::Site::TGI::Observers;

# A white-list of GSC modules which can be used on the /usr/bin/perl interpreter
my @lims_whitelist = (
    'GSC::Clone' => [
        ['Genome/Site/TGI/objects-load.t', 13],
        ['Genome/Site/TGI/use-gscapp-in-modules-fails.t', 17],
    ],
    'App::DB::TableRow::Iterator' => [
        ['Genome/Site/TGI/Synchronize/Classes/SangerRun.pm', 58],
    ],
    'GSC::Setup::CaptureSet' => [
        ['Genome/Site/TGI/CaptureSet.pm', 110],
    ],
);

# re-structured for use in the autoloader callback:
# $lims_whitelist{$class}{$file}{$line} = 1
our %lims_whitelist;
while (@lims_whitelist) {
    my $class = shift @lims_whitelist;
    my $locs = shift @lims_whitelist;
    for my $loc (@$locs) {
        my ($file,$line) = @$loc;
        $lims_whitelist{$class}{$file}{$line} = 1;
    }
}


# this callback will load the GSCApp module, and initialize the app to work with GSC modules
my $initialized = '';
our $checks = 0;
use Data::Dumper;
my $callback = sub {
    my ($pkg, $method, $class) = @_;
    $checks++;
    # print "ck @_: $initialized\n";

    return if $initialized eq 'complete' or $initialized eq 'in progress';
    return unless substr($pkg,0,5) eq 'GSC::' or substr($pkg,0,5) eq 'App::';

    if ($^X !~ /\/gsc\/bin/) {  # when not using the LIMS interpreter, only certain GSC:: modules can be used...
        my ($calling_pkg,$calling_file,$calling_line) = caller(1);
        if ($calling_file !~ /^\// and $calling_file ne '-e') {
            $calling_file =~ s/^\.\///;
            $calling_file = $ENV{PWD} . '/' . $calling_file;
        }
        $calling_file =~ s/^.*?(Genome\/.*)/$1/;
        if ($calling_file eq '-e') {
            warn "LIMS package $pkg from -e: importing GSCApp...\n";
        }
        elsif ($lims_whitelist{$pkg}{$calling_file}{$calling_line}) {
            warn "LIMS package $pkg from ($calling_file : $calling_line) is on the whitelist, bringing in GSCApp ...\n";
        }
        else {
            Carp::confess "LIMS package $pkg is not on the whitelist for file $calling_file, line $calling_line.  Update the whitelist near " . __FILE__ . " line " . __LINE__ . "!";
        }
    }

    # load and initialize GSCApp the first time something GSC:: or App:: is used.
    # since App::Init configures its own dynamic loader we dont' do anything
    # afterward, but we do need to wrap its configuration the first time to prevent conflicts

    warn "using internal LIMS modules...";

    if ($GSCApp::{BEGIN}) {
        # We've already done "use GSCApp" somewhere, and it was not
        # done before "use Genome".  Just bail.
        $initialized = 'error';
        Carp::confess("Some code in the Genome tree has a 'use GSCApp' in it.  Please remove this.");
    }

    if ($initialized eq 'error') {
        # the above happened earlier, and apparently the app did not exit
        Carp::confess("Cannot work with $pkg.  Some code in the Genome tree has a 'use GSCApp' in it.  Please remove this.");
    }

    $initialized = 'in progress';

    # remove the placeholder from above so we can actually load this module
    delete $INC{"App/Init.pm"};
    delete $INC{"GSCApp.pm"};

    require GSCApp;
    GSCApp->import();

    # ensure our access to the GSC schema is rw, and that our special env variables match up
    unless (App::Init->initialized) {
        App::DB->db_access_level('rw');
    }
    _sync_env();

    # GSCApp removes our overrides to can/isa for Class::Autoloader.  Tell it to put them back.
    App::Init->_restore_isa_can_hooks();

    # call the init process to prepare the object cache if it needs being created.
    App->init;

    $initialized = 'complete';

    return $class->can($method);
};

if ($GSCApp::{BEGIN}) {
    # GSCApp is was used first.

    # configure Genome & UR to follow its configuration.
    _sync_env();

    # GSCApp removes our overrides to can/isa for Class::Autoloader.  Tell it to put them back.
    App::Init->_restore_isa_can_hooks();
}
else {
    # No code has used GSCApp yet.
    Class::Autouse->sugar($callback);

    # The following ensures that, if someone uses GSCApp directly later, instead
    # of using the GSC classes directly, the callback will catch it and raise an error.
    # Since App::Init messes with UNIVERSAL::{can,isa} directly we need to
    # wrap the actual use of this module and restore those methods.
    $INC{"App/Init.pm"} ||= 'virtual';
    $INC{"GSCApp.pm"} ||= 'virtual';
}

sub _sync_env {
    if (App::DB::TableRow->use_dummy_autogenerated_ids || UR::DataSource->use_dummy_autogenerated_ids) {
        unless (App::Init->initialized) {
            App::DB::TableRow->use_dummy_autogenerated_ids(1);
        }
        UR::DataSource->use_dummy_autogenerated_ids(1);
    }
    if (App::DBI->no_commit || UR::DBI->no_commit) {
        unless (App::Init->initialized) {
            App::DBI->no_commit(1);
        }
        UR::DBI->no_commit(1);
    }
}


1;

=pod

=head1 NAME

Genome::Site::TGI - internal configuration for The Genome Institute at Washington University

=head1 DESCRIPTION

Configures the Genome Modeling system to work on the internal network at
The Genome Institute at Washington University

This module ensures that GSCApp and related modules are avialable to the running application.

It is currently a goal that GSCApp need not be used by this module, and that individual
modules under it provide transparent wrappers for WUTGI-specific infrastructure.

=head1 BUGS

For defects with any software in the genome namespace,
contact software@genome.wustl.edu.

=head1 SEE ALSO

B<Genome>, B<Genome::Config>, B<Genome::Site>

=cut



