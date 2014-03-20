package Genome::Site::TGI;
use strict;
use warnings;

# do this first so we get usage metrics even if something crashes below
use Genome::Site::TGI::UsageLog;

BEGIN {
    my @shell_vars = qw(
        GENOME_DB
        GENOME_DB_QUERY_PAUSE
        GENOME_DISK_GROUP_ALIGNMENTS
        GENOME_DISK_GROUP_ARCHIVE
        GENOME_DISK_GROUP_DEV
        GENOME_DISK_GROUP_MODELS
        GENOME_DISK_GROUP_REFERENCES
        GENOME_DISK_GROUP_RESEARCH
        GENOME_DISK_GROUP_TRASH
        GENOME_DS_DGIDB_AUTH
        GENOME_DS_DGIDB_LOGIN
        GENOME_DS_DGIDB_OWNER
        GENOME_DS_DGIDB_SERVER
        GENOME_DS_DGIDB_TYPE
        GENOME_DS_DWRAC_AUTH
        GENOME_DS_DWRAC_LOGIN
        GENOME_DS_DWRAC_OWNER
        GENOME_DS_DWRAC_SERVER
        GENOME_DS_DWRAC_TYPE
        GENOME_DS_GMSCHEMA_AUTH
        GENOME_DS_GMSCHEMA_LOGIN
        GENOME_DS_GMSCHEMA_OWNER
        GENOME_DS_GMSCHEMA_SERVER
        GENOME_DS_GMSCHEMA_TYPE
        GENOME_DS_OLTP_AUTH
        GENOME_DS_OLTP_LOGIN
        GENOME_DS_OLTP_OWNER
        GENOME_DS_OLTP_SERVER
        GENOME_DS_OLTP_TYPE
        GENOME_EMAIL_ANNOTATION
        GENOME_EMAIL_DOMAIN
        GENOME_EMAIL_ILLUMINA_BWA
        GENOME_EMAIL_NOREPLY
        GENOME_EMAIL_PIPELINE
        GENOME_EMAIL_PIPELINE_NOISY
        GENOME_EMAIL_SMTP_SERVER
        GENOME_EMAIL_TEST
        GENOME_EMAIL_VIROME_SCREENING
        GENOME_FS_LOCAL_NETWORK_CACHE
        GENOME_HOME
        GENOME_LOCK_DIR
        GENOME_LSF_QUEUE_ALIGNMENT_DEFAULT
        GENOME_LSF_QUEUE_ALIGNMENT_PROD
        GENOME_LSF_QUEUE_ASSEMBLY
        GENOME_LSF_QUEUE_BIGMEM
        GENOME_LSF_QUEUE_BUILD_WORKER
        GENOME_LSF_QUEUE_BUILD_WORKER_ALT
        GENOME_LSF_QUEUE_BUILD_WORKFLOW
        GENOME_LSF_QUEUE_DV2_WORKER
        GENOME_LSF_QUEUE_DV2_WORKFLOW
        GENOME_LSF_QUEUE_SHORT
        GENOME_ROOT
        GENOME_SW
        GENOME_SW_IGNORE
        GENOME_SW_LEGACY_JAVA
        GENOME_SYS_GROUP
        GENOME_SYS_ID
        GENOME_SYS_SERVICES_FILES_URL
        GENOME_SYS_SERVICES_MEMCACHE
        GENOME_SYS_SERVICES_SEARCH_URL
        GENOME_SYS_SERVICES_SOLR
        GENOME_SYS_SERVICES_WEB_VIEW_URL
        GENOME_SYS_SERVICES_WIKI_URL
        GENOME_SYS_UMASK
        GENOME_TEST_INPUTS
        GENOME_TEST_URL
    );

    if ($ENV{GENOME_DEV_MODE}) {
        $ENV{GENOME_SYS_SERVICES_MEMCACHE} = 'apipe-dev.gsc.wustl.edu:11211';
        $ENV{GENOME_SYS_SERVICES_SOLR} = 'http://solr-dev:8080/solr';
    }

    if ($ENV{UR_DBI_NO_COMMIT}) {
        delete $ENV{GENOME_STATSD_HOST};
        delete $ENV{GENOME_STATSD_PORT};
        delete $ENV{GENOME_DB_PAUSE};
    } else {
        # these get pushed on because if Genome::Site::TGI gets called by a
        # subshell during tests then it will fail because the parent deleted
        # these environment variables
        push @shell_vars, qw(
            GENOME_DB_PAUSE
            GENOME_STATSD_HOST
            GENOME_STATSD_PORT
        );
    }

    my @unset_shell_vars = grep { !defined($ENV{$_}) } @shell_vars;
    if (@unset_shell_vars) {
        die q(ERROR: You need to start a new shell so that environment variables are setup for Genome.pm.  The following environment variables are not set: ) . join(', ', @unset_shell_vars);
    }

    if (!$ENV{GENOME_SYS_ID} || $ENV{GENOME_SYS_ID} ne 'GMS1') {
        die q(ERROR: At TGI we expect the GENOME_SYS_ID to be 'GMS1'.  Other sites should not use the Genome::Site::TGI module.)
    }
};

# this conflicts with all sorts of Finishing/Finfo stuff
# ironicall it is used by Pcap stuff
BEGIN { $INC{"UNIVERSAL/can.pm"} = 'no' };
BEGIN { $INC{"UNIVERSAL/isa.pm"} = 'no' };

# ensure nothing loads the old Genome::Config module
BEGIN { $INC{"Genome/Config.pm"} = 'no' };

BEGIN { $INC{ "UR/Time.pm"} = "no" };

# a unique ID for each program execution.  Used got logging saves to the database
$ENV{GENOME_EXECUTION_ID} = UR::Object::Type->autogenerate_new_object_id_uuid();

# If running either the genome or gmt command, log the arguments, user,
# host, etc in a log file
if ($0 =~ /(?:gmt|genome)(?:5\.12\.1)?$/ and not `grep log_command $0`) {
    require Genome::Site::TGI::Security;
    Genome::Site::TGI::Security::log_command(@ARGV);
}

# this keeps available parts of the UR pre-0.01 API we still use
use UR::ObjectV001removed;

# we removed UR::Time, but lots of things still depend on it
# this brings back UR::Time as a namespace, but only for legacy things
use Genome::Site::TGI::LegacyTime;

# bring in the regular Genome::Sys, then extend
use Genome::Sys;
use Genome::Site::TGI::Extension::Sys;      # extensions to Genome::Sys

use Genome::Site::TGI::Extension::Logger;

# the old Genome::Config is all deprecated
# the core stuff about looking up your host config is now in Genome::Site
use Genome::Site::TGI::LegacyConfig;

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
    'GSC::PSEParam' => [
        ['Genome/Model/Tools/Lims/ApipeBridge/InstrumentDataStatus.pm', 109],
    ],
    'GSC::PSE' => [
        ['Genome/Model/Tools/Lims/ApipeBridge/FixPidfaParamsForBase.pm', 89],
        ['Genome/Model/Tools/Lims/ApipeBridge/FixPidfaParamsForBase.pm', 118],
        ['Genome/Model/Tools/Lims/ApipeBridge/FixPidfaParamsForGenotype.pm', 25],
        ['Genome/Model/Tools/Lims/ApipeBridge/InstrumentDataStatus.pm', 115],
    ],
    'GSC::Sequence::Genome' => [
        ['Genome/Model/Tools/Snp/GetDbsnps.pm', 282],
    ],
    'GSC::Clone' => [
        ['Genome/Site/TGI/objects-load.t', 13],
        ['Genome/Site/TGI/use-gscapp-in-modules-fails.t', 17],
    ],
    'GSC::Genotyping::External' => [
        ['Genome/Model/Tools/Lims/ApipeBridge/UpdateGenomeGenotypes', 48],
    ],
    'GSC::Genotyping::Internal::Illumina' => [
        ['Genome/Model/Tools/Lims/ApipeBridge/UpdateGenomeGenotypes', 48],
    ],
    'App::DB::TableRow::Iterator' => [
        ['Genome/Model/Tools/Lims/ImportSangerRuns.pm', 105],
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



