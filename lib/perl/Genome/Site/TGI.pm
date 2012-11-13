package Genome::Site::TGI;
use strict;
use warnings;

BEGIN {
    if ($ENV{GENOME_DEV_MODE}) {
        $ENV{GENOME_SYS_SERVICES_MEMCACHE} ||= 'apipe-dev.gsc.wustl.edu:11211';
        $ENV{GENOME_SYS_SERVICES_SOLR} ||= 'http://solr-dev:8080/solr';
    }
    else {
        $ENV{GENOME_SYS_SERVICES_MEMCACHE} ||= 'imp-apipe.gsc.wustl.edu:11211';
        $ENV{GENOME_SYS_SERVICES_SOLR} ||= 'http://solr:8080/solr';
    }
	
}

BEGIN {
    if (!defined $ENV{GENOME_DB_SKIP_POSTGRES}) {
        $ENV{GENOME_DB_SKIP_POSTGRES} ||= '/gsc/scripts/opt/genome/run/skip_postgres_sync_new';
    }
}

BEGIN {
	if (defined $ENV{GENOME_QUERY_POSTGRES}) {
		no warnings;
		use UR::Context;
		use UR::DataSource::Pg;
		use Workflow;
		*UR::Context::resolve_data_sources_for_class_meta_and_rule_genome_filtered = \&UR::Context::resolve_data_sources_for_class_meta_and_rule;;
		*UR::Context::resolve_data_sources_for_class_meta_and_rule = \&Genome::Site::TGI::resolve_data_sources_for_class_meta_and_rule;
	
		*UR::Object::Type::table_name_filtered = \&UR::Object::Type::table_name;
		*UR::Object::Type::table_name = \&Genome::Site::TGI::table_name_patch;
		use warnings;
	}
}

sub undo_table_name_patch {
    no warnings;
    *UR::Object::Type::table_name = \&UR::Object::Type::table_name_filtered;
    use warnings;
}

sub redo_table_name_patch {
    no warnings;
    *UR::Object::Type::table_name = \&Genome::Site::TGI::table_name_patch;
    use warnings;
}

sub table_name_patch {
        my $self = shift;
        my $is_generating_id = ((caller(1))[3] eq 'UR::DataSource::RDBMS::autogenerate_new_object_id_for_class_name_and_rule');

        if (@_ || $is_generating_id) {
                return $self->table_name_filtered(@_);
        } else {
                my $table_name = $self->table_name_filtered;
                my $mapped_table_name = Genome::DataSource::Main->postgres_table_name_for_oracle_table(lc($table_name));
#    || Workflow::DataSource::InstanceSchemaPostgres->postgres_table_name_for_oracle_table(lc($table_name));
        
                if ($mapped_table_name) {
                        $table_name =$mapped_table_name;
                } 
                return $table_name;
        }

}

sub resolve_data_sources_for_class_meta_and_rule {
        my $context = shift; 
        my $data_source = $context->resolve_data_sources_for_class_meta_and_rule_genome_filtered(@_);

        if ($data_source) {
                my $caller = (caller(1))[3];
                return $data_source if ($caller eq 'UR::Object::Type::autogenerate_new_object_id');
                if ($data_source->isa('Genome::DataSource::GMSchema')) {
                        $data_source = Genome::DataSource::PGTest->get();
#                } elsif ($data_source->isa('Workflow::DataSource::InstanceSchema')) {
#                        $data_source = Workflow::DataSource::InstanceSchemaPostgres->get();
                }
        }
                
        return $data_source;
}

# configure local statsd server
BEGIN {
    unless ($ENV{UR_DBI_NO_COMMIT}) {
        $ENV{GENOME_STATSD_HOST} ||= 'apipe-statsd.gsc.wustl.edu';
        $ENV{GENOME_STATSD_PORT} ||= 8125;
    }
};

# this conflicts with all sorts of Finishing/Finfo stuff
# ironicall it is used by Pcap stuff
BEGIN { $INC{"UNIVERSAL/can.pm"} = 'no' };
BEGIN { $INC{"UNIVERSAL/isa.pm"} = 'no' };

# ensure nothing loads the old Genome::Config module
BEGIN { $INC{"Genome/Config.pm"} = 'no' };

BEGIN { $INC{ "UR/Time.pm"} = "no" };

# set our internal paths for data and software
$ENV{GENOME_DB} ||= '/gsc/scripts/opt/genome/db';
$ENV{GENOME_SW} ||= '/gsc/pkg/bio';
$ENV{GENOME_LOCK_DIR} ||= '/gsc/var/lock';

# testsuite data
$ENV{GENOME_TEST_INPUTS} ||= '/gsc/var/cache/testsuite/data';
$ENV{GENOME_TEST_TEMP} ||= '/gsc/var/cache/testsuite/running_testsuites';
$ENV{GENOME_TEST_URL} ||= 'https://gscweb.gsc.wustl.edu/gscmnt/gc4096/info/test_suite_data/';

# configure file that signals that database updates should be paused
if (!$ENV{UR_DBI_NO_COMMIT}) {
    $ENV{GENOME_DB_PAUSE} ||= $ENV{GENOME_LOCK_DIR} . '/database/pause_updates';
}

# configure our local ensembl db
$ENV{GENOME_DB_ENSEMBL_DEFAULT_IMPORTED_ANNOTATION_BUILD} ||= '122704720';
$ENV{GENOME_DB_ENSEMBL_HOST} ||= 'mysql1';
$ENV{GENOME_DB_ENSEMBL_USER} ||= 'mse';
$ENV{GENOME_DB_ENSEMBL_PORT} ||= '3306';

# default nomenclature for new instrument data and sample attributes
$ENV{GENOME_NOMENCLATURE_DEFAULT} ||= 'WUGC';

# Log directory
$ENV{GENOME_LOG_DIR} ||= '/gsc/var/log/genome';

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



