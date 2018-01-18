#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

# FIXME Tests to cover:
# Allocation - all allocation in builds are not tested
# Reports - limited report testing

use strict;
use warnings;

use above 'Genome';
use Genome::Test::Factory::Build;
require Genome::Test::Factory::DiskAllocation;
use Genome::Test::Factory::Model::ReferenceAlignment;

use Data::Dumper 'Dumper';
use Test::More;

# Create temporary test subclasses
class Genome::ProcessingProfile::Test {
    is => 'Genome::ProcessingProfile',
};

class Genome::InstrumentData::Test {
    is => 'Genome::InstrumentData',
};

class Genome::Model::Test {
    is => 'Genome::ModelDeprecated',
};
sub Genome::Model::Test::_execute_build { return 1 };
sub Genome::Model::Test::files_ignored_by_build_diff { return 'meh'; }

class Genome::Model::Build::Test {
    is => 'Genome::Model::Build',
    has_optional => [
        metric1 => { is_metric => 1, },
   ],
};
ok(Genome::Model::Build::Test->does('Genome::Role::ObjectWithAllocations'), 'Test build class does ObjectWithAllocations');

my $build_meta = Genome::Model::Build::Test->__meta__;
ok($build_meta, 'build meta') or die;
my $metric1_property = $build_meta->property_meta_for_name('metric1');
my %expected_metric_property_names = (
    via => 'metrics',
    where => [ name => 'metric1', ],
    to => 'value',
    is_optional => 1,
    is_delegated => 1,
    is_mutable => 1,
);
for my $name ( keys %expected_metric_property_names ) {
    is_deeply($metric1_property->{$name}, $expected_metric_property_names{$name}, "metric property $name is correct");
}

# Create sample, library, instrument data, processing profile
my $sample = Genome::Sample->create(
    name => 'test sample'
);

my $library = Genome::Library->create(
    sample_id => $sample->id,
    name => 'test library',
);

my $inst_data = Genome::InstrumentData::Test->create(
    sequencing_platform => 'test',
    library_id => $library->id,
);

my $pp = Genome::ProcessingProfile::Test->create(
    name => 'test pp',
);

# Create test model and make sure inputs and such are correct
my $model = Genome::Model::Test->create(
    subject_id => $sample->id,
    subject_class_name => $sample->class,
    processing_profile_id => $pp->id,
    name => 'test model',
);
ok($model, 'created test model') or die;
isa_ok($model, 'Genome::Model::Test', 'model subclass automagically generated');
$model->add_instrument_data(value => $inst_data);
my @model_inst_data = $model->instrument_data;
ok(@model_inst_data, 'Added instrument data to model');
my @model_inputs = $model->inputs;
is(scalar(@model_inputs), 1, 'Correct number of model inputs');

# Create test build
my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $build = Genome::Model::Build->create(
    model_id => $model->id,
    data_directory => $tmpdir,
);
ok($build, 'Created build') or die;
isa_ok($build, 'Genome::Model::Build::Test', 'build subclass automagically generated');
is(Genome::Model::Build::Test->model_class, $model->class, 'build model_class');
ok($build->the_master_event, 'master event created');
is($build->model->id, $model->id, 'indirect model accessor');
is($build->status, 'New', 'status set to New');

# Make sure inputs were copied to build
my @build_inputs = $build->inputs;
is(scalar(@build_inputs), 1, 'Correct number of build inputs');
my @build_inst_data = $build->instrument_data;
is($build->instrument_data_count, 1, 'Instrument data count');
is_deeply(\@build_inst_data, \@model_inst_data, 'Build instrument data matches model instrument data');

# TODO Check that model is updated via check_for_updates and that _initialize_build is called on pp
my @validate_for_start_tags = ( UR::Object::Tag->create(properties => ['foo'], desc => 'test tag') );
no warnings;
*Genome::Model::Build::Test::validate_for_start = sub { return @validate_for_start_tags };
*Genome::Report::Email::send_report = sub{ 
    my ($class, %params) = @_;
    _test_expected_report_params(\%params); 
    return 1;
};
use warnings;

my $disk_allocation = Genome::Disk::Allocation->create(
    kilobytes_requested => 1,
    owner_id => $build->id,
    owner_class_name => $build->class,
    disk_group_name => Genome::Config::get('disk_group_models'),
    allocation_path => 'temp/Model/Build.t_test/' . $build->id,
);
is($build->disk_allocation, $disk_allocation, 'disk_allocation');
is_deeply([$build->disk_allocations], [$disk_allocation], 'disk_allocations');

ok(!$build->start, 'build failed to start, as expected');
my @notes = Genome::MiscNote->get(
    subject => $build,
    header_text => 'Unstartable',
);
ok(@notes, 'got unstartable notes from build');
ok(scalar(@notes) == 1, 'got exactly one unstartable note back');
ok($notes[-1]->body_text =~ /could not be validated for start/, 'error message matches expected');
is($build->status, 'Unstartable', 'build status set to unstartable');

# Now remove the build validation fail tags so we can start ok
@validate_for_start_tags = ();
$model->build_requested(1);
my $guard = Genome::Config::set_env('workflow_builder_backend', 'inline');
ok($build->start(), 'build started!');
is($build->status, 'Succeeded', 'build completed successfully');
ok(!$model->current_running_build_id, 'Current running build id set to undef in success');
ok(!$model->build_needed, 'This succeeded build satisfies the model');
ok($build->data_directory, 'data directory resolved');
ok($build->software_revision, 'software revision set on build');
is($model->build_requested, 0, 'build requested unset on model');
undef($guard);

# FAIL
is($build->the_master_event->user_name('apipe-builder'), 'apipe-builder', 'set user name to apipe-builder');
isnt($disk_allocation->kilobytes_requested, 2000, 'disk_allocation reallocated');
_test_fail($build);
is($build->status, 'Failed', 'Status is Failed');

# ABANDON
$disk_allocation->kilobytes_requested(2000); # reset
ok($build->abandon, 'Abandon');
is($build->status, 'Abandoned', 'Status is Abandoned');
isnt($model->last_complete_build_id, $build->id, 'Model last complete build is not this build\'s id in abandon');

UR::Context->commit(); # triggers the allocation purge
is($disk_allocation->status, 'purged', 'disk_allocation purged');

# Init, fail and succeed a abandoned build
ok(!$build->initialize, 'Failed to initialize an abandoned build');
ok(!$build->fail, 'Failed to fail an abandoned build');
ok(!$build->success, 'Failed to success an abandoned build');

# DIFF
my @files_ignored_by_build_diff = $build->files_ignored_by_diff;

is_deeply(\@files_ignored_by_build_diff, ['meh'], 'files_ignored_by_diff');
my $build2 = Genome::Model::Build::Test->create(
    model => $model,
    data_directory => $tmpdir, # TODO actually test file diffs?
);
$build->add_metric(name => 'one', value => 1);
$build->add_metric(name => 'two', value => 2);
$build2->add_metric(name => 'two', value => 'not 2');
$build->add_metric(name => 'three', value => 3);
$build2->add_metric(name => 'three', value => 3);
$build2->add_metric(name => 'four', value => 4);
my %diffs = $build->compare_output($build2->id);
my %expected_diffs = (
    'one' => qr|no build metric with name one found for build \-?[[:xdigit:]]+|,
    'two' => qr|metric two has value 2 for build \-?[[:xdigit:]]+ and value not 2 for build \-?[[:xdigit:]]+|,
    'four' => qr|no build metric with name four found for build \-?[[:xdigit:]]+|,
);
for my $diff ( keys %expected_diffs ) {
    like($diffs{$diff}, $expected_diffs{$diff}, "diff message for $diff is correct");
}

# Delete build, verify disk allocation deleted
ok($build->delete, 'Deleted build');
UR::Context->commit; # this fails
isa_ok($disk_allocation, 'UR::DeletedRef');

test_diff_vcf();

done_testing();

####

sub _test_fail {
    #Ensure that when builds fail, they add notes, equivlent to the failure email's stage, step, and error
    my $build = shift;

    my $error = Genome::Model::Build::Error->create(
        build_event_id => $build->build_event->id,
        stage_event_id => $build->build_event->id,
        stage => 'all stages',
        step_event_id => $build->build_event->id,
        step => 'main',
        error => 'Testing error',
    );
    isa_ok($error, 'Genome::Model::Build::Error', 'error');

    $build->fail($error);

    my $failed_step_note = $build->notes(header_text => 'Failed Step');
    isa_ok($failed_step_note, 'Genome::MiscNote', 'failed_step_note');
    is($failed_step_note->body_text, $error->step, 'failed step matches error step');

    my $failed_stage_note = $build->notes(header_text => 'Failed Stage');
    isa_ok($failed_stage_note, 'Genome::MiscNote', 'failed_stage_note');
    is($failed_stage_note->body_text, $error->stage, 'failed stage matches error stage');

    my $failed_error_note = $build->notes(header_text => 'Failed Error');
    isa_ok($failed_error_note, 'Genome::MiscNote', 'failed_error_note');
    is($failed_error_note->body_text, $error->error, 'failed error matches error error');

    # delete fail notes
    for ( $failed_step_note, $failed_error_note, $failed_stage_note) { $_->delete; }

    return 1;
}

sub test_diff_vcf {
    my $test_data_dir = Genome::Config::get('test_inputs') . '/Genome-Model-Build';
    my $input_dir = join('/', $test_data_dir, 'input');

    my $control_file = join('/', $input_dir, 'indels.vcf');
    ok(-s $control_file, 'control_file has size') || return;

    my $similar_file = join('/', $input_dir, 'indels_similar.vcf');
    ok(-s $similar_file, 'similar_file has size') || return;
    is(Genome::Model::Build::diff_vcf($control_file, $similar_file), 1, 'similar_file matches control_file');

    my $different_file = join('/', $input_dir, 'indels_modified.vcf');
    ok(-s $different_file, 'different_file has size') || return;
    ok(!Genome::Model::Build::diff_vcf($control_file, $different_file), 'different_file does not match control_file');
}

sub _test_expected_report_params {
    my ($report_params) = @_;
    my $user = Genome::Sys::User->get(username => $build->the_master_event->user_name);
    my %expected_params = (
        to => $user->email,
        from => Genome::Config::get('email_pipeline'),
        replyto => Genome::Config::get('email_noreply'),
    );
    my %got_params = map { $_ => $report_params->{$_} } keys %expected_params; 
    die 'No report params!' if not %got_params;
    is_deeply(
        \%got_params,
        \%expected_params,
        "Sent report params",
    );
    return 1;
}
