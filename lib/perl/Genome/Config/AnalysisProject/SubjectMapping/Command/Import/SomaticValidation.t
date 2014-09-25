#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
};

use Test::More;
use Genome::Test::Factory::Sample;

use above 'Genome';

my $number_of_mappings = 2;
my $class = 'Genome::Config::AnalysisProject::SubjectMapping::Command::Import::SomaticValidation';
use_ok($class);

my $test_file = _test_file();
my $analysis_project = Genome::Config::AnalysisProject->__define__(name => 'test proj');

my $cmd = $class->create(
    analysis_project => $analysis_project,
    file_path => $test_file,
);

isa_ok($cmd, $class);

my $res = $cmd->execute();

ok($res, 'command ran successfully');

is($res, $number_of_mappings,
    "we expected to create $number_of_mappings subject mappings and did");

my @mappings = $analysis_project->subject_mappings;
is(scalar(@mappings), $number_of_mappings,
    'we associated the correct number of pairings with the AnalysisProject');

done_testing();

sub _test_file {
    my ($fh, $path) = Genome::Sys->create_temp_file();

    for(1..$number_of_mappings) {
        if ($_ % 2 == 0) {
            $fh->printf("%s\t%s\t%s\t%s\t%s\n",
                Genome::Test::Factory::Sample->setup_object()->id,
                Genome::Test::Factory::Sample->setup_object()->id,
                Genome::Model::Tools::DetectVariants2::Result::Manual->__define__()->id,
                Genome::Model::Tools::DetectVariants2::Result::Manual->__define__()->id,
                Genome::Model::Tools::DetectVariants2::Result::Manual->__define__()->id
            );
        } else {
            $fh->printf("%s\t%s\t%s\t%s\t%s\n",
                Genome::Test::Factory::Sample->setup_object()->name,
                Genome::Test::Factory::Sample->setup_object()->name,
                Genome::Model::Tools::DetectVariants2::Result::Manual->__define__()->id,
                Genome::Model::Tools::DetectVariants2::Result::Manual->__define__()->id,
                Genome::Model::Tools::DetectVariants2::Result::Manual->__define__()->id
            );
        }
    }

    $fh->close();
    return $path;
}
