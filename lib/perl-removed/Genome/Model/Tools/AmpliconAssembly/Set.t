#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::MockObject;
use Test::More;
require File::Temp;
require File::Path;

use_ok('Genome::Model::Tools::AmpliconAssembly::Set') or die;

my $base_test_dir = $ENV{GENOME_TEST_INPUTS} . '';
my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
    
my $amplicon_assembly = Genome::Model::Tools::AmpliconAssembly::Set->get(
    directory => $base_test_dir.'/Genome-Model/AmpliconAssembly/build',
);
ok($amplicon_assembly, 'get amplicon assembly');

my $create_amplicon_assembly = Genome::Model::Tools::AmpliconAssembly::Set->create(
    directory => $tmp_dir,
);
ok($create_amplicon_assembly, 'create');
unlink $create_amplicon_assembly->_properties_file;
$create_amplicon_assembly->delete;

my %invalid_params = (
    sequencing_center => 'washu',
    sequencing_platform => '373',
);
for my $invalid_attr ( keys %invalid_params ) {
    ok(!Genome::Model::Tools::AmpliconAssembly::Set->create(
            directory => $tmp_dir,
            $invalid_attr => $invalid_params{$invalid_attr},
        ),
        "failed as expected - create w/ $invalid_attr\: ".$invalid_params{$invalid_attr},
    );
}

# amplicons
my $amplicons = $amplicon_assembly->get_amplicons;
is_deeply(
    [ map { $_->name } @$amplicons ],
    [qw/ HMPB-aad13a05 HMPB-aad13e12 HMPB-aad15e03 HMPB-aad16a01 HMPB-aad16c10 /],
    'Got 5 amplicons',
);
# reads for amplicon
my @reads = $amplicon_assembly->get_all_amplicons_reads_for_read_name(
    ($amplicons->[0]->reads)[0],
);
is_deeply(\@reads, $amplicons->[0]->get_reads, 'Got all amplicons reads for read name');

done_testing();
exit;

