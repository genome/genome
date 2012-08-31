#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;

use_ok('Genome::Model::Tools::Soap::DaccDownload') or die;

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Soap-DaccDownload';
ok(-d $data_dir, "Data dir exists");

my %available_files;
no warnings;
*Genome::Model::Tools::Dacc::Download::available_files_and_sizes = sub{ return %available_files; };
*Genome::Model::Tools::Dacc::Download::execute = sub{ return 1; };
*File::Copy::move = sub{ return 1; };
use warnings;

note('fail - no files');
my $import = Genome::Model::Tools::Soap::DaccDownload->create(
    version => 'dacc',
    import_location => '/WholeMetagenomic/03-Assembly/PGA/SRS078176_LANL/',
    output_dir_and_file_prefix => "$data_dir/SRS000000_LANL",
);
ok($import, 'create');
$import->dump_status_messages(1);
ok(!$import->execute, 'execute');

note('fail - some files missing');
$available_files{'SRS078176_LANL.agp'} = 6793140;
$import = Genome::Model::Tools::Soap::DaccDownload->create(
    version => 'dacc',
    import_location => '/WholeMetagenomic/03-Assembly/PGA/SRS078176_LANL/',
    output_dir_and_file_prefix => "$data_dir/SRS000000_LANL",
);
ok($import, 'create');
$import->dump_status_messages(1);
ok(!$import->execute, 'execute');

note('fail - extra scafSeq file and not a PGA');
%available_files = (
    'SRS078176_LANL.readInGap' => 0,
    'SRS078176_PGA.scaffolds.fa' => 0,
    'SRS078176_LANL.scaffolds.fa' => 0,
    'SRS078176_PGA.contigs.fa' => 0,
    'SRS078176_LANL.contig' => 0,
    'SRS078176_LANL.readOnContig' => 0,
    'SRS078176_LANL.contigs.fa' => 0,
    'SRS078176_LANL.agp' => 0,
    'SRS078176_LANL.scafSeq' => 0,
);
$available_files{'SRS078176_Baylor.agp'} = 6793140;
$import = Genome::Model::Tools::Soap::DaccDownload->create(
    version => 'dacc',
    import_location => '/WholeMetagenomic/03-Assembly/PGA/SRS078176_LANL/',
    output_dir_and_file_prefix => "$data_dir/SRS000000_LANL",
);
ok($import, 'create');
$import->dump_status_messages(1);
ok(!$import->execute, 'execute');

note('success');
delete $available_files{'SRS078176_Baylor.agp'};
$import = Genome::Model::Tools::Soap::DaccDownload->create(
    version => 'dacc',
    import_location => '/WholeMetagenomic/03-Assembly/PGA/SRS078176_LANL/',
    output_dir_and_file_prefix => "$data_dir/SRS000000_LANL",
);
ok($import, 'create');
$import->dump_status_messages(1);
ok($import->execute, 'execute');

done_testing();
exit;

