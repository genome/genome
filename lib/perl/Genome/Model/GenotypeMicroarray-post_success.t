use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::InstrumentData::Imported;
use Genome::Test::Factory::Library;
use Genome::Test::Factory::Sample;

my $genotype_sample = Genome::Test::Factory::Sample->setup_object(
);

my $genotype_library => Genome::Test::Factory::Library->setup_object(
   sample_id => $genotype_sample, 
);

my $genotype_data = Genome::Test::Factory::InstrumentData::Imported->setup_object(
    library => $genotype_library,
);

my $qc_sample = Genome::Test::Factory::Sample->setup_object(
   extraction_type => 'genomic_dna',
   default_genotype_data_id => , $genotype_data->id,
);

my $qc_library = Genome::Test::Factory::Library->setup_object(
   sample_id => $qc_sample, 
);

my $qc_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object(
    library => $qc_library,
);


#make qc_model
my $qc_model = Genome::Model::ReferenceAlignment->create(
    processing_profile => Genome::ProcessingProfile->get('2653572'),
    reference_sequence_build_id => '106942997',
    auto_assign_inst_data => '0',
    build_requested => '0',
    roi_track_name => 'tiled_region',
    subject => $qc_sample,
    instrument_data => [$qc_data],
);

#make genotype microarray model
my $genotype_microarray_model = Genome::Model::GenotypeMicroarray->create(
    processing_profile_id => '2166945',
    # reference_sequence_build_id => '106942997',
    dbsnp_build_id => 127786607,
    subject => $genotype_sample,
    instrument_data => [$genotype_data],
);


#call ->success on genotype_model
my $tmp_dir = Genome::Sys->create_temp_directory;
my $gmb = Genome::Model::Build->create(model_id => $genotype_microarray_model->id, data_directory => $tmp_dir );
$gmb->success();

#check $qc_model->build_requested
ok($qc_model->build_requested, 'QC model has build_requested');

done_testing();
