package Genome::Model::Build::MetagenomicComposition16s::TestBuildFactory;

use strict;
use warnings;

use Genome;

use File::Temp;
use Test::More;

my $id = -9000;
my %entities;

sub library {
    return $entities{library} if $entities{library};
    my $sample = sample();
    $entities{library} = Genome::Library->__define__(
        id => --$id,
        name => '__TEST_LIBRARY__',
        sample_id => $sample->id,
    );
    die 'Failed to create library!' if not $entities{library};
    return $entities{library};
}

sub sample {
    return $entities{sample} if $entities{sample};
    $entities{sample} = Genome::Sample->__define__(
        id => --$id,
        name => '__TEST_SAMPLE__',
    );
    die 'Failed to create sample!' if not $entities{sample};
    return $entities{sample};
}

sub build_with_example_build {
    my ($sequencing_platform) = @_;

    return ( $entities{build}, $entities{example_build} ) if $entities{build};

    $entities{pp} = Genome::ProcessingProfile->__define__(
        id => --$id,
        type_name => 'metagenomic composition 16s',
        amplicon_processor => 'filter by-min-length --length 200',
        chimera_detector => 'chimera-slayer',
        chimera_detector_params => "--nastier-params '-num_top_hits 10' --chimera-slayer-params '-windowSize 50 -printCSalignments -windowStep 5'",
        classifier => 'rdp2-2',
        classifier_params => '-training-set 6 -format hmp_fix_ranks -version 2x2',
        sequencing_platform => $sequencing_platform,
    ) if not $entities{pp};

    my $sample = sample();
    $entities{model} = Genome::Model::MetagenomicComposition16s->__define__(
        id => --$id,
        name => '__TEST_MC16S_MODEL__',
        processing_profile => $entities{pp},
        subject_name => $sample->name,
        subject_type => 'sample_name'
    ) if not $entities{model};
    die 'Failed to create MC16s model!' if not $entities{model};

    $entities{build} = Genome::Model::Build::MetagenomicComposition16s->__define__(
        id => --$id,
        model => $entities{model},
        data_directory => File::Temp::tempdir(CLEANUP => 1),
    );
    die 'Failed to create MC16s model!' if not $entities{build};
    $entities{build}->create_subdirectories;

    my $instrument_data = instrument_data_454();
    $entities{build}->add_instrument_data($instrument_data) or die 'Failed to add instrument data to build!';

    $entities{example_build} = Genome::Model::Build->__define__(
        model=> $entities{model},
        id => --$id,
    ) if not $entities{example_build};
    die 'Failed to create MC16s model!' if not $entities{example_build};

    my %example_data_directories = (
        454 => $ENV{GENOME_TEST_INPUTS} . '/Genome-Model/MetagenomicComposition16s454/build_v5.2chimeras', # start w/ 2 chimeras
    );
    $entities{example_build}->data_directory( $example_data_directories{$sequencing_platform} ) or die 'Failed to get example data directory!';

    return ($entities{build}, $entities{example_build});
}

sub build_with_example_build_for_454 {
    my ($build, $example_build) = build_with_example_build('454');
    return ($build, $example_build);
}

sub instrument_data_454 {
    my $library = library();

    my $inst_data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model/MetagenomicComposition16s454/inst_data';
    $entities{instrument_data_454} = Genome::InstrumentData::454->create(
        id => --$id,
        run_name => 'R_2010_01_09_11_08_12_FLX08080418_Administrator_100737113',
        region_number => 3,
        total_reads => 20,
        library => $library,
        sequencing_platform => '454',
        archive_path =>  $inst_data_dir.'/archive.tgz',
    ) if not $entities{instrument_data_454};
    die 'Failed to create instrument data 454!' if not $entities{instrument_data_454};

    die 'archive_path' if not -s $entities{instrument_data_454}->attributes(attribute_label => 'archive_path')->attribute_value;

    return $entities{instrument_data_454};
}

1;

