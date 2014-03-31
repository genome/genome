package Genome::Model::Build::MetagenomicComposition16s::TestBuildFactory;

use strict;
use warnings;

use Genome;

use Date::Format;
use File::Temp;
use Test::More;

my $id = -9000;
my %entities;
my %pp_params = (
    454 => {
        amplicon_processor => 'filter by-min-length --length 200',
        chimera_detector => 'chimera-slayer',
        chimera_detector_params => "--nastier-params '-num_top_hits 10' --chimera-slayer-params '-windowSize 50 -printCSalignments -windowStep 5'",
        classifier => 'rdp2-2',
        classifier_params => '-training-set 6 -format hmp_fix_ranks -version 2x2',
    },
    sanger => {
        amplicon_processor => 'filter by-min-length --length 1150',
        assembler => 'phred_phrap',
        assembler_params => '-vector_bound 0 -trim_qual 0',
        classifier => 'rdp2-1',
        classifier_params => '-training-set broad -format hmp_fix_ranks -version 2x1',
    }
);

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

    return ( $entities{'build_'.$sequencing_platform}, $entities{'example_build_'.$sequencing_platform} ) if $entities{'build_'.$sequencing_platform};

    my $pp = Genome::ProcessingProfile->__define__(
        id => --$id,
        type_name => 'metagenomic composition 16s',
        %{$pp_params{$sequencing_platform}},
    );
    die 'Failed to create MC16s processing profile!' if not $pp;

    my $sample = sample();
    my $model = Genome::Model::MetagenomicComposition16s->__define__(
        id => --$id,
        name => '__TEST_MC16S_MODEL__',
        processing_profile => $pp,
        subject_name => $sample->name,
        subject_type => 'sample_name'
    );
    die 'Failed to create MC16s model!' if not $model;

    $entities{'build_'.$sequencing_platform} = Genome::Model::Build::MetagenomicComposition16s->create(
        id => --$id,
        model => $model,
        data_directory => File::Temp::tempdir(CLEANUP => 1),
    );
    die 'Failed to create MC16s model!' if not $entities{'build_'.$sequencing_platform};
    $entities{'build_'.$sequencing_platform}->create_subdirectories;
    $entities{'build_'.$sequencing_platform}->the_master_event->event_status('Succeeded');
    my $time = time();
    my $date_template = UR::Context->date_template;
    my $timestamp = Date::Format::time2str($date_template, $time);
    $entities{'build_'.$sequencing_platform}->the_master_event->date_completed($timestamp);

    $entities{'example_build_'.$sequencing_platform} = Genome::Model::Build->create(
        model=> $model,
        id => --$id,
    );
    die 'Failed to create MC16s model!' if not $entities{'example_build_'.$sequencing_platform};

    my %example_data_directories = (
        454 => $ENV{GENOME_TEST_INPUTS} . '/Genome-Model/MetagenomicComposition16s454/build_v5.2chimeras', # start w/ 2 chimeras
        sanger => $ENV{GENOME_TEST_INPUTS} . '/Genome-Model/MetagenomicComposition16sSanger/build_v3',
    );
    $entities{'example_build_'.$sequencing_platform}->data_directory( $example_data_directories{$sequencing_platform} ) or die 'Failed to get example data directory!';
    $entities{'example_build_'.$sequencing_platform}->the_master_event->event_status('Succeeded');
    my $past_timestamp = Date::Format::time2str($date_template, $time - 60);
    $entities{'example_build_'.$sequencing_platform}->the_master_event->date_completed($past_timestamp);

    my $instrument_data_method = 'instrument_data_'.$sequencing_platform;
    my $instrument_data = __PACKAGE__->$instrument_data_method;
    $model->add_instrument_data($instrument_data) or die 'Failed to add instrument data to model!';
    $entities{'build_'.$sequencing_platform}->add_instrument_data($instrument_data) or die 'Failed to add instrument data to build!';
    $entities{'example_build_'.$sequencing_platform}->add_instrument_data($instrument_data) or die 'Failed to add instrument data to example build!';

    return ($entities{'build_'.$sequencing_platform}, $entities{'example_build_'.$sequencing_platform});
}

sub build_with_example_build_for_454 {
    my ($build, $example_build) = build_with_example_build('454');
    return ($build, $example_build);
}

sub instrument_data_454 {
    my $library = library();

    my $inst_data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model/MetagenomicComposition16s454/inst_data';
    $entities{instrument_data_454} = Genome::InstrumentData::Solexa->create(# don't look
        id => --$id,
        run_name => 'R_2010_01_09_11_08_12_FLX08080418_Administrator_100737113',
        region_number => 3,
        total_reads => 20,
        read_count => 20,
        clusters => 20,
        library => $library,
        sequencing_platform => '454',
        archive_path =>  $inst_data_dir.'/archive.tgz',
        analysis_software_version => 'CASAVA-1.8',
    ) if not $entities{instrument_data_454};
    die 'Failed to create instrument data 454!' if not $entities{instrument_data_454};

    die 'archive_path 454!' if not -s $entities{instrument_data_454}->attributes(attribute_label => 'archive_path')->attribute_value;

    return $entities{instrument_data_454};
}

sub build_with_example_build_for_sanger {
    my ($build, $example_build) = build_with_example_build('sanger');
    return ($build, $example_build);
}

sub instrument_data_sanger {
    $entities{instrument_data_sanger} = Genome::InstrumentData::Sanger->create(
        id => '01jan00.101amaa',
        library => $entities{library},
    );
    die 'Failed to create instrument data sanger!' if not $entities{instrument_data_sanger};

    no warnings qw/ once redefine /;
    *Genome::InstrumentData::Sanger::full_path = sub{ $ENV{GENOME_TEST_INPUTS} . '/Genome-Model/MetagenomicComposition16sSanger/inst_data/'.$entities{instrument_data_sanger}->id; };
    use warnings;

    die 'full_path sanger!' if not -d $entities{instrument_data_sanger}->full_path;

    return $entities{instrument_data_sanger};
}

1;

