package Genome::Model::DeNovoAssembly::Build::Test;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
use Test::More;

sub base_test_dir {
    return $ENV{GENOME_TEST_INPUTS} . '/Genome-Model/DeNovoAssembly';
}

my %entities;
sub instrument_data_for_assembler {
    my ($class, $assembler) = @_;

    Carp::confess('No assembler to get instrument data!') if not $assembler;
    Carp::confess('Unsupported assembler to get instrument data!') if $assembler ne 'soap de-novo-assemble';

    my $key = $assembler.'_inst_data';
    return $entities{$key} if $entities{$key};

    if ( not $entities{sample} ) {
        my $taxon = Genome::Taxon->create(
            name => 'Escherichia coli TEST',
            domain => 'Bacteria',
            estimated_genome_size => 4500000,
            species_latin_name => 'Escherichia coli',
            strain_name => 'TEST',
        );
        ok($taxon, 'taxon') or die;
        $entities{taxon} = $taxon;

        my $source = Genome::Individual->create(
            id => -123,
            name => 'TEST-000',
        );
        ok($source, 'source') or die;

        my $sample = Genome::Sample->create(
            id => -1234,
            name => 'TEST-000',
            source_id => $source->id,
        );
        ok($sample, 'sample') or die;
        $entities{sample} = $sample;

        $entities{library} = Genome::Library->create(
            id => -12345,
            name => $sample->name.'-testlibs',
            sample_id => $sample->id,
            library_insert_size => 260,
        );
        ok($entities{library}, 'library') or die;
    }

    my $archive_path = base_test_dir().'/inst_data/-7777/archive.tgz';
    ok(-s $archive_path, 'inst data archive path') or die;
    $entities{$key} = Genome::InstrumentData::Solexa->create(
        id => -7777,
        sequencing_platform => 'solexa',
        read_length => 100,
        subset_name => '8-CGATGT',
        index_sequence => 'CGATGT',
        run_name => 'XXXXXX/8-CGATGT',
        run_type => 'Paired',
        flow_cell_id => 'XXXXXX',
        lane => 8,
        library => $entities{library},
        archive_path => $archive_path,
        clusters => 15000,
        fwd_clusters => 15000,
        rev_clusters => 15000,
        analysis_software_version => 'GAPipeline-0.3.0',
    );
    ok($entities{$key}, 'instrument data');
    ok($entities{$key}->is_paired_end, 'inst data is paired');
    ok(-s $entities{$key}->archive_path, 'inst data archive path');

    return $entities{$key};
}

sub processing_profile_for_assembler {
    my ($class, $assembler) = @_;

    Carp::confess('No assembler to get processing profile!') if not $assembler;

    my $key = $assembler.'_pp';
    return $entities{$key} if $entities{$key};

    my %assembler_and_pp_params = (
        'soap de-novo-assemble' => {
            name => 'De Novo Assembly Soap PGA Test',
            assembler_name => 'soap de-novo-assemble',
            assembler_version => '1.04',
            assembler_params => '-kmer_size 31 -resolve_repeats -kmer_frequency_cutoff 1',
            read_processor => 'DEFAULT (trim bwa-style -trim-qual-level 10 | filter by-length -filter-length 35 | rename illumina-to-pcap, coverage 10X)' ,
            post_assemble => 'standard-outputs -min_contig_length 10',
        },
    );

    Carp::confess('No processing profile params for assembler! '.$assembler) if not $assembler_and_pp_params{$assembler};

    $entities{$key} = Genome::ProcessingProfile::DeNovoAssembly->create($assembler_and_pp_params{$assembler});
    ok($entities{$key}, 'pp') or die;

    return $entities{$key};
}

sub model_for_assembler {
    my ($class, $assembler) = @_;

    Carp::confess('No assembler to get model!') if not $assembler;

    my $key = $assembler.'_model';
    return $entities{$key} if $entities{$key};

    my $pp = $class->processing_profile_for_assembler($assembler); # dies on error
    my $instrument_data = $class->instrument_data_for_assembler($assembler); # dies on error
    $entities{$key} = Genome::Model::DeNovoAssembly->create(
        processing_profile => $pp,
        subject_name => $entities{taxon}->name,
        subject_type => 'species_name',
        center_name => 'WUGC',
    );
    ok($entities{$key}, 'model') or die;
    ok($entities{$key}->add_instrument_data($instrument_data), 'add inst data to model') or die;

    return $entities{$key};
}

sub build_for_assembler { 
    my ($class, $assembler) = @_;

    Carp::confess('No assembler to get model!') if not $assembler;

    my $model = $class->model_for_assembler($assembler);
    my $tmpdir_template = "/DeNovoAssembly-Soap.t-XXXXXXXX";
    my $tmpdir = File::Temp::tempdir($tmpdir_template, CLEANUP => 1, TMPDIR => 1);
    ok(-d $tmpdir, 'temp dir: '.$tmpdir);

    my $build = Genome::Model::Build::DeNovoAssembly->create(
        model => $model,
        data_directory => $tmpdir,
    );
    ok($build->instrument_data, 'build instrument data');

    return $build;
}

sub example_build {
    my ($class, $assembler) = @_;

    Carp::confess('No assembler to get example build!') if not $assembler;

    my $example_dir = base_test_dir().'/soap_v17';
    ok(-d $example_dir, 'example dir') or die;

    my $model = $class->model_for_assembler($assembler);
    my $example_build = Genome::Model::Build->create(
        model => $model,
        data_directory => $example_dir,
    );
    ok($example_build, 'create example build');

    return $example_build;
}

1;

