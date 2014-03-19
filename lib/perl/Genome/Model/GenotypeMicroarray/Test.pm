package Genome::Model::GenotypeMicroarray::Test;

use strict;
use warnings;

require File::Temp;

sub testdir {
    return $ENV{GENOME_TEST_INPUTS} . '/GenotypeMicroarray/v4';
}

my %cache;
sub reference_sequence_build {
    return $cache{reference_sequence_build} if $cache{reference_sequence_build};

    $cache{reference_sequence_build} = Genome::Model::Build::ReferenceSequence->__define__(
        name => '__TEST_REF__',
    );
    die 'Failed to define reference sequence build' if not $cache{reference_sequence_build};

    my %pos_seq = (    # Alleles
        752566 => 'A', # AG 
        752721 => 'G', # AG
        798959 => 'C', # AG SNV
        800007 => 'T', # TC
        838555 => 'C', # CC
        846808 => 'C', # CC
        854250 => 'A', # AA
        861808 => 'G', # GG DUP
        873558 => 'G', # TT SNV
        882033 => 'G', # GG
    );

    no warnings;
    *Genome::Model::Build::ReferenceSequence::sequence = sub{ return $pos_seq{$_[2]}; };
    *Genome::Model::Build::ReferenceSequence::chromosome_array_ref = sub{ return [qw/ 1 /]; };
    use warnings;

    return $cache{reference_sequence_build};
}

sub variation_list_model {
    return $cache{variation_list_model} if $cache{variation_list_model};

    $cache{variation_list_model} = Genome::Model::ImportedVariationList->__define__(
        name => '__TEST_DBSNP__',
        subject => Genome::Taxon->__define__(name => '__TEST_HUMAN__'),
    );

    return $cache{variation_list_model};
}

sub variation_list_build {
    return $cache{variation_list_build} if $cache{variation_list_build};

    my $dbsnp_file = testdir().'/dbsnp/dbsnp.132';
    my $fl = Genome::Model::Tools::DetectVariants2::Result::Manual->__define__(
        description => '__TEST__DBSNP132__',
        username => 'apipe-tester',
        file_content_hash => 'c746fb7b7a88712d27cf71f8262dd6e8',
        output_dir => testdir().'/dbsnp',
    );
    $fl->lookup_hash($fl->calculate_lookup_hash());

    $cache{variation_list_build} = Genome::Model::Build::ImportedVariationList->__define__(
        model => variation_list_model(),
        reference => reference_sequence_build(),
        snv_result => $fl,
        version => 'vT',
    );
    die 'Failed to define variation list build' if not $cache{variation_list_build};

    return $cache{variation_list_build};
}

sub sample {
    return $cache{sample} if $cache{sample};
    instrument_data();
    return $cache{sample};
}

sub instrument_data {
    return $cache{instrument_data} if $cache{instrument_data};

    $cache{sample} = Genome::Sample->__define__(
        id => 2879594813,
        name => '__TEST_SAMPLE__',
    );
    die 'Failed to define sample' if not $cache{sample};

    my $library = Genome::Library->__define__(
        name => $cache{sample}->name.'-microarraylib',
        sample => $cache{sample},
    );

    $cache{instrument_data} = Genome::InstrumentData::Imported->__define__(
        library => $library,
        import_format => 'genotype file',
        sequencing_platform => 'infinium',
        import_source_name => 'WUGC',
    );
    die 'Failed to define instrument data' if not $cache{instrument_data};
    $cache{instrument_data}->add_attribute(attribute_label => 'genotype_file', attribute_value => testdir().'/instdata/snpreport/genotypes.tsv');
    die 'Failed to set genotype file on instdata!' if not -s $cache{instrument_data}->genotype_file;

    $cache{sample}->default_genotype_data_id($cache{instrument_data}->id);
    die 'Failed to set default genotype data id on sample!' if $cache{sample}->default_genotype_data_id ne $cache{instrument_data}->id;

    return $cache{instrument_data};
}

sub processing_profile {
    return $cache{pp} if $cache{pp};

    $cache{pp} = Genome::ProcessingProfile::GenotypeMicroarray->__define__(
        name => '__TEST_PP__ wugc',
    );

    return $cache{pp};
}

sub model {
    return $cache{model} if $cache{model};

    my $instrument_data = instrument_data();
    my $sample = $instrument_data->sample;
    $cache{model} = Genome::Model::GenotypeMicroarray->create(
        name => 'Test Genotype Microarray pp',
        processing_profile => processing_profile(),
        subject => $sample,
        subject_id => $sample->id,
        subject_class_name => $sample->class,
        dbsnp_build => variation_list_build(),
    );
    $cache{model}->add_instrument_data($instrument_data);

    return $cache{model};
}

sub example_build {
    return $cache{example_build} if $cache{example_build};

    my $model = model();
    $cache{example_build} = Genome::Model::Build::GenotypeMicroarray->create(
        model => $model,
        data_directory => testdir().'/build',
    );

    return $cache{example_build};
}

sub example_legacy_build {
    return $cache{example_legacy_build} if $cache{example_legacy_build};

    my $model = model();
    $cache{example_legacy_build} = Genome::Model::Build::GenotypeMicroarray->create(
        model => $model,
        data_directory => testdir().'/build-legacy',
    );

    return $cache{example_legacy_build};
}

sub build {

    my $model = model();
    my $tempdir = File::Temp::tempdir(CLEANUP => 1);
    my $build = Genome::Model::Build::GenotypeMicroarray->create(
        model => $model,
        data_directory => $tempdir,
    );

    return $build;
}

1;

