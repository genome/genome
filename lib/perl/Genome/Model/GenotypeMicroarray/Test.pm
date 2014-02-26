package Genome::Model::GenotypeMicroarray::Test;

use strict;
use warnings;


sub testdir {
    return $ENV{GENOME_TEST_INPUTS} . '/GenotypeMicroarray/v2';
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

sub variation_list_build { # FIXME upgrade to VCF
    return $cache{variation_list_build} if $cache{variation_list_build};

    my $dbsnp_file = testdir().'/dbsnp/dbsnp.132';
    my $fl = Genome::Model::Tools::DetectVariants2::Result::Manual->__define__(
        description => '__TEST__DBSNP132__',
        username => 'apipe-tester',
        file_content_hash => 'c746fb7b7a88712d27cf71f8262dd6e8',
        output_dir => testdir().'/dbsnp',
    );
    $fl->lookup_hash($fl->calculate_lookup_hash());

    my $variation_list_model = Genome::Model::ImportedVariationList->__define__(
        name => '__TEST_DBSNP__',
        subject => Genome::Taxon->__define__(name => '__TEST_HUMAN__'),
    );

    $cache{variation_list_build} = Genome::Model::Build::ImportedVariationList->__define__(
        model => $variation_list_model,
        reference => reference_sequence_build(),
        snv_result => $fl,
        version => 'TEST',
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
    );
    die 'Failed to define instrument data' if not $cache{instrument_data};
    $cache{instrument_data}->add_attribute(attribute_label => 'genotype_file', attribute_value => testdir().'/instdata/snpreport/genotypes.tsv');
    die 'Failed to set genotype file on instdata!' if not -s $cache{instrument_data}->genotype_file;

    $cache{sample}->default_genotype_data_id($cache{instrument_data}->id);
    die 'Failed to set default genotype data id on sample!' if $cache{sample}->default_genotype_data_id ne $cache{instrument_data}->id;

    return $cache{instrument_data};
}

sub model {
    return $cache{model} if $cache{model};

    my $instrument_data = instrument_data();
    my $sample = $instrument_data->sample;
    my $pp = Genome::ProcessingProfile::GenotypeMicroarray->__define__(name => '__TEST_PP__');
    $cache{model} = Genome::Model::GenotypeMicroarray->create(
        name => 'Test Genotype Microarray pp',
        processing_profile => $pp,
        subject => $sample,
        subject_id => $sample->id,
        subject_class_name => $sample->class,
        reference_sequence_build => reference_sequence_build(),
        dbsnp_build => variation_list_build(),
    );
    $cache{model}->add_instrument_data($instrument_data);

    return $cache{model};
}

sub example_build {
    return $cache{example_build} if $cache{example_build};

    my $model = model();
    $cache{example_build} = Genome::Model::Build::GenotypeMicroarray->__define__(
        model => $model,
        reference_sequence_build => $model->reference_sequence_build,
        dbsnp_build => $model->dbsnp_build,
        data_directory => testdir().'/build',
    );

    return $cache{example_build};
}

sub example_legacy_build {
    return $cache{example_legacy_build} if $cache{example_legacy_build};

    my $model = model();
    $cache{example_legacy_build} = Genome::Model::Build::GenotypeMicroarray->__define__(
        model => $model,
        reference_sequence_build => $model->reference_sequence_build,
        dbsnp_build => $model->dbsnp_build,
        data_directory => testdir().'/build-legacy',
    );
    $cache{example_legacy_build}->add_instrument_data( $model->instrument_data );

    return $cache{example_legacy_build};
}

sub build {

    my $model = model();
    my $tempdir = File::Temp::tempdir(CLEANUP => 1);
    my $build = Genome::Model::Build::GenotypeMicroarray->create(
        model => $model,
        reference => $model->reference,
        dbsnp_build => $model->dbsnp_build,
        data_directory => $tempdir,
    );
    ok($build, 'create genotype microarray build');
    $build->add_instrument_data( $model->instrument_data );

    return $build;
}

sub expected_genotypes {
    return $cache{expected_genotypes} if $cache{expected_genotypes};
    $cache{expected_genotypes} = [
        {
            'id' => 'rs3094315',
            'chromosome' => '1',
            'position' => '752566',
            'allele1' => 'A',
            'allele2' => 'G',
            'cnv_confidence' => 'NA',
            'cnv_value' => '2.0',
            'gc_score' => '0.8931',
            'log_r_ratio' => '-0.3639',
        },
        {
            'id' => 'rs3131972',
            'chromosome' => '1',
            'position' => '752721',
            'allele1' => 'A',
            'allele2' => 'G',
            'cnv_confidence' => 'NA',
            'cnv_value' => '2.0',
            'gc_score' => '0.9256',
            'log_r_ratio' => '-0.0539',
        },
        {
            'id' => 'rs11240777',
            'chromosome' => '1',
            'position' => '798959',
            'allele1' => 'A',
            'allele2' => 'G',
            'cnv_confidence' => 'NA',
            'cnv_value' => '2.0',
            'gc_score' => '0.8729',
            'log_r_ratio' => '-0.0192',
        },
        {
            'id' => 'rs6681049',
            'chromosome' => '1',
            'position' => '800007',
            'allele1' => 'T',
            'allele2' => 'C',
            'cnv_confidence' => 'NA',
            'cnv_value' => '2.0',
            'gc_score' => '0.7156',
            'log_r_ratio' => '0.2960',
        },
        {
            'id' => 'rs4970383',
            'chromosome' => '1',
            'position' => '838555',
            'allele1' => 'C',
            'allele2' => 'C',
            'cnv_confidence' => 'NA',
            'cnv_value' => '2.0',
            'gc_score' => '0.8749',
            'log_r_ratio' => '0.4694',
        },
        {
            'id' => 'rs4475691',
            'chromosome' => '1',
            'position' => '846808',
            'allele1' => 'C',
            'allele2' => 'C',
            'cnv_confidence' => 'NA',
            'cnv_value' => '2.0',
            'gc_score' => '0.8480',
            'log_r_ratio' => '-0.0174',
        },
        {
            'id' => 'rs7537756',
            'chromosome' => '1',
            'position' => '854250',
            'allele1' => 'A',
            'allele2' => 'A',
            'cnv_confidence' => 'NA',
            'cnv_value' => '2.0',
            'gc_score' => '0.8670',
            'log_r_ratio' => '0.0389',
            'sample_id' => '2879594813',
        },
        {
            'id' => 'rs1110052',
            'chromosome' => '1',
            'position' => '873558',
            'allele1' => 'T',
            'allele2' => 'T',
            'cnv_confidence' => 'NA',
            'cnv_value' => '2.0',
            'gc_score' => '0.7787',
            'log_r_ratio' => '0.1487',
        },
        {
            'id' => 'rs2272756',
            'chromosome' => '1',
            'position' => '882033',
            'allele1' => 'G',
            'allele2' => 'G',
            'cnv_confidence' => 'NA',
            'cnv_value' => '2.0',
            'gc_score' => '0.8677',
            'log_r_ratio' => '-0.0801',
        },
    ];

    my $refseq = reference_sequence_build();
    my $sam = sample();
    for my $genotype ( @{$cache{expected_genotypes}} ) {
        $genotype->{sample_id} = $sam->id;
        #$genotype->{sample_name} = $sam->name;
        $genotype->{identifiers} = [$genotype->{id}];
        $genotype->{chrom} = $genotype->{chromosome};
        $genotype->{reference} = $refseq->sequence($genotype->{chromosome}, $genotype->{position}, $genotype->{position});
        $genotype->{reference_allele} = $genotype->{reference};
        $genotype->{alleles} = $genotype->{allele1}.$genotype->{allele2};
        $genotype->{alternate_alleles} = Genome::Model::GenotypeMicroarray::GenotypeFile::ReadTsvAndAnnotate->_alts_for_genotype($genotype);
        $genotype->{quality} = '.';
        $genotype->{_filter} = [];
        $genotype->{_format} = [Genome::Model::GenotypeMicroarray::GenotypeFile::ReadTsvAndAnnotate->_format_for_genotype($genotype)];
        $genotype->{info_fields} = Genome::Model::GenotypeMicroarray::GenotypeFile::ReadTsvAndAnnotate->_info_hash_for_genotype($genotype);
        $genotype->{_sample_data} = [ [] ], #FIXME
    }

    return $cache{expected_genotypes};
}

1;

