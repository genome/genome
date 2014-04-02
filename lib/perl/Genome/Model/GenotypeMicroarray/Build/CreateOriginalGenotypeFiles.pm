package Genome::Model::GenotypeMicroarray::Build::CreateOriginalGenotypeFiles;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::Build::CreateOriginalGenotypeFiles {
    is => 'Command::V2',
    has_input_output => {
        build => { is => 'Genome::Model::Build::GenotypeMicroarray', },
    },
    has_optional_transient => {
        alleles => { is => 'Hash', },
        genotypes_output => { is => 'Number', },
    },
};

sub help_brief {
    return 'Create original genotype files as VCF and TSV';
}

sub help_detail {
    return <<HELP;
HELP
}

sub execute {
    my $self = shift;
    $self->debug_message('Create original genotype files...');

    my $build = $self->build;
    my $genotype_file = $build->instrument_data->genotype_file;
    if ( not -s $genotype_file ) {
        $self->error_message('Instrument data genotype file does not exist! '.$genotype_file);
        return;
    }

    my $reader = Genome::Model::GenotypeMicroarray::GenotypeFile::ReaderFactory->build_reader(
        source => $build,
    );
    return if not $reader;

    my $original_genotype_vcf_file = $build->original_genotype_vcf_file_path;
    $self->debug_message('Original genotype VCF: '.$original_genotype_vcf_file);
    my $vcf_writer = Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory->build_writer(
        header => $reader->header,
        string => $original_genotype_vcf_file,
    );
    return if not $vcf_writer;

    my $original_genotype_file = $build->original_genotype_file_path;
    $self->debug_message('Original genotype TSV: '.$original_genotype_file);
    my $tsv_writer = Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory->build_writer(
        header => $reader->header,
        string => $original_genotype_file.':format=csv:separator=tab',
    );
    return if not $tsv_writer;

    my %alleles;
    my $genotypes_output = 0;
    GENOTYPE: while ( my $genotype = $reader->read ) {
        $vcf_writer->write($genotype);
        $tsv_writer->write($genotype);
        $alleles{ $genotype->sample_field(0, 'ALLELES') }++;
        $genotypes_output++;
    }
    $self->alleles(\%alleles);
    $self->genotypes_output($genotypes_output);
    $self->debug_message("Genotypes output: $genotypes_output");

    if ( not -e $original_genotype_vcf_file ) {
        $self->error_message('Executed extract command to create original genotype VCF file and genotypes were output, but file is gone!');
        return;
    }

    if ( not -e $original_genotype_file ) {
        $self->error_message('Executed command to create original genotype file and genotypes were output, but file is gone!');
        return;
    }

    my @alleles = grep { $_ ne '--' } keys %alleles;
    if ( not @alleles ) {
        $self->error_message('Executed command to create original genotype file, but there are no alleles!');

        return;
    }

    $self->debug_message('Create original genotype files...done');
    return 1;
}

1;

