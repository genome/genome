package Genome::Model::GenotypeMicroarray::Build::CreateOriginalGenotypeFiles;

use strict;
use warnings;

use Genome;

require File::Basename;

class Genome::Model::GenotypeMicroarray::Build::CreateOriginalGenotypeFiles {
    is => 'Command::V2',
    has => {
        build => {
            is => 'Genome::Model::Build::GenotypeMicroarray',
            doc => 'The genotype build to use.',
        },
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
        $self->error_message('Instrument dat genotype file does not exist! '.$genotype_file);
        return;
    }

    my $reader = Genome::Model::GenotypeMicroarray::GenotypeFile::ReaderFactory->build_reader(
        source => $build,
    );
    return if not $reader;

    my $vcf_writer = Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory->build_writer(
        header => $reader->header,
        string => $build->original_genotype_vcf_file_path,
    );
    return if not $vcf_writer;

    my $tsv_writer = Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory->build_writer(
        header => $reader->header,
        string => $build->original_genotype_file_path.':format=csv:separator=tab',
    );
    return if not $tsv_writer;

    my %alleles;
    my $genotypes_output = 0;
    GENOTYPE: while ( my $genotype = $reader->next ) {
        $vcf_writer->write($genotype);
        $tsv_writer->write($genotype);
        $alleles{ $genotype->sample_field(0, 'ALLELES') }++;
        $genotypes_output++;
    }
    $self->alleles(\%alleles);
    $self->genotypes_output($genotypes_output);
    $self->debug_message("Genotypes output: $genotypes_output");

    $self->debug_message('Create original genotype files...done');
    return 1;
}

1;

