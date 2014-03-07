package Genome::Model::GenotypeMicroarray::GenotypeFile::ReaderFactory;

use strict;
use warnings;

use Genome;

use Genome::File::Vcf::Reader;
use Genome::Model::GenotypeMicroarray::GenotypeFile::ReadUnannotatedCsv;

class Genome::Model::GenotypeMicroarray::GenotypeFile::ReaderFactory { 
    is => 'UR::Singleton',
};

sub build_reader {
    my ($class, %params) = @_;

    Carp::confess('Nothing given to build reader!') if not %params;
    my $source = delete $params{source};
    Carp::confess('No source given to build reader!') if not $source;

    if ( $source->isa('Genome::InstrumentData') ) {
        my $variation_list_build = delete $params{variation_list_build};
        Carp::confess('No variation list build given to build reader!') if not $variation_list_build;
        return $class->_build_reader_for_instrument_data($source, $variation_list_build);
    }
    elsif ( $source->isa('Genome::Model::Build::GenotypeMicroarray') ) {
        return $class->_build_reader_for_build($source);
    }
    else {
        Carp::confess('Do not know how to build genotype file reader for source! '.$source->__display_name__);
    }

}

sub _build_reader_for_instrument_data {
    my ($class, $instrument_data, $variation_list_build) = @_;

    my $genotype_file = eval{ $instrument_data->attributes(attribute_label => 'genotype_file')->attribute_value; };
    if ( not $genotype_file or not -s $genotype_file ) {
        $class->error_message('No genotype file for instrument data! '.$instrument_data->id);
        return;
    }

    my $snp_id_mapping = Genome::InstrumentData::Microarray->get_snpid_hash_for_variant_list(
        $instrument_data, $variation_list_build
    );

    my $reader = Genome::Model::GenotypeMicroarray::GenotypeFile::ReadUnannotatedCsv->create(
        input => $genotype_file,
        variation_list_build => $variation_list_build,
        snp_id_mapping => $snp_id_mapping,
    );

    return $reader;
}

sub _build_reader_for_build {
    my ($class, $build) = @_;

    # VCF
    my $genotype_file = $build->original_genotype_vcf_file_path;
    if ( -s $genotype_file ) {
        return Genome::File::Vcf::Reader->new($genotype_file);
    }

    # Use inst data
    my $instrument_data = $build->instrument_data;
    if ( not $instrument_data ) {
        $class->error_message('No instrument data for genotype build! '.$build->id);
        return;
    }

    my $variation_list_build = $build->dbsnp_build;
    if ( not $variation_list_build ) {
        $class->error_message('No variation list build for genotype build! '.$build->id);
        return;
    }

    return $class->_build_reader_for_instrument_data($instrument_data, $variation_list_build);
}

1;

