package Genome::Model::GenotypeMicroarray::GenotypeFile::ReaderFactory;

use strict;
use warnings;

use Genome;

use Genome::File::Vcf::Reader;
use Genome::Model::GenotypeMicroarray::GenotypeFile::FromInstDataWithAnnotationReader;

class Genome::Model::GenotypeMicroarray::GenotypeFile::ReaderFactory { 
    is => 'UR::Singleton',
};

my $current_sample_name;

sub build_reader {
    my ($class, %params) = @_;

    Carp::confess('Nothing given to build reader!') if not %params;
    my $source = delete $params{source};
    Carp::confess('No source given to build reader!') if not $source;

    my $reader;
    if ( $source->isa('Genome::InstrumentData') ) {
        my $variation_list_build = delete $params{variation_list_build};
        Carp::confess('No variation list build given to build reader!') if not $variation_list_build;
        $reader = $class->_build_reader_for_instrument_data($source, $variation_list_build);
    }
    elsif ( $source->isa('Genome::Model::Build::GenotypeMicroarray') ) {
        $reader = $class->_build_reader_for_build($source);
    }
    else {
        Carp::confess('Do not know how to build genotype file reader for source! '.$source->__display_name__);
    }
    return if not $reader;

    my $annotate_header = $class->_annotate_header($reader);
    return if not $annotate_header;

    return $reader;
}

sub _build_reader_for_instrument_data {
    my ($class, $instrument_data, $variation_list_build) = @_;

    my $reader = Genome::Model::GenotypeMicroarray::GenotypeFile::FromInstDataWithAnnotationReader->create(
        instrument_data => $instrument_data,
        variation_list_build => $variation_list_build,
    );

    $current_sample_name = $instrument_data->sample->name;

    return $reader;
}

sub _build_reader_for_build {
    my ($class, $build) = @_;

    # VCF
    my $genotype_file = $build->original_genotype_vcf_file_path;
    if ( -s $genotype_file ) {
        $current_sample_name = $build->subject->name;
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

sub _annotate_header {
    my ($class, $reader) = @_;

    Carp::confess('No reader given to annotate header!') if not $reader;

    my $header = $reader->header;
    Carp::confess('No header found on reader to annotate!') if not $header;

    # Add sample name
    if ( not grep { $_ eq $current_sample_name } $header->sample_names ) {
        $header->add_sample_name($current_sample_name);
    }

    # Add sample format fields
    my @header_format_type_ids;
    if ( $header->format_types ) {
        @header_format_type_ids = keys %{$header->format_types};
    }
    for my $format_type ( Genome::Model::GenotypeMicroarray->format_types ) {
        next if grep { $_ eq $format_type->{id} } @header_format_type_ids;
        $header->add_format_str('<ID='.$format_type->{id}.$format_type->{header});
    }

    $reader->header($header);

    return 1;
}

1;

