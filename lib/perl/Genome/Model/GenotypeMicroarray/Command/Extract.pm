package Genome::Model::GenotypeMicroarray::Command::Extract;

use strict;
use warnings;

use Genome;

require File::Basename;

class Genome::Model::GenotypeMicroarray::Command::Extract {
    is => 'Command::V2',
    has_optional => {
        # SOURCE
        model => {
            is => 'Genome::Model::GenotypeMicroarray',
            doc => 'The genotype model to work with. This will get the most recent succeeded build.',
        },
        build => {
            is => 'Genome::Model::Build::GenotypeMicroarray',
            doc => 'The genotype build to use.',
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'The genotype instrument data to work with.',
        },
        sample => {
            is => 'Genome::Sample',
            doc => 'The sample instrument data to work with.',
        },
        sample_type_priority => {
            is => 'Text',
            is_many => 1,
            default_value => [qw/ default internal external /],
            valid_values => [qw/ default internal external /],
            doc => 'Priority of the sample type to favor when getting microarray instrument data.',
        },
        variation_list_build => { # req for sample and instdata
            is => 'Genome::Model::Build::ImportedVariationList',
            doc => 'Imported variation list build. Give id from command line. Commonly used:
                     ID          REFERENCE                   VERSION
                     106227442   dbSNP-NCBI-human-build36    130
                     106375969   dbSNP-g1k-human-build37     132',
        },
        # OUTPUT
        output => {
            is => 'Text',
            doc => <<DOC
Output config. Give give key value (key=value) pairs separated by colons (:).
Defaults to writing all fields with headers as CSV w/ tabs to STDOUT.
Options:

 file           File to write to. Use '-' for STDOUT.
 format         Format to output. Valid vlaues: csv [character separated values], vcf

CSV Options:
 separator  The separator to use. Default is 'TAB', which writes values separated by tabs.
 headers    Include headers. Use 1 to include, 0 to exclude.
 fields     Write these fields. Default writes all.
             Fields: 
              id
              chromosome
              position
              sample_id
              log_r_ratio
              gc_score
              cnv_value
              cnv_confidence
              alleles
              allele1
              allele2 

VCF Options:
 None.

Examples:
             
 Write to FILE as comma separated values, with headers and all available fields:

  output=FILE:separator=,

 Write to FILE as TAB separated values, without headers and only chromosome, position
  and alleles fields:

  output=FILE:headers=chromosome,position,alleles:print_headers=0

DOC
            ,
        },
        filters => {
            is => 'Text',
            is_many => 1,
            doc => "Filter genotypes. Give name and parameters, if required. Filters:\n gc_scrore => filter by min gc score (Ex: gc_score:min=0.7)\n invalid_iscan_ids => list of invalid iscan snvs compiled by Nate",
        },
    },
    has_optional_transient => {
        alleles => { is => 'Hash', },
        metrics => { is => 'Hash', },
        source => { is => 'Text', },
        source_type => { is => 'Text', },
    },
    has_calculated => {
        genotypes_input => { is => 'Number', calculate => q( return $self->metrics->{input}; ), },
        genotypes_filtered => {
            is => 'Number', calculate => q( return $self->metrics->{output} - $self->metrics->{input}; ), 
        },
        genotypes_output => { is => 'Number', calculate => q( return $self->metrics->{output}; ), },
    },
};

sub help_brief {
    return 'extract genotype data from a build or inst data';
}

sub help_detail {
    return <<HELP;
HELP
}

sub execute {
    my $self = shift;
    $self->status_message('Extract genotytpes...');

    $DB::single=1;
    my $resolve_source = $self->resolve_source;
    return if not $resolve_source;

    my $filters = $self->_create_filters;
    return if not $filters;

    my $reader = Genome::Model::GenotypeMicroarray::GenotypeFile::ReaderFactory->build_reader($self->source, $self->variation_list_build);
    return if not $reader;

    my $writer = Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory->build_writer($self->output.':sample_name='.$self->sample->name);
    return if not $writer;

    my %alleles;
    my %metrics = ( input => 0, filtered => 0, output => 0, );
    GENOTYPE: while ( my $genotype = $reader->next ) {
        $metrics{input}++;
        for my $filter ( @$filters ) {
            next GENOTYPE if not $filter->filter($genotype);
        }
        $metrics{output}++;
        $alleles{ $genotype->{alleles} }++;
        $writer->write($genotype);
    }
    $self->metrics(\%metrics);
    $self->alleles(\%alleles);
    for my $name ( map { 'genotypes_'.$_ } (qw/ input filtered output /) ) {
        $self->status_message(ucfirst(join(' ', split('_', $name))).": ".$self->$name);
    }

    $self->debug_message('Extract genotytpes...done');
    return 1;
}

sub resolve_source {
    my $self = shift;

    return 1 if $self->source;

    my @sources = grep { $self->$_ } (qw/ build model sample instrument_data /);
    if ( not @sources ) {
        $self->error_message('No source given! Can be build, model, instrument data or sample.');
        return;
    }

    my $source_method = '_resolve_source_for_'.$sources[0];
    my $resolve_source = $self->$source_method;
    return if not $resolve_source;

    return 1;
}

sub _resolve_source_for_build {
    my $self = shift;

    $self->source($self->build);
    $self->source_type('build');
    $self->sample( $self->build->subject );

    return 1;
}

sub _resolve_source_for_model {
    my $self = shift;

    my $build = $self->model->last_succeeded_build;
    if ( not $build ) {
        $self->instrument_data( $self->model->instrument_data );
        return $self->_resolve_source_for_instrument_data;
    }

    $self->source($build);
    $self->source_type('build');
    $self->sample( $build->subject );

    return 1;
}

sub _resolve_source_for_sample {
    my ($self) = @_;

    # Need variation list build
    my $variation_list_build = $self->variation_list_build;
    if ( not $variation_list_build ) {
        $self->error_message('Variation list build is required to get genotypes for a sample!');
        return;
    }

    # Get InstData
    my $sample = $self->sample;
    my $library = Genome::Library->get(name => $sample->name.'-microarraylib', sample => $sample);
        #'import_source_name in' => ( $self->use_external ) ? [qw/ BGI bgi Broad broad CSHL cshl external /] : [qw/ wugsc wugc wutgi tgi /],
    my @instrument_data = Genome::InstrumentData->get(library => $library);
    if ( not @instrument_data ) {
        $self->error_message('No microarray instrument data for sample!' .$sample->__display_name__);
        return;
    }

    # Restrict by priority
    my @filtered_instrument_data;
    for my $priority ( $self->sample_type_priority ) {
        for my $instrument_data ( @instrument_data ) {
            my $verification_method = '_is_instrument_data_'.$priority;
            push @filtered_instrument_data, $instrument_data if $self->$verification_method($instrument_data);
        }
    }

    if ( not @filtered_instrument_data ) {
        $self->error_message('No instrument data found matches the indicated priorities!');
        return;
    }

    # Maybe there is a build already
    my $build;
    for my $instrument_data ( @filtered_instrument_data ) {
        $build = $self->_last_succeeded_build_from_model_for_instrument_data($instrument_data);
        last if $build;
    }

    if ( $build ) {
        $self->source($build);
        $self->source_type('build');
        $self->sample( $build->subject);
    }
    else {
        $self->source($filtered_instrument_data[$#filtered_instrument_data]);
        $self->source_type('instrument_data');
        $self->sample( $filtered_instrument_data[$#filtered_instrument_data]->sample );
    }

    return 1;
}

sub _resolve_source_for_instrument_data {
    my $self = shift;

    my $variation_list_build = $self->variation_list_build;
    if ( not $variation_list_build ) {
        $self->error_message('Variation list build is required to get genotypes for an instrument data!');
        return;
    }

    # Maybe there is a build already
    my $instrument_data = $self->instrument_data;
    my $build = $self->_last_succeeded_build_from_model_for_instrument_data($instrument_data);
    if ( $build ) {
        $self->source($build);
        $self->source_type('build');
        $self->sample( $build->subject);
    }
    else {
        $self->source($instrument_data);
        $self->source_type('instrument_data');
        $self->sample( $instrument_data->sample );
    }

    return 1;
}

sub _is_instrument_data_default {
    my ($self, $instrument_data) = @_;

    return 1 if $instrument_data->id eq $instrument_data->sample->default_genotype_data_id;
    return;
}

my @internal_source_names = (qw/ wugsc wugc wutgi tgi /);
sub _is_instrument_data_internal {
    my ($self, $instrument_data) = @_;

    for my $internal_source_name ( @internal_source_names ) {
        return 1 if $instrument_data->import_source_name =~ /^$internal_source_name$/i;
    }

    return;
}

sub _is_instrument_data_external {
    my ($self, $instrument_data) = @_;

    for my $internal_source_name ( @internal_source_names ) {
        return 1 if $instrument_data->import_source_name !~ /^$internal_source_name$/i;
    }

    return;
}

sub _last_succeeded_build_from_model_for_instrument_data {
    my ($self, $instrument_data) = @_;

    my $variation_list_build = $self->variation_list_build;
    my @builds = sort { $b->date_completed cmp $a->date_completed } Genome::Model::Build::GenotypeMicroarray->get(
        instrument_data => [ $instrument_data ],
        reference_sequence_build => $variation_list_build->reference,
        dbsnp_build => $variation_list_build,
        status => 'Succeeded',
    );

    return $builds[0];
}
#<>#

#< FILTERS >#
sub _create_filters {
    my $self = shift;

    my @filters;
    for my $filter_string ( $self->filters ) {
        $self->debug_message('Filter: '.$filter_string);
        my ($name, $config) = split(':', $filter_string, 2);
        my %params;
        %params = map { split('=') } split(':', $config) if $config;
        my $filter_class = 'Genome::Model::GenotypeMicroarray::Filter::By'.Genome::Utility::Text::string_to_camel_case($name);
        my $filter = $filter_class->create(%params);
        if ( not $filter ) {
            $self->error_message("Failed to create fitler for $filter_string");
            return;
        }
        push @filters, $filter;
    }

    return \@filters;
}
#<>#

#< OPEN WRITER >#
sub _open_writer {
    my $self  = shift;
    $self->status_message('Open writer...');

    my $output = $self->output;
    $self->status_message('Output file: '.$output);
    my $output_format = $self->output_format;
    $self->status_message('Output format: '.$output_format);
    my $writer;
    if ( $output_format eq 'vcf' ) {
        die "VCF not implmented!";
    }
    else {
        my %writer_params = ( output => $output );
        $writer_params{ headers} = [ $self->fields ] if $output_format =~ /\-headers$/;
        $writer = Genome::Utility::IO::Writer->create(%writer_params);
    }

    $self->status_message('Open writer...done');
    return $writer;
}
#<>#

1;

