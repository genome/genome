package Genome::Model::GenotypeMicroarray::Command::ExtractToCsv;

use strict;
use warnings;

use Genome;

require File::Basename;
use Genome::Model::GenotypeMicroarray::GenotypeFile::ReaderForBuildOriginalTsv

class Genome::Model::GenotypeMicroarray::Command::ExtractToCsv {
    is => 'Command::V2',
    has => {
        build => {
            is => 'Genome::Model::Build::GenotypeMicroarray',
            doc => 'The genotype build to use.',
        },
    },
    has_optional => {
        output => {
            is => 'Text',
            default_value => '-',
            doc => 'Output file. Default is STDOUT.',
        },
        fields => {
            is => 'Text',
            is_many => 1,
            valid_values => [qw/ chromosome position alleles reference id sample_name log_r_ratio gc_score cnv_value cnv_confidence allele1 allele2 / ],
            default_value => [qw/ chromosome position alleles reference id sample_name log_r_ratio gc_score cnv_value cnv_confidence allele1 allele2 / ],
            doc => 'The fields to include in the output CSV output.',
        },
        separator => {
            is => 'Text',
            default_value => "\t",
            doc => 'Separator to use. Default is tabs.'
        },
        headers => {
            is => 'Boolean',
            default_value => 1,
            doc => 'Include headers in the output.',
        },
        filters => {
            is => 'Text',
            is_many => 1,
            doc => "Filter genotypes. Give name and parameters, if required. Filters:\n gc_scrore => filter by min gc score (Ex: gc_score:min=0.7)\n invalid_iscan_ids => list of invalid iscan snvs compiled by Nate",
        },
    },
    has_optional_transient => {
        metrics => { is => 'Hash', default_value => { input => 0, output => 0, }, },
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
    $self->debug_message('Extract genotypes to CSV...');

    my $filters = $self->_create_filters;
    return if not $filters;

    my $build = $self->build;
    my $reader = Genome::Model::GenotypeMicroarray::GenotypeFile::ReaderForBuildOriginalTsv->create(
        build => $build,
    );
    return if not $reader;

    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->output,
        headers => [$self->fields],
        separator => $self->separator,
        in_place_of_null_value => 'NA',
        print_headers => $self->headers,
        ignore_extra_columns => 1,
    );
    return if not $writer;

    my $write_sub = sub { $writer->write($_[0]); };
    if ( grep { $_ eq 'reference' } $self->fields and not grep { $_ eq 'reference' } @{$reader->{headers}} ) {
        # Add the refence base if requested to legacy build files
        my $reference_sequence_build = $build->reference_sequence_build;
        $write_sub = sub{
            $_[0]->{reference} = $reference_sequence_build->sequence(
                $_[0]->{chromosome}, $_[0]->{position}, $_[0]->{position},
            );
            $writer->write($_[0]);
        };
    }

    my %metrics = ( input => 0, filtered => 0, output => 0, );
    GENOTYPE: while ( my $genotype = $reader->read ) {
        $metrics{input}++;
        for my $filter ( @$filters ) {
            next GENOTYPE if not $filter->filter($genotype);
        }
        $metrics{output}++;
        $write_sub->($genotype);
    }
    $self->metrics(\%metrics);
    for my $name ( map { 'genotypes_'.$_ } (qw/ input filtered output /) ) {
        $self->debug_message(ucfirst(join(' ', split('_', $name))).": ".$self->$name);
    }

    $self->debug_message('Extract genotypes to CSV...done');
    return 1;
}

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
            $self->error_message("Failed to create filter for $filter_string");
            return;
        }
        push @filters, $filter;
    }

    return \@filters;
}



1;

