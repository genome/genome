package Genome::Model::GenotypeMicroarray::GenotypeFile::ReaderForOldBuildFile;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::GenotypeFile::ReaderForOldBuildFile { 
    is => 'Genome::Utility::IO::SeparatedValueReader',
    has => [
        separator => { default => "\t", },
        total => { is => 'Number', is_transient => 1, default => 0, },
        filtered => { is => 'Number', calculate => q| return $self->total - $self->passed; |, },
        passed => { is => 'Number', is_transient => 1, default => 0, },
    ],
    has_optional => [
        filters => {
            is => 'Text',
            is_many => 1,
            doc => "Filter genotypes. Give name and parameters, if required. Filters:\n gc_scrore => filter by min gc score (Ex: gc_score:min=0.7)\n invalid_iscan_ids => list of invalid iscan snvs compiled by Nate",
        },
        _filters => { is => 'Array', is_transient => 1 },
    ]
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $filters = $self->_create_filters;
    return if not $filters;

    return $self;
}

sub _create_filters {
    my $self = shift;
    $self->debug_message('Create filters...');

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
    $self->_filters(\@filters);

    $self->debug_message('Create '.@filters.' filter(s)...OK');
    return 1;
}

sub read {
    my $self = shift;

    GENOTYPE: while ( my $genotype = $self->next ) {
        return if not $genotype;
        $self->total( $self->total + 1 );
        $genotype->{id} = $genotype->{snp_name};
        for my $filter ( @{$self->_filters} ) {
            next GENOTYPE if not $filter->filter($genotype);
        }
        $self->passed( $self->passed + 1 );
        $genotype->{alleles} = $genotype->{allele1}.$genotype->{allele2};
        return $genotype;
    }

    return;
}

1;

