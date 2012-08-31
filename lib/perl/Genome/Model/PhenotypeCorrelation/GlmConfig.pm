package Genome::Model::PhenotypeCorrelation::GlmConfig;

use Genome;
use Carp qw/confess/;

use strict;
use warnings;

my @HEADER_FIELDS = qw( analysis_type clinical_data_trait_name variant/gene_name covariates memo );
my @SHORT_NAMES = qw( type attr_name var_gene_name covariates memo );

class Genome::Model::PhenotypeCorrelation::GlmConfig {
    is => "UR::Object",
    has_transient_optional => [
        _attribute_hash => {
            is => "HASH",
            doc => "Hash of hashrefs containing attribute information for the GLM model, keyed on attribute name",
            default_value => {},
        },
        _attribute_order => {
            is => "ARRAY",
            doc => "Original ordering of attributes.",
            default_value => [],
        },
    ],
};

sub from_file {
    my ($class, $file) = @_;
    return $class->from_filehandle(Genome::Sys->open_file_for_reading($file));
}

sub from_filehandle {
    my ($class, $fh) = @_;
    my $header = <$fh>;
    chomp $header;
    my $obj = Genome::Model::PhenotypeCorrelation::GlmConfig->create;
    my $expected_header = join("\t", @HEADER_FIELDS);
    if ($header ne $expected_header) {
        confess "Invalid header in glm config file:\n$header\nEXPECTED:\n$expected_header";
    }

    my $attributes = [];
    while (my $line = <$fh>) {
        $obj->_add_attribute_string($line)
    }

    return $obj;
}

sub from_string_arrayref {
    my ($class, $lines) = @_;
    my $obj = Genome::Model::PhenotypeCorrelation::GlmConfig->create;
    for my $line (@$lines) {
        $obj->_add_attribute_string($line);
    }
    return $obj;
}

sub to_string {
    my $self = shift;
    my $str =  join("\t", @HEADER_FIELDS) . "\n"
        . join("\n", map {$self->_attribute_to_line($_)} $self->attributes) . "\n";
    return $str;
}

sub to_file {
    my ($self, $file) = @_;
    return $self->to_filehandle(Genome::Sys->open_file_for_writing($file));
}

sub to_filehandle {
    my ($self, $fh) = @_;
    $fh->write($self->to_string);
}

sub to_arrayref {
    my $self = shift;
    return [map {$self->_attribute_to_line($_)} $self->attributes];
}

sub categorical_attributes {
    my $self = shift;
    return grep {$_->{type} eq 'B'} $self->attributes;
}

sub quantitative_attributes {
    my $self = shift;
    return grep {$_->{type} eq 'Q'} $self->attributes;
}

sub attributes {
    my ($self, @query) = @_;
    @query = @{$self->_attribute_order} unless @query;
    return @{$self->_attribute_hash}{@query};
}

sub _add_attribute_string {
    my ($self, $line) = @_;
    chomp $line;
    my @fields = split("\t", $line);
    my $missing = scalar @HEADER_FIELDS - scalar @fields;
    if ($missing < 0 || $missing > 1) {
        confess "Invalid entry in glm config file: $line";
    }
    push(@fields, ""x$missing) if $missing;
    
    my %entry;
    @entry{@SHORT_NAMES} = @fields;
    $entry{covariates} = [grep {$_ ne "NA"} split(/\s*\+\s*/, $entry{covariates})];
    push(@{$self->_attribute_order}, $entry{attr_name});
    $self->_attribute_hash->{$entry{attr_name}} = \%entry;
}

sub _attribute_to_line {
    my ($self, $attribute) = @_;
    my $str .= join("\t", @{$attribute}{qw(type attr_name var_gene_name)});
    my $cov = "NA";
    $cov = join("+", @{$attribute->{covariates}}) if scalar @{$attribute->{covariates}};
    my $memo = $attribute->{memo} || "";
    $str .= "\t$cov\t$memo";
    return $str;
}

1;
