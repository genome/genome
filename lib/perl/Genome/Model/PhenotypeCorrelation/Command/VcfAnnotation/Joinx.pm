package Genome::Model::PhenotypeCorrelation::Command::VcfAnnotation::Joinx;

use Carp qw/confess/;
use Data::Dumper;
use Genome;

use strict;
use warnings;

my $target_package = "Genome::Model::Tools::Joinx::VcfAnnotateMulti";

class Genome::Model::PhenotypeCorrelation::Command::VcfAnnotation::Joinx {
    is => "Genome::Model::PhenotypeCorrelation::Command::VcfAnnotation::Base",
    has_input => [
        params => {
            is => "ARRAY",
            doc => "Parsed parameters from strategy",
        },
    ]
};

sub create {
    my $class = shift;
    return $class->SUPER::create(@_);
}

sub _parse_info_fields {
    return split(";", shift);
}

sub _resolve_annotation_source {
    my $self = shift;
    my $raw_spec = shift;

    if (exists $raw_spec->{source_name} and exists $raw_spec->{source_version}) {
        my $source = delete $raw_spec->{source_name};
        my $version = delete $raw_spec->{source_version};

        my @models = Genome::Model::ImportedVariationList->get(source_name => $source);
        my @builds =
            grep {$_->version eq $version}
            grep {$_->status eq "Succeeded"}
            map {$_->builds} @models;

        if (!@builds) {
            confess "No ImportedVariationList builds found for "
                . "source_name: $source, version: $version";
        }
        if (@builds != 1) {
            confess "Multiple ImportedVariationList builds found for "
                . " source_name: $source, version: $version. Build ids:\n\t"
                . join("\n\t", map {$_->id} @builds);
        }

        return $builds[0]->snvs_vcf;
    }
    elsif (!exists $raw_spec->{source_build}) {
        confess "Must specify either (source_name, version_name), or source_build: " . Dumper($raw_spec);
    }

    return Genome::Model::Build->get(delete $raw_spec->{source_build})->snvs_vcf;
}

sub _parse_spec {
    my ($self, $raw_spec) = @_;
    my $annotation_file = $self->_resolve_annotation_source($raw_spec);
    confess "Unable to resolve annotation source from spec: " . Dumper($raw_spec) unless $annotation_file;
    my %params = (
        annotation_file => $annotation_file,
    );

    if (exists $raw_spec->{info_fields}) {
        $params{info_fields} = [_parse_info_fields(delete $raw_spec->{info_fields})];
    }

    if (exists $raw_spec->{identifiers}) {
        $params{identifiers} = delete $raw_spec->{identifiers}
    }

    if (%$raw_spec) {
        confess "Unknown params in annotation spec: " . Dumper($raw_spec);
    }

    my $spec = Genome::Model::Tools::Joinx::VcfAnnotationSpec->create(%params);
    return $spec;
}

sub command_class_and_params {
    my $self = shift;

    my @specs = map {$self->_parse_spec($_)} @{$self->params};
    my $input_vcf = $self->input_file;

    return ($target_package,
        input_file => $input_vcf,
        annotation_specs => \@specs,
        output_file => $self->output_file
        );
}

sub _build_command {
    my $self = shift;
    my ($cls, %params) = $self->command_class_and_params;
    return $cls->create(%params);
}

sub execute {
    my $self = shift;
    my $command = $self->_build_command;
    return $command->execute;
}

1;
