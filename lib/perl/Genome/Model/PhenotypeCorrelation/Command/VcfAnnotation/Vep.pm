package Genome::Model::PhenotypeCorrelation::Command::VcfAnnotation::Vep;

use Data::Dumper;
use File::Basename qw/basename/;
use Genome;
use Carp qw/confess/;
use File::Path qw/mkpath/;

use strict;
use warnings;

my %VEP_ENSEMBL_VERSIONS = (
    # vep version -> ensembl version
    "2_7" => "69_37n_v3",
);

class Genome::Model::PhenotypeCorrelation::Command::VcfAnnotation::Vep {
    is => "Genome::Model::PhenotypeCorrelation::Command::VcfAnnotation::Base",
    has_input => [
        params => {
            is => "HASH",
            doc => "Parsed parameters from strategy",
        },

    ]
};

sub _resolve_annotation_build {
    my ($self, $vep_version) = @_;

    confess "Unable to determine correct ensembl version for vep version ($vep_version)"
        unless exists $VEP_ENSEMBL_VERSIONS{$vep_version};

    my $ensembl_version = $VEP_ENSEMBL_VERSIONS{$vep_version};

    my $species = $self->species_name;

    my @models = Genome::Model::ImportedAnnotation->get(
        annotation_source => "ensembl",
        species_name => $species
        );

    $self->status_message("Annotation models:\n\t" . join("\n\t", map {$_->__display_name__} @models));
    my @builds = grep {defined} map {$_->build_by_version($ensembl_version)} @models;
    if (!@builds) {
        confess "Failed to resolve annotation build for species $species, ensembl version $ensembl_version";
    }

    $self->status_message("Annotation builds:\n\t" . join("\n\t", map {$_->name} @builds));
    if (@builds > 1) {
        confess "Multiple annotation builds found species $species, ensembl version $ensembl_version";
    }

    $self->status_message(sprintf("Chose annotation build %s for %s, ensembl version %s",
        $builds[0]->name, $species, $ensembl_version));
    return $builds[0];
}

sub execute {
    my $self = shift;
    my %input_params = %{$self->params};
    my $condel = delete $input_params{condel};
    my $version = delete $input_params{version} || "2_7";
    my $build = $self->_resolve_annotation_build($version);
    my %params = (
        %input_params,
        ensembl_annotation_build_id => $build->id,
        input_file => $self->input_file,
        output_file => $self->output_file,
        format => "vcf", # input format
        quiet => 1, # shhhh
        vcf => 1, # sets output format to vcf.
        version => $version,
        );

    if ($condel) {
        if (exists $params{"plugins"}) {
            confess "The options 'condel' and 'plugins' may not be used together";
        }
        push @{$params{"plugins"}}, "Condel,PLUGIN_DIR,$condel,2";
    }

    $self->status_message("Running vep with parameters: " . Dumper(\%params));
    my $vep_command = Genome::Db::Ensembl::Command::Vep->create(%params);

    return $vep_command->execute();
}

1;
