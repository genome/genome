#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 7;
use above 'Genome';
use Genome::Model::ReferenceAlignment;
use Genome::Model::Build::ReferenceAlignment;

sub create_ref_align_model {
    my $model = {};
    bless $model, 'Genome::Model::ReferenceAlignment';
    return $model;
}

sub create_ref_align_build {
    my $build = {};
    bless $build, 'Genome::Model::Build::ReferenceAlignment';
    my $geno_build = create_default_gold_snp_build();
    $build->{gold_snp_build} = $geno_build;
    return $build;
}

sub Genome::Model::Build::ReferenceAlignment::genotype_microarray_build_id {
    my $self = shift;
    return $self->{genotype_microarray_build_id};
}

sub Genome::Model::Build::GenotypeMicroarray::id {
    my $self = shift;
    return $self->{id};
}

sub Genome::Model::Build::GenotypeMicroarray::class {
    my $self = shift;
    return $self->{class};
}

sub Genome::Model::Build::ReferenceAlignment::gold_snp_build {
    my $self = shift;
    return $self->{gold_snp_build};
}

sub create_default_gold_snp_build {
    my $build = {id => 2, class => 'Genome::Model::Build::GenotypeMicroarray'};
    bless $build, 'Genome::Model::Build::GenotypeMicroarray';
    return $build;
}

sub Genome::Model::Build::ReferenceAlignment::genotype_microarray_build {
    my $self = shift;
    my $genotype_microarray_build;
    if ($self->genotype_microarray_build_id) {
        $genotype_microarray_build = {id => $self->genotype_microarray_build_id};
        bless $genotype_microarray_build, 'Genome::Model::Build::GenotypeMicroarray';
    }
    return $genotype_microarray_build;
}

sub Genome::Model::Build::GenotypeMicroarray::formatted_genotype_file_path {
    my $self = shift;
    return 'file ' . $self->id;
}

sub Genome::Model::ReferenceAlignment::add_input {
    my $self = shift;
    my %params = @_;
    $self->{genotype_microarray_build_id} = $params{value_id};
}

{
    my $model = create_ref_align_model();
    my $build = create_ref_align_build();
    $build->{genotype_microarray_build_id} = undef;
    $build->{genotype_microarray_build} = undef;
    $build->{gold_snp_build} = undef;
    is($build->genotype_microarray_build_id, undef, 'genotype_microarray_build_id is not set');
    is($build->gold_snp_build, undef, 'gold_snp_build is not set');
    is($build->gold_snp_path, undef, 'got undef');
}

{
    my $model = create_ref_align_model();
    my $build = create_ref_align_build();
    $build->{genotype_microarray_build_id} = undef;
    is($build->genotype_microarray_build_id, undef, 'genotype_microarray_build_id is not set');
    is($build->gold_snp_path, 'file 2', 'got file 2');
}

{
    my $model = create_ref_align_model();
    my $build = create_ref_align_build();
    my $geno_build = $build->gold_snp_build;
    $geno_build->{id} = 1;
    $build->{gold_snp_build} = $geno_build;
    $build->{genotype_microarray_build_id} = '1';
    is($build->genotype_microarray_build_id, '1', 'genotype_microarray_build_id is set');
    is($build->gold_snp_path, 'file 1', 'got file 1');
}
