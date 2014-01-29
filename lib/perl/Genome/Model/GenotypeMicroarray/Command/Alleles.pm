package Genome::Model::GenotypeMicroarray::Command::Alleles;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::Command::Alleles {
    is => 'Command::V2',
    has_optional => [
        model => {
            is => 'Genome::Model',
            doc => 'The genotype model to work with. This will get the most recent succeeded build.',
        },
        build => {
            is => 'Genome::Model::Build',
            doc => 'The genotype build to use.',
        },
    ],
};

sub help_brief { return 'Show the alleles counts in the original genotype file.'; }
sub hlep_detail { return help_brief(); }

sub execute {
    my $self = shift;
    $self->debug_message('Genotype alleles...');

    my $build = $self->_resolve_build;
    return if not $build;

    my $original_genotype_file_path = $build->original_genotype_file_path;
    if ( not $original_genotype_file_path or not -s $original_genotype_file_path ) {
        $self->error_message('Original genotype file does not exist! '.$original_genotype_file_path);
        return;
    }

    my $reader = Genome::Model::GenotypeMicroarray::OriginalGenotypeReader->create(
        input => $build->original_genotype_file_path,
    );
    if ( not $reader ) {
        $self->error_message('Failed to create original genotype file reader.');
        return;
    }

    my %alleles;
    while ( my $genotype = $reader->read ) {
        $alleles{ $genotype->{alleles} }++;
    }

    $self->debug_message('Alleles and instances:');
    for my $alleles ( sort keys %alleles ) {
        print join(' ', $alleles, $alleles{$alleles})."\n";
    }
    $self->debug_message('Total genotypes: '.$reader->total);
    
    $self->debug_message('Done');
    return 1;
};

sub _resolve_build {
    my $self = shift;

    if ( $self->build ) {
        return $self->build;
    }

    my $model = $self->model;
    if ( not $model ) {
        $self->error_message('No model or build given!');
        return;
    }

    my $build = $model->last_succeeded_build;
    if ( not $build ) {
        $build = $model->latest_build;
        if ( not $build ) {
            $self->error_message('No last succeeded or latest build for model! '.$model->id);
            return;
        }
    }

    return $self->build($build);
}

1;

