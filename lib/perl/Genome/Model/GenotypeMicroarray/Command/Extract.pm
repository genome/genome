package Genome::Model::GenotypeMicroarray::Command::Extract;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::Command::Extract {
    is => 'Command::V2',
    has_optional => [
        output => {
            is => 'Text',
            default_value => '-',
            doc => 'The output. Defaults to STDOUT.',
        },
        fields => {
            is => 'Text',
            is_many => 1,
            default_value => [qw/ chromosome position alleles /],
            valid_values => [qw/ chromosome position alleles id sample_id log_r_ratio gc_score cnv_value cnv_confidence allele1 allele2 /],
            doc => 'The fields to output in the genotype file.',
        },
        separator => {
            is => 'Text',
            default_value => 'tab',
            doc => 'Field separator of the output. Use "tab" for tab delineated.',
        },
        model => {
            is => 'Genome::Model',
            doc => 'The genotype model to work with. This will get the most recent succeeded build.',
        },
        build => {
            is => 'Genome::Model::Build',
            doc => 'The genotype build to use.',
        },
        filters => {
            is => 'Text',
            is_many => 1,
            doc => "Filter genotypes. Give name and parameters, if required. Filters:\n gc_scrore => filter by min gc score (Ex: gc_score:min=0.7)\n invalid_iscan_ids => list of invalid iscan snvs compiled by Nate",
        },
    ],
};

sub help_brief {
    return 'extract genotype data from a build';
}

sub help_detail {
    return <<HELP;
HELP
}

sub execute {
    my $self = shift;
    $self->debug_message('Extract genotytpes from build...');

    my $build = $self->_resolve_build;
    return if not $build;

    my $reader = $self->_open_oringinal_genotype_reader;
    return if not $reader;

    my $output_fh = $self->_open_output;
    return if not $output_fh;

    my $sep = ( $self->separator eq 'tab' ? "\t" : $self->separator );
    my @fields = $self->fields;
    my ($total, $pass) = (qw/ 0 0 /);
    GENOTYPE: while ( my $genotype = $reader->read ) {
        $output_fh->print( join($sep, map { defined $genotype->{$_} ? $genotype->{$_} : 'NA' } @fields)."\n" );
    }
    $output_fh->flush;

    $self->debug_message('Wrote '.$reader->passed.' of '.$reader->total.' genotypes. Output genotypes...OK');
    $self->debug_message('Done');
    return 1;
}

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

    my $last_succeeded_build = $model->last_succeeded_build;
    if ( not $last_succeeded_build ) {
        $self->error_message('No last succeeded build for model! '.$model->id);
        return;
    }

    return $self->build($last_succeeded_build);
}

sub _open_oringinal_genotype_reader {
    my $self = shift;
    $self->debug_message('Open original genotype file...');

    my $original_genotype_file = $self->build->original_genotype_file_path;
    $self->debug_message('Original genotype file: '.$original_genotype_file);
    if ( not -s $original_genotype_file ) {
        $self->error_message('Original genotype file does not exist!');
        return;
    }

    my $reader = Genome::Model::GenotypeMicroarray::OriginalGenotypeReader->create(
        input => $self->build->original_genotype_file_path,
    );
    if ( not $reader ) {
        $self->error_message('Failed to create original genotype file reader.');
        return;
    }
    $self->debug_message('Found headers in genotype file: '. join(', ', $reader->headers));

    $self->debug_message('Open original genotype reader...OK');
    return $reader;
}

sub _open_output {
    my $self = shift;

    $self->debug_message('Open output file...');

    my $output = $self->output;
    unlink $output if -e $output;
    $self->debug_message('Output file: '.$output);
    my $output_fh = eval{ Genome::Sys->open_file_for_writing($output); };
    if ( not $output_fh ) {
        $self->error_message("Failed to open output file ($output): $@");
        return;
    }

    $self->debug_message('Open output file...OK');
    return $output_fh;
}


1;

