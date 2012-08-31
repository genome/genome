package Genome::Model::GenotypeMicroarray::Command::CompareToReference;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::Command::CompareToReference {
    is => 'Genome::Command::Base',
    doc => 'Compare reference bases in a genotype file of a given genotype microarray build to a reference sequence',
    has_input => [
        build_id => {
            is => 'Number',
            doc => 'ID of genotype microarray build',
        },
        build => {
            is => 'Genome::Model::Build::GenotypeMicroarray',
            id_by => 'build_id',
            doc => 'GenotypeMicroarray build',
            shell_args_position => 1,
        },
        reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            id_by => 'reference_id',
            doc => 'The reference sequence build id to use',
            shell_args_position => 2,
        },
        reference_id => {
            is => 'Integer',
        },
        max_base_count => {
            is => 'Number',
            doc => 'Maximum number of bases to compare',
            default => 10000,
        },
    ],
};

sub execute {
    my $self = shift;
    my $input_file = $self->build->formatted_genotype_file_path;
    my $fh = new IO::File("<$input_file") or die "Failed to open genotype file $input_file for reading.";
    my $hits = 0;
    my $total = 0;
    print "Comparing up to " . $self->max_base_count . " reference bases...\n";
    while (<$fh>) {
        next unless /ref\tref\tref\tref/;
        my ($chr, $start, $stop, $call, @extra) = split("\t");
        my $allele = $self->reference->sequence($chr, $start, $stop);
        ++$hits if $call eq $allele;
        last if (++$total == $self->max_base_count);    
    }
    $fh->close();

    print "Checked $total reference bases from GenotypeMicroarray build " . $self->build->__display_name__ .
        " against reference " . $self->reference->__display_name__ ."\n";
    printf "$hits / $total matches, %.02f%% match\n", (100.0*$hits/$total);
    return 1;
}

1;
