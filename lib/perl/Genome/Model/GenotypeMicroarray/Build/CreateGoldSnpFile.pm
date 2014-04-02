package Genome::Model::GenotypeMicroarray::Build::CreateGoldSnpFile;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::Build::CreateGoldSnpFile {
    is => 'Command::V2',
    has_input_output => {
        build => { is => 'Genome::Model::Build::GenotypeMicroarray', },
    },
};

sub execute {
    my $self = shift;
    $self->debug_message("Create gold SNP array file...");

    # TODO bdericks: I'm guessing that second genotype file is supposed to be the replicate. It should be changed
    # to be the actual replicate when we know how to figure it out.
    # abrummet: This is the only place in the tree where this Command is used.  I've stripped out the second input
    # file to fix a bug where it would not read from the "second" file when switching chromosomes and the next position
    # is numerically higher than the last position
    my $build = $self->build;
    my $genotype_file = $build->genotype_file_path;
    $self->debug_message("Genotype file: ".$genotype_file);
    my $snp_array_file = $build->formatted_genotype_file_path;
    $self->debug_message("Gold SNP file: ".$snp_array_file);

    my $gold_snp = Genome::Model::GenotypeMicroarray::Command::CreateGoldSnpFileFromGenotypes->create(
        genotype_file => $genotype_file,
        output_file => $snp_array_file,
        reference_sequence_build => $build->reference_sequence_build, 
    );
    if ( not $gold_snp ) {
        $self->error_message("Cannot create gold snp tool.");
        return;
    }

    $gold_snp->dump_status_messages(1);

    if ( not $gold_snp->execute ) {
        $self->error_message("Cannot execute gold snp tool");
        return;
    }

    if ( not -s $snp_array_file ) {
        $self->error_message("Executed gold snp tool, but snp array file ($snp_array_file) does not exist");
        return;
    }

    $self->debug_message("Create gold SNP array file...OK");
    return 1;
}

1;

