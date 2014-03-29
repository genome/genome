package Genome::Model::GenotypeMicroarray::Build::CreateGoldSnpBedFile;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::Build::CreateGoldSnpBedFile {
    is => 'Command::V2',
    has => {
        build => {
            is => 'Genome::Model::Build::GenotypeMicroarray',
            doc => 'Build to operate on.',
        },
    },
};

sub execute {
    my $self = shift;
    $self->debug_message('Create gold snp bed file...');

    my $build = $self->build;
    my $snp_array_file = $build->formatted_genotype_file_path;
    $self->debug_message('Formatted genotpye file: '.$snp_array_file);
    my $snvs_bed = $build->snvs_bed;
    $self->debug_message('Gold snp bed file: '.$snvs_bed);
    my $gold_snp_bed = Genome::Model::GenotypeMicroarray::Command::CreateGoldSnpBed->create(
        input_file => $snp_array_file,
        output_file => $snvs_bed,
        reference => $build->reference_sequence_build,
    );
    if ( not $gold_snp_bed ) {
        $self->error_message('Failed to create gold snp bed tool!');
        return;
    }

    $gold_snp_bed->dump_status_messages(1);

    unless ($gold_snp_bed->execute) {
        $self->error_message("Could not generate gold snp bed file at $snvs_bed from snp array file $snp_array_file");
        return;
    }

    if ( not -s $snvs_bed ) {
        $self->error_message("Executed 'create gold snp bed', but snvs bed file ($snvs_bed) does not exist");
        return;
    }

    $self->debug_message("Create gold snp bed file...OK");
    return 1;
}

1;

