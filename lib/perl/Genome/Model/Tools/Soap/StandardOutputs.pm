package Genome::Model::Tools::Soap::StandardOutputs;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Model::Tools::Soap::StandardOutputs {
    is => 'Genome::Model::Tools::Soap::Base',
    has => [
	assembly_directory => {
	    is => 'Text',
	    doc => 'Input assembly directory',
	},
        min_contig_length => {
            is => 'Number',
            doc => 'Minimum contig size to output',
            is_optional => 1,
            default_value => 50,
        }
    ],
};

sub help_brief {
    'Tool to create default post assembly output files, including contigs.bases, supercontigs.fasta, supercontigs.agp and stats.txt files';
}

sub help_detail {
    return <<"EOS"
gmt soap standard-outputs --assembly-directory /gscmnt/111/soap_assembly --min-contig-length 200
EOS
}

sub execute {
    my $self = shift;

    my %params = (
        assembly_directory => $self->assembly_directory,
        min_contig_length => $self->min_contig_length,
    );

    #create contigs.bases files
    $self->debug_message("Creating contigs fasta file");
    my $contigs = Genome::Model::Tools::Soap::CreateContigsBasesFile->create( %params );
    unless ($contigs->execute) {
        $self->error_message("Failed to successfully execute creating contigs fasta file");
        return;
    }
    $self->debug_message("Finished creating contigs fasta file");

    #create supercontigs fasta file
    $self->debug_message("Creating supercontigs fasta file");
    my $supercontigs = Genome::Model::Tools::Soap::CreateSupercontigsFastaFile->create( %params );
    unless ($supercontigs->execute) {
        $self->error_message("Failed to successfully execute creating scaffolds fasta file");
        return;
    }
    $self->debug_message("Finished creating scaffolds fasta file");

    #create supercontigs agp file
    $self->debug_message("Creating supercontigs agp file");
    my $agp = Genome::Model::Tools::Soap::CreateSupercontigsAgpFile->create( %params );
    unless ($agp->execute) {
        $self->error_message("Failed to successfully execute creating agp file");
        return;
    }
    $self->debug_message("Finished creating agp file");
    
    return 1;
}

1;
