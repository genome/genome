package Genome::Model::Tools::Velvet::StandardOutputs;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Model::Tools::Velvet::StandardOutputs {
    is => 'Genome::Model::Tools::Velvet::Base',
    has => [
	assembly_directory => {
	    is => 'Text',
	    doc => 'Directory where the assembly is located',
	},
        min_contig_length => {
            is => 'Number',
            doc => 'Minimum contig length to consider for post assembly files',
        }
    ],
};


sub help_brief {
    "Tool to create default post assembly output files for velvet assembly";
}

sub help_detail {
    return <<"EOS"
gmt velvet standard-outputs --assembly-directory /gscmnt/111/velvet_assembly
EOS
}

sub execute {
    my $self = shift;

    unless ( -d $self->assembly_directory ) {
	$self->error_message("Failed to find assembly directory: ".$self->assembly_directory);
	return;
    }

    unless ( $self->create_edit_dir ) {
	$self->error_message("Assembly edit_dir does not exist and could not create one");
	return;
    }

    my %params = (
        assembly_directory => $self->assembly_directory,
        min_contig_length => $self->min_contig_length,
    );

    #create gap.txt file
    $self->status_message("Creating gap.txt file");
    my $gap = Genome::Model::Tools::Velvet::CreateGapFile->create( %params );
    unless ($gap->execute) {
        $self->error_message("Execute failed to to create gap.txt file");
        return;
    }
    $self->status_message("Completed creating gap.txt file");


    #create contigs.bases and contigs.quals files
    $self->status_message("Creating contigs.bases and contigs.quals files");
    my $contigs = Genome::Model::Tools::Velvet::CreateContigsFiles->create ( %params );
    unless ($contigs->execute) {
	$self->error_message("Failed to execute creating contigs.bases and quals files");
	return;
    }
    $self->status_message("Completed creating contigs.bases and contigs.qual files");
    

    #create reads.placed and readinfo.txt files
    $self->status_message("Creating reads.placed and readinfo files");
    my $reads = Genome::Model::Tools::Velvet::CreateReadsFiles->create( %params );
    unless ($reads->execute) {
	$self->error_message("Failed to execute creating reads files");
	return;
    }
    $self->status_message("Completed creating reads.placed and readinfo files");


    #create reads.unplaced and reads.unplaced.fasta files
    $self->status_message("Creating reads.unplaced and reads.unplaced.fasta files");
    my $unplaced = Genome::Model::Tools::Velvet::CreateUnplacedReadsFiles->create( %params );
    unless ($unplaced->execute) {
	$self->error_message("Failed to execute creating reads.unplaced files");
	return;
    }
    $self->status_message("Completed creating reads.unplaced and reads.unplaced.fasta files");


    #create supercontigs.fasta and supercontigs.agp file
    $self->status_message("Creating supercontigs fasta and agp files");
    my $supercontigs = Genome::Model::Tools::Velvet::CreateSupercontigsFiles->create( %params );
    unless ($supercontigs->execute) {
	$self->error_message("Failed execute creating of supercontigs files");
	return;
    }
    $self->status_message("Completed creating supercontigs.fasta and agp files");

    return 1;
}

1;
