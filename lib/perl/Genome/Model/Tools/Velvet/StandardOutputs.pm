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
        },
        sequencing_platform => {
            is => 'Text',
            is_optional => 1,
            is_many => 1,
            doc => 'Technology used to sequence data, solexa, 454, sanger, etc',
        },
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
    $self->debug_message("Creating gap.txt file");
    my $gap = Genome::Model::Tools::Velvet::CreateGapFile->create( %params );
    unless ($gap->execute) {
        $self->error_message("Execute failed to to create gap.txt file");
        return;
    }
    $self->debug_message("Completed creating gap.txt file");


    #create contigs.bases and contigs.quals files
    $self->debug_message("Creating contigs.bases and contigs.quals files");
    my $contigs = Genome::Model::Tools::Velvet::CreateContigsFiles->create ( %params );
    unless ($contigs->execute) {
	$self->error_message("Failed to execute creating contigs.bases and quals files");
	return;
    }
    $self->debug_message("Completed creating contigs.bases and contigs.qual files");
    

    #create reads.placed and readinfo.txt files
    $self->debug_message("Creating reads.placed and readinfo files");
    my $reads = Genome::Model::Tools::Velvet::CreateReadsFiles->create( %params );
    unless ($reads->execute) {
	$self->error_message("Failed to execute creating reads files");
	return;
    }
    $self->debug_message("Completed creating reads.placed and readinfo files");


    #create reads.unplaced and reads.unplaced.fasta files
    $self->debug_message("Creating reads.unplaced and reads.unplaced.fasta files");
    my $unplaced = Genome::Model::Tools::Velvet::CreateUnplacedReadsFiles->create( %params );
    unless ($unplaced->execute) {
	$self->error_message("Failed to execute creating reads.unplaced files");
	return;
    }
    $self->debug_message("Completed creating reads.unplaced and reads.unplaced.fasta files");


    #create supercontigs.fasta and supercontigs.agp file
    $self->debug_message("Creating supercontigs fasta and agp files");
    my $supercontigs = Genome::Model::Tools::Velvet::CreateSupercontigsFiles->create( %params );
    unless ($supercontigs->execute) {
	$self->error_message("Failed execute creating of supercontigs files");
	return;
    }
    $self->debug_message("Completed creating supercontigs.fasta and agp files");


    # create contigs.cmt file
    $self->debug_message("Creating contigs.cmt file");
    my $assembler_version = ( $self->version )
        ? $self->version
        : 'Unknown' ;
    my @sequencing_platform = ( $self->sequencing_platform )
        ? $self->sequencing_platform
        : 'unknown' ;

    my $cmt_file = Genome::Model::Tools::Velvet::CmtFile->create(
        %params,
        version => $assembler_version,
        sequencing_technologies => \@sequencing_platform,
    );
    if( not $cmt_file->execute ) {
        $self->error_message("Failed to execute creating of contigs.cmt file");
        return;
    }
    $self->debug_message("Completed creating contigs.cmt file");

    return 1;
}

1;
