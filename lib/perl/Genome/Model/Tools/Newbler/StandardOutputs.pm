package Genome::Model::Tools::Newbler::StandardOutputs;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Newbler::StandardOutputs {
    is => 'Genome::Model::Tools::Newbler',
    has => [
        assembly_directory => {
            is => 'Text',
            doc => 'Path to assembly',
        },
        min_contig_length => {
            is => 'Number',
            doc => 'Minimum contig length to export',
        },
        default_gap_size => {
            is => 'Number',
            doc => 'Gap size to assign when newbler does not assign one',
            is_optional => 1,
            default_value => 10,
        },
    ],
};

sub help_brief {
    'Tools to create standard output files for newbler assemblies';
}

sub help_detail {
    return <<"EOS"
gmt newbler standard-outputs --assembly-directory /gscmnt/111/newbler_assembly --min-contig-length 200 --default-gap-size 10
EOS
}

sub execute {
    my $self = shift;
    
    my %params = (
        assembly_directory => $self->assembly_directory,
        min_contig_length  => $self->min_contig_length,
    );
    $params{default_gap_size} = $self->default_gap_size if $self->default_gap_size;

    #create consed/edit_dir if not there
    unless ( -d $self->consed_edit_dir ) {
        $self->create_consed_dir;
    }

    #pcap style ace file
    my $ec_ace = Genome::Model::Tools::Newbler::ToPcapAce->create( %params );
    if ( not $ec_ace->execute ) {
        $self->error_message("Failed to execute newbler to-pcap-ace");
        return;
    }

    #contigs.bases and contigs.quals files
    my $ec_contigs = Genome::Model::Tools::Newbler::CreateContigsFiles->create( %params );
    if ( not $ec_contigs->execute ) {
        $self->error_message("Failed to execute newbler create-contigs-files");
        return;
    }

    #readinfo.txt and reads.placed files
    my $placed = Genome::Model::Tools::Newbler::CreatePlacedReadsFiles->create( %params );
    if ( not $placed->execute ) {
        $self->error_message("Failed to execute to create placed reads files");
        return;
    }

    #reads.unplaced and reads.unplaced.fasta files
    my $unplaced = Genome::Model::Tools::Newbler::CreateUnplacedReadsFiles->create( %params );
    if ( not $unplaced->execute ) {
        $self->error_message("Failed to execute to create unplaced reads files");
        return;
    }

    #supercontigs.fasta and supercontigs.agp files
    my $ec_sctgs = Genome::Model::Tools::Newbler::CreateSupercontigsFiles->create( %params );
    if ( not $ec_sctgs->execute ) {
        $self->error_message("Failed to execute newbler create-supercontigs-files");
        return;
    }

    return 1;
}

1;
