package Genome::Model::Tools::Assembly::CreateOutputFiles::Msi;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Assembly::CreateOutputFiles::Msi {
    is => 'Genome::Model::Tools::Assembly::CreateOutputFiles',
    has => [
        acefile => {
            is => 'Text',
            doc => 'Ace file to get fasta and qual from',
        },
        assembly_directory => {
            is => 'Text',
            doc => 'Main Assembly directory, not assembly edit_dir',
        },
	assembler => {
	    is => 'Text',
	    doc => 'Assembler used to create the assembly',
	    valid_values => ['pcap', 'newbler', 'Velvet'],
	},
    ],
    has_optional_transient => [
	_contigs_bases_file      => { is => 'Text', doc => 'Assembly contigs.bases file'},
	_contigs_quals_file      => { is => 'Text', doc => 'Assembly contigs.quals file'},
	_gap_file                => { is => 'Text', doc => 'Assembly gap.txt file'},
	_reads_placed_file       => { is => 'Text', doc => 'Assembly reads.placed file'},
	_read_info_file          => { is => 'Text', doc => 'Assembly readinfo.txt file'},
	_supercontigs_agp_file   => { is => 'Text', doc => 'Assembly supercontigs.agp file'},
	_supercontigs_fasta_file => { is => 'Text', doc => 'Assembly supercontigs.fasta file'},
	_stats_file              => { is => 'Text', doc => 'Assembly stats.txt file'},
    ],
};

sub help_brief {
    'Tool to create assembly output files for msi assemblies';
}

sub help_detail {
"Tool to create assembly output files from modified newbler, pcap or
velvet assemblies.  It will need msi.gap.txt file to get gap info
for scaffold contigs.  If get size is not provided, it will use 100 bp
which is default value for unknown gap sizes for submission to ncbi"
}

sub execute {
    my $self = shift;

    #resolve acefile location
    my $acefile = $self->acefile;
    unless (-s $acefile) {
	$acefile = $self->assembly_directory."/edit_dir/$acefile";
	unless (-s $acefile) {
	    $self->error_message("Can't find ace file, neither $acefile nor ".$self->assembly_directory."/edit_dir/$acefile exists");
	    return;
	}
    }

    #set in/out put file names values
    unless ($self->_set_transients()) {
	$self->error_message("Failed to define and set assembly files");
	return;
    }

    #create sorted msi.contigs.bases and msi.contigs.qual files
    $self->debug_message("Creating msi.contigs.bases and msi.contigs.quals files");
    my $contigs = Genome::Model::Tools::Assembly::CreateOutputFiles::ContigsFromAce->create (
	acefile => $acefile,
	fasta_out => $self->_contigs_bases_file,
	qual_out => $self->_contigs_quals_file,
	directory => $self->assembly_directory,
	);
    unless ($contigs->execute) {
	$self->error_message("Failed to create contigs bases and qual files from acefile");
	return;
    }
    sleep(1);

    #create msi.readinfo.txt file
    $self->debug_message("Creating msi.readinfo.txt file");
    my $readinfo = Genome::Model::Tools::Assembly::CreateOutputFiles::ReadInfo->create (
	acefile => $acefile,
	output_file => $self->_read_info_file,
	directory => $self->assembly_directory,
	);
    unless ($readinfo->execute) {
	$self->error_message("Failed to create msi.readinfo.txt file");
	return;
    }

    #create msi.reads.placed file
    $self->debug_message("Creating msi.reads.placed file");
    my $reads_placed = Genome::Model::Tools::Assembly::CreateOutputFiles::ReadsPlaced->create (
	read_info_file => $self->_read_info_file,
	gap_file => $self->_gap_file,
	contigs_bases_file => $self->_contigs_bases_file,
	output_file => $self->_reads_placed_file,
	directory => $self->assembly_directory,
	);
    unless ($reads_placed->execute) {
	$self->error_message("Failed to create msi.reads.placed file");
	return;
    }
    sleep(1);

    #create supercontigs.agp file
    $self->debug_message("Creating msi.supercontigs.agp file");
    my $agp = Genome::Model::Tools::Assembly::CreateOutputFiles::SupercontigsAgp->create (
	contigs_bases_file => $self->_contigs_bases_file,
	gap_file => $self->_gap_file,
	output_file => $self->_supercontigs_agp_file,
	directory => $self->assembly_directory,
	);
    unless ($agp->execute) {
	$self->error_message("Failed to create msi.supercontigs.agp file");
	return;
    }

    #create msi.supercontigs.fasta
    $self->debug_message("Creating msi.supercontigs.fasta file");
    my $s_fa = Genome::Model::Tools::Assembly::CreateOutputFiles::SupercontigsFasta->create (
	contigs_bases_file => $self->_contigs_bases_file,
	gap_file => $self->_gap_file,
	output_file => $self->_supercontigs_fasta_file,
	directory => $self->assembly_directory,
	);
    unless ($s_fa->execute) {
	$self->error_message("Failed to create msi.supercontigs.fasta file");
	return;
    }

    #create msi.stats.txt file
    $self->debug_message("Creating msi.stats.txt file");
    my $class = 'Genome::Model::Tools::Assembly::Stats::'. ucfirst $self->assembler;
    my $stats = $class->create (
	assembly_directory => $self->assembly_directory,
	out_file => $self->_stats_file,
	msi_assembly => 1,
	);
    unless ($stats->execute) {
	$self->error_message("Failed to create msi.stats.txt file");
	return;
    }
    return 1;
}

sub _set_transients {
    my $self = shift;
    
    my $edit_dir = $self->assembly_directory.'/edit_dir';

    $self->_contigs_bases_file("$edit_dir/msi.contigs.bases");
    $self->_contigs_quals_file("$edit_dir/msi.contigs.quals");
    $self->_gap_file("$edit_dir/msi.gap.txt");
    $self->_reads_placed_file("$edit_dir/msi.reads.placed");
    $self->_read_info_file("$edit_dir/msi.readinfo.txt");
    $self->_supercontigs_agp_file("$edit_dir/msi.supercontigs.agp");
    $self->_supercontigs_fasta_file("$edit_dir/msi.supercontigs.fasta");
    $self->_stats_file("$edit_dir/msi.stats.txt");
    
    return 1;
}

1;
