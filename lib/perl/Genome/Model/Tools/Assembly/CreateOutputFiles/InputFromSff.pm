package Genome::Model::Tools::Assembly::CreateOutputFiles::InputFromSff;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Assembly::CreateOutputFiles::InputFromSff {
    is => 'Genome::Model::Tools::Assembly::CreateOutputFiles',
    has => [
	bin_path => {
	    is => 'Text',
	    doc => 'Path to sffinfo program to use',
	    is_optional => 1,
	},
	directory => {
	    is => 'Text',
	    doc => 'Main assembly directory above edit_dir',
	},
    ],
};

sub help_brief {
    'Tool to create fasta and qual files from sff files';
}

sub help_detail {
    "Tool to create fasta and qual files from sff files";
}

sub execute {
    my $self = shift;

    #needs to run from 64 bit machine
    my $archos = `uname -a`;
    unless ($archos =~ /64/) {
	$self->error_message("Create input files from sff must be run from 64 bit machine");
	return;
    }

    #validate directory
    if ($self->directory =~ /edit_dir/) {
	#for consistancy with other tools
	$self->error_message("Directory should be main assembly directory and not edit_dir");
	return;
    }
    if (! -d $self->directory.'/edit_dir') {
	my $edit_dir = Genome::Sys->create_directory($self->directory.'/edit_dir');
	unless (-d $edit_dir) {
	    $self->error_message("Failed to create assembly edit_dir");
	    return;
	}
    }

    #figure out what version of sffinfo to use
    my $bin_path = ($self->bin_path) ? $self->bin_path : $ENV{GENOME_SW} . '/454/DataProcessing-2.6';
    my $sff_info = $bin_path.'/bin/sffinfo';
    unless (-s $sff_info) {
	$self->error_message("Failed to find sffinfo in bin: ".$bin_path.'/bin');
	return;
    }

    #get sff files from sff file
    #TODO - maybe make sff dir an input just so it's absolutely clear that it has to be there
    my $sff_dir;
    if (-d $self->directory.'/sff') {
	$sff_dir = $self->directory.'/sff';
    }
    elsif (-d $self->directory.'/../sff') {
	$sff_dir = $self->directory.'/../sff';
    }
    else {
	$self->error_message("Failed to find sff_dir to create input fasta and qual files in either "
			     .$self->directory.'/sff'.' or '.$self->directory.'/../sff');
	return;
    }
    my @sff_files = glob("$sff_dir/*sff");

    unless (@sff_files) {
	$self->error_message("Failed to find any sff files in sff dir: $sff_dir");
	return;
    }
    
    #create fasta and qual from sff files using sffinfo
    foreach my $file (@sff_files) {
	my $file_name = File::Basename::basename ($file);
	my $fasta = $self->directory."/edit_dir/$file_name";
	my $qual = $self->directory."/edit_dir/$file_name";

	$fasta =~ s/\.sff$/\.fasta/;
	$qual =~ s/\.sff$/\.fasta\.qual/;

	#skip if fasta and qual already exists
	if (-s $fasta.'.gz' and -s $qual.'.gz') {
	    $self->status_message("Input files for $file already exists so skipping:\n\t$fasta".'.gz'."\n\t$qual".'.gz');
	    next;
	}

	my $cmd = "$sff_info -q $file > $qual";
	if (system($cmd)) {
	    $self->error_message("Failed to create qual file from sff file: $file\nCommand: $cmd");
	    return;
	}
	
	$cmd = "$sff_info -s $file > $fasta";
	if (system($cmd)) {
	    $self->error_message("Failed to create fasta file from sff file: $file\nCommand: $cmd");
	    return
	}

	#zip afsta and qual files
	if (system("gzip $fasta $qual")) {
	    $self->error_message("Failed to zip fasta and qual files");
	    return;
	}
    }

    return 1;
}

1;
