package Genome::Model::Build::DeNovoAssembly::Newbler;

use strict;
use warnings;

use Genome;

require Carp;
use Data::Dumper 'Dumper';
use Regexp::Common;

class Genome::Model::Build::DeNovoAssembly::Newbler {
    is => 'Genome::Model::Build::DeNovoAssembly',
};

#< Methods for new fastq processing pipeline >#
sub read_processor_output_files_for_instrument_data {
    my ( $self, $inst_data ) = @_;
    my $file_name = $inst_data->id.'-input.fastq';
    return $self->data_directory."/$file_name";
}

sub fastq_input_files {
    my $self = shift;
    my @fastq_files = glob( $self->data_directory."/*input.fastq" );
    unless ( @fastq_files ) {
        $self->debug_message( "Did not find any input fastq files in build data directory" );
        return;
    }
    #return join (',', map{$_} @fastq_files );
    return @fastq_files;
}

#< Files / Dirs >#
sub assembly_directory {
    return $_[0]->data_directory.'/assembly';
}

sub sff_directory {
    return $_[0]->data_directory.'/sff';
}

sub existing_assembler_input_files {
    my $self = shift;
    return grep { -s $_ } $self->fastq_input_files;
}

sub stats_file {
    return $_[0]->data_directory.'/consed/edit_dir/stats.txt';
}

sub all_contigs_fasta_file {
    return $_[0]->data_directory.'/454AllContigs.fna';
}

sub all_contigs_qual_file {
    return $_[0]->data_directory.'/454AllContigs.qual';
}

sub all_contigs_ace_file {
    return $_[0]->data_directory.'/consed/edit_dir/454Contigs.ace.1';
}

sub scaffolds_fasta_file {
    return $_[0]->data_directory.'/454Scaffolds.fna';
}

sub scaffolds_qual_file {
    return $_[0]->data_directory.'/454Scaffolds.qual';
}

sub scaffolds_agp_file {
    return $_[0]->data_directory.'/454Scaffolds.txt';
}

sub newbler_metrics_file {
    return $_[0]->data_directory.'/454NewblerMetrics.txt';
}


sub input_fastas {
    my $self = shift;
    my @files;
    foreach (glob($self->data_directory."/*fasta.gz")) {
	#make sure qual file exists for the fasta
	my ($qual_file) = $_ =~ s/\.gz$/\.qual\.gz/;
	next unless -s $qual_file;
	push @files, $_;
    }
    
    return @files;
}

sub fasta_file {
    my $self = shift;
    # FIXME
    my @instrument_data = $self->instrument_data;
    #SINGULAR FOR NOW .. NEED TO GET IT TO WORK FOR MULTIPLE INPUTS
    my $fasta = $instrument_data[0]->fasta_file;
    unless ($fasta) {
	$self->error_message("Instrument data does not have a fasta file");
	return;
    }
    #COPY THIS FASTA TO BUILD INPUT_DATA_DIRECTORY
    File::Copy::copy ($fasta, $self->input_data_directory);
    #RENAME THIS FASTA TO SFF_NAME.FA ??
    #RETURN TO PREPARE-INSTRUMENT DATA FOR CLEANING
    return 1;
}

sub input_data_directory {
    my $self = shift;
    mkdir $self->data_directory.'/input_data' unless
	-d $self->data_directory.'/input_data';
    return $self->data_directory.'/input_data';
}
#<>#

#<ASSEMBLE>#
sub assembler_params {
    my $self = shift;

    my %params = $self->processing_profile->assembler_params_as_hash;
    $params{version} = $self->processing_profile->assembler_version;
    $params{output_directory} = $self->data_directory;
    $params{input_files} = [ $self->fastq_input_files ] if $self->fastq_input_files;
    $params{import_location} = $self->model->import_location if $self->model->import_location;

    return %params;
}

#< Metrics >#
#run base methods
#<>#

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/Genome/Model/Build/DeNovoAssembly/Newbler.pm $
#$Id: Newbler.pm 59774 2010-06-09 23:46:44Z ebelter $
