package Genome::Model::MetagenomicShotgun::Command::MapxAlignmentToAlignedFastq;

use strict;
use warnings;
use Genome;

class Genome::Model::MetagenomicShotgun::Command::MapxAlignmentToAlignedFastq{
    is => 'Command::V2',
    has =>[
        mapx_alignment => {
            is => 'File',
            doc => 'mapx alignment.txt file containing aligned reads'
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'instrument data that produced the alignment, used to extract the original fastq lines',
        },
        output_directory => {
            is => 'Path',
            doc => 'directory to put the output in, output is in format <output_directory>/<instrument_data_id>/s_<lane>[_(1|2)]_sequence.txt',
        }
    ],
};

sub execute{
    my $self = shift;

    my $instrument_data = $self->instrument_data;

    unless (-e $self->mapx_alignment){
        die $self->error_message("Mapx alignment file ".$self->mapx_alignment."doesn't exist!");
    }
    
    my $lane = $instrument_data->lane;
    unless ($lane){
        die $self->error_message("Couldn't get lane info from instrument data ".$instrument_data->id);
    }
    
    my $output_dir = $self->output_directory.'/'.$instrument_data->id;
    if (glob($output_dir.'/*')){
        die $self->error_message("Files already exist in specified output directory, aborting.");
    }

    unless ( -d $output_dir){
        unless (File::Path::make_path($output_dir)){
            die $self->error_message("Failed to create output dir $output_dir")
        }
    }

    $self->debug_message("Collecting read names for extracting from instrument data");
    
    my $sorted_reads = IO::File->new("awk ' \$0 !~ /^#/ {print \$3}' ".$self->mapx_alignment." | sort -u |");
    unless ($sorted_reads){
        die $self->error_message("can't get file handle for sorted reads");
    }

    # dump fastqs
    $self->debug_message("Getting fastq files from source instrument data");
    my @fastq_files = $instrument_data->dump_sanger_fastq_files;
    unless (@fastq_files){
        die $self->error_message("no fastq files for instrument data");
    }
    my $reader = Genome::Model::Tools::Sx::Reader->create(config => [ map { $_.':type=sanger' } @fastq_files ]);
    if ( not $reader ) {
        die $self->error_message("Failed to create SX reader for fastqs: @fastq_files");
    }

    # open outputs
    my @config;
    if ( $instrument_data->is_paired_end ) {
        @config = ( 
            $self->output_directory.'/'.$instrument_data->id.'/s_'.$lane.'_1_sequence.txt:type=sanger:name=fwd',
            $self->output_directory.'/'.$instrument_data->id.'/s_'.$lane.'_2_sequence.txt:type=sanger:name=rev',
            $self->output_directory.'/'.$instrument_data->id.'/s_'.$lane.'_sequence.txt:type=sanger:name=sing',
        );
    }
    else {
        @config = ( $self->output_directory.'/'.$instrument_data->id.'/s_'.$lane.'_sequence.txt:type=sanger' );
    }
    my $writer = Genome::Model::Tools::Sx::Writer->create(config => \@config); 
    if ( not $writer ) {
        die $self->error_message('Failed to create SX writer for files in output directory: '.$self->output_directory);
    }

    $self->debug_message("Extracting aligned reads from source instrument data");

    my $cached_read = $sorted_reads->getline;
    while ($cached_read){
        chomp $cached_read;
        my $read = $sorted_reads->getline;
        my $seqs;
        do { 
            $seqs = $reader->read;
        } until not $seqs or grep { $_->{id} eq $cached_read } @$seqs;
        if ( not $seqs ) {
            die $self->error_message("Aligned read ($cached_read) not present in fastq files!");
        }
        my $cached_read_base = $cached_read;
        $cached_read_base =~ s/\/[12]$//;
        if ($read and $read =~/$cached_read_base/){
            $writer->write($seqs);
            $cached_read = $sorted_reads->getline;
        }else{
            $writer->write([ grep {$_->{id} eq $cached_read} @$seqs ]);
            $cached_read = $read;
        }
    }
    $self->debug_message("Finished extraction");
    return 1;
}

1;

