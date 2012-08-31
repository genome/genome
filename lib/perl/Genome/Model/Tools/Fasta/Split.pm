package Genome::Model::Tools::Fasta::Split;

use strict;
use warnings;

use Genome;

use FASTAParse;

class Genome::Model::Tools::Fasta::Split {
    is => 'Command',
    has => [
        fasta_file => { is => 'Text', },
    ],
    has_optional => {
        _split_fasta_files => { },
        min_sequence => {
            is => 'Number',
            doc => 'The minimum number of base pair sequence to include in each file.',
        },
        number_of_files => {
            is => 'Number',
            doc => 'The number of files to return.',
        },
        output_directory => {
            is => 'Text',
            doc => 'The directory to output the split fasta files',
        }
    }
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self->min_sequence || $self->number_of_files) {
        die('Must provide min_sequence or number_of_files');
    }
    return $self;
}

sub execute {
    my $self = shift;

    my @fasta_files;

    my $file_counter = 1;
    my ($fasta_basename,$fasta_dirname) = File::Basename::fileparse($self->fasta_file);
    unless (Genome::Sys->validate_directory_for_read_write_access($fasta_dirname)) {
        $self->error_message('Failed to validate directory '. $fasta_dirname ." for read/write access:  $!");
        die($self->error_message);
    }
    unless ($self->output_directory) {
        $self->output_directory($fasta_dirname);
    }
    my $output_file = $self->output_directory .'/'. $fasta_basename .'_'. $file_counter;

    my $output_fh = Genome::Sys->open_file_for_writing($output_file);
    unless ($output_fh) {
        $self->error_message('Failed to open output file '. $output_file);
        die($self->error_message);
    }
    push @fasta_files, $output_file;
    my $fasta_reader = Genome::Sys->open_file_for_reading($self->fasta_file);
    unless ($fasta_reader) {
        $self->error_message('Failed to open fasta file '. $self->fasta_file);
        die($self->error_message);
    }
    local $/ = "\n>";
    if ($self->number_of_files) {
        my $seq_len;
        while (<$fasta_reader>) {
            if ($_) {
                chomp;
                if ($_ =~ /^>/) { $_ =~ s/\>//g }
                my $myFASTA = FASTAParse->new();
                $myFASTA->load_FASTA( fasta => '>' . $_ );
                $seq_len += length( $myFASTA->sequence() );
            }
        }
        $fasta_reader->seek(0,0);
        $self->min_sequence( int( ( $seq_len / $self->number_of_files ) ) );
    }
    my $total_seq = 0;
    while (<$fasta_reader>) {
        if ($_) {
            chomp;
            if ($_ =~ /^>/) { $_ =~ s/\>//g }
            my $myFASTA = FASTAParse->new();
            $myFASTA->load_FASTA( fasta => '>' . $_ );
            my $seqlen = length( $myFASTA->sequence() );
            $total_seq += $seqlen;
            if ($total_seq > ($self->min_sequence)) {
                $file_counter++;
                $output_fh->close;
                $output_file = $self->output_directory .'/'. $fasta_basename .'_'. $file_counter;
                $output_fh = Genome::Sys->open_file_for_writing($output_file);
                unless ($output_fh) {
                    $self->error_message('Failed to open fasta output file '. $output_file);
                    return;
                }
                push @fasta_files, $output_file;
                $total_seq = 0;
            }
            print $output_fh '>'. $myFASTA->id() ."\n". $myFASTA->sequence() . "\n";
        }
    }
    $fasta_reader->close;
    $output_fh->close;

    if ($self->number_of_files && $self->number_of_files != scalar(@fasta_files)) {
        $self->warning_message('Expected to return '. $self->number_of_files .' but returning '. scalar(@fasta_files) .' fasta files.');
    }
    $self->_split_fasta_files(\@fasta_files);
    return 1;
}

1;
