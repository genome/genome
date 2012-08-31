package Genome::Model::Tools::Kmer::DuplicationRate;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Kmer::DuplicationRate {
    is => ['Genome::Model::Tools::Kmer'],
    has => [
        fastq_files => {
            is => 'Text',
            doc => 'All the fastq files to plot complexity of seperated by commas',
        },
        output_file => {
            is => 'Text',
            doc => 'The path to write occratio output to',
        },
    ],
};

sub execute {
    my $self = shift;

    my @fastq_files = split(/,/,$self->fastq_files);
    my @fasta_files;
    my @fastq_basenames;
    my $max_read_length;
    my $tmp_dir = Genome::Sys->create_temp_directory();
    for my $fastq_file (@fastq_files) {
        my ($fastq_basename,$fastq_dirname,$fastq_suffix) = File::Basename::fileparse($fastq_file,qw/\.fastq \.fq \.txt/);
        unless ($fastq_basename) {
            die('Failed to parse fastq file path '. $fastq_file);
        }
        my $fasta_file = $tmp_dir .'/'. $fastq_basename .'.fa';
        push @fastq_basenames, $fastq_basename;
        unless (-e $fasta_file) {
            unless (Genome::Model::Tools::Fastq::ToFasta->execute(
                fastq_file => $fastq_file,
                fasta_file => $fasta_file,
            )) {
                die('Failed to convert fastq_file '. $fastq_file .' to fasta file '. $fasta_file);
            }
        }
        my $fasta_fh = Genome::Sys->open_file_for_reading($fasta_file);
        while (my $line = $fasta_fh->getline){
            chomp($line);
            if ($line =~ /^>/) { next; }
            my $read_length = length($line);
            if (!defined($max_read_length) || ($read_length > $max_read_length)) {
                $max_read_length = $read_length;
            }
        }
        push @fasta_files, $fasta_file;
    }
    my $index_name = $tmp_dir .'/suffixerator';
    my $log_file = $index_name .'.log';
    unless (Genome::Model::Tools::Kmer::Suffixerator->execute(
        fasta_files => \@fasta_files,
        index_name => $index_name,
        log_file => $log_file,
    )) {
        die("Failed to run suffixerator on fasta files:\n". join("\n", @fasta_files) );
    }
    for my $fasta_file (@fasta_files) {
        unlink $fasta_file || die('Failed to remove fasta file '. $fasta_file);
    }
    #TODO: Make scan a variable that will do the right thing for really large index files
    unless (Genome::Model::Tools::Kmer::OccurrenceRatio->execute(
        index_name => $index_name,
        minimum_mer_size => 1,
        maximum_mer_size => $max_read_length,
        scan => 1,
        output_file => $self->output_file,
        output_type => 'nonuniquemulti relative',
    )) {
        die('Failed to generate occurence ratio file '. $self->output_file .' from index '. $index_name);
    }
    return 1;
}
