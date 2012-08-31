package Genome::Model::Tools::Fastq::Dust;

use strict;
use warnings;
use Genome;
use File::Basename;

class Genome::Model::Tools::Fastq::Dust{
    is => 'Command',
    has =>[
        fastq_file => {
            is => 'Path',
            doc => 'Fastq to dust',
        },
        output_file => {
            is => 'Text',
            doc => 'Path to write output fastq file',
        },
    ],
};

sub execute{
    my $self = shift;

    my $fastq_file = $self->fastq_file;
    unless (-e $fastq_file){
        die $self->error_message("Input fastq file($fastq_file) to dust doesn't exist!");
    }
    unless (-s $fastq_file){
        die $self->error_message("Input fastq file($fastq_file) to dust has zero size!");
    }

    my $fastq_basename = basename $fastq_file;
    my $temp_dir = Genome::Sys->create_temp_directory;
    my $temp_base = "$temp_dir/$fastq_basename";

    my $fasta_file = $temp_base.".FASTA";
    my $qual_file = $temp_base.".QUAL";

    my $dusted_file = $fasta_file.".DUSTED";
    

    $self->status_message("Running dust on $fastq_file");

    my $fastq_fh  = Genome::Sys->open_file_for_reading($fastq_file);
    unless ($fastq_fh) {
        $self->error_message('Failed to open fastq file ' . $fastq_file . ": $!");
        return;
    }
    binmode $fastq_fh, ":utf8";

    my $fasta_output_fh = Genome::Sys->open_file_for_writing($fasta_file);
    unless ($fasta_output_fh) {
        $self->error_message('Failed to open output file ' . $fasta_file . ": $!");
        return;
    }
    binmode $fasta_output_fh, ":utf8";

    my $qual_output_fh = Genome::Sys->open_file_for_writing($qual_file);
    unless ($qual_output_fh) {
        $self->error_message('Failed to open output file ' . $qual_file . ": $!");
        return;
    }
    binmode $qual_output_fh, ":utf8";

    while (my $header = $fastq_fh->getline) {
        my $seq  = $fastq_fh->getline;
        my $sep  = $fastq_fh->getline;
        my $qual = $fastq_fh->getline;

        unless (substr($header,0,1) eq '@') {
            die "Unexpected header in fastq! $header";
        }
        substr($header,0,1) = '>';

        $fasta_output_fh->print($header, $seq);
        $qual_output_fh->print($sep, $qual);
    }
    $fasta_output_fh->close();
    $fastq_fh->close;
    $qual_output_fh->close();

    unless (-s $fasta_file){
        die $self->error_message("fasta file($fasta_file) separated for dusting has no size!");
    }
    unless (-s $qual_file){
        die $self->error_message("qual file($qual_file) separated from fastq has no size!");
    }

    #2. run dust command
    my $cmd = "dust $fasta_file > $dusted_file";
    my $rv = Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$fasta_file],
        output_files => [$dusted_file],
    );

    unless ($rv){
        die $self->error_message("failed to execute cmd: $cmd : $!");
    }

    #3. re-produce fastq 

    my $dusted_fh  = Genome::Sys->open_file_for_reading($dusted_file);
    unless ($dusted_fh) {
        $self->error_message('Failed to open fastq file ' . $dusted_file . ": $!");
        return;
    }
    binmode $dusted_fh, ":utf8";

    my $qual_input_fh = Genome::Sys->open_file_for_reading($qual_file);
    unless ($qual_input_fh) {
        $self->error_message('Failed to open input file ' . $qual_file . ": $!");
        return;
    }
    binmode $qual_input_fh, ":utf8";

    my $out_fh = Genome::Sys->open_file_for_writing($self->output_file);
    unless ($out_fh) {
        $self->error_message('Failed to open output file ' . $self->output_file . ": $!");
        return;
    }
    binmode $out_fh, ":utf8";

    # since dusting wraps sequences, may have to read multiple lines to reconstruct sequence
    # pull header then concat lines until next header encountered
    my ($header, $seq, $sep, $qual);
    while (my $line = $dusted_fh->getline) {
        if ($line=~/^>.*/) { #found a header 
            # this only grabs the header on the first sequence
            # other sequences in the file will have their header pre-caught below
            # confusing :(
            $header = $line;
        }
        else {
            chomp($seq .= $line);
            #$seq .= $line;
        }

        while ($line = $dusted_fh->getline) { #accumulate lines for read, until next header encountered 
            if ($line=~/^>.*/) { #found a new header - read has been accumulated 
                last;
            }
            else {
                chomp($seq .= $line);
                #$seq .= $line;
            }
        }

        $sep = $qual_input_fh->getline;
        $qual = $qual_input_fh->getline;

        unless (substr($header,0,1) eq '>') {
            die "Unexpected fasta header: $header";
        }
        substr($header,0,1) = '@';
        $out_fh->print("$header$seq\n$sep$qual");

        #reset
        $seq = '';
        $header = $line;
    }
    return 1;
}

1;
