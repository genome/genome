package Genome::Model::Tools::Picard::ReplaceSequenceNames;

use strict;
use warnings;

use Genome;
use Bio::DB::Sam;
use File::Basename;

class Genome::Model::Tools::Picard::ReplaceSequenceNames {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_bam_file => {
            is  => 'String',
            doc => 'BAM file from which SAMRecords will be read.  The sequence names will be replaced with those from HEADER file',
        },
        output_bam_file => {
            is  => 'String',
            doc => 'SAMFileHeader from HEADER file will be written to this file, followed by SAMRecords from INPUT file.  File type is determined by suffix.',
        },
        header_file => {
            is => 'String',
            doc => 'SAM/BAM file from which SAMFileHeader will be read.  The sequence names from this file will be used to replace the sequence names in the input SAM/BAM file',
        },
    ],
};

sub help_brief {
    'Tool to replace the sequence names in a BAM file using BioSamtools and Picard';
}

sub help_detail {
    return <<EOS
    Tool to replace the sequence names in a BAM file using BioSamtools and Picard.  
EOS
}

sub execute {
    my $self = shift;

    my ($basename,$dirname,$suffix) = File::Basename::fileparse($self->header_file,qw/\.sam \.bam/);
    my $header;
    if ($suffix eq '.sam') {
        my $tam = Bio::DB::Tam->open($self->header_file);
        unless ($tam) {
            die('Failed to open sequence dictionary header SAM file: '. $self->header_file);
        }
        $header = $tam->header_read;
        unless ($header) {
            die('Failed to read TAM header!');
        }
    } elsif ($suffix eq '.bam') {
        my $bam = Bio::DB::Bam->open($self->header_file);
        unless ($bam) {
            die('Failed to open sequence dictionary header BAM file: '. $self->header_file);
        }
        $header = $bam->header();
        unless ($header) {
            die('Failed to read BAM header!');
        }
    } elsif ($suffix eq '.fai') {
        die('Please Fix.  Add support for FASTA index files.');
    }

    # Open the input BAM file and read the header
    my $input_bam = Bio::DB::Bam->open($self->input_bam_file);
    unless ($input_bam) {
        die('Failed to read input BAM file: '. $self->input_bam_file);
    }
    my $input_header = $input_bam->header();
    unless ($input_header) {
        die('Failed to read input BAM file header!');
    }
    
    # Check the total number of targets.  If there are more in the header than the input, this process will not work
    my $n_targets = $header->n_targets;
    my $input_n_targets = $input_header->n_targets;
    unless ($input_n_targets <= $n_targets) {
        die('There are more targets in the input BAM file than the sequence dictionary header file!');
    }

    my %seq_id_map;
    # Assume that if a chromosome is identical length, it's the same chromosome
    for (my $i = 0; $i < $input_n_targets; $i++) {
        my $input_seq_id = $input_header->target_name->[$i];
        my $input_length = $input_header->target_len->[$i];
        for ( my $j = 0; $j < $n_targets; $j++) {
            my $seq_id = $header->target_name->[$j];
            my $length = $header->target_len->[$j];
            if ($length == $input_length) {
                $seq_id_map{$input_seq_id} = $seq_id;
            }
        }
        unless (defined($seq_id_map{$input_seq_id})) {
            die('Failed to find seq id mapping for input chromosome '. $input_seq_id);
        }
    }
    
    # Parse the header text of the input BAM file and replace with the chromosome id from the header BAM 
    # Only look at BAM header lines that are part of the sequence dictionary
    my $header_text = $input_header->text;
    my @header_lines = split("\n",$header_text);
    for my $input_seq_id (keys %seq_id_map) {
        my $input_pattern = 'SN:'.$input_seq_id;
        my $output_pattern = 'SN:'.$seq_id_map{$input_seq_id};
        for my $header_line (@header_lines) {
            if ($header_line =~ /^\@SQ/ && $header_line =~ /$input_pattern/) {
                $header_line =~ s/$input_pattern/$output_pattern/;
            }
        }
    }

    #TODO: Add status message with old to new mapping
    
    # A new SAM header that can be used to replace the input BAM header
    my ($tmp_sam_fh,$tmp_sam_file) = Genome::Sys->create_temp_file();
    # Print out all the lines from the input BAM file now that the sequence names have been replaced above
    for my $header_line (@header_lines) {
        print $tmp_sam_fh $header_line ."\n";
    }
    $tmp_sam_fh->close;

    # TODO: Add status message with old to new seqeuence dictionary and header

    my $jar_path = $self->picard_path .'/ReplaceSamHeader.jar';
    unless (-e $jar_path) {
        die('Failed to find '. $jar_path .'!  This command may not be available in version '. $self->use_version);
    }
    my $input_file = $self->input_bam_file;
    my $header_file = $tmp_sam_file;
    my $output_file = $self->output_bam_file;
    my $sort_cmd = $jar_path .' net.sf.picard.sam.ReplaceSamHeader OUTPUT='. $output_file .' INPUT='. $input_file .' HEADER='. $header_file;
    $self->run_java_vm(
        cmd => $sort_cmd,
        input_files => [$input_file,$header_file],
        output_files => [$output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
