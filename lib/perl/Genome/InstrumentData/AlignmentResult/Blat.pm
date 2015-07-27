package Genome::InstrumentData::AlignmentResult::Blat;

use strict;
use warnings;
use File::Basename;

use Genome;

class Genome::InstrumentData::AlignmentResult::Blat {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => { value => 'Blat', is_param=>1 },
    ],
    has_transient_optional => [
         _Blat_sam_cmd => { is=>'Text' }
    ]
};

sub required_arch_os { 'x86_64' }

sub required_rusage { 
    "-R 'select[tmp>90000 && mem>12000] span[hosts=1] rusage[tmp=90000, mem=12000]' -M 12000000 -n 6";
}

sub _run_aligner {
    my $self = shift;
    my @input_pathnames = @_;
    
    my $tmp_dir = $self->temp_scratch_directory;

    my $tmp_blastn_file = $tmp_dir . "/" . "aligned.blastn"; 
    my $tmp_sam_file = $tmp_dir . "/" . "aligned.sam";

    my $staging_sam_file = $tmp_dir . "/" . "all_sequences.sam";
    
    my $reference_build = $self->reference_build;
    my $ref_path = $reference_build->full_consensus_path('fa');

    my $blat_path = Genome::Model::Tools::Blat->path_for_blat_version($self->aligner_version);

    my $aligner_params = $self->decomposed_aligner_params;

    foreach my $cur_input (@input_pathnames) {

        my $cur_input_fa = $cur_input.".fa";
        my $fastq_fh = IO::File->new($cur_input);
        my $fasta_fh = IO::File->new(">$cur_input_fa");
        my $line_type;
        while (<$fastq_fh>){
            $line_type = $. % 4;
            if ($line_type == 1) {
                $fasta_fh->print(">",substr($_,1));
            } elsif ($line_type == 2) {
                $fasta_fh->print($_);
            }
        }
        $fasta_fh->close();
        $fastq_fh->close();

        my $align_cmd = $blat_path . sprintf(' %s %s %s %s',
            $ref_path, $cur_input_fa, $aligner_params, $tmp_blastn_file);
        Genome::Sys->shellcmd(
            cmd             => $align_cmd,
            input_files     => [ $ref_path, $cur_input_fa ],
            output_files    => [ $tmp_blastn_file ],
            skip_if_output_is_present => 0,
        );

        my $convert_cmd = "/usr/bin/blast2sam.pl" . sprintf(' %s > %s',
            $tmp_blastn_file, $tmp_sam_file);

        Genome::Sys->shellcmd(
            cmd             => $convert_cmd,
            input_files     => [ $tmp_blastn_file ],
            output_files    => [ $tmp_sam_file ],
            skip_if_output_is_present => 0,
        );

        die "Failed to process sam command line, error_message is ".$self->error_message unless $self->_filter_sam_output($tmp_sam_file, $staging_sam_file);
    }

    return 1;
}

sub _filter_sam_output {
    my ($self, $cur_sam_file, $all_sequences_sam_file) = @_;

    my $cur_sam_fh = IO::File->new( $cur_sam_file );
    if ( !$cur_sam_fh ) {
        $self->error_message("Error opening current sam file for reading $!");
        return;
    }
    $self->debug_message("Opened $cur_sam_file");

    my $all_seq_fh = IO::File->new(">>$all_sequences_sam_file");
    if ( !$all_seq_fh ) {
        $self->error_message("Error opening all seq sam file for writing $!");
        return;
    }
    $self->debug_message("Opened $all_sequences_sam_file");
    
    while (<$cur_sam_fh>) {
        #write out the aligned map, excluding the default header- all lines starting with @.
        $all_seq_fh->print($_) unless $_ =~ /^@/;
    }

    $cur_sam_fh->close;
    $all_seq_fh->close;
    return 1;
}

sub decomposed_aligner_params {
    my $self = shift;
    my $params = $self->aligner_params || "-out=blast";
    
    return ($params);
}

sub aligner_params_for_sam_header {
    my $self = shift;
    
    my $params = $self->decomposed_aligner_params;
    
    return "blat $params";
}

1;

