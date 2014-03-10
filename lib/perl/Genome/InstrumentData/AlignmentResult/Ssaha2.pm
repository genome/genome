package Genome::InstrumentData::AlignmentResult::Ssaha2;

use strict;
use warnings;
use File::Basename;

use Genome;

class Genome::InstrumentData::AlignmentResult::Ssaha2 {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => { value => 'ssaha2', is_param=>1 },
    ],
    has_transient => [
        static_params => { is => 'String', is_optional => 1 }
    ]
};

sub required_arch_os { 'x86_64' }

sub required_rusage { 
    "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>90000 && mem>24000] rusage[tmp=90000, mem=24000]' -M 24000000";
}

sub _run_aligner {
    my $self = shift;

    my $input_pathnames = join ' ', @_;
    my $aligner_params = $self->aligner_params;
    
    # collect filepaths
    my $ssaha_path = Genome::Model::Tools::Ssaha2->path_for_ssaha2_version($self->aligner_version);
    my $ref_index = $self->get_reference_sequence_index->full_consensus_path('fa.ssaha2.body');
    $ref_index =~ s/\.body$//;;  # all we want to pass is all_sequences.ssaha2
    
    unless (defined $ref_index && -s $ref_index) {
      $self->error_message("Ssaha2 index file either doesn't exist or is empty at $ref_index or " . $self->reference_build->data_directory);  
      die;
    };
    
    my $output_file = $self->temp_scratch_directory . "/all_sequences.sam";
    my $log_file = $self->temp_staging_directory . "/aligner.log";

    # construct the command (using hacky temp-file to append)
    $self->static_params('-best 1 -udiff 1 -align 0 -output sam_soft');
    $aligner_params = $aligner_params || '';
    if ( @_ > 1 && not $aligner_params =~ /-pair/ ){
        my ($lower,$upper) = $self->_derive_insert_size_bounds(600, 100);
        $self->static_params($self->static_params . " -pair $lower,$upper");
    }
    my $static_params = $self->static_params;

    my $cmd = "$ssaha_path $aligner_params $static_params -outfile $output_file.tmp -save $ref_index $input_pathnames >>$log_file && cat $output_file.tmp >>$output_file";

    Genome::Sys->shellcmd(
        cmd          => $cmd,
        input_files  => \@_,
        output_files => [ $output_file, $log_file ],
        skip_if_output_is_present => 0,
    );

    # 
    # Ssaha doesn't output unaligned reads by default. There may be a way to do so, or it might take comparing input / output.
    #
    unless (-s $output_file){
        $self->error_message('The sam output file is missing or empty.');
        return 0;
    }
    $self->debug_message('SSAHA2 alignment finished.');
    return 1;
}

sub aligner_params_for_sam_header {
    my $self = shift;
    return 'ssaha2 ' . $self->aligner_params . ' ' . $self->static_params;
}

sub fillmd_for_sam { return 1; } 

sub _check_read_count {
    my $self = shift;
    $self->warning_message("SSAHA2 does not support outputting unaligned reads.  Cannot verify read count against original fasta.");
    return 1;
}


sub prepare_reference_sequence_index {
    my $class = shift;
    my $refindex = shift;

    my $aligner_params = $refindex->aligner_params;
    my $reference_build_params = ""; 
    if ($aligner_params =~ m/rtype (.*?)/)  {
        $reference_build_params .= " -rtype $1";
    }
    if ($aligner_params =~ m/kmer\s*(.*?)(\s.*|)$/) {
        $reference_build_params .= " -kmer $1";
    }
    if ($aligner_params =~ m/skip\s*(.*?)(\s.*|)$/) {
        $reference_build_params .= " -skip $1";
    }

    my $staging_dir = $refindex->temp_staging_directory;
    my $staged_fasta_file = sprintf("%s/all_sequences.fa", $staging_dir);
    my $target_fasta_file = readlink($staged_fasta_file);
    unless ($target_fasta_file) {
        $class->error_message("Cant resolve the target of the fasta symlink dropped in the staging directory.");
        return;
    }
    unless (symlink($target_fasta_file, $staged_fasta_file.'.ssaha2')) {
        $class->error_message("Couldn't make an all_sequences.fa.ssaha2 symlink to the raw fasta file");
        return;
    }
        
    

    my $ss_path = Genome::Model::Tools::Ssaha2->path_for_ssaha2_version($refindex->aligner_version);
    my $ss_cmd = sprintf('%sBuild %s -save %s.ssaha2 %s', $ss_path, $reference_build_params, $staged_fasta_file, $staged_fasta_file);
    my $rv = Genome::Sys->shellcmd(
        cmd => $ss_cmd,
    );

    unless($rv) {
        $class->error_message('ssaha2 indexing failed');
        return;
    }

    return 1;
}

sub aligner_params_required_for_index {
    return 1;
}


