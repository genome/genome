package Genome::InstrumentData::AlignmentResult::Smalt;

use strict;
use warnings;
use File::Basename;

use Genome;

class Genome::InstrumentData::AlignmentResult::Smalt {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => { value => 'smalt', is_param=>1 },
    ],
    has_transient => [
        static_params => { is => 'String', is_optional => 1 }
    ]
};

sub required_arch_os { 'x86_64' }

sub required_rusage { 
    "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>90000 && mem>12000] span[hosts=1] rusage[tmp=90000, mem=12000]' -M 12000000";
}

sub _run_aligner {
    my $self = shift;

    my $input_pathnames = join ' ', @_;
    my $aligner_params = $self->aligner_params || '';
    
    # collect filepaths
    my $smalt_path = Genome::Model::Tools::Smalt->path_for_smalt_version($self->aligner_version);
    my $ref_index_base = substr($self->get_reference_sequence_index->full_consensus_path('fa.smi'),0,-4);
    
    unless (defined $ref_index_base && -s "$ref_index_base.smi" && -s "$ref_index_base.sma") {
      $self->error_message("Smalt index files either don't exist or are empty at $ref_index_base or " . $self->reference_build->data_directory);  
      return 0;
    };
    
    my $output_file = $self->temp_scratch_directory . "/all_sequences.sam";
    my $log_file = $self->temp_staging_directory . "/aligner.log";

    # construct the command
    $self->static_params('-f samsoft');
    my $static_params = $self->static_params;

    my $cmd = "$smalt_path map $aligner_params $static_params -o $output_file.tmp $ref_index_base $input_pathnames >>$log_file && cat $output_file.tmp >>$output_file";

    Genome::Sys->shellcmd(
        cmd          => $cmd,
        input_files  => \@_,
        output_files => [ $output_file, $log_file ],
        skip_if_output_is_present => 0,
    );

    unless (-s $output_file){
        $self->error_message('The sam output file is missing or empty.');
        return 0;
    }
    $self->debug_message('Smalt alignment finished.');
    return 1;
}

sub aligner_params_for_sam_header {
    my $self = shift;
    return 'smalt map' . $self->aligner_params . ' ' . $self->static_params;
}

sub fillmd_for_sam { return 1; } 

sub _check_read_count {return 1;}

sub prepare_reference_sequence_index {
    my $class = shift;
    my $refindex = shift;

    my $staging_dir = $refindex->temp_staging_directory;
    my $staged_fasta_file = sprintf("%s/all_sequences.fa", $staging_dir);

    my $actual_fasta_file = $staged_fasta_file;

    if (-l $staged_fasta_file) {
        $class->debug_message(sprintf("Following symlink for fasta file %s", $staged_fasta_file));
        $actual_fasta_file = readlink($staged_fasta_file);
        unless($actual_fasta_file) {
            $class->error_message("Can't read target of symlink $staged_fasta_file");
            return;
        } 
    }
    my $smalt_path = Genome::Model::Tools::Smalt->path_for_smalt_version($refindex->aligner_version);

    my $index_cmd = sprintf("%s index -k 13 -s 6 %s/all_sequences.fa %s", $smalt_path, $staging_dir, $staged_fasta_file);

    Genome::Sys->shellcmd(cmd=>$index_cmd, output_files => [$staging_dir.'/all_sequences.fa.sma', $staging_dir.'/all_sequences.fa.smi']);

    # if it got to here then we have valid files 

    $class->debug_message("Successfully ran smalt index");

    return 1;

}
