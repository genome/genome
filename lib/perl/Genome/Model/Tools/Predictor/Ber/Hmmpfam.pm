package Genome::Model::Tools::Predictor::Ber::Hmmpfam;

use strict;
use warnings;
use Carp;
use Genome;

class Genome::Model::Tools::Predictor::Ber::Hmmpfam {
    is  => ['Command::V2'],
    has => [
        output_directory => {
            is => 'DirectoryPath',
            is_input => 1,
            doc => 'Directory in which raw and parsed output from this predictor should go',
        },
        input_fasta_file => {
            is => 'FilePath',
            is_input => 1,
            doc => 'File containing assembly sequence (typically fasta) to be used as input to predictor',
        },
        ber_source_path => {
            is => 'DirectoryPath',
            is_input => 1,
            doc => 'Directory in which ber source app lives',
        },
        htab_file => { 
            is => 'FilePath',  
            is_optional => 1,
            is_output => 1,
            doc => 'File path of the output htab file',
        },
        lsf_queue => {
            is_param => 1,
            is_input => 1,
            default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        },
        lsf_resource => {
            is_param => 1,
            is_input => 1,
            default => "-M 4000000 -R 'select[type==LINUX64 && mem4000] rusage[mem=4000]'",
        },
    ],
};

sub execute { 
    my $self = shift;

    my $config_dir = $self->output_directory.'/';
    my $seq_in = Bio::SeqIO->new(-file => $self->input_fasta_file, -format => 'Fasta')
        or croak "failed to open: ".$self->input_fasta_file;

    my $seq  = $seq_in->next_seq();
    my $gene = $seq->primary_id;
    $seq_in->close;


    my $hmmpfam_db = $self->ber_source_path.'/data/ALL_LIB.HMM';
    my $hmmpfam_source = '/gsc/bin/hmmpfam';

    my $hmmpfam_output = $self->hmmpfam_file_name($config_dir, $gene);
    my $htab_output = $self->htab_file_name($config_dir,$gene);

    my $hmmpfam_cmd = join(' ', $hmmpfam_source,$hmmpfam_db, $self->input_fasta_file,
                           '>', $hmmpfam_output
                       );
    #execute hmmpfam_cmd job
    my $rv = Genome::Sys->shellcmd(
        cmd => $hmmpfam_cmd,
    );
    unless ($rv) {
        die "Failed to execute hmmpfam command: $hmmpfam_cmd!";
    }
    return 0 unless(-e $hmmpfam_output and -s $hmmpfam_output);

    my $hmm2htab_source = $self->ber_source_path.'/src/hmmToHtab.pl';
    my $hmm2htab_source_info = $self->ber_source_path.'/src/hmm_info.txt';
    my $htab_cmd = 'genome-perl '.$hmm2htab_source. ' -H '.$hmm2htab_source_info.' < '.$hmmpfam_output.
        ' > '.$htab_output; 

    #execute hmm2htab command
    $rv = Genome::Sys->shellcmd(
        cmd => $htab_cmd,
    );
    unless ($rv) {
        die "Failed to execute hmmToHtab command: $htab_cmd!";
    }

    return $self->validate_output($hmmpfam_output, $htab_output) ? $self->htab_file($htab_output) : 0;
}

sub hmmpfam_file_name {
    my ($self, $dir, $gene) = @_;
    return $dir.'hmm/'.$gene.'.hmmpfam';
}

sub htab_file_name {
    my ($self, $dir, $gene) = @_;
    return $dir.'hmm/'.$gene.'.hmmpfam.htab';
}


sub no_domain_hits {
    my $self = shift;
    my $file = shift;
    
    #if there are no hits found in a zero htab file, there is a problem, otherwise it is OK
    my @no_hits = grep {$_ =~ /no hits above thresholds/} Genome::Sys->read_file($file);
    return (@no_hits==3) ? 1 : 0;
}

sub validate_output {
    my ($self, $hmmpfam_output, $htab_output) = @_;

    return 1 if(-e $htab_output and -s $htab_output);
    if (-e $htab_output and -z $htab_output) {
        return $self->no_domain_hits($hmmpfam_output) ? 1 : 0;
    }
}

1;
