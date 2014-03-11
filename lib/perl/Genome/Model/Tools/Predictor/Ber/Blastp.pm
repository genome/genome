package Genome::Model::Tools::Predictor::Ber::Blastp;

use strict;
use warnings;
use Carp;
use Genome;

class Genome::Model::Tools::Predictor::Ber::Blastp {
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
        btab_file => { 
            is => 'FilePath',  
            is_optional => 1,
            is_output => 1,
            doc => 'File path of the output btab file',
        },
        lsf_queue => {
            is_param => 1,
            is_input => 1,
            default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        },
        lsf_resource => {
            is_param => 1,
            is_input => 1,
            default => "-M 4000000 -R 'select[type==LINUX64 && mem>4000] rusage[mem=4000]'",
        },
      ],
};

sub execute { 

    my $self = shift;
    my $config_dir = $self->output_directory.'/';
    my $blastp_db = $self->ber_source_path.'/data/panda/AllGroup/AllGroup.niaa';
    my $blastp_source = '/gsc/scripts/bin/blastp';

    my $seq_in = Bio::SeqIO->new(-file => $self->input_fasta_file, -format => 'Fasta')
        or croak "failed to open: ".$self->input_fasta_file;

    my $seq  = $seq_in->next_seq();
    my $gene = $seq->primary_id;
    $seq_in->close;

    my $blastp_output = $self->blastp_file_name($config_dir, $gene);
    my $btab_output = $self->btab_file_name($config_dir, $gene);
    my $blastp_cmd = join( ' ', $blastp_source,$blastp_db,$self->input_fasta_file,'>',
                           $blastp_output
                       );

    my $rv = Genome::Sys->shellcmd(
        cmd => $blastp_cmd,
    );
    unless ($rv) {
        die "Failed to execute blastp command: $blastp_cmd!";
    }

    return 0 unless(-e $blastp_output and -s $blastp_output);
    my $blast2btab_source = $self->ber_source_path.'/src/wu-blast2btab.pl';

    my $btab_cmd = 'genome-perl '.$blast2btab_source. ' --input '.$blastp_output.' --output '
        .$btab_output; 
    $rv = Genome::Sys->shellcmd(
        cmd => $btab_cmd,
    );
    unless ($rv) {
        die "Failed to execute blast2btab command: $btab_cmd!";
    }
    
    return $self->validate_output($blastp_output, $btab_output) ? $self->btab_file($btab_output) : 0;
}

sub blastp_file_name {
    my ($self, $dir, $gene) = @_;
    return $dir.'ber/'.$gene.'.nr';
}

sub btab_file_name {
    my ($self, $dir, $gene) = @_;
    return $dir.'ber/'.$gene.'.nr.btab';
}

sub no_blast_hits {
    my $self = shift;
    my $file = shift;

    #if there are no hits in a zero btab file, there is a problem, otherwise it is OK
    my @no_hits = grep {$_ =~ /\*\*\* NONE \*\*\*/} Genome::Sys->read_file($file);
    return (@no_hits==1) ? 1 : 0;
}

sub validate_output {
    my ($self, $blastp_output, $btab_output) = @_;
    
    return 1 if(-e $btab_output and -s $btab_output);
    if (-e $btab_output and -z $btab_output) {
        return $self->no_blast_hits($blastp_output) ? 1 : 0;
    }
}

1;
