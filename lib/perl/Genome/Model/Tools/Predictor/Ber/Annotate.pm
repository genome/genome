package Genome::Model::Tools::Predictor::Ber::Annotate;

use strict;
use warnings;
use Carp;
use DateTime;

use Genome;

class Genome::Model::Tools::Predictor::Ber::Annotate {
    is  => ['Command::V2'],
    has => [
        output_directory => {
            is => 'DirectoryPath',
            is_input => 1,
            doc => 'Directory in which raw and parsed output from this predictor should go',
        },
        ber_source_path => {
            is => 'DirectoryPath',
            is_input => 1,
            doc => 'Directory in which ber source app lives',
        },
         locus_id => {
            is => 'Text',
            is_input => 1,
            doc => 'locus id of annotation',
        },
        converge_result => {
            is => 'Text',
            is_optional => 1,
            is_input => 1,
            doc => '',
        },
        gram_stain => {
            is => 'Text',
            is_input => 1,
            default_value => 'negative',
            valid_values => ['positive', 'negative'],
            doc => 'Gram stain of species on which prediction is to be run, if relevant',
        },
        output_file => { 
            is => 'FilePath',  
            is_optional => 1,
            is_output => 1,
            doc => 'File path of the output dat file from BER',
        },
        lsf_queue => {
            is_param => 1,
            default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        },
    ],
};

sub execute { 

    my $self = shift;
    
    # original command /gsc/bin/perl autoAnnotate.dbi -blast `which blastp` -D $1 -dbd $dbd -F $gstring -nosignalp -now -o $cwd/../out/$label -nointerpolate -trace 5 -verbose -nodebug -nothreads > $outfile 2>&1
    my $blastp_source = '/gsc/scripts/bin/blastp';
    my $gscbin_perl = '/gsc/bin/perl';
    my $ber_source  = $self->ber_source_path.'/src/autoAnnotate.dbi';
    my $today = DateTime->now(time_zone => "America/Chicago");
    my $gram_stain = $self->gram_stain eq 'negative' ? 0 : 1;
    my $ber_output_file_key  = $self->output_directory.'/'.$self->locus_id.'-'.$today;
    my $ber_output_log       = $ber_output_file_key.'.output';
    my $options = '-blast '.$blastp_source.' -D '.$self->locus_id.
        ' -dbd SQLite -F -gram '.$gram_stain.' -nosignalp -now -o '.
            $ber_output_file_key.' -nointerpolate -trace 5 -verbose -nodebug -nothreads';
    my $cmd = join(' ', $gscbin_perl, $ber_source, $options, '>&'.$ber_output_log);
    #execute ber
    my $rv = Genome::Sys->shellcmd(
        cmd => $cmd,
    );
    unless ($rv) {
        die "Failed to execute blast2btab command: $cmd!";
    }

    my $ber_dat_file = $ber_output_file_key.'.dat';
    unless (-e $ber_dat_file) {
        die "failed to create BER output dat file\n";
    }

    $self->output_file($ber_dat_file);
    return 1;
}




1;
