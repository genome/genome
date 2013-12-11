package Genome::Model::Tools::ViromeScreening;

use strict;
use warnings;

use Genome;
use Command;
use Workflow::Simple;
use Data::Dumper;
use File::Basename;
use Mail::Sender;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is => 'Command',
    has => [
        fasta_file => {
            doc => 'file of reads to be checked for contamination',
            is => 'String',
            is_input => 1,
        },
        barcode_file => { 
            doc => 'list of samples for screening',
            is => 'String',
            is_input => 1,
        },
        dir => {
            doc => 'directory of inputs',
            is => 'String',
            is_optional => 1,
            default => $ENV{"PWD"},
        },
        logfile => {
            doc => 'output file for monitoring progress of pipeline',
            is => 'String',
            is_optional => 1,
            default => "logfile.txt",
        },
        human_db => {
            doc => 'human blast db',
            is => 'String',
            is_optional => 1,
            default => '/gscmnt/sata835/info/medseq/virome/blast_db/human_genomic/2009_07_09.humna_genomic',
        },
        nt_db => {
            doc => 'nt sequence blast db',
            is => 'String',
            is_input => 1,
            is_optional => 1,
            default => '/gscmnt/sata835/info/medseq/virome/blast_db/nt/nt',
        },
        virus_db => {
            doc => 'Virus sequence blast db',
            is => 'String',
            is_input => 1,
            is_optional => 1,
            default => '/gscmnt/sata835/info/medseq/virome/blast_db/viral/viral.genomic.fna',
        },
        taxonomy_db => {
            doc => 'taxonomy db',
            is => 'String',
            is_input => 1,
            is_optional => 1,
            default => '/gscmnt/sata835/info/medseq/virome/taxonomy_db_2010_03_16',
        },
    ],
);

sub help_brief {
    "Runs virome screening workflow";
}

sub help_detail {
    'Runs the virome screening pipeline, using ViromeEvent modules.  Takes directory path, fasta, sample log, and logfile';
}

sub execute {
    my $self = shift;
    unlink($self->logfile) if (-e $self->logfile);

    $self->_log_dbs_used;
    
    $self->_validate_dbs;

    my $rv = $self->_run_workflow;

    if ( $rv ) {
        $self->_send_succeeded_mail;
        return 1;
    } else {
        $self->_send_failed_mail;
        return;
    }

    return 1;
}

sub _run_workflow {
    my $self = shift;

    my $rv = run_workflow_lsf(
	$ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ViromeScreening/virome-screening6.xml',
	'fasta_file'  => $self->fasta_file,
	'barcode_file'=> $self->barcode_file,
	'dir'         => $self->dir,
	'logfile'     => $self->logfile,
        'human_db'    => $self->human_db,
        'nt_db'       => $self->nt_db,
        'virus_db'    => $self->virus_db,
        'taxonomy_db' => $self->taxonomy_db,
    );

    return if not $rv;

    return 1;
}

sub _send_failed_mail {
    my $self = shift;

    my $log_file = $self->logfile;
    my $reason = `tail -5 $log_file`;

    my $msg = "Virome screening completed for run:\n".
              "\tfasta file: ".$self->fasta_file."\n".
              "\tbarcode file: ".$self->barcode_file."\n".
              "\tlog file: ".$self->logfile."\n".
              "\tdir: ".$self->dir."\n\n".
              "Reason from logfile:\n\n$reason\n";
                  
    my $mail_dest = Genome::Config->user_email;

    my $sender = Mail::Sender->new({
        smtp => 'gscsmtp.wustl.edu',
        from => $ENV{GENOME_EMAIL_VIROME_SCREENING},
        replyto => $ENV{GENOME_EMAIL_VIROME_SCREENING},
    });
    $sender->MailMsg({
        to => $mail_dest,
        subject => "Virome Screening Failed",
        msg     => $msg,
    });
}

sub _send_succeeded_mail {
    my $self = shift;

    my $msg = "Virome screening completed for run:\n".
              "\tfasta file: ".$self->fasta_file."\n".
              "\tbarcode file: ".$self->barcode_file."\n".
              "\tlog file: ".$self->logfile."\n".
              "\tdir: ".$self->dir."\n";

    my $mail_dest = Genome::Config->user_email;
    my $sender = Mail::Sender->new({
        smtp => $ENV{GENOME_EMAIL_SMTP_SERVER},
        from => $ENV{GENOME_EMAIL_VIROME_SCREENING},
        replyto => $ENV{GENOME_EMAIL_VIROME_SCREENING},
    });
    $sender->MailMsg({
        to => $mail_dest,
        subject => "Virome Screen completed",
        msg     => $msg,
    });
}

sub _log_dbs_used {
    my $self = shift;
    unlink $self->logfile;
    my $fh = Genome::Sys->open_file_for_writing( $self->logfile );
    for my $name ( qw/ human nt virus taxonomy/ ) {
        my $db = $name.'_db';
        $fh->printf("%-15s%40s\n", uc $name.' db:', $self->$db);
    }
    $fh->close;
}

sub _validate_dbs {
    my $self = shift;
    for my $db_name ( qw/ human nt virus / ) {
        my $db_param_name = $db_name.'_db';
        if( not -s $self->$db_param_name ) {
            die 'Failed to find db file: '.$self->$db_param_name."\n";
        }
    }
}

1;
