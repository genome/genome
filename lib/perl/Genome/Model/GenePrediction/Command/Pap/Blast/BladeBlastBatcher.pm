package Genome::Model::GenePrediction::Command::Pap::Blast::BladeBlastBatcher;

# bdericks: This module was originally written by Todd Wylie and named 
# BladeBlastBatcher_lowmemory_LSF. As of August 27, 2010, the original
# module can be found at /gscmnt/233/analysis/sequence_analysis/lib. I 
# wasn't able to find it in version control anywhere. I've made some
# changes to turn this into a command object and put in the PAP namespace

use strict;
use warnings;

use PAP;
use PP::LSF;
use PP::JobScheduler;
use Benchmark;
use Carp 'confess';

class Genome::Model::GenePrediction::Command::Pap::Blast::BladeBlastBatcher {
    is => 'Command::V2',
    has => [
        query_fasta_path => {
            is => 'Path',
            doc => 'Path to the original query FASTA file',
        },
        subject_fasta_path => {
            is => 'Path',
            doc => 'Path to the subject FASTA file',
        },
        lsf_queue => {
            is => 'Text',
            doc => 'The name of the LSF queue which the jobs will be submitted to',
        },
        lsf_job_limit => {
            is => 'Number',
            doc => 'The number of jobs to split the session into',
        },
        output_directory => {
            is => 'Path',
            doc => 'Output directory where all output will be directed',
        },
        blast_name => {
            is => 'Text',
            doc => 'Name of the BLAST type (e.g. tblastn, blastn, blastx)',
        },
    ],
    has_optional => [
        reports => {
            is => 'ARRAY',
            doc => 'Array of results produced by this module',
        },
        blast_params => {
            is => 'Text',
            default => "",
            doc => 'A string reserved for BLAST parameters',
        },
        lsf_resource => {
            is => 'Text',
            default => 'select[mem>1024] rusage[mem=1024]',
            doc => 'Resource string for each scheduled LSF job',
        },
        lsf_max_memory => {
            is => 'Number',
            default => '1024000',
            doc => 'Maximum allowable memory usage for each scheduled LSF job',
        },
        fasta_chunk_size => {
            is => 'Number',
            default => 20,
            doc => 'Maximum number of sequences allowed in each fasta chunk',
        },
    ],
};

sub help_detail {
    return <<EOS
PURPOSE:
This module was written to help developers who wish to smash a large
FASTA file into smaller sections and then run BLAST on them via the blade
center (qsub). Other scripts aid in doing this, but this module was
designed to accomplish the job inside of another, longer perl
application. Also, routines are supplied that help determine if the jobs
are finished on the blades, so the parent perl script may continue.
EOS
}

# This routine (based on incoming values) will cut a given query FASTA file
# into XX number of smaller FASTA files, then run the chosen BLAST program
# on them by farming the jobs out to the blades. BLAST reports are written
# to a specified output directory. When all jobs are finished on the
# blades, this module sends the calling application back a list of report
# file paths.
sub execute {
    my $self = shift;
    my @reports;
    my @chunks;
    my @jobs;
    my $sequence_counter = 0;
    my $fasta_chunk;

    $self->debug_message("Splitting up sequences of " . $self->query_fasta_path . " into separate files in " . 
        $self->output_directory . " containing no more than " . $self->fasta_chunk_size . 
        " sequences each and executing in batches of " . $self->lsf_job_limit);

    my $query_fasta = Bio::SeqIO->new(
        -file => $self->query_fasta_path,
        -format => 'fasta',
    );

    my $chunk_counter = 0;
    while (my $seq = $query_fasta->next_seq()) {
        # Create a new fasta chunk!
        if (not defined $fasta_chunk or $sequence_counter >= $self->fasta_chunk_size) {

            # Run and delete chunks (if our maximum number of jobs limit has been reached)
            if ((defined $fasta_chunk) and ((scalar @chunks) % $self->lsf_job_limit == 0)) {
                $self->_run_lsf_jobs(\@jobs); # This blocks and waits for the jobs to finish

                # Keep the number of fasta chunks to a reasonable level. This process can end up scheduling
                # thousands of jobs, and having a directory full of thousands of files is not cool.
                for my $chunk (@chunks) {
                    unlink $chunk;
                }
                @chunks = ();
            }

            # Create new fasta chunk, output report path, and LSF job object
            my $chunk_file_name = $self->output_directory . "/CHUNK-$chunk_counter.fasta";
            $fasta_chunk = Bio::SeqIO->new(
                -file => ">$chunk_file_name",
                -format => 'fasta',
            );
            $sequence_counter = 0;
            $chunk_counter++;
            push @chunks, $fasta_chunk;

            my $report = $chunk_file_name . ".blast.report";
            push @reports, $report;

            my $cmd;
            if ($self->blast_name eq "wu-blastall") {
                $cmd = $self->blast_name . " -d " . $self->subject_fasta_path . 
                    " -i $chunk_file_name " . $self->blast_params . " > $report";
            } 
            else {
                $cmd = join(" ", $self->blast_name, $self->subject_fasta_path, $chunk_file_name, $self->blast_params) . 
                    " > $report";
            }
            my $job = PP::LSF->create(
                'q'           => $self->lsf_queue,
                'o'           => $self->output_directory,
                'e'           => $self->output_directory,
                'R'           => "'" . $self->lsf_resource . "'",
                'M'           => $self->lsf_max_memory,
                'command'     => "\"" . $cmd . "\"",
            );
            push @jobs, $job;
        }

        $fasta_chunk->write_seq($seq);
        $sequence_counter++;
    }

    # Catch any jobs that didn't get executed above
    $self->_run_lsf_jobs(\@jobs); 

    # Clean up any remaining fasta chunks
    for my $chunk (@chunks) {
        unlink $chunk;
    }

    $self->reports(\@reports);
    $self->debug_message("All jobs are done and fasta chunks have been cleaned up. Reports can be found at " . $self->output_directory . ".");
    return 1;
}

sub _run_lsf_jobs {
    my ($self, $jobs) = @_;

    my $num_jobs = scalar @$jobs;

    # Jobs are run and throttled by this scheduler
    my $scheduler = new PP::JobScheduler(
        pp_type => 'lsf',
        job_list => $jobs,
        default_max => $self->lsf_job_limit,
        refresh_interval => 60,
    );
    $scheduler->start();

    $self->debug_message("Starting a batch of $num_jobs jobs!");
    my $batch_start = Benchmark->new;

    # Wait for jobs to finish
    while(1) {
        sleep 60;
        # bdericks: It seems that this doesn't get called by default when 
        # calling running jobs, which can lead to an infinite loop...
        $scheduler->_update_running_jobs();  
        my @running_jobs = @{$scheduler->running_jobs()};
        unless (@running_jobs) {
            my $batch_stop = Benchmark->new;
            my $batch_time = timediff($batch_stop, $batch_start);
            $self->debug_message("Batch complete, took " . timestr($batch_time, 'noc'));
            last;
        }
    }

    $jobs = (); # Pretty sure that the job scheduler does this (shifts stuff off the job list)
    return 1;
}
# ---------------------------------------------------------------------------
# D O C U M E N T A T I O N
# ---------------------------------------------------------------------------

=head1 NAME

I<BladeBlastBatcher.pm>

=head1 VERSION

version 1.0 [April 2005]

=head1 DESCRIPTION

This module was written to help developers who wish to smash a large FASTA file into smaller sections and then run BLAST on them via the blade center (qsub). Other scripts aid in doing this, but this module was designed to accomplish the job inside of another, longer perl application. Also, routines are supplied that help determine if the jobs are finished on the blades, so the parent perl script may continue.

=head1 SYNOPSIS

To perform a simple call from a perl script:

  my @reports = &Genome::Model::GenePrediction::Command::Pap::Blast::BladeBlastBatcher->execute(
                          query        => $query,
                          subject      => $revised_subject,
                          queue        => "compbio\@qblade",
                          load         => $blades,
                          outdir       => $SessionOutDir,
                          blast_name   => "blastx",
                          blast_params => "-b 10000",
                          mailoptions  => "abe",
                          mailto       => "twylie\@watson.wustl.edu"
                                               );


CONFIGURATION:

=over 5

=item 1

query:        The path to the original query FASTA file.				 

=item 2

subject:      The path to the subject FASTA file.					 

=item 3

queue:        The name of the blade queue which the jobs will be submitted to.	 

=item 4

load:         The number of jobs to split the session into.			 

=item 5

outdir:       Output directory where all output will be directed.			 

=item 6

blast_name:   Name of the BLAST type (e.g. tblastn, blastn, blastx).		 

=item 7

blast_params: A string reserved for BLAST parameters [OPTIONAL].		 

=item 8

mailoptions:  Qsub switch reserved to provide mailing options [OPTIONAL].		 

=item 9

mailto:       Qsub switch reserved for mail address [OPTIONAL].                    

=back

NOTE: The following params are optional: blast_params, mailoptions, mailto. The mailoptions & mailto params must both be present to work correctly.

=head1 AUTHOR

 Todd Wylie
 CompBio
 Genome Sequencing Center
 Washington University
 School of Medicine
 4444 Forest Park Boulevard
 St. Louis, MO 63108

CONTACT: twylie@watson.wustl.edu

=head1 LIMITATION/BUGS

=over 5

=item 1

Mail function has not been tested properly as of Tue Apr 26 16:05:53 CDT 2005. This may be a qsub glitch--but should be resolved at a later date (time permitting).

=item 2

Only "blastx" has been tested for v1.0 release. Other blast types should be tested.

=back

Please contact the author if any complications are encountered.

=head1 COPYRIGHT

Copyright (C) 2005 by Todd Wylie and Washington University School of Medicine Genome Sequencing Center.

=head1 NOTE:

This software was written using the latest version of GNU Emacs, the extensible, real-time text editor. Please see http://www.gnu.org/software/emacs/ for more information and download sources.

=cut
