package Genome::Model::Tools::Predictor::Keggscan::BladeBlastBatcher;

# bdericks: This module was originally written by Todd Wylie and named 
# BladeBlastBatcher_lowmemory_LSF. As of August 27, 2010, the original
# module can be found at /gscmnt/233/analysis/sequence_analysis/lib. I 
# wasn't able to find it in version control anywhere. I've made some
# changes to turn this into a command object and put in the PAP namespace

use strict;
use warnings;

use Genome;
use PP::LSF;
use PP::JobScheduler;
use Benchmark;
use Carp 'confess';

class Genome::Model::Tools::Predictor::Keggscan::BladeBlastBatcher {
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

    $self->status_message("Splitting up sequences of " . $self->query_fasta_path . " into separate files in " . 
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
    $self->status_message("All jobs are done and fasta chunks have been cleaned up. Reports can be found at " . $self->output_directory . ".");
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

    $self->status_message("Starting a batch of $num_jobs jobs!");
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
            $self->status_message("Batch complete, took " . timestr($batch_time, 'noc'));
            last;
        }
    }

    $jobs = (); # Pretty sure that the job scheduler does this (shifts stuff off the job list)
    return 1;
}

1;
