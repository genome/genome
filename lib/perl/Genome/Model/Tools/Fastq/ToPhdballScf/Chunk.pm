package Genome::Model::Tools::Fastq::ToPhdballScf::Chunk;

use strict;
use warnings;

use PP;
use Genome;
use File::Temp;
use File::Path;
use File::Basename;
use Bio::SeqIO;
use Bio::Seq::Quality;


class Genome::Model::Tools::Fastq::ToPhdballScf::Chunk {
    is  => 'Genome::Model::Tools::Fastq::ToPhdballScf',
    has => [
        chunk_size => {
            is  => 'Integer',
            doc => 'number of each chunk fastq',
        },
    ],
};


sub help_brief {
    "Parallel run fastq-to-phdball-scf as chunk through LSF queue to handle the big fastq file"
}


sub help_detail {
    return <<EOS

EOS
}

sub execute {
    my $self = shift;

    my $ball_file = $self->ball_file;
    my $ball_dir  = dirname($self->ball_file);

    my $cmd = 'gmt fastq to-phdball-scf';

    for my $property (qw(time scf_dir base_fix solexa_fastq)) {
        if ($self->$property) {
            my $opt_name = $property;
            $opt_name =~ s/_/-/g;
            $cmd .= " --$opt_name";

            unless ($property =~ /^(base_fix|solexa_fastq)$/) {
                my $prop_val = $self->$property;
                $prop_val = '"'.$prop_val.'"' if $property eq 'time';
                $cmd .= ' '.$prop_val;
            }
        }
    }
    my $fq_split_cmd = Genome::Model::Tools::Fastq::Split->create(
        fastq_file       => $self->fastq_file,
        split_size       => $self->chunk_size,
        output_directory => $ball_dir,
    );
    unless ($fq_split_cmd->execute) {
        die $self->error_message("Failed to split fastq file.");
    }

    my %jobs;
    my @ball_files;
    my $fq_split_dir = $fq_split_cmd->output_directory;
    my $bl_ct = 0;

    for my $fq_split_file ($fq_split_cmd->split_files) {
        unless (-s $fq_split_file) {
            $self->error_message("fastq split file: $fq_split_file not existing");
            return;
        }
        $bl_ct++;
        my $ball_chunk_file = $ball_dir.'/phd.ball.'.$bl_ct;
        push @ball_files, $ball_chunk_file;

        my $job = $self->_lsf_job($cmd, $fq_split_file, $ball_chunk_file);
        $jobs{$bl_ct} = $job;
    }

    map{$_->start()}values %jobs;

    my %run;
    map{$run{$_}++}keys %jobs;

    while (%run) {
        sleep 120;
        for my $id (sort{$a<=>$b}keys %run) {
            my $job = $jobs{$id};

            if ($job->has_ended) {
                delete $run{$id};
                if ($job->is_successful) {
                    $self->status_message("Job $id done");
                }
                elsif ($job->has_exited) {
                    $self->warning_message("Job $id : ".$job->command." exited");
                }
                else {
                    $self->error_message("Job $id : ".$job->command." ended with neither DONE nor EXIT");
                }
            }
        }
    }
    my $files = join ' ', @ball_files;
    system "cat $files > $ball_file";
    map{unlink $_}@ball_files;
    rmtree $fq_split_dir;

    return 1;
}


sub _lsf_job {
    my ($self, $command, $fq_chunk_file, $ball_chunk_file) = @_;

    $command .= ' --ball-file '.$ball_chunk_file;
    $command .= ' --fastq-file '.$fq_chunk_file;

    return PP->create(
        pp_type => 'lsf',
        command => $command,
        q       => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        J       => basename($ball_chunk_file),
    );
}

1;
