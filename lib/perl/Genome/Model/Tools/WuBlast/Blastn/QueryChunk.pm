package Genome::Model::Tools::WuBlast::Blastn::QueryChunk;

use strict;
use warnings;

use PP;
use Genome;
use File::Path;
use File::Basename;

class Genome::Model::Tools::WuBlast::Blastn::QueryChunk {
    is  => 'Genome::Model::Tools::WuBlast::Blastn',
    has => [
        chunk_size => {
            is  => 'Integer',
            doc => 'Number of query fasta for each chunk',
        },
    ],
    has_optional => [
        chunk_dir  => {
            is  => 'String',
            doc => 'Root directory of temp dir to hold chunk fasta files, default is dir of query_file',
        },
    ],
};


sub help_brief {
    'Divide query fasta into chunk and run Wu-Blast on each chunk' 
}


sub help_detail {  
    return <<EOS
EOS
}

    
sub execute {
    my $self = shift;

    unless ($self->params) {
        $self->error_message("Must have --params option to run WuBlast_querychunk");
        return;
    }

    my $fa_chunk = Genome::Model::Tools::Fasta::Chunk->create(
        fasta_file => $self->query_file,
        chunk_size => $self->chunk_size,
    );
    my $chunk_dir = $fa_chunk->chunk_dir;
    my $out = $fa_chunk->execute;

    unless ($out) {
        $self->error_message("Fasta_chunk fails");
        return;
    }

    my $out_file = $self->output_file;
    my $out_dir  = dirname $out_file;
    my $cmd = 'gmt wu-blast blastn --database '.$self->database.' --params "'.$self->params.'"';
    
    my @blast_outs;
    my %jobs;
    my $ct = 0;
    
    for my $chunk_fa_file (@{$out->[0]}) {
        my $chunk_name = basename $chunk_fa_file;
        my ($count) = $chunk_name =~ /(\d+)\.fasta$/; #hard coded chunk fasta file name Chunkxxx.fasta from Fasta::Chunk
        my $blast_out  = $out_dir."/$chunk_name.blast";
        my $command = $cmd . " --query-file $chunk_fa_file --output-file $blast_out";

        my $pp = PP->create(
            pp_type => 'lsf',
            command => $command,
            q       => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
            J       => $chunk_name,
        );    
        
        $jobs{$count} = $pp;
        push @blast_outs, $blast_out;
    }
    
    my %run;
    map{$run{$_}++}keys %jobs;
    map{$_->start}values %jobs;

    while (%run) {
        sleep 10;
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
    my $files = join ' ', @blast_outs;
    system "cat $files > $out_file";
        
    map{unlink $_}@blast_outs;
    rmtree $chunk_dir;
    
    return 1;
}
        
1;

#$HeadURL$
#$Id$
