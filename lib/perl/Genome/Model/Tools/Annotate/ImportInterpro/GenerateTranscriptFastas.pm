package Genome::Model::Tools::Annotate::ImportInterpro::GenerateTranscriptFastas;

use strict;
use warnings;
use Genome;
use Benchmark qw(:all) ;
use Bio::Seq;
use Bio::SeqIO;
use File::Temp;

class Genome::Model::Tools::Annotate::ImportInterpro::GenerateTranscriptFastas{
    is => 'Genome::Model::Tools::Annotate',
    has => [
        build => {
            is => 'Genome::Model::Build',
            is_input => 1,
            is_optional => 0,
        },
        chunk_size => {
            is => 'Number',
            is_input => 1,
            is_optional => 1,
            default => 25000,
            doc => 'Number of sequences per fasta. Cannot exceed 50000. Defaults to 25000',
        },
        benchmark => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            is_input => 1,
            doc => 'if set, run times are displayed as status messages after certain steps are completed (x, y, z, etc)',
        },
        scratch_dir => { 
            is => 'Path',
            is_input => 1,
            doc => 'files for fasta generation, iprscan output, etc. are written to this directory'
        },
    ]
};

#TODO: Write me
sub help_synopsis {
    return <<EOS
TODO
EOS
}

#TODO: Write me
sub help_detail{
    return <<EOS
TODO
EOS
}

sub execute{
    my $self = shift;

    #Generate .fasta files frome the build to be submitted wtih iprscan 
    my $pre_fasta_generation = Benchmark->new;

    my $build = $self->build;
    my $transcript_iterator = $build->transcript_iterator;
    die "Could not get iterator" unless $transcript_iterator;
    my $chunk_size = $self->chunk_size;
    die "Could not get chunk-size $chunk_size" unless $chunk_size; 
    die "chunk-size of $chunk_size is invalid.  Must be between 1 and 50000" if($chunk_size > 50000 or $chunk_size < 1);
    my $scratch_dir = $self->scratch_dir;
    die "Could not get tmp directory $scratch_dir" unless $scratch_dir; #TODO: Sanity check this

    # db disconnect to avoid Oracle failures killing our long running stuff
    if (Genome::DataSource::GMSchema->has_default_handle) {
        $self->debug_message("Disconnecting GMSchema default handle.");
        Genome::DataSource::GMSchema->disconnect_default_dbh();
    }

    my %fastas; 
    my ($fasta_temp, $fasta, $fasta_writer);
    my $transcript_counter = 0; 
    while (my $transcript = $transcript_iterator->next){
        if ($transcript_counter >= $chunk_size or not defined $fasta_temp){
            if (defined $fasta_temp){
                $fastas{$fasta} = $fasta_temp;
            }
            $fasta_temp = File::Temp->new(UNLINK => 0, 
                                          DIR => $scratch_dir,
                                          TEMPLATE => 'import-interpro_fasta_XXXXXX');
            $fasta = $fasta_temp->filename; 
            $fasta_writer = Bio::SeqIO->new(-file => ">$fasta", -format => 'fasta');
            die "Could not get fasta writer" unless $fasta_writer;
            $transcript_counter = 0;
        }
        my $protein = $transcript->protein;
        next unless $protein;
        my $amino_acid_seq = $protein->amino_acid_seq;
        my $bio_seq = Bio::Seq->new(-display_id => $transcript->transcript_name,
                                    -seq => $amino_acid_seq);
        $fasta_writer->write_seq($bio_seq);
        $transcript_counter++; 
    }

    unless (exists $fastas{$fasta}){
        $fastas{$fasta} = $fasta_temp;
    }

    my $post_fasta_generation = Benchmark->new;
    my $fasta_generation_time = timediff($post_fasta_generation, $pre_fasta_generation);
    $self->status_message('.fasta generation: ' . timestr($fasta_generation_time, 'noc')) if $self->benchmark;
    return 1;
}
1;
