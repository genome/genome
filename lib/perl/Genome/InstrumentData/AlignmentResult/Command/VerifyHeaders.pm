package Genome::InstrumentData::AlignmentResult::Command::VerifyHeaders;
use warnings;
use strict;
use Genome;

class Genome::InstrumentData::AlignmentResult::Command::VerifyHeaders {
    is => 'Genome::Command::Base',
    #has_many_input => [
    #    alignments => { 
    #        is => 'Genome::InstrumentData::AlignmentResult', 
    #        id_by => 'alignment_id',
    #        shell_args_position => 1, 
    #        doc => 'the alignment result, specified by properties'
    #    },
    #],
    has => [
        alignment_spec => { shell_args_position => 1 },
    ],
    doc => 'Generate a set of bitmask data for a given specific reference sequence.'
};

sub help_synopsis {
    my $class = shift;
    return <<EOS;
genome instrument-data alignment-result verify-headers 123456 
genome instrument-data alignment-result verify-headers reference_build_id=56789
EOS
}

sub help_detail {
    my $class = shift;
    return <<'EOS';
Verify that the headers in the BAM file match those in the reference sequence.
EOS
}

sub execute {
    my $self = shift;


    my $spec = $self->alignment_spec;
    my @p = eval $spec;
    my @a = sort { $b->id cmp $a->id } Genome::InstrumentData::AlignmentResult->get(@p);

    

    unless (Genome::Config->arch_os() =~ /64/) {
        $self->error_message("Run me on a 64-bit machine.");
        return;
    }
    
    #my @a = sort { $b->id <=> $a->id } $self->alignments;
    $self->status_message(sprintf("Found %s alignments...",scalar(@a)));
    
    my %reference_build_header_count;
    my @bad;
    for my $a (@a) {
        $self->status_message(sprintf("checking alignment %s...",$a->__display_name__));
        
        my $reference_build_id = $a->reference_build_id;
        unless (exists $reference_build_header_count{$reference_build_id}) {
            my $build = $a->reference_build;
            my $data_directory = $build->data_directory;
            my $path = $data_directory . '/all_sequences.fa.fai';
            my $cnt = `wc -l $path`;
            die $! unless $cnt;
            chomp $cnt;
            $cnt =~ s/ .*//;
            $self->status_message("expected count for $data_directory is $cnt\n");
            $reference_build_header_count{$reference_build_id} = $cnt;
        }
        my $expected_count = $reference_build_header_count{$reference_build_id};

        my $dir = $a->output_dir;
        my $bam = $dir . '/all_sequences.bam';
        my $cnt = `samtools view -H $bam | grep '^\@SQ' | wc -l `;
        chomp $cnt;
        $cnt =~ s/ .*//;

        my $msg;
        if ($cnt == $expected_count) {
            $msg = "OK";
        }
        else {
            push @bad, $a;
            $msg = "BAD";
        }

        $self->status_message(
            sprintf(
                " alignment %s has header size %s with expected %s %s",
                $a->__display_name__,
                $cnt,
                $expected_count,
                $msg,
            )
        );
    }

    $self->status_message(sprintf("Found %s bad alignments.",scalar(@bad)));
    for my $bad (@bad) {
        print $bad->id,"\n";
    }
    
    return 1;
}

