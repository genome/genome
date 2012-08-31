package Genome::Model::Tools::BioSamtools::SimulateRnaSeqReads;

use strict;
use warnings;

use Genome;
use Math::Random;

class Genome::Model::Tools::BioSamtools::SimulateRnaSeqReads {
    is => ['Genome::Model::Tools::BioSamtools'],
    has => [
        reference_fasta_file => {
            is => 'Text',
            doc => 'The reference fasta(not fai) indexed by samtools.',
        },
        annotation_gtf_file => {
            is => 'Text',
            doc => 'The transcript annotation GTF file.',
        },
        fastq_file => {
            is =>'Text',
            doc => 'The read 1 fastq output file',
        },
    ],
    has_optional => [
        read_length => {
            is => 'Integer',
            default_value => 100,
        },
        paired_end => {
            is => 'Boolean',
            default_value => 1,
        },
        fastq2_file => {
            is => 'Text',
            doc => 'The read 2 fastq output file',
        },
        mean_insert_size => {
            is => 'Boolean',
            default_value => '250',
        },
        mean_insert_size_sd => {
            is => 'Boolean',
            default_value => 20,
        },
        feature_types => {
            is => 'Text',
            default_value => 'exon',
        },
        approximate_depth => {
            is => 'Integer',
            default_value => 10,
        },
        # TODO: Error rate
        # The problem is simulating an Illumina error model with errors concentrated at the 3' end of the read

        # TODO: Adapter sequence/frequency
        # Again, the big issue is concentrating the adapter to a single end 5' or 3' without making it static

        # TODO: SNP frequency
        # Or, feed a dbsnp file?

        # TODO: indel frequency
        
        # TODO: introduce exon skipping or fusions...
        # This could work now if the GTF represented such mutations/editing events.
        # For the GTF to work, the transcript_id must match the two sides of the fusion
    ],
};

sub execute {
    my $self = shift;

    my $gtf_reader = Genome::Utility::IO::GffReader->create(
        input => $self->annotation_gtf_file,
    );
    unless ($gtf_reader) {
        die('Failed to load GTF file: '. $self->annotation_gtf_file);
    }
    my $fai = Bio::DB::Sam::Fai->load($self->reference_fasta_file);
    unless ($fai) {
        die('Failed to load fai index for FASTA file: '. $self->reference_fasta_file);
    }
    my $fastq_1_fh = Genome::Sys->open_file_for_writing($self->fastq_file);
    unless ($fastq_1_fh) {
        die('Failed to open fastq file for read1: '. $self->fastq_file);
    }
    my $fastq_2_fh;
    if ($self->paired_end) {
        $fastq_2_fh = Genome::Sys->open_file_for_writing($self->fastq2_file);
        unless ($fastq_2_fh) {
            die('Failed to open fastq file for read2: '. $self->fastq2_file);
        }
    }

    my @feature_types = split(',',$self->feature_types);
    my %transcripts;
    while (my $data = $gtf_reader->next_with_attributes_hash_ref) {
        my $keep = 0;
        for my $feature_type (@feature_types) {
            if ($data->{type} eq $feature_type) {
                $keep = 1;
                last;
            }
        }
        unless ($keep) { next; }
        my $attributes = delete($data->{attributes_hash_ref});
        push @{$transcripts{$attributes->{transcript_id}}{features}}, $data;
    }
    # Build transcript models
    for my $transcript (keys %transcripts) {
        my @transcript_dna;
        #TODO: Check for overlaps.  This should be biologically impossible, right?  unless CDS or UTR combined with exons...
        for my $feature (sort { $a->{start} <=> $b->{start} } @{$transcripts{$transcript}{features}} ) {
            my $id = $feature->{chr} .':'. $feature->{start} .'-'. $feature->{end};
            my $dna_string = $fai->fetch($id);
            $feature->{dna} = $dna_string;
            $feature->{length} = $feature->{end} - $feature->{start} + 1;
            if ( !defined($transcripts{$transcript}{ref_start}) || ($transcripts{$transcript}{ref_start} > $feature->{start}) ) {
                $transcripts{$transcript}{ref_start} = $feature->{start};
            }
            if ( !defined($transcripts{$transcript}{ref_end}) || ($transcripts{$transcript}{ref_end} < $feature->{end}) ) {
                $transcripts{$transcript}{ref_end} = $feature->{end};
            }
            push @transcript_dna, $dna_string;
        }
        my $dna_string = join('',@transcript_dna);
        $transcripts{$transcript}{dna} = $dna_string;
    }

    # Divide the read length by the approximate depth to get the offset used spacing reads out over the length of the reference ROI
    # Also, divide by two since we will make one pass on the forward strand and one pass on the reverse
    my $offset = int($self->read_length / $self->approximate_depth) * 2;

    my $quality_string = '';
    for (1 .. $self->read_length) {
        # TODO: hard coded sanger 30 quality value(make phred value a parameter and resolve quality value conversion)
        $quality_string .= '?';
    }

    #Make reads from transcript dna(includes splice junctions)
    for my $transcript (keys %transcripts) {
        my $positive_strand = $transcripts{$transcript}{dna};
        my $reverse_strand = reverse $positive_strand;
        $reverse_strand =~ tr/ACGTacgt/TGCAtgca/;
        my $length = length($positive_strand);
        my $strand_symbol = 0;
        # TODO: It might be worth adding some logic to get an even mix of R1<->R2 and R2<->R1 fragments from both strands
        # This might require some chalk talk with lab folk
        for my $strand ( ($positive_strand, $reverse_strand) ) {
            my $start_site_number = ( $length / $offset );
            my @insert_sizes;
            if ($self->paired_end) {
                # Divide by two for paired end
                $start_site_number /= 2;

                # NOTE: These are not integer values... yet... see below
                @insert_sizes = random_normal($start_site_number,$self->mean_insert_size,$self->mean_insert_size_sd);
            }
            my @start_sites = random_uniform_integer($start_site_number,0,($length - 1));
            for (my $i = 0; $i < scalar(@start_sites); $i++) {
                my $start_site = $start_sites[$i];

                my $read_1_offset = $start_site;
                my $read_1_end = $start_site + $self->read_length;
                if ($read_1_end > $length) { next; }

                my $read_1_dna = substr($strand,$read_1_offset,$self->read_length);
                # NOTE: Don't print to FASTQ just yet in case we have paired_end reads and read2 falls off the defined transcript
                my $fragment_start = ($read_1_offset + 1);
                my $fragment_end = $read_1_end;
                my $insert_size = 0;
                if ($self->paired_end) {
                    # NOTE:  now this we've got an integer insert size
                    $insert_size = int($insert_sizes[$i]);
                    my $read_2_end = $start_site + $insert_size;
                    if ($read_2_end > $length) { next; }
                    $fragment_end = $read_2_end;
                    my $read_2_offset =  $read_2_end - $self->read_length;
                    my $read_2_dna = substr($strand,$read_2_offset,$self->read_length);

                    my $revcom = reverse $read_2_dna;
                    $revcom =~ tr/ACGTacgt/TGCAtgca/;

                    print $fastq_2_fh '@'. $transcript .':'. $length .':'. $fragment_start .':'. $insert_size . ':'. $strand_symbol .'/2' ."\n";
                    print $fastq_2_fh $revcom ."\n";
                    print $fastq_2_fh '+' ."\n";
                    print $fastq_2_fh $quality_string ."\n";
                }
                # TODO: potentially resolve genomic coordinate start stop and use in read name....
                # This might require using each exon fragment rather than a squashed transcript sequence
                # That would make things much more difficult... I think...
                print $fastq_1_fh '@'. $transcript .':'. $length .':'. $fragment_start .':'. $insert_size .':'. $strand_symbol .'/1' ."\n";
                print $fastq_1_fh $read_1_dna ."\n";
                print $fastq_1_fh '+' ."\n";
                print $fastq_1_fh $quality_string ."\n";
            }
            $strand_symbol++;
        }
    }
    $fastq_1_fh->close;
    $fastq_2_fh->close;
    return 1;
}
