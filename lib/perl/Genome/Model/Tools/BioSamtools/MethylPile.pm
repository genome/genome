package Genome::Model::Tools::BioSamtools::MethylPile;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BioSamtools::MethylPile {
    is => ['Genome::Model::Tools::BioSamtools'],
    has => [
        roi_file_path => {
            doc => 'The expected CpG sites.  Only C nucleotides will be evaluated.',
        },
        bam_file_path => {
            doc => 'The aligned reads in BAM format',
        },
        reference_fasta => {
            doc => 'The FASTA reference genome.',
        },
        print_region => {
            doc => 'This prints the entire region that overlaps each ROI.  When false, only the exact ROI position will be printed.',
            is => 'Boolean',
            default_value => 1,
        },
    ],
};

# Output tsv format(all positons in 1-based coordinates
# [1] ROI Name
# [2] ROI Id ($name:$start-$end)
# [3] Pileup Start
# [4] Pileup End
# [5] Pileup Position
# [6] Count Cs
# [7] Count Ts
# [8] Count Gs
# [9] Count As
# [10] Count Ns

sub execute {
    my $self = shift;

    # Load the ROI file
    my $roi_file = Genome::Model::Tools::RefCov::ROI::Bed->create(
        file => $self->roi_file_path,
    );

    # Load the BAM file
    my $bam = Genome::Model::Tools::RefCov::Bam->create(bam_file =>$self->bam_file_path);
    unless ($bam) {
        die('Failed to load BAM file: '. $self->bam_file_path);
    }

    my $bio_db_bam = $bam->bio_db_bam;
    my $bio_db_index = $bam->bio_db_index;

    my $header = $bio_db_bam->header;
    my $tid_names = $header->target_name;
    my $tid_lengths = $header->target_len;

    # Load the reference genome
    my $fai = Bio::DB::Sam::Fai->load($self->reference_fasta);
    unless ($fai) {
        die('Failed to load fai for: '. $self->reference_fasta);
    }

    # The subroutine to count C/T nucleotides in the query at all reference C positions
    my $callback = sub {
        my ($tid,$pos,$pileups,$callback_data) = @_;

        my ($ct_count_hash_ref,$region,$pileup) = @{$callback_data};

        # set boundary positions of the pileup "cluster"
        if (!defined($pileup->{start}) && !defined($pileup->{end})) {
            $pileup->{start} = $pos;
            $pileup->{end} = $pos;
        } elsif ($pos < $pileup->{start}) {
            $pileup->{start} = $pos;
        } elsif ($pos > $pileup->{end}) {
            $pileup->{end} = $pos;
        }

        my $chr = $tid_names->[$tid];
        my $ref_base = $fai->fetch($chr .':'. ($pos+1) .'-'. ($pos+1) );
        # hash the number of C/T query nucleotides at each reference C position
        if ($ref_base eq 'C') { 
            for my $pileup (@$pileups) {
                if ($pileup->indel || $pileup->is_refskip) { next; }
                my $b = $pileup->alignment;
                my $qpos = $pileup->qpos;
                my $qbase = substr($b->qseq,$qpos,1);
                $ct_count_hash_ref->{$pos}->{$qbase}++;
            }
        }
    };

    # Loop over each ROI 
    while (my $region = $roi_file->next_region) {
        my %nuc_counts;
        my %pileup;
        my $tid = $bam->tid_for_chr($region->{chrom});
        $bio_db_index->pileup($bio_db_bam,$tid,($region->{start} -1),$region->{end},$callback,[\%nuc_counts,$region,\%pileup]);
        # For all reference C positions, print a count of C/T nucleotides
        if ($self->print_region) {
            if (defined($pileup{start}) && defined($pileup{end})) {
                for my $pos ($pileup{start} .. $pileup{end}) {
                    if ($nuc_counts{$pos}) {
                        print $region->{name} ."\t". $region->{id} ."\t". ($pileup{start} +1) ."\t". ($pileup{end} + 1) ."\t". ($pos+1)
                            ."\t". ($nuc_counts{$pos}->{C} || 0) ."\t". ($nuc_counts{$pos}->{T} || 0)
                                ."\t". ($nuc_counts{$pos}->{G} || 0) ."\t". ($nuc_counts{$pos}->{A} || 0)
                                    ."\t". ($nuc_counts{$pos}->{N} || 0)
                                        ."\n";
                    }
                }
            } else {
                # NO COVERAGE
            }
        } else {
        # Only print the C/T nucleotides at the exact CpG site
            my $pos = ($region->{start} - 1);
            if ($nuc_counts{$pos}) {
                print $region->{name} ."\t". $region->{id} ."\t". $region->{start} ."\t". $region->{end} ."\t". $region->{start}
                    ."\t". ($nuc_counts{$pos}->{C} || 0) ."\t". ($nuc_counts{$pos}->{T} || 0)
                        ."\t". ($nuc_counts{$pos}->{G} || 0) ."\t". ($nuc_counts{$pos}->{A} || 0)
                            ."\t". ($nuc_counts{$pos}->{N} || 0)
                                ."\n";
            } else {
                # NO COVERAGE
            }
        }
    }

    return 1;
}


1;
