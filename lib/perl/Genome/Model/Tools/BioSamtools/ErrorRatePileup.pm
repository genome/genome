package Genome::Model::Tools::BioSamtools::ErrorRatePileup;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BioSamtools::ErrorRatePileup {
    is => ['Genome::Model::Tools::BioSamtools'],
    has => [
        bam_file => {},
        reference_fasta => {},
    ],
    has_optional => [
        output_file => {},
    ],
};

sub execute {
    my $self = shift;
    $DB::single=1;
    my $fai = Bio::DB::Sam::Fai->load($self->reference_fasta);
    unless ($fai) {
        die('Failed to load fai for: '. $self->reference_fasta);
    }

    my $bam = Bio::DB::Bam->open($self->bam_file);
    unless ($bam) {
        die('Failed to load BAM file: '. $self->bam_file);
    }

    my $index  = Bio::DB::Bam->index($self->bam_file);
    unless ($index) {
        die('Failed to load BAM index for: '. $self->bam_file);
    }

    my $header = $bam->header;
    my $tid_names = $header->target_name;
    my $tid_lengths = $header->target_len;

    my %read_counts;
    my %del_pileups;
    my $callback = sub {
        my ($tid,$pos,$pileups,$callback_data) = @_;
        my $chr = $tid_names->[$tid];
        my $ref_base = $fai->fetch($chr .':'. ($pos+1) .'-'. ($pos+1) );
        for my $pileup (@$pileups) {
            my $b = $pileup->alignment;
            my $qname = $b->qname;
            my $flag = $b->flag;
            # Unmapped
            if ($flag & 4) {
                next;
            }
            my $read_end = 0;
            if ($flag & 1) {
                if ($flag & 64) {
                    $read_end = 1;
                } elsif ($flag & 128) {
                    $read_end = 2;
                }
            }
            $qname .= '/'. $read_end;

            my $qpos = $pileup->qpos;
            $read_counts{$read_end}{total}->[$qpos]++;

            my $qbase = substr($b->qseq,$qpos,1);

            my $indel = $pileup->indel;
            if ($indel) {
                my $indel_size = abs($indel);
                if ($indel > 0) {
                    for (1 .. $indel_size) {
                        my $ipos = ($qpos + $_);
                        $read_counts{$read_end}{insertion}->[$ipos]++;
                        $read_counts{$read_end}{total}->[$ipos]++;
                    }
                } else {
                    $read_counts{$read_end}{deletion}->[$qpos] += $indel_size;
                    for (1 .. $indel_size) {
                        # This is relative to the reference position
                        my $dpos = ($pos + $_);
                        $del_pileups{$qname}{$dpos} = 1;
                    }
                }
            } elsif ($del_pileups{$qname}{$pos}) {
                my $del_pileup = delete($del_pileups{$qname}{$pos});
                $read_counts{$read_end}{total}->[$qpos]--;
                #$DB::single = 1;
                #print $pos ."\t". $ref_base ."\t". $qname ."\t". $del_pileup ."\t". $del_pileup ."\t". $pileup->qpos ."\t". $qbase ."\t". $pileup->indel ."\n";
                next;
            }
            if ($qbase =~ /[nN]/) {
                $read_counts{$read_end}{ambiguous}->[$qpos]++;
            } elsif ($ref_base ne $qbase) {
                $read_counts{$read_end}{mismatch}->[$qpos]++;
            } else {
                $read_counts{$read_end}{match}->[$qpos]++;
            }
        }
    };
    for (my $i = 0; $i < scalar(@{$tid_lengths}); $i++) {
        my $end = $tid_lengths->[$i];
        $index->pileup($bam,$i,'0',$end,$callback);
    }
    my @headers = qw/read_end position total match error error_rate mismatch mismatch_rate ambiguous ambiguous_rate insertion insertion_rate deletion deletion_rate/;
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        separator => "\t",
        headers => \@headers,
        output => $self->output_file,
    );
    unless ($writer) {
        die('Failed to create output writer!');
    }
    $DB::single=1;
    for my $read_end (sort keys %read_counts) {
        my @positions = @{$read_counts{$read_end}{total}};
        my $sum_total;
        my $sum_match;
        my $sum_mismatch;
        my $sum_ambiguous;
        my $sum_insertion;
        my $sum_deletion;
        my $sum_error;
        for (my $i = 0; $i < scalar(@positions); $i++) {
            my $position_count = $positions[$i];
            if (!$position_count) {
                my %data = (
                    read_end => $read_end,
                    position => $i,
                    total => 0,
                    match => 0,
                    error => 0,
                    error_rate => 0,
                    mismatch => 0,
                    mismatch_rate => 0,
                    ambiguous => 0,
                    ambiguous_rate => 0,
                    insertion => 0,
                    insertion_rate => 0,
                    deletion => 0,
                    deletion_rate => 0,
                );
                $writer->write_one(\%data);
                next;
            }
            my $match = $read_counts{$read_end}{match}->[$i] || 0;
            my $mismatch = $read_counts{$read_end}{mismatch}->[$i] || 0;
            my $ambiguous = $read_counts{$read_end}{ambiguous}->[$i] || 0;
            my $insertion = $read_counts{$read_end}{insertion}->[$i] || 0;
            my $deletion = $read_counts{$read_end}{deletion}->[$i] || 0;
            
            my $total = $match + $mismatch + $ambiguous + $insertion + $deletion;
            my $error = $mismatch + $ambiguous + $insertion + $deletion;
            my $error_rate = $error / $total;
            my $mismatch_rate = $mismatch / $total;
            my $ambiguous_rate = $ambiguous / $total;
            my $insertion_rate = $insertion / $total;
            my $deletion_rate = $deletion / $total;
            my %data = (
                read_end => $read_end,
                position => $i,
                total => $total,
                match => $match,
                error => $error,
                error_rate => $error_rate,
                mismatch => $mismatch,
                mismatch_rate => $mismatch_rate,
                ambiguous => $ambiguous,
                ambiguous_rate => $ambiguous_rate,
                insertion => $insertion,
                insertion_rate => $insertion_rate,
                deletion => $deletion,
                deletion_rate => $deletion_rate,
            );
            $writer->write_one(\%data);
            $sum_total += $total;
            $sum_match += $match;
            $sum_mismatch += $mismatch;
            $sum_ambiguous += $ambiguous;
            $sum_insertion += $insertion;
            $sum_deletion += $deletion;
            $sum_error += $error;
        }
        my $error_rate = $sum_error / $sum_total;
        my $mismatch_rate = $sum_mismatch / $sum_total;
        my $ambiguous_rate = $sum_ambiguous / $sum_total;
        my $insertion_rate = $sum_insertion / $sum_total;
        my $deletion_rate = $sum_deletion / $sum_total;
        my %data = (
            read_end => $read_end,
            position => 'SUM',
            total => $sum_total,
            match => $sum_match,
            error => $sum_error,
            error_rate => $error_rate,
            mismatch => $sum_mismatch,
            mismatch_rate => $mismatch_rate,
            ambiguous => $sum_ambiguous,
            ambiguous_rate => $ambiguous_rate,
            insertion => $sum_insertion,
            insertion_rate => $insertion_rate,
            deletion => $sum_deletion,
            deletion_rate => $deletion_rate,
        );
        $writer->write_one(\%data);
    }
    return 1;
}
