package Genome::Model::Tools::Bed::Convert::Snv::MutectToBed;

use strict;
use warnings;

use Genome;
use Genome::Info::IUB;

class Genome::Model::Tools::Bed::Convert::Snv::MutectToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Snv'],
    has_param => [
        limit_variants_to => { 
            is => 'Text', 
            valid_values => ['hq','lq'], 
            is_optional => 1,
            doc => 'set to "hq" or "lq" to only get variants which pass or fail filter, respectively',
        },
    ],
};

sub help_brief {
    "Tools to convert Mutect SNV format to BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert snv mutect-to-bed --source snps_all_sequences --output snps_all_sequences.bed
EOS
}

sub help_detail {                           
    return <<EOS
    This is a small tool to take SNV calls in raw Mutect format and convert them to a common BED format (using the first five columns).
EOS
}

=cut 
## muTector v1.0.47986
contig  position    context ref_allele  alt_allele  tumor_name  normal_name score   dbsnp_site  covered power   tumor_power normal_power    total_pairs improper_pairs  map_Q0_reads    t_lod_fstar tumor_f contaminant_fraction    contaminant_lod t_ref_count t_alt_count t_ref_sum   t_alt_sum   t_ref_max_mapq  t_alt_max_mapq  t_ins_count t_del_count normal_best_gt  init_n_lod  n_ref_count n_alt_count n_ref_sum   n_alt_sum   judgement
MT  73  GGTxTGC A   G   H_KA-452198-0912806 H_KA-452198-1227537 0   NOVEL   COVERED 1   1   1   1982    708 0   1887.387697 1   0.02    40.615705   0   495 0   16516   0   60  0   0   GG  -3537.592164    1   977 22  30754   REJECT
MT  150 CATxCTA C   T   H_KA-452198-0912806 H_KA-452198-1227537 0   NOVEL   COVERED 1   1   1   1996    488 0   922.690548  1   0.02    40.469904   0   238 0   8093    0   60  0   0   TT  -3609.916026    3   970 84  31577   REJECT
=cut

sub process_source {
    my $self = shift;
    
    my $input_fh = $self->_input_fh;

    my $judgement_string;
    if(defined $self->limit_variants_to) {
        if($self->limit_variants_to eq 'hq') {
            $judgement_string = 'KEEP';
        }
        else {
            $judgement_string = 'REJECT';
        }
    }

    
    my $fn;
    my @headers;
    while(my $line = <$input_fh>) {
        next if $line =~ /^#/;
        chomp $line;
        my @fields = split("\t", $line);
#        die "Non-comment line without 33 columns found.\n" unless @fields == 34;
        my %entry;
        if($fields[0] eq 'contig') {
            @headers = @fields;
            next;
        }
        else {
            @entry{@headers} = @fields;
        }
        $self->write_bed_line(
            $entry{contig},
            $entry{position} - 1,
            $entry{position},
            $entry{ref_allele},
            Genome::Info::IUB::iub_for_alleles($entry{ref_allele}, $entry{alt_allele}),
            int($entry{t_lod_fstar} + 0.5),
            $entry{t_ref_count} + $entry{t_alt_count},
        ) if(!defined $judgement_string || $entry{judgement} eq $judgement_string);
    }
    
    return 1;
}

1;

