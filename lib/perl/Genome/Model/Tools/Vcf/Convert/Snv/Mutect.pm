package Genome::Model::Tools::Vcf::Convert::Snv::Mutect;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;
use File::Basename;

class Genome::Model::Tools::Vcf::Convert::Snv::Mutect {
    is  => 'Genome::Model::Tools::Vcf::Convert::Base',
    doc => 'Generate a VCF file from raw Mutect output',
    has_optional => [
       _headers => {
           is => 'ArrayRef',
       },
   ],

};

sub help_synopsis {
    <<'HELP';
    Generate a VCF file from mutect snv output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the snvs.
HELP
}

sub source {
    my $self = shift;
    return "muTect";
}

=cut
## muTector v1.0.47986
contig  position    context ref_allele  alt_allele  tumor_name  normal_name score   dbsnp_site  covered power   tumor_power normal_power    total_pairs improper_pairs  map_Q0_reads    t_lod_fstar tumor_f contaminant_fraction    contaminant_lod t_ref_count t_alt_count t_ref_sum   t_alt_sum   t_ref_max_mapq  t_alt_max_mapq  t_ins_count t_del_count normal_best_gt  init_n_lod  n_ref_count n_alt_count n_ref_sum   n_alt_sum   judgement
MT  73  GGTxTGC A   G   H_KA-452198-0912806 H_KA-452198-1227537 0   NOVEL   COVERED 1   1   1   1982    708 0   1887.387697 1   0.02    40.615705   0   495 0   16516   0   60  0   0   GG  -3537.592164    1   977 22  30754   REJECT
MT  150 CATxCTA C   T   H_KA-452198-0912806 H_KA-452198-1227537 0   NOVEL   COVERED 1   1   1   1996    488 0   922.690548  1   0.02    40.469904   0   238 0   8093    0   60  0   0   TT  -3609.916026    3   970 84  31577   REJECT
=cut

sub parse_line {
    my ($self, $line) = @_;
    return if $line =~ /^##/; # no vcf header here. It would be better to embed the version in the header, but maybe that should be part of the harness.
    
    chomp $line;
    my @fields = split("\t", $line);
    my %entry;
    if($fields[0] eq 'contig') {
        $self->_headers(\@fields);
        return;
    }
    else {
        @entry{@{$self->_headers}} = @fields;
    }

    my ($ref, $alt) = @entry{qw(ref_allele alt_allele)};

    #add the ref and alt alleles' positions in the allele array to the GT field
    my $tumor_gt = '0/1';   # all are heterozygous b/c mutect doesn't care
    my $normal_gt = '0/0';  # muTect calculates based on the assumption that the normal is reference

    # allele depth
    my $normal_dp =  $entry{n_alt_count} + $entry{n_ref_count};
    my $tumor_dp =  $entry{t_alt_count} + $entry{t_ref_count};
    my $normal_ad = join(",",$entry{n_ref_count}, $entry{n_alt_count});
    my $tumor_ad = join(",", $entry{t_ref_count}, $entry{t_alt_count});
    # fraction of reads supporting alt
    my $normal_fa = $normal_dp ? sprintf("%0.02f",$entry{n_alt_count} / $normal_dp) : 0;
    my $tumor_fa = $tumor_dp ? sprintf("%0.02f",$entry{t_alt_count} / $tumor_dp) : 0;
    # somatic status
    my $normal_ss = ".";
    my $tumor_ss  = 2;

    # Placeholder for later adjustment
    my $dbsnp_id = ".";
    my $qual = "."; # Can also be $tumor_vaq
    my $filter = $entry{judgement} eq 'REJECT' ? 'REJECT' : 'PASS';
    my $format = "GT:DP:AD:FA:SS:TLOD";
    my $info = ".";
    my $tumor_sample_string = join (":", ($tumor_gt, $tumor_dp, $tumor_ad, $tumor_fa, $tumor_ss, $entry{t_lod_fstar}));
    my $normal_sample_string = join (":", ($normal_gt, $normal_dp, $normal_ad, $normal_fa, $normal_ss, "."));

    my $vcf_line = join("\t", $entry{contig}, $entry{position}, $dbsnp_id, $ref, $alt, $qual, $filter, $info, $format, $normal_sample_string, $tumor_sample_string);

    return $vcf_line;
}


sub get_format_meta {
    my $self = shift;

    # Get all of the base FORMAT lines
    my @tags = $self->SUPER::get_format_meta;

    my $tlod = {MetaType => "FORMAT", ID => "TLOD",    Number => ".", Type => "Float", Description => "Log of (likelihood tumor event is real / likelihood event is sequencing error)"};
    
    return (@tags, $tlod,);
}

sub get_filter_meta {
    my $self = shift;

    ##FILTER=<ID=REJECT,Description="Rejected as a confident somatic mutation">
    my $filter = {MetaType => "FILTER", ID => "REJECT", Description => "Rejected as a confident somatic mutation by MuTect"};
    return ($filter);

}


