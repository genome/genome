package Genome::Model::Tools::Vcf;

use strict;
use warnings;

use Genome;

#This is the variable to change if you wish to change the version of all vcf files being created
# 4: because of a bug in VcfFilter.pm. When filtered twice, FT fields were being wiped out rather than propagated from "FILTER".
# 5: Varscan now rounds tumor vaq to the nearest integer value so it agrees with the header type field
# 6: some TCGA-compliance format, fix VcfFilter bug to mis-treat some samtools mpileup indel
# 7: more TCGA-compliance format, add TCGA format output of snv and indel to streka tool, add fix to Varscan Somatic snv vcf
# 8: TCGA compliant vcf headers
# 9: When combining vcfs in DV2, keep the original per-detector sample columns. We will now have one column per sample and detector plus a per-sample consensus column.
#10: Change the description of FORMAT "FT" to be TCGA-compliant
#11: Joinx now creates new ##SAMPLE columns when using the -D option
#12: Samtools indel -> vcf conversion didn't handle ins/del or del/del calls properly
#13: gmt vcf snv varscan changes, adding BQ and AD values for the reference allele, add FT to FORMAT subfield and populate it in sample colums for strelka detector vcf output, as well as a few other corrections.
#14: We produce vcfs during indel combination operations now
#15: 14 had undefined sample names in the Vcf software results
#17: Indel Normalization
#18: Indel Normalization using joinx 1.7 (since 1.6 is bad)
#19: Sort normalized indels (numeric)
#20: Switch to natural sorting of normalized indels and change null GT from '.' to './.'
#21: Tabix index all vcf results

my $VCF_VERSION = "21";

class Genome::Model::Tools::Vcf {
    is => ['Command'],
    has => [
        vcf_version => {
            is => 'Text',
            default => $VCF_VERSION,
        },
    ],
};

sub get_vcf_version {
    return $VCF_VERSION;
}

sub help_brief {
    "Tools and scripts to create and manipulate VCF files."
}

sub help_detail {
    return <<EOS
Tools and scripts to create and manipulate VCF files.
EOS
}

1;
