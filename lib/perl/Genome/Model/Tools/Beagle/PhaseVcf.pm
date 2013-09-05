package Genome::Model::Tools::Beagle::PhaseVcf;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Beagle::PhaseVcf {
    is => 'Genome::Model::Tools::Beagle',
    has_input => [
    vcf_file => {
        is => 'Text',
        doc => 'unphased VCF input file',
    },
    output_file => {
        is => 'Text',
        doc => 'The output file prefix',
    },
    chromosome => {
        is => 'Text',
        doc => 'The chromosome to phase',
    },
    ],
    has_optional_input => [
    ped_file => {
        is => 'Text',
        doc => 'Pedigree file in PLINK ped format. See http://pngu.mgh.harvard.edu/~purcell/plink/ for more information',
    },
    impute => {
        is => 'Boolean',
        default => 0,
        doc => 'Specifies wheather variants that are present in the reference panel but absent in your data will be imputed',
    },
    ref_pop => {
        is => 'Text',
        doc => 'Specifies the reference population to use for imputation (AFR = African, AMR = American, ASN = Asian, EUR = European)',
    },
    exclude_samples => {
        is => 'Text',
        doc => 'specifies a file containing non-reference samples to be excluded from analysis and output files',
    },
    exclude_markers => {
        is => 'Text',
        doc => 'Specifies a file containing markers to be excluded from analysis and output files',
    },
    ],
};

sub help_brief {
    "Phase a VCF file containing family or unrelated cohorts using Beagle v4 with or without a reference panel"
}

sub help_synopsis {
    my $self = shift;
    "gmt beagle phasevcf --vcf-file=a.vcf --output-file=phased --chromosome=1 [--ped-file=a.ped] [--impute=0] [--ref-pop=AFR/AMR/ASN/EUR] [--exclude-samples=sample.txt] [--exclude-markers=positions.txt]"
}

sub flags {
    my $self = shift;
    my @flags;
    push(@flags, "ped=" . $self->ped_file) if defined $self->ped_file;
    push(@flags, "excludesamples=" . $self->exclude_samples) if defined ($self->exclude_samples);
    push (@flags, "excludemarkers=" . $self->exclude_markers) if defined ($self->exclude_markers);

    return @flags;
}

sub execute {
    my $self = shift;
    my $input = $self->vcf_file;
    my $output = $self->output_file;
    my $chrom = $self->chromosome;
    my($ref_file,$ref_pop);
    my $flags = join(" ", $self->flags);
    my $path = $self->path_for_version($self->version);
    my $cmd = "java -Xmx2G -jar $path gt=$input out=$output chrom=$chrom $flags";

    if($self->impute){
        $cmd .= " impute=true";
        $ref_pop = $self->ref_pop;
        $ref_file = "/gscmnt/ams1161/info/model_data/kmeltzst/reference_files/beagle_files/ref_vcf/".$ref_pop.".chr".$chrom.".1kg.ref.phase1_release_v3.20101123.vcf.gz";
        $cmd .= " ref=$ref_file";
    }
    else {
        $cmd .=" impute=false";
    }

    Genome::Sys->shellcmd (cmd=>$cmd);
    return 1;
}
1;
