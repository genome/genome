package Genome::Model::Tools::FastTier::TierVcf;

use strict;
use warnings;

use Genome;     
use File::Basename;
use File::Copy;

class Genome::Model::Tools::FastTier::TierVcf {
    is => 'Command',
    has => [
        vcf_file => {
            type => 'Text',
            doc => "VCF formatted file, gzipped, containing variants to be tiered",
            is_input => 1,
        },
        tier_file_location => {
            type => 'Text',
            is_input => 1,
            doc => 'Use this to point to a directory containing tiers.bed.gz',
        },
        vcf_output_file => {
            is => 'Text',
            doc => 'Use this only when running on VCF files, it will be the location of your output',
        },
    ],
};

sub sub_command_sort_position { 15 }

sub help_brief {
    "This tool adds TIER info tags to VCF files."
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt fast-tier tier-vcf
EOS
}

sub execute {
    my $self=shift;

    unless($self->vcf_file){
        die $self->error_message("You must specify either variant-bed-file or variant-vcf-file.");
    }

    my ($input, $output) = (0,0);

    my $vcf = $self->vcf_file;
    my $tiering_file = $self->tier_file_location . "/tiers.bed.gz";
    unless(-s $tiering_file){
        die $self->error_message("Could not locate tiering file at: ".$tiering_file);
    }
    my $gz = ($vcf =~ m/gz$/);

    unless($gz){
        die $self->error_message("variant-vcf-file must be gzipped, and named *.gz");
    }
    my $output_file = $self->vcf_output_file;
    if(-e $output_file){
        die $self->error_message("Found a file in the location of the proposed output: ".$output_file);
    }
    my $info_record = "key=INFO,ID=TIER,Number=1,Type=Integer,Description='Location of variant by tier'";
    my $tier_cmd = "zcat ".$vcf." | vcf-annotate -a ".$tiering_file." -d ".$info_record." -c CHROM,FROM,TO,INFO/TIER | bgzip -c > ".$output_file;

    unless(Genome::Sys->shellcmd( cmd => $tier_cmd)){
        die $self->error_message("vcf-annotate command did not return correctly.");
    }

    return 1;
}
