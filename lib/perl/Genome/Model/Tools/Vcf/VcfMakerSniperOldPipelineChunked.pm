package Genome::Model::Tools::Vcf::VcfMakerSniperOldPipelineChunked;

use strict;
use warnings;
use Genome;
use File::stat;
#use Time::localtime;
use IO::File;
use File::Basename;
use Getopt::Long;
use FileHandle;
#use POSIX qw(log10);
use List::MoreUtils qw(firstidx);
use List::MoreUtils qw(uniq);
use Data::Dumper;


class Genome::Model::Tools::Vcf::VcfMakerSniperOldPipelineChunked {
    is => 'Command',
    has => [
    output_file => {
        is => 'Text',
        is_output => 1,
        doc => "List of mutations in Vcf format",
    },

    tumor_bam_file => {
        is => 'Text',
        doc => "Tumor sample bam file (for header)" ,
        is_optional => 0,
        is_input => 1},

    normal_bam_file => {
        is => 'Text',
        doc => "Normal sample bam file (for header)" ,
        is_optional => 0,
        is_input => 1},

    tumor_snp_file => {
        is => 'Text',
        doc => "all snps from the tumor" ,
        is_optional => 0,
        is_input => 1},

    normal_snp_file => {
        is => 'Text',
        doc => "all snps from the normal" ,
        is_optional => 0,
        is_input => 1},

    file_source => {
        is => 'Text',
        doc => "source of the bam files",
        is_optional => 1,
        is_input => 1,
        default =>"dbGap" },

    build_dir => {
        is => 'Text',
        doc => "Build directory",
        is_optional => 0,
        is_input => 1},

    dbsnp_file => {
        is => 'Text',
        doc => "dbsnp File - if specified, will label dbSNP sites" ,
        is_optional => 1,
        is_input => 1,
        default => ""},

    individual_id => {
        is => 'Text',
        doc => "Individual ID",
        is_optional => 0,
        is_input => 1},


    center => {
        is => 'Text',
        doc => "Genome center name (WUSTL, Broad, Baylor)" ,
        is_optional => 1,
        default => "WUSTL",
        is_input => 1},

    filterOut => {
        is => 'Text',
        doc => "file containing SNVs to label as filtered out in format: filterName:description:file,filterName2,description2,file2,..." ,
        is_optional => 1,
        default => "",
        is_input => 1},


    filterPass => {
        is => 'Text',
        doc => "file containing SNVs that pass filters (all SNVs not contained in the file will be marked filtered out) in format:  filterName:description:file,filterName2,description2:file2,... " ,
        is_optional => 1,
        default => "",
        is_input => 1},


    genome_build => {
        is => 'Text',
        doc => "Reference genome build" ,
        is_optional => 1,
        default => "36",
        is_input => 1},

    ],
};


sub help_brief {                            # keep this to just a few words <---
    "Generate Vcf File without memory issues by running sequentially"
}


sub help_synopsis {
    <<'HELP';
Generate a VCF File
HELP
}

sub help_detail {                  # this is what the user will see with the longer version of help. <---
    <<'HELP';
Parses the relevant files and creates a VCF containing all the SNVs. This includes those that fail filters (noted in the FILTER field). This wraps gmt vcf vcf-make-sniper-old-pipeline
HELP
}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
    my $self = shift;

    #first do chromosome 1 with the header
    unless(Genome::Model::Tools::Vcf::VcfMakerSniperOldPipeline->execute( build_dir => $self->build_dir, normal_snp_file => $self->normal_snp_file, tumor_snp_file => $self->tumor_snp_file, normal_bam_file => $self->normal_bam_file, tumor_bam_file => $self->tumor_bam_file, output_file => $self->output_file, individual_id => $self->individual_id, chrom => 1, dbsnp_file => $self->dbsnp_file )) {
        $self->error_message("Failed to run chromosome 1");
        return;
    }

    for my $chromosome (2..22,'X','Y','MT') {
        unless(Genome::Model::Tools::Vcf::VcfMakerSniperOldPipeline->execute( build_dir => $self->build_dir, normal_snp_file => $self->normal_snp_file, tumor_snp_file => $self->tumor_snp_file, normal_bam_file => $self->normal_bam_file, tumor_bam_file => $self->tumor_bam_file, output_file => $self->output_file, individual_id => $self->individual_id, chrom => $chromosome, dbsnp_file => $self->dbsnp_file, skip_header => 1)) {
            $self->error_message("Failed to run chromosome $chromosome");
            return;
        }
    }
    return 1;
}
1;
