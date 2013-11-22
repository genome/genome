package Genome::Model::Tools::Validation::LongIndelsPartTwo;

use warnings;
use strict;
use Genome;
use IO::File;
use File::Basename;

class Genome::Model::Tools::Validation::LongIndelsPartTwo {
    is => 'Command',
    has_input => [
        output_dir => {
            is => 'String',
            doc => 'directory for output files from LongIndelsPartOne.pm; must contain a file "contigs.fa"',
        },
        tumor_val_model_copy_id => {
            is => 'Number',
            doc => 'model ID for copy of tumor validation model made in LongIndelsPartOne.pm',
        },
        normal_val_model_copy_id => {
            is => 'Number',
            doc => 'model ID for copy of normal validation model made in LongIndelsPartOne.pm',
        },
    ],
    has_optional_input => [
        contigs_file => {
            is => 'String',
            doc => 'contigs.fa file output from LongIndelsPartOne.pm',
            is_optional => 1,
        },
        tier_file_location => {
            is => 'String',
            doc => 'tiering file location to be used by gmt fast-tier fast-tier.',
            is_optional => 1,
            example_values => ['/gscmnt/ams1100/info/model_data/2771411739/build102550711/annotation_data/tiering_bed_files_v3'],
        },
    ],
    doc => 'Final steps in the validation of 3bp indels.',
};

sub help_detail {
    return <<"EOS"
    This tool performs the last steps (#'s 6-8) of the 3bp indel validation process outlined on this wiki page: $ENV{GENOME_SYS_SERVICES_WIKI_URL}Medical_Genomics/Nimblegen_Solid_Phase_Capture_Validation/Analysis#.3E3bp_Indels. It also then creates an annotate-able file using an adaptor 'gmt annotate adaptor indel-contig'. It also creates a bed file from this adapted list, and then uses fast-tier to tier the final calls. Lastly, the tool prints some details for possible manual review tickets as well, so be sure to SAVE THE STDOUT.
EOS
}

sub execute {

    my $self = shift;

    #parse input params
    my $output_dir = $self->output_dir;
    my $tumor_val_model_copy_id = $self->tumor_val_model_copy_id;
    my $normal_val_model_copy_id = $self->normal_val_model_copy_id;
    my $tier_file_location = $self->tier_file_location;

    #look for contigs.fa file
    my $contigs_file = $self->contigs_file;
    unless (defined $contigs_file) {
        $contigs_file = $output_dir . "/contigs.fa";
    }
    unless (-s $contigs_file) {
        $self->error_message("Contigs file $contigs_file does not exist with size > 0.");
        return;
    }

    #get builds and bam files
    my $tumor_val_model_copy = Genome::Model->get($tumor_val_model_copy_id) or die "Could not find tumor model with id $tumor_val_model_copy_id.\n";
    my $normal_val_model_copy = Genome::Model->get($normal_val_model_copy_id) or die "Could not find normal model with id $normal_val_model_copy_id.\n";
    my $tumor_build = $tumor_val_model_copy->last_succeeded_build or die "Could not find last succeeded build from tumor model $tumor_val_model_copy_id.\n";
    my $normal_build = $normal_val_model_copy->last_succeeded_build or die "Could not find last succeeded build from normal model $normal_val_model_copy_id.\n";
    my $tumor_bam = $tumor_build->whole_rmdup_bam_file or die "Cannot find tumor .bam.\n";
    my $normal_bam = $normal_build->whole_rmdup_bam_file or die "Cannot find normal .bam.\n";

    my $rv = Genome::Model::Tools::Validation::LongIndelsParseRemapped->execute(
        contigs_file => $contigs_file,
        tumor_bam => $tumor_bam,
        normal_bam => $normal_bam,
        output_dir => $output_dir,
        tier_file_location => $tier_file_location,
    );

    return $rv;
}

1;
