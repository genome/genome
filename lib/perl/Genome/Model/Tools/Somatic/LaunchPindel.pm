package Genome::Model::Tools::Somatic::LaunchPindel;

use warnings;
use strict;

use Genome;
use Carp;
use IO::File;
use Genome::Info::IUB;

class Genome::Model::Tools::Somatic::LaunchPindel{
    is => 'Command',
    has => [
        tumor_build => {
            is  => 'String',
            is_input  => 1,
            is_optional=>1,
            doc => 'tumor build id to get the bam from',
        },
        normal_build => {
            is  => 'String',
            is_input  => 1,
            is_optional => 1,
            doc => 'normal build id to get the bam from',
        },
        tumor_bam => {
            is  => 'String',
            is_input  => 1,
            is_optional => 1,
            doc => 'use this to directly specify tumor bam if you do not have a build',
        },
        normal_bam => {
            is  => 'String',
            is_input  => 1,
            is_optional => 1,
            doc => 'use this to directly specify normal bam if you do not have a build',
        },
        reference_sequence_build => {
            is => 'Text',
            is_optional => 0,
            example_values => ['101947881'],
            doc => 'reference sequence build id, e.g. 101947881 for NCBI human build36',
        },
        output_dir => {
            is => 'Text',
            is_input => 1,
            is_output => 1,
            doc => "place to put results"
        },
    ],
};

sub help_brief {
    "Separate LOH calls from non-LOH calls",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
 gmt somatic launch-pindel --normal-bam=/gscmnt/sata921/info/medseq/indel_validation_bams/true_positive_normal_validation.bam  --tumor-bam=/gscmnt/sata921/info/medseq/indel_validation_bams/true_positive_normal_validation.bam --output-dir=/gscmnt/sata921/info/medseq/testing_launch_pindel
EOS
}

sub help_detail {                           
    return <<EOS 
This runs pindel version .4 with the Fisher's Exact test germline filter at the conclusion.  It is implemented inside the DetectVariants API so this wrapper makes it easier to use.  This new version of pindel allows more reads into consideration which should generate more events and also filter out more germline events.  The new filter does a FET between the distribution of reads mapping with an insertion/deletion vs normally mapped reads at the site in tumor and normal.  If the p-value is .15 or below that they are different, we keep the site under consideration.  The idea here is that GATK will resolve similar indel distribution sites more accurately than pindel has been able to do.  In addition, pindel recovers some reads that are unmapped previously and we attempt to pass any of these  sites through the filter.  That 'override' condition is: Tumor reads at site < 10 and pindel_reads > indel mapping reads in tumor
EOS
}

sub execute {
    my $self = shift;
    my $normal_bam = $self->normal_bam;
    my $tumor_bam = $self->tumor_bam;
    if($self->normal_build && $self->tumor_build) { 
        my $normal_build = Genome::Model::Build::ReferenceAlignment->get($self->normal_build);
        my $tumor_build = Genome::Model::Build::ReferenceAlignment->get($self->tumor_build);
         $normal_bam = $normal_build->whole_rmdup_bam_file;
         $tumor_bam = $tumor_build->whole_rmdup_bam_file;
    }
    Genome::Sys->validate_file_for_reading($normal_bam); 
    Genome::Sys->validate_file_for_reading($tumor_bam);

    my $output = $self->output_dir;
    Genome::Sys->create_directory($output);  #this seems to be a no op if it exists
    
        my $email_address = $ENV{'LOGNAME'} . "\@genome.wustl.edu";
        $self->status_message("sending a completion mail to: $email_address");
        my $include;
        if($INC[0] !~ m/noarch/) {
            my $newlib  = $INC[0];
            $include = "-I $newlib";
            $self->status_message("using $newlib include on bsub");
        }
        my $reference_build_id=$self->reference_sequence_build;
        print `bsub -u $email_address -N -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKFLOW} "perl $include -S gmt detect-variants2 dispatcher --aligned-reads-input $tumor_bam    --control-aligned-reads-input $normal_bam --reference-build-id $reference_build_id --output-directory $output --indel-detection-strategy 'pindel 0.4 filtered by pindel-somatic-calls v1 then pindel-read-support v1'"`;
        return 1;
    }
1;
