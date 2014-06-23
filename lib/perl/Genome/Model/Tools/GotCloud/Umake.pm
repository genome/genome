package Genome::Model::Tools::GotCloud::Umake;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::GotCloud::Umake {
    is => 'Genome::Command::Base',
    has => [
        model_group => {
            is => 'Genome::ModelGroup',
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
        },
        output_directory => {
            is => 'path',
        },
        pedigree_file => {
            is => 'path',
            is_optional => 1,
        },
        vcf_file => {
            is => 'path',
        },
        roi_list => {
            is => 'Genome::FeatureList',
            is_optional => 1,
        },
        target_track => {
            is => 'text',
            is_optional => 1,
        },
        wingspan => {
            is => 'integer',
            is_optional => 1,
        },
        chromosomes => {
            is => 'text',
            is_many => 1,
            valid_values => [qw (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)]
        },
    ],
};

sub execute {
    my $self = shift;
    $self->validate;
    my $config_file = $self->generate_config_file;
    my $gotcloud_path = $self->gotcloud_path;
    Genome::Sys->shellcmd(cmd =>"$gotcloud_path snpcall --conf $config_file");
    Genome::Sys->shellcmd(cmd =>"$gotcloud_path ldrefine --conf $config_file");
    return 1;
}

sub validate {
    my $self = shift;
    Genome::Sys->create_directory($self->output_directory);
    if(defined $self->pedigree_file){
        Genome::Sys->validate_file_for_reading($self->pedigree_file);
    }
    Genome::Sys->validate_file_for_reading($self->vcf_file);
    if(defined $self->wingspan && !defined $self->roi_list){
        die $self->error_message("you must define an roi_list if you define a wingspan!");
    }
    if(defined $self->roi_list && !defined $self->wingspan){
        die $self->error_message("you must define an roi_list if you define a wingspan!");
    }
    if(defined $self->target_track && !defined $self->roi_list){
        die $self->error_message("you must define an roi_list if you define a target_track!");
    }
    if(grep {$_ eq "X"} $self->chromosomes && !defined $self->pedigree_file){
        die $self->error_message("you must give a pedigree file to run chromosome X!");
    }

    return 1;
}

sub build_config_hash {
    my $self = shift;
    my %config;
    $config{"AS"}=$self->reference_build->name;
    $config{"REF_DIR"}=$self->reference_build->data_directory;
    $config{"REF"}=$self->reference_build->full_consensus_path("fa");
    $config{"DBSNP_VCF"}="/gscuser/kmeltzst/gscmnt/reference_files/dbsnp_135.b37.vcf.gz";
    $config{"HM3_VCF"}="/gscuser/kmeltzst/gscmnt/reference_files/hapmap_3.3.b37.sites.vcf.gz";
    $config{"OMNI_VCF"}="/gscuser/kmeltzst/gscmnt/reference_files/1000G_omni2.5.b37.sites.PASS.vcf.gz";
    $config{"INDEL_PREFIX"}="/gscuser/kmeltzst/gscmnt/reference_files/1kg.pilot_release.merged.indels.sites.hg19";
    #FIXME:don't make me hardcoded
    $config{"OUT_DIR"}=$self->output_directory;
    $config{"BAM_INDEX"}=$self->generate_bam_index;
    $config{"PED_INDEX"}=$self->pedigree_file;
    $config{"RUN_INDEX"}="TRUE";
    $config{"RUN_PILEUP"}="TRUE";
    $config{"RUN_SPLIT"}="TRUE";
    $config{"RUN_BEAGLE"}="TRUE";
    $config{"RUN_SUBSET"}="TRUE";
    $config{"RUN_THUNDER"}="TRUE";
    $config{"RUN_EXTRACT"}="TRUE";
    $config{"VCF_EXTRACT"}=$self->vcf_file;
    $config{"WRITE_TARGET_LOCI"}="TRUE";
    $config{"UNIFORM_TARGET_BED"}=$self->roi_list->get_target_track_only($self->target_track);
    #FIXME: add subroutine get_target_track_only to Genome::FeatureList and make sure it dies if target track doesn't exist
    $config{"OFFSET_OFF_TARGET"}=$self->wingspan;
    $config{"TARGET_DIR"}="target";

    return %config;
}

sub generate_bam_index{
    my $self = shift;
    my($outfile_fh,$outfile_path) = Genome::Sys->create_temp_file;
    my @models = $self->model_group->models;
    for my $model (@models){
        my $model_name = $model->name;
        my $sample_name = $model->subject->name;
        my $data_dir = $model->last_succeeded_build->data_directory;
        my @bam_files = glob($data_dir."/alignments/tumor/*.bam");
        my $pop;
        if($sample_name =~m/H_OS/) {
            $pop = "METSIM";
        }
        if($sample_name =~m/H_OY/){
            $pop = "FINRISK";
        }
        #FIXME:don't make this hardcoded--make the pop be a class definition
        for my $file (@bam_files){
            $outfile_fh->print("$sample_name\t$pop\t$file\n");
        }
    }
    return $outfile_path;
}

sub generate_config_file{
    my $self = shift;
    my %config = $self->build_config_hash;
    my($config_fh,$config_path) = Genome::Sys->create_temp_file;
    for my $key (keys %config){
        $config_fh->print("$key=".$config{$key}."\n");
    }

    return $config_path;
}
