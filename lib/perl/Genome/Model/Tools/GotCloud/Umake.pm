package Genome::Model::Tools::GotCloud::Umake;

use strict;
use warnings;

use Genome;
use Memoize;

class Genome::Model::Tools::GotCloud::Umake {
    is => ['Command::V2', 'Genome::Model::Tools::GotCloud'],
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
    my $thunder_config_file = $self->generate_config_file($self->build_thunder_config_hash, File::Spec->join($self->output_directory,"thunder.conf"));
    my $beagle_config_file = $self->generate_config_file($self->build_beagle_config_hash, File::Spec->join($self->output_directory,"beagle.conf"));
    my $vcf_extract_config_file = $self->generate_config_file($self->build_vcf_extract_config_hash, File::Spec->join($self->output_directory,"vcf_extract.conf"));
    my $gotcloud_path = $self->gotcloud_path;
    Genome::Sys->shellcmd(cmd =>"$gotcloud_path vc --conf $vcf_extract_config_file --numjobs 1");
    Genome::Sys->shellcmd(cmd =>"$gotcloud_path vc --conf $beagle_config_file --numjobs 1");
    Genome::Sys->shellcmd(cmd =>"$gotcloud_path vc --conf $thunder_config_file --numjobs 1");
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

sub build_thunder_config_hash {
    my $self= shift;
    my $config;
    $config->{"RUN_THUNDER"}="TRUE";
    return $self->build_config_hash($config);
}

sub build_beagle_config_hash {
    my $self = shift;
    my $config;
    $config->{"RUN_BEAGLE"}="TRUE";
    $config->{"RUN_SUBSET"}="TRUE";
    return $self->build_config_hash($config);
}

sub build_vcf_extract_config_hash {
    my $self = shift;
    my $config;
    $config->{"RUN_PILEUP"}="TRUE";
    $config->{"RUN_SPLIT"}="TRUE";
    $config->{"RUN_EXTRACT"}="TRUE";
    $config->{"VCF_EXTRACT"}=$self->vcf_file;
    return $self->build_config_hash($config);
}

sub single_track_roi {
    my $self = shift;
    my $uniform_target_bed = $self->output_directory."/roi.bed";
    Genome::Sys->copy_file($self->roi_list->get_target_track_only($self->target_track),$uniform_target_bed);
    return $uniform_target_bed;
}

Memoize::memoize('single_track_roi');

sub build_config_hash {
    my $self = shift;
    my $config = shift;
    $config->{"AS"}=$self->reference_build->name;
    $config->{"REF_DIR"}=$self->reference_build->data_directory;
    $config->{"REF"}=$self->reference_build->full_consensus_path("fa");
    $config->{"DBSNP_VCF"}="/gscuser/kmeltzst/gscmnt/reference_files/gotcloud.ref/dbsnp_135.b37.vcf.gz";
    $config->{"HM3_VCF"}="/gscuser/kmeltzst/gscmnt/reference_files/gotcloud.ref/hapmap_3.3.b37.sites.vcf.gz";
    $config->{"OMNI_VCF"}="/gscuser/kmeltzst/gscmnt/reference_files/gotcloud.ref/1000G_omni2.5.b37.sites.PASS.vcf.gz";
    $config->{"INDEL_PREFIX"}="/gscuser/kmeltzst/gscmnt/reference_files/gotcloud.ref/1kg.pilot_release.merged.indels.sites.hg19";
    #FIXME:don't make me hardcoded
    $config->{"OUT_DIR"}=$self->output_directory;
    $config->{"BAM_INDEX"}=$self->generate_bam_index;
    $config->{"PED_INDEX"}=$self->pedigree_file;
    $config->{"WRITE_TARGET_LOCI"}="TRUE";
    $config->{"UNIFORM_TARGET_BED"}=$self->single_track_roi;
    $config->{"OFFSET_OFF_TARGET"}=$self->wingspan;
    $config->{"TARGET_DIR"}="target";
    $config->{"CHRS"}=join " ",$self->chromosomes;
    return $config;
}

sub generate_bam_index{
    my $self = shift;
    my $outdir = $self->output_directory;
    my $outfile = $outdir."/bam_index.txt";
    my $outfile_fh = Genome::Sys->open_file_for_writing($outfile);
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
    return $outfile;
}

Memoize::memoize('generate_bam_index');

sub generate_config_file{
    my $self = shift;
    my $config = shift;
    my $config_file = shift;
    my $outdir = $self->output_directory;
    my $config_fh = Genome::Sys->open_file_for_writing($config_file);
    for my $key (keys %$config){
        $config_fh->print("$key=".$config->{$key}."\n");
    }

    return $config_file;
}
