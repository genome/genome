package Genome::ModelGroup::Command::FamilyBasedSubmission;

use strict;
use warnings;
use File::Basename;
use Genome;
use Workflow;
use Workflow::Simple;

use List::MoreUtils qw/ uniq /;
class Genome::ModelGroup::Command::FamilyBasedSubmission {
    is => 'Genome::Command::Base',
    has => [
    model_group => {
        is => 'Genome::ModelGroup',
        is_optional=>0,
        doc => 'this is the model group you wish to QC',
    },
    sample_mapping_file => {
        is=> 'String',
        doc => 'this command will generate this list of bam file names & samples',
    }, 
    bam_list_file => {
        is=> 'String',
        doc => 'this command will generate this list of bam file names, newline separated suitable for gxfer to submit bams',
    },
    md5_list_file => {
        is=> 'String',
        doc => 'this will generate a list of bam md5 file names, newline separated',
    },
    vcf_list_file =>{
        is=> 'String',
        doc => 'this will generate a list of vcf file names, newline separated',
    },
    vcf_dir=> {
        is=>'String',
        doc=>"directory to copy vcfs into so that we can change their names from the default pipeline snvs.vcf.gz",

    },
    ],
    has_optional => {
    },
    doc => 'run fastIBD' 
};

sub help_synopsis {
    return <<EOS
genome model-group relationship-qc --model-group=1745 --output-dir=/foo/bar/

EOS
}


sub execute {
    my $self=shift;
    my $model_group  = $self->model_group;
    my @models = $model_group->models;
    my $samp_map = IO::File->new(">".$self->sample_mapping_file);
    unless ($samp_map) {
        $self->error_message("Failed to open sample mapping file for writing ". $self->sample_mapping_file);
        return;
    }
    my $bam_list = IO::File->new(">".$self->bam_list_file);
    unless ($bam_list) {
        $self->error_message("Failed to open bam list file for writing ". $self->bam_list_file);
        return;
    }

    my $md5_list = IO::File->new(">".$self->md5_list_file);
    unless ($md5_list) {
        $self->error_message("Failed to open md5 list file for writing ". $self->md5_list_file);
        return;
    }
    my $vcf_list = IO::File->new(">".$self->vcf_list_file);
    unless ($vcf_list) {
        $self->error_message("Failed to open vcf list file for writing ". $self->vcf_list_file);
        return;
    }

    my @builds = map {$_->last_succeeded_build} @models;
    for my $build (@builds) {
        my @bams;
        my @results = $build->results;
        for my $result (@results) {
            if ($result->class eq "Genome::InstrumentData::AlignmentResult::Merged") {
                push @bams, $result->bam_file;
            }
        }
        my @base_bam_name = map {basename($_)} @bams;
        my $roi_name = $build->model->roi_list ? $build->model->roi_list->name : 'N/A';
        my $refbuild_name = $build->model->reference_sequence_build->name ? $build->model->reference_sequence_build->name : 'N/A';
        my @md5_name = map {$_ . ".md5"}@bams;
        my $subject_name = $build->model->subject_name;
        ###cleft lip hack  FIXME when this gets generalized
        if($subject_name =~m/model group/) {
            $subject_name =~ s/population group for model group //;
            $subject_name =~ s/ /_/g;
        }
        
        print $samp_map join ("\t", $subject_name, $refbuild_name, $roi_name, @base_bam_name);
        print $samp_map "\n";
        my $num_bams = scalar(@bams);
        $bam_list->print(join("\n", @bams) . "\n");
        $md5_list->print(join("\n", @md5_name) . "\n");
        my $new_vcf_file =  $self->vcf_dir  . "/$subject_name.vcf.gz";
        my $vcf_file = $build->data_directory . "/variants/snvs.vcf.gz";
        `cp $vcf_file $new_vcf_file`;
        $vcf_list->print("$new_vcf_file\n");
    }
    $samp_map->close;
    $bam_list->close;
    $md5_list->close;
    $vcf_list->close;

}



1;
