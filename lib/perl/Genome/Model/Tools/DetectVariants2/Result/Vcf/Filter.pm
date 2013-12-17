package Genome::Model::Tools::DetectVariants2::Result::Vcf::Filter;

use strict;
use warnings;

use Genome;
use File::Copy;
use Sys::Hostname;

class Genome::Model::Tools::DetectVariants2::Result::Vcf::Filter {
    is  => ['Genome::Model::Tools::DetectVariants2::Result::Vcf'],
    has_param => [
        incoming_vcf_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Vcf',
            doc => 'This is the vcf-result of the detector or filter being run on',
        },
        filter_name => {
            is => 'Text',
            doc => 'Name of the filter run to generate this vcf',
        },
        filter_version => {
            is => 'Text',
            doc => 'Version of the filter',
            is_optional => 1,
        },
        filter_params => {
            is => 'Text',
            doc => 'Params for the filter',
            is_optional => 1,
        },
        previous_filter_strategy => {
            is => 'Text',
            doc => 'Name version and params for the previous filter, if there is one.',
            is_optional => 1,
        },
    ],
    has_optional => [
        # This is not saved in the db but is used to create the vcf file
        filter_description => {
            is => 'Text',
            doc => 'Description of the filter applied to this data',
            default => 'Filter variants',
            is_transient => 1,
        },
    ],
};

sub _generate_vcf {
    my $self = shift;
    my $retval=1;
    my $detector = $self->input->detector_name;
    my $path = $self->input_directory;

    for my $variant_type ("snvs","indels"){
        my $variant_file = $path."/".$variant_type.".hq";
        unless( -e $variant_file ){
            next;
        }
        if(Genome::Model::Tools::DetectVariants2::Result::Vcf->conversion_class_name($detector,$variant_type)){
            $self->status_message("Generating Vcf");
            $retval &&= $self->_run_vcf_converter($variant_type);
        } else {
            $self->status_message("Skipping ".$variant_type." vcf generation for ".$detector);
        }
    }

    return $retval;
}

sub _run_vcf_converter {
    my $self = shift;
    my $type = shift;
    my $dirname = $self->output_dir;
    unless($dirname){
        die $self->error_message("Could not get dirname!");
    }
    my $output_file = $dirname . '/'.$type.'.vcf.gz';

    my $incoming_vcf = $self->incoming_vcf_result->output_dir."/".$type.".vcf.gz";
    unless(-s $incoming_vcf){
        $self->status_message("Skipping VCF generation, no vcf in the previous result: $incoming_vcf");
        return 0;
    }

    my $hq_filter_file = $self->input_directory."/".$type.".hq.bed";

    my $filter_name = $self->filter_name;
    my @names = split /\:\:/,$filter_name;
    $filter_name = $names[-1];

    my $filter_description = $self->filter_description;

    my %params = (
        output_file => $output_file,
        vcf_file => $incoming_vcf,
        filter_file => $hq_filter_file,
        filter_keep => 1,
        filter_name => $filter_name,
        filter_description => $filter_description,
        bed_input => 1,
    );

    #Make sure the input vcf does not have multiple indels in ALT column
    my $vcf_filter_cmd = Genome::Model::Tools::Vcf::VcfFilter->create(%params);
    unless($vcf_filter_cmd->execute){
        die $self->error_message("Could not complete call to gmt vcf vcf-filter!");
    }
    return 1;
}

sub _needs_symlinks_followed_when_syncing { 
    return 0;
}

sub _working_dir_prefix {
    return "detector_vcf_results";
}

sub resolve_allocation_disk_group_name { 
    $ENV{GENOME_DISK_GROUP_MODELS};
}

sub allocation_subdir_prefix {
    return "detector_vcf_results";
}

sub _combine_variants {
    die "overload this function to do work";
}

sub estimated_kb_usage {
    return 10_000_000;
}

sub _staging_disk_usage {
    return 10_000_000;
}

sub _add_as_user_of_inputs {
    my $self = shift;

    my $prev_vcf_result = $self->incoming_vcf_result;
    $prev_vcf_result->add_user(user => $self, label => 'uses');

    my $input = $self->input;

    return (
        $input->add_user(user => $self, label => 'uses')
    );
}

1;
