package Genome::Model::Tools::DetectVariants2::Result::Vcf::Combine;

use strict;
use warnings;

use Genome;
use File::Copy;
use Sys::Hostname;

class Genome::Model::Tools::DetectVariants2::Result::Vcf::Combine {
    is  => ['Genome::Model::Tools::DetectVariants2::Result::Vcf'],
    has_param => [
        input_a_id => {
            is => 'Text',
            doc => 'ID of the first incoming software result',
        },
        input_b_id => {
            is => 'Text',
            doc => 'ID of the second incoming software result',
        },
        variant_type => {
            is => 'Text',
            valid_values => ['snvs','indels'],
            doc => 'type of variants being combined',
        },
        joinx_version => {
            is => 'Text',
            doc => 'Version of joinx to use for the combination',
        },
        #This isn't set on combine results--they use the samples of their inputs
        aligned_reads_sample => {
            is => 'Text',
            is_optional => 1,
        },
        incoming_vcf_result_a_id => {
            is => 'Number',
        },
        incoming_vcf_result_b_id => {
            is => 'Number',
        },
    ],
    has => [
        incoming_vcf_result_a => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Vcf',
            id_by => 'incoming_vcf_result_a_id',
            doc => 'This is the vcf-result of the first detector or filter being run on',
        },
        incoming_vcf_result_b => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Vcf',
            id_by => 'incoming_vcf_result_b_id',
            doc => 'This is the vcf-result of the second detector or filter being run on',
        },
    ],
};

sub _generate_vcf {
    my $self = shift;
    my $retval=1;
    my $path = $self->input_directory;

    for my $variant_type ("snvs", "indels"){
        $self->_run_vcf_converter($variant_type);
    }

    return $retval;
}

sub _run_vcf_converter {
    my $self = shift;
    my $type = shift;

    my $input = $self->input;
    my $dirname = $self->output_dir;
    unless($dirname){
        die $self->error_message("Could not get dirname!");
    }
    my $output_file = $dirname . '/'.$type.'.vcf.gz';

    my @input_files;
    my @labeled_input_files;
    for my $vcf_result ($self->incoming_vcf_result_a, $self->incoming_vcf_result_b) {
        my $input_vcf = $vcf_result->output_dir."/".$type.".vcf.gz";
        unless(-s $input_vcf){
            $self->status_message("Skipping VCF generation for type $type, no vcf in the previous result: $input_vcf");
            return 0;
        }

        # Using a label (-D in joinx) preserves that input column.
        # We want to preserve the original input columns (via joinx) for each detector result
        # but not combination operations (because that would be redundant as they would have already been preserved)
        if ($vcf_result->isa("Genome::Model::Tools::DetectVariants2::Result::Vcf::Combine")) {
            push @input_files, $input_vcf;
        } else {
            my $detector_class = $vcf_result->input->detector_name;
            my @class_path = split "::", $detector_class;
            my $detector_name = $class_path[-1];
            unless ($detector_name) {
                die $self->error_message("Could not get a detector name from detector class $detector_class from software result id" . $vcf_result->id);
            }
            my $tag = "-[$detector_name]";
            $input_vcf .= "=$tag";
            push @labeled_input_files, $input_vcf;
        }
    }

    my %params = ( 
        input_files => \@input_files,
        labeled_input_files => \@labeled_input_files,
        output_file => $output_file,
        merge_samples => 1,
        clear_filters => 1,
        use_bgzip => 1,
        use_version => $self->joinx_version,
    );

    # If we are doing an intersection, set the ratio filter to mark things as filtered where they do not agree
    if ($input->class =~ m/Intersect/) {
        $params{ratio_filter} = "1.0,IntersectionFailure,Variant callers do not agree on this position";
        $params{sample_priority} = "filtered";
    } else {
        $params{sample_priority} = "unfiltered";
    }

    my $merge_cmd = Genome::Model::Tools::Joinx::VcfMerge->create(%params);

    unless($merge_cmd->execute){
        die $self->error_message("Could not complete call to gmt vcf vcf-filter!");
    }
    return 1;
}

sub _validate_input {
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

    for my $prev_vcf_result ($self->incoming_vcf_result_a,$self->incoming_vcf_result_b){
        $prev_vcf_result->add_user(user => $self, label => 'uses');
    }   

    my $input = $self->input;

    return (
        $input->add_user(user => $self, label => 'uses')
    );
}

1;
