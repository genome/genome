package Genome::Model::Tools::DetectVariants2::Result::Vcf::Detector;

use strict;
use warnings;

use Genome;
use File::Copy;
use Sys::Hostname;

class Genome::Model::Tools::DetectVariants2::Result::Vcf::Detector {
    is  => ['Genome::Model::Tools::DetectVariants2::Result::Vcf'],
};

sub _generate_vcf {
    my $self = shift;
    my $path = $self->input_directory;
    my $retval = 1;
    my %files;
  
    my $detector_class;
    if($self->input->can("detector_name")){
        $detector_class = $self->input->detector_name;
    } else {
        die $self->error_message("Could not call detector_name on input!");
    }
    my $convertable = 0;
    for my $variant_type ("snvs","indels"){
        my $variant_file = $path."/".$variant_type.".hq";
        unless(-e $variant_file){
            next;
        }

        my $vcf_module = $self->conversion_class_name($detector_class,$variant_type);
        if($vcf_module){
            $convertable++;
            $self->status_message("Generating Vcf");
            $self->status_message("executing $vcf_module on file $variant_file");
            $retval &&= $self->_run_vcf_converter($vcf_module, $variant_file,$variant_type);
            
        } else {
            $self->status_message("Couldn't find a working vcf converter for $detector_class $variant_type");
        }
    }
    unless($convertable>0){
        die $self->error_message("Was not able to convert any outputs to VCF");
    }

    return $retval;
}

sub _run_vcf_converter {
    my $self = shift;
    my $converter = shift;
    my $input_file = shift;
    my $type = shift;
    my $detector_result = $self->input;
    my $dirname = $self->output_dir;
    unless($dirname){
        die $self->error_message("Could not get dirname!");
    }
    my $output_file = $dirname . '/'.$type.'.vcf.gz';
    my $aligned_reads_sample = $self->aligned_reads_sample;
    my $control_aligned_reads_sample = $self->control_aligned_reads_sample;
    my $reference_sequence_build = Genome::Model::Build->get($detector_result->reference_build_id);
    my %params = (
        input_file => $input_file,
        output_file => $output_file, 
        aligned_reads_sample => $aligned_reads_sample,
        sequencing_center => 'WUSTL',
        reference_sequence_build => $reference_sequence_build,
    );
    $params{control_aligned_reads_sample} = $control_aligned_reads_sample if $control_aligned_reads_sample;

    my $command = $converter->create(%params); 
    unless($command->execute) {
        $self->error_message('Failed to convert ' . $input_file . ' to the standard format.');
        return;
    }

    return 1;
}

sub _gather_params_for_get_or_create {
    my $class = shift;

    my $bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);

    my %params = $bx->params_list;
    my %is_input;
    my %is_param;
    my $class_object = $class->__meta__;
    for my $key ($class->property_names) {
        my $meta = $class_object->property_meta_for_name($key);
        if ($meta->{is_input} && exists $params{$key}) {
            $is_input{$key} = $params{$key};
        } elsif ($meta->{is_param} && exists $params{$key}) {
            $is_param{$key} = $params{$key};
        }
    }

    my $inputs_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_input);
    my $params_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_param);

    my %software_result_params = (
        params_id => $params_bx->id,
        inputs_id => $inputs_bx->id,
        subclass_name => $class,
    );

    return {
        software_result_params => \%software_result_params,
        subclass => $class,
        inputs => \%is_input,
        params => \%is_param,
    };
}

sub _needs_symlinks_followed_when_syncing { 
    return 0;
}

sub _working_dir_prefix {
    return "detector_vcf_results";
}

sub resolve_allocation_disk_group_name { 
    return "info_genome_models";
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

    my $input = $self->input;

    return (
        $input->add_user(user => $self, label => 'uses')
    );
}

1;
