package Genome::Model::Tools::DetectVariants2::Result::Vcf;

use strict;
use warnings;

use Genome;
use File::Copy;
use Sys::Hostname;

my $VCF_VERSION = Genome::Model::Tools::Vcf->get_vcf_version;

class Genome::Model::Tools::DetectVariants2::Result::Vcf {
    is  => ['Genome::Model::Tools::DetectVariants2::Result::Base'],
    is_abstract => 1,
    has => [
        _disk_allocation => {
            is => 'Genome::Disk::Allocation',
            is_optional => 1,
            is_many => 1,
            reverse_as => 'owner'
        },
    ],
    has_input => [
        input_id => {
            is => 'Text',
        },

    ],
    has_param => [
        vcf_version => {
            is => 'Text',
            default => $VCF_VERSION,
        },
        aligned_reads_sample => {
            is => 'Text',
            doc => 'sample name of the aligned_reads input',
        },
        control_aligned_reads_sample => {
            is => 'Text',
            doc => 'sample name of the control_aligned_reads input',
            is_optional => 1,
        },
    ],
    has_optional => [
        input => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            id_by => 'input_id',
        },
        input_directory => {
            is => 'Text',
            via => 'input',
            to => 'output_dir',
        },
    ],
};

sub _variant_type { die 'override _variant_type' };

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    unless ($self){
        die $self->error_message("Failed to create software-result!");
    }

    unless($self->_prepare_output_directory){
        die $self->error_message("Could not prepare output directory!");
    }

    unless($self->_validate_input){
        die $self->error_message("Could not validate input");
    }

    unless($self->_remove_existing_vcf){
        die $self->error_message("Could not remove existing vcf..");
    }

    unless($self->_generate_vcf){
        die $self->error_message("Could not generate vcf..");
    }

    unless($self->_reallocate_disk_allocation) {
        #this is suboptimal, but no need to die as our data should still be good and complete
        $self->warning_message('Failed to reallocate disk allocation.');
    }

    unless($self->_add_as_user_of_inputs) {
        die $self->error_message('Failed to add self as user of inputs.');
    }

    return $self;
}

sub conversion_class_name {
    my $self = shift;
    my $detector_class = shift;
    my $variant_type = shift;

    my @words = split "::",$detector_class;
    my $detector = $words[-1];

    my $vcf_module_base = 'Genome::Model::Tools::Vcf::Convert';

    my $vt = $variant_type;
    $vt =~ s/s$//;
    my $vcf_module = join('::', $vcf_module_base, ucfirst($vt) , $detector);
    eval {
        $vcf_module->__meta__;
    };
    if($@){
        $self->status_message("Did not find vcf converter for: ".$vcf_module);
        return undef;
    }
    return $vcf_module;
}

sub get_vcf {
    my $self = shift;
    my $type = shift;
    my $file = $self->output_dir."/".$type.".vcf.gz";
    unless(-e $file){
        die $self->status_message("Cannot find a vcf in this result for ".$type." variants at: ".$file);
        return;
    }
    return $file;
}

sub _generate_vcf {
    my $self = shift;
    die $self->error_message("Overload _generate_vcf in the subclass!");
}

sub _run_vcf_converter {
    my $self = shift;
    die $self->error_message("Overload _run_vcf_converter in the subclass!");
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

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $hostname = hostname;
    my $user = $ENV{'USER'};
    my $base_dir = sprintf("detect-variants--%s-%s-%s-%s", $hostname, $user, $$, $self->id);
    return join('/', 'build_merged_alignments', $base_dir);
};

sub _reallocate_disk_allocation {
    my $self = shift;
    my $allocation = $self->_disk_allocation;
    $self->status_message('Resizing the disk allocation...');
    my $rv = eval { $allocation->reallocate };
    my $error = $@;
    if ($rv != 1) {
        my $warning_message = 'Failed to reallocate disk allocation (' . $allocation->__display_name__ . ').';
        $warning_message   .= " Error: '$error'." if $error;
        $self->warning_message($warning_message);
    }
    return $rv;
}


sub _validate_input {
    my $self = shift;

    my $input_dir = $self->input_directory;
    unless (Genome::Sys->check_for_path_existence($input_dir)) {
        $self->error_message("input_directory input $input_dir does not exist");
        return;
    }

    return 1;
}

sub get_module_name_from_class_name {
    my $self = shift;
    my $class = shift;
    my @words = split('::', $class);
    my $retval = 1;

    unless(scalar(@words) > 2 and $words[0] eq 'Genome') {
        die('Could not determine proper class-name automatically.');
    }
    return $words[-1];
}

sub line_count {
    my $self = shift;
    my $input = shift;
    unless( -e $input ) {
        die $self->error_message("Could not locate file for line count: $input");
    }
    my $result = `wc -l $input`; 
    my ($answer)  = split /\s/,$result;
    return $answer
}

sub _add_as_user_of_inputs {
    my $self = shift;

    my $input = $self->input;

    return (
        $input->add_user(user => $self, label => 'uses')
    );
}

sub _remove_existing_vcf {
    my $self = shift;
    
    my $path = $self->input_directory;    

    my @existing_vcfs;
    for my $type ("snvs","indels"){
        my $file = $path."/".$type.".vcf.gz";
        push @existing_vcfs, $file if -e $file;
    }
    $self->status_message("Found no existing vcfs.") unless @existing_vcfs;

    map { $self->unlink_existing_vcf($_); } @existing_vcfs;

    return 1;
}

sub unlink_existing_vcf {
    my $self = shift;
    my $file = shift;
    if( not -l $file){
        unless(-e $file){
            next;
        }
        $self->status_message("Removing existing vcf at: ".$file);
        unless(unlink ( $file )){
            die $self->error_message("Could not unlink existing vcf at: ".$file);
        }
    } else {
        my @vcf_results = Genome::Model::Tools::DetectVariants2::Result::Vcf->get(input_id => $self->input_id);
        my $vcf_version = Genome::Model::Tools::Vcf->get_vcf_version;
        if(@vcf_results > 0){
            for my $existing_result (@vcf_results) {
                if($self->compare_vcf_versions($existing_result->vcf_version,$vcf_version)){
                    die $self->error_message("Found an existing vcf result with a greater vcf_version (".$existing_result->vcf_version.") than the one I wish to make (".$vcf_version.").");
                }
            }
        }
        unless(unlink $file){
            die $self->error_message("Could not unlink vcf link at: ".$file);
        }
        $self->status_message("Unlinked vcf symlink.");
    }
    return 1;
}


sub compare_vcf_versions {
    my $self = shift;
    my $vcf_a = shift;
    my $vcf_b = shift;

    my ($vcf_a_is_int) = ( $vcf_a =~ m/^\d+$/ );
    my ($vcf_b_is_int) = ( $vcf_b =~ m/^\d+$/ );
    
    if($vcf_a_is_int && $vcf_b_is_int){
        return ($vcf_a > $vcf_b);
    } elsif( $vcf_a_is_int) {
        return 1;
    } elsif( $vcf_b_is_int) {
        return 0;
    } else {
        return ($vcf_a gt $vcf_b);
    }
}

1;
