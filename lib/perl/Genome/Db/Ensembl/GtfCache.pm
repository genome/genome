package Genome::Db::Ensembl::GtfCache;

use strict;
use warnings;
use Genome;
use Sys::Hostname;

my ($VEP_DIR) = Cwd::abs_path(__FILE__) =~ /(.*)\//;
my $VEP_SCRIPT_PATH = $VEP_DIR . "/Command/Vep.d/gtf2vep";

class Genome::Db::Ensembl::GtfCache {
    is => "Genome::SoftwareResult::Stageable",
    has => [
        version => {
            is => 'Text',
            doc => 'Version of ensembl db',
            is_param => 1,
        },
        species => {
            is => 'Text',
            doc => 'Species contained in cache',
            is_param => 1,
        },
        gtf_content_hash => {
            is => 'Text',
            doc => 'MD5 hash of the gtf file this cache is based on',
            is_input => 1,
        },
        reference_build_id => {
            is => "String",
            doc => 'Id of Reference build that this annotation is based on',
            is_param => 1,
        },
    ],
    has_optional_metric => [
        gtf_file_path => {
            is => 'Text',
            doc => 'Path to the original file used to create this result',
        },
        vep_version => {
            is => 'Text',
            doc => 'Version of vep to use',
        },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    $self->_prepare_staging_directory;
    $self->gtf_content_hash(Genome::Sys->md5sum($self->gtf_file_path));
    $self->lookup_hash($self->calculate_lookup_hash); #reset after modifying file_content_hash

    my $script_path = $VEP_SCRIPT_PATH.$self->vep_version.".pl";
    my $reference_build = Genome::Model::Build->get($self->reference_build_id);
    my $reference_build_temp_path = Genome::Sys->create_temp_directory;
    #Create a symlink to the reference sequence because the script creates a .index file next to the fasta file and
    #we don't want to modify the reference build directory
    Genome::Sys->create_symlink($reference_build->full_consensus_path("fa"), $reference_build_temp_path."/reference.fa");
    my $cmd = $script_path." -i ".$self->gtf_file_path." -f $reference_build_temp_path/reference.fa -d ".$self->version." -s ".$self->species." --dir ".$self->temp_staging_directory."/";
    my %params = (
        cmd=>$cmd,
        skip_if_output_is_present => 0,
        input_files => [$self->gtf_file_path],
    );

    my $annotation_api = Genome::Db::Ensembl::Api->get_or_create(version => $self->version);
    unless ($annotation_api) {
        $self->error_message("Couldn't get ensembl api for version ".$self->version);
        return;
    }

    $annotation_api->prepend_api_path_and_execute(
        %params
    );

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;
    return $self;
}

sub _modify_params_for_lookup_hash {
    my ($class, $params_ref) = @_;
    my $original_file_path = delete $params_ref->{'gtf_file_path'};
    my $vep_version = delete $params_ref->{'vep_version'};
    my $specified_checksum = $params_ref->{'gtf_content_hash'};
    $params_ref->{'gtf_content_hash'} = $class->_calculate_and_compare_md5_hashes(
        $original_file_path, $specified_checksum);
}

sub _calculate_and_compare_md5_hashes {
    my ($class, $original_file_path, $specified_checksum) = @_;

    my $checksum = $specified_checksum;
    if (defined($original_file_path) and -e $original_file_path) {
        $checksum = Genome::Sys->md5sum($original_file_path);
        if (defined($specified_checksum) and $specified_checksum ne $checksum) {
            die $class->error_message(
                'file_content_hash does not match md5sum output for original_file_path');
        }
    }

    return $checksum;
}

sub _gather_params_for_get_or_create {
    my $class = shift;
    my $bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);

    if($bx->specifies_value_for('gtf_file_path')) {
        my $original_file_path = $bx->value_for('gtf_file_path');
        my $specified_checksum;
        if ($bx->specifies_value_for('gtf_content_hash')) {
            $specified_checksum = $bx->value_for('gtf_content_hash');
        }
        $bx = $bx->add_filter('gtf_content_hash',
            $class->_calculate_and_compare_md5_hashes(
                $original_file_path, $specified_checksum));
    }

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

    #my $inputs_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_input);
    #my $params_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_param);

    my %software_result_params = (
        #params_id => $params_bx->id,
        #inputs_id => $inputs_bx->id,
        subclass_name => $class,
    );
    return {
        software_result_params => \%software_result_params,
        subclass => $class,
        inputs => \%is_input,
        params => \%is_param,
    };
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $hostname = hostname;

    my $user = $ENV{'USER'};
    my $base_dir = sprintf("ensemblgtfcache-%s-%s-%s-%s",           $hostname,       $user, $$, $self->id);
    my $directory = join('/', 'build_merged_alignments',$self->id,$base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

1;

