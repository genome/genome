package Genome::VariantReporting::Suite::Vep::RunResult;

use strict;
use warnings FATAL => 'all';
use Genome;
use Sys::Hostname;
use IPC::Run qw(run);
use File::Basename qw(dirname);
use JSON;
use Memoize qw();

my $_JSON_CODEC = new JSON->allow_nonref;

class Genome::VariantReporting::Suite::Vep::RunResult {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Result',
    has_input => [
        ensembl_version => {
            is => 'String',
        },
        custom_annotation_tags => {
            is => 'String',
            is_many => 1,
            is_optional => 1,
        },
        feature_list_ids => {
            is => 'Text',
            doc => 'A json-encoded hash keyed on INFO TAG with values of FeatureList IDs',
        },
        species => {
            is => 'String',
        },
        reference_fasta_lookup => {is => 'Path'},
    ],
    has_param => [
        plugins => {is => 'String',
                    is_many => 1},
        plugins_version => {is => 'String',},
        joinx_version => {is => 'String',},
    ],
    has_transient_optional => [
        reference_fasta => {is => 'Path'},
    ],
};

sub decoded_feature_list_ids {
    my $self = shift;

    return $_JSON_CODEC->decode($self->feature_list_ids);
}

my $BUFFER_SIZE = '5000';

sub output_filename_base {
    return 'vep.vcf';
}

sub output_filename {
    my $self = shift;
    return $self->output_filename_base.'.gz';
}

sub _run {
    my $self = shift;

    $self->_sort_input_vcf;
    $self->_strip_input_vcf;
    $self->_split_alternate_alleles;
    $self->_annotate;
    $self->_sort_annotated_vcf;
    $self->_merge_annotations;
    $self->_fix_feature_list_descriptions;
    $self->_zip;

    return;
}

sub _sort_input_vcf {
    my $self = shift;

    Genome::Model::Tools::Joinx::Sort->execute(
        input_files => [$self->input_vcf],
        use_version => $self->joinx_version,
        output_file => $self->sorted_input_vcf,
    );
}

sub sorted_input_vcf {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, 'sorted_input.vcf');
}

sub _strip_input_vcf {
    my $self = shift;

    Genome::Sys->shellcmd(
        cmd => sprintf('cat %s | cut -f -8 > %s',
            $self->sorted_input_vcf, $self->stripped_input_vcf
        ),
        input_files => [$self->sorted_input_vcf],
        output_files => [$self->stripped_input_vcf],
        set_pipefail => 0, # cut can close the pipe early
        keep_dbh_connection_open => 1, # this should run fast
    );
}

sub stripped_input_vcf {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, 'stripped_input.vcf');
}

sub _split_alternate_alleles {
    my $self = shift;

    Genome::VariantReporting::Suite::Vep::SplitAlternateAlleles->execute(
        input_file => $self->stripped_input_vcf,
        output_file => $self->split_vcf,
    );
    unlink($self->stripped_input_vcf);
}

sub split_vcf {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, 'split.vcf');
}

sub _annotate {
    my $self = shift;

    die "Failed to run VEP" unless Genome::Db::Ensembl::Command::Run::Vep->execute(
        input_file => $self->split_vcf,
        output_file => $self->vep_output_file,
        $self->vep_params,
    )->result;
    unlink $self->split_vcf;
}

sub _get_file_path_for_feature_list {
    my ($self, $id) = @_;
    my $feature_list = Genome::FeatureList->get($id);
    return $feature_list->get_tabix_and_gzipped_bed_file,
}
Memoize::memoize("_get_file_path_for_feature_list", LIST_CACHE => 'MERGE');

sub custom_annotation_inputs {
    my $self = shift;

    my $result = [];
    for my $tag ($self->custom_annotation_tags) {
        my $id = $self->decoded_feature_list_ids->{$tag};
        push @$result, join("@",
            $self->_get_file_path_for_feature_list($id),
            $tag,
            "bed",
            "overlap",
            "0",
        );
    }
    return $result;
}

sub vep_params {
    my $self = shift;

    my %param_hash = $self->_add_condel_to_plugins($self->param_hash);

    my %params = (
        %param_hash,
        fasta => $self->reference_fasta,
        ensembl_version => $self->ensembl_version,
        custom => $self->custom_annotation_inputs,
        format => "vcf",
        vcf => 1,
        quiet => 0,
        hgvs => 1,
        pick => 1,
        polyphen => 'b',
        sift => 'b',
        regulatory => 1,
        canonical => 1,
        terms => "SO",
        buffer_size => $BUFFER_SIZE,
    );
    delete $params{variant_type};
    delete $params{test_name};
    delete $params{joinx_version};

    return %params;
}

sub vep_output_file {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, 'vep.vcf');
}

sub _sort_annotated_vcf {
    my $self = shift;

    Genome::Model::Tools::Joinx::Sort->execute(
        input_files => [$self->vep_output_file],
        use_version => $self->joinx_version,
        output_file => $self->sorted_vep_output,
    );
    unlink $self->vep_output_file;
}

sub sorted_vep_output {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, 'vep_sorted.vcf');
}

sub _merge_annotations {
    my $self = shift;

    Genome::Model::Tools::Joinx::VcfMerge->execute(
        input_files => [$self->sorted_input_vcf, $self->sorted_vep_output],
        output_file => $self->final_vcf_file,
        use_version => $self->joinx_version,
        merge_strategy_file => $self->joinx_merge_strategy_file,
        exact_pos => 1,
    );
    unlink($self->sorted_input_vcf);
    unlink($self->sorted_vep_output);
}

sub joinx_merge_strategy_file {
    my $self = shift;
    return File::Spec->join(dirname(__FILE__), 'joinx_merge.strategy');
}

sub final_vcf_file {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, $self->output_filename_base);
}

sub _fix_feature_list_descriptions {
    my $self = shift;
    my @substitutions;
    for my $tag ($self->custom_annotation_tags) {
        my $id = $self->decoded_feature_list_ids->{$tag};
        my $feature_list = Genome::FeatureList->get($id);
        push @substitutions, "-e ".join("|", "s", $self->_get_file_path_for_feature_list($id),
            $feature_list->name)."|";
    }
    run("sed", "-i", @substitutions, $self->final_vcf_file)
}

sub _zip {
    my $self = shift;

    # deletes $self->final_vcf_file and creates $self->final_output_file
    run(['bgzip', $self->final_vcf_file]);
}

sub final_output_file {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, $self->output_filename);
}

sub _add_condel_to_plugins {
    my $self = shift;
    my %params = @_;

    my $condel = 'Condel@PLUGIN_DIR@b@2';
    my $condel_is_set = 0;
    for my $plugin (@{$params{plugins}}) {
        if ($plugin eq $condel) {
            $condel_is_set = 1;
            last;
        }
    }
    unless ($condel_is_set) {
        push @{$params{plugins}}, $condel;
    }
    return %params;
}
