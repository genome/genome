package Genome::Model::ReferenceAlignment::Command::CreateMetrics::DbSnpConcordance;

use strict;
use warnings;

use File::Basename;
use Genome;

class Genome::Model::ReferenceAlignment::Command::CreateMetrics::DbSnpConcordance {
    is => 'Genome::Command::Base',
    has => [
        build => {
            doc => 'The build for which to compute dbSNP concordance',
            is => 'Genome::Model::Build::ReferenceAlignment',
            id_by => 'build_id',
            shell_args_position => 1,
        },
        build_id => {
            is => 'Integer',
            is_input => 1,
        },
    ],
    has_optional => [
        output_dir => {
            doc => "Override the default output directory",
            is => 'File',
            is_input => 1,
        },
        _snvs_bed => {
            is => 'File',
            doc => "instance variable: path to build's snv bed file",
        },
        _filtered_snvs_bed => {
            is => 'File',
            doc => "instance variable: path to build's filtered snv bed file",
        },
        _dbsnp_file => {
            is => 'File',
            doc => "instance variable: path to dbsnp build's snv bed file",
         },
    ],
    doc => "Compute dbSNP concordance for a build and store the resulting metrics the the database",
};

sub _verify_build_and_set_paths {
    my ($self, $build) = @_;

    my $bname = $build->__display_name__;
    my $dbsnp_build = $build->model->dbsnp_build;
    if (!defined $dbsnp_build) {
        die "No dbsnp_build property found on build $bname.";
    }

    my $build_rsb = $build->model->reference_sequence_build;
    my $dbsnp_rsb = $dbsnp_build->model->reference;
    if (!defined $dbsnp_rsb) {
        die "DbSnp build " . $dbsnp_build->__display_name__ . " does not define a reference sequence!";
    }

    if (!$build_rsb->is_compatible_with($dbsnp_rsb)) {
        die "Build $bname has reference sequence " . $build_rsb->__display_name__ .
            " which is incompatible with " .  $dbsnp_rsb->__display_name__ . " specified by " .
            $dbsnp_build->__display_name__;
    }

    $self->_dbsnp_file($dbsnp_build->snvs_bed());
    if (!defined $self->_dbsnp_file()) {
        die "Failed to get dbsnp file from dbsnp build " . $dbsnp_build->__display_name__;
    }

    for my $type ("snvs_bed", "filtered_snvs_bed") {
        if (!$build->can($type)) {
            die "Don't know how to find snv bed file for build $bname.";
        }
        my $snv_file = $build->$type("v1");
        if (!defined $snv_file || ! -f $snv_file) {
            die "No suitable snv bed file found for build $bname [$type].";
        }
        my $propname = "_$type";
        $self->$propname($snv_file);
    }
}

sub _gen_concordance {
    my ($self, $bed_file, $dbsnp_file, $output_path) = @_;

    my $joinx_cmd = Genome::Model::Tools::Joinx::SnvConcordanceByQuality->create(
        input_file_a => $bed_file,
        input_file_b => $dbsnp_file,
        output_file  => $output_path,
    );
    $joinx_cmd->execute() or die "joinx failed!";
}

sub execute {
    my $self = shift;

    eval {
        $self->_verify_build_and_set_paths($self->build);

        my $out_filt = $self->build->dbsnp_file_filtered;
        my $out_unfilt = $self->build->dbsnp_file_unfiltered;
        if ($self->output_dir) {
            $out_filt = join('/', $self->output_dir, basename($out_filt));
            $out_unfilt = join('/', $self->output_dir, basename($out_unfilt));
        }
        $self->_gen_concordance($self->_snvs_bed, $self->_dbsnp_file, $out_unfilt);
        $self->_gen_concordance($self->_filtered_snvs_bed, $self->_dbsnp_file, $out_filt);
    };
    if ($@) {
        $self->error_message($@);
        return;
    }

    return 1;
}

1;
