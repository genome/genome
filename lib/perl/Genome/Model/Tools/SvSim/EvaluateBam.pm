package Genome::Model::Tools::SvSim::EvaluateBam;

use strict;
use warnings;

use Data::Dumper;
use YAML qw(LoadFile DumpFile);
use File::Slurp qw(write_file);
use Genome;

my $BD_SEQUENCES = [1..22, qw(X Y MT CTX)];

class Genome::Model::Tools::SvSim::EvaluateBam {
    is => ["Genome::Model::Tools::SvSim::Worklet", "Command::V2"],
    has_input => [
        build_id => {
            doc => "Build to evaluate. Use this OR bam_path",
            is_optional => 1,
        },

        bam_path => {
            doc => "Full path to bam file to analyze. Use this OR build_id",
            is_optional => 1,
        },

        breakdancer_version => {
            doc => "Version of breakdancer to use",
            default_value => "1.4.2",
        },

        output_directory => {
            doc => "Where to store results",
        }
    ],
};

sub execute {
    my $self = shift;

    my ($dag, %inputs) = $self->generate_workflow;

    my $workflow_xml_path = $self->subpath("workflow.xml");
    my $workflow_inputs_path = $self->subpath("workflow-inputs.yaml");
    DumpFile($workflow_inputs_path, \%inputs);
    write_file($workflow_xml_path, $dag->get_xml);

    my $rv = $dag->execute(%inputs);
    die "Workflow failed" unless $rv;
    $self->status_message("Return value: " . Dumper($rv));

    return 1;
}


sub subpath {
    my ($self, @components) = @_;
    return join("/", $self->output_directory, @components);
}

sub _get_bam {
    my $self = shift;

    my $bam_path = $self->bam_path;
    if (defined $bam_path) {
        die "Bam file $bam_path does not exist or is empty" unless -s $bam_path;
        return $bam_path;
    }

    my $build_id = $self->build_id;
    if (!defined $build_id) {
        die "Specify one of --build-id or --bam-path";
    }

    my $build = Genome::Model::Build->get($build_id);
    die "Unable to find build $build_id" unless $build;

    my $alignment = $build->merged_alignment_result;
    die "No merged alignment result for build $build_id" unless $alignment;

    $bam_path = $alignment->bam_path;
    die "Bam file $bam_path not found for build $build_id" unless -s $bam_path;

    $self->status_message("Analyzing bam file $bam_path");
    return $bam_path;
}


sub generate_workflow {
    my ($self) = @_;
    my $bam_path = $self->_get_bam;

    my $sv_dir = $self->subpath("sv");
    my $analysis_dir = $self->subpath("analysis");
    Genome::Sys->create_directory($sv_dir);
    Genome::Sys->create_directory($analysis_dir);

    my $bam_config_path = "$sv_dir/bam_config";
    my $breakdancer_output_directory = "$sv_dir/breakdancer";
    my $merged_breakdancer_results_path = "$sv_dir/breakdancer/svs.txt";

    my @sv_types = qw(ins del inv);
    my %bd_bed_paths = map {
        sprintf("breakdancer_%s_bed_path", $_) => "$analysis_dir/bd-$_.bed"
        } @sv_types;

    $bd_bed_paths{breakdancer_ctx_bedpe_path} = "$analysis_dir/bd-ctx.bedpe";

    my %inputs = (
        %bd_bed_paths,
        input_bams => [$bam_path],
        merged_breakdancer_results_path => $merged_breakdancer_results_path,
        bam_config_path => $bam_config_path,
        breakdancer_output_directory => $breakdancer_output_directory,
        breakdancer_sequences => $BD_SEQUENCES,
        );


    my $log_dir = $self->subpath("logs");
    Genome::Sys->create_directory($log_dir);

    my $model = $self->create_model(
        name => "sv evaluation",
        log_dir => $log_dir
        );

    my $bam_config_op = $self->add_command(
        name => "Generate Bam Config",
        command => "Genome::Model::Tools::SvSim::CreateBamConfig"
        );

    $self->connect_inputs($bam_config_op,
        ["bam_config_path" => "output_file"],
        ["input_bams" => "bam_files"]
        );

    $self->connect_outputs($bam_config_op,
        ["result" => "bam_config_result"]);

    my $breakdancer_op = $self->add_command(
        name => "Breakdancer",
        command => "Genome::Model::Tools::SvSim::RunBreakdancer"
        );

    $breakdancer_op->parallel_by("sequence_name");
    $self->connect_inputs($breakdancer_op,
        ["breakdancer_output_directory" => "output_directory"],
        ["breakdancer_sequences" => "sequence_name"],
        );

    $self->create_links($bam_config_op, $breakdancer_op,
        ["output_file" => "config_path"],
        );

    $self->connect_outputs($breakdancer_op,
        ["result" => "breakdancer_result"]
        );

    my $merge_breakdancer_op = $self->add_command(
        name => "Merge breakdancer results",
        command => "Genome::Model::Tools::SvSim::MergeBreakdancerResults",
        );

    $self->create_links($breakdancer_op, $merge_breakdancer_op,
        ["output_path_svs" => "input_files"],
        );

    $self->connect_inputs($merge_breakdancer_op,
        ["merged_breakdancer_results_path" => "output_file"]
        );

    $self->connect_outputs($merge_breakdancer_op,
        ["result" => "merge_bd_result"]
        );

    my $bd2bed_op = $self->add_command(
        name => "Breakdancer to bed",
        command => "Genome::Model::Tools::SvSim::BreakdancerToBed",
        );

    $self->create_links($merge_breakdancer_op, $bd2bed_op,
        ["output_file" => "input_file"]
        );

    $self->connect_inputs($bd2bed_op,
        ["breakdancer_del_bed_path" => "deletion_output_bed"],
        ["breakdancer_ins_bed_path" => "insertion_output_bed"],
        ["breakdancer_inv_bed_path" => "inversion_output_bed"],
        ["breakdancer_ctx_bedpe_path" => "ctx_output_bedpe"],
        );

    $self->connect_outputs($bd2bed_op,
        ["result" => "bd2bed_result"]
        );

    return $model, %inputs;
}


1;
