package Genome::Model::Tools::Vcf::Convert::Indel::Pindel;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;
use File::Basename;

class Genome::Model::Tools::Vcf::Convert::Indel::Pindel {
    is =>  'Genome::Model::Tools::Vcf::Convert::Base' ,
    doc => 'Generate a VCF file from pindel output',
    has => [
        _refseq => {
            is => 'Text',
            calculate_from => ['reference_sequence_input'],
            calculate => q| $reference_sequence_input |,
        },
        _unsorted_output_file => {
            is => 'Text',
            calculate_from => ['output_file'],
            calculate => q/"$output_file.unsorted"/,
        },
    ],
};

sub help_synopsis {
    <<'HELP';
    Generate a VCF file from pindel indel output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the indels.
HELP
}

sub source {
    my $self = shift;
    return "Pindel";
}

sub execute {
    my $self = shift;
    my $input_id       = $self->input_id;
    my $output         = $self->_unsorted_output_file;
    my $pindel_raw     = $self->input_file;
    my $refbuild_id    = $self->reference_sequence_build->id;
    my $target_sample  = $self->aligned_reads_sample;
    my $control_sample = $self->control_aligned_reads_sample;
    my ($output_directory) = dirname($self->output_file);
    $self->debug_message("Output Directory for pindel vcf creation will be: ".$output_directory);

    my $pindel2vcf_path;
    if ($input_id) {
        my $sr = Genome::SoftwareResult->get($input_id);
        my $pindel_version  = $sr->detector_version;
        $pindel2vcf_path = Genome::Model::Tools::Pindel::RunPindel->path_for_pindel2vcf_version($pindel_version);
        die 'Failed to get pindel2vcf_path for pindel version: '. $pindel_version unless $pindel2vcf_path;
    }
    else {
        my $default_pindel_version = Genome::Model::Tools::Pindel::RunPindel->default_pindel_version;
        $pindel2vcf_path = Genome::Model::Tools::Pindel::RunPindel->path_for_pindel2vcf_version($default_pindel_version);
        $self->debug_message("pindel2vcf tool: $pindel2vcf_path used for default pindel version: $default_pindel_version");
    }

    my @prop_names = qw(
        tool_path
        pindel_raw_output 
        output_file
        reference_build_id
        aligned_reads_sample
    );
    push @prop_names, 'control_aligned_reads_sample' if $control_sample;

    my %inputs = (
        tool_path            => $pindel2vcf_path,
        pindel_raw_output    => $pindel_raw,
        output_file          => $output,
        reference_build_id   => $refbuild_id,
        aligned_reads_sample => $target_sample,
    );
    $inputs{control_aligned_reads_sample} = $control_sample if $control_sample;

    $self->debug_message("VCF conversion output will be at: ".$output);
    
    my $workflow = Genome::WorkflowBuilder::DAG->create(
        name => 'Multi-Vcf Merge',
    );

    my $pindel2vcf = Genome::WorkflowBuilder::Command->create(
        name => "Pindel2Vcf",
        command => "Genome::Model::Tools::Pindel::RunPindel2Vcf",
    );
    $workflow->add_operation($pindel2vcf);
    $workflow->recursively_set_log_dir($output_directory);

    for my $prop_name (@prop_names) {
        $workflow->connect_input(
            input_property   => $prop_name,
            destination => $pindel2vcf,
            destination_property  => $prop_name,
        );
    }

    $workflow->connect_output(
        source  => $pindel2vcf,
        source_property   => "output_file",
        output_property  => "output",
    );

    $self->debug_message("Now launching the vcf-merge workflow.");
    my $result = $workflow->execute(inputs => \%inputs);

    unless($result){
        die $self->error_message("Workflow did not return correctly.");
    }
    my $unzipped_output = Genome::Sys->create_temp_file_path;
    Genome::Model::Tools::Joinx::Sort->execute(
        input_files => [$output],
        output_file => $unzipped_output,
    );

    Genome::Sys->gzip_file($unzipped_output, $self->output_file);

    return 1;
}

