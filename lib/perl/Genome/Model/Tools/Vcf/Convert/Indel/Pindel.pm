package Genome::Model::Tools::Vcf::Convert::Indel::Pindel;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;
use Workflow;
use Workflow::Simple;
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
    my $output         = $self->output_file;
    my $pindel_raw     = $self->input_file;
    my $refbuild_id    = $self->reference_sequence_build->id;
    my $target_sample  = $self->aligned_reads_sample;
    my $control_sample = $self->control_aligned_reads_sample;

    my ($output_directory) = dirname($output);
    $self->debug_message("Output Directory for pindel vcf creation will be: ".$output_directory);

    my @prop_names = qw(
        pindel_raw_output 
        output_file
        reference_build_id
        aligned_reads_sample
    );
    push @prop_names, 'control_aligned_reads_sample' if $control_sample;

    my %inputs;
    $inputs{pindel_raw_output} = $pindel_raw;
    $inputs{output_file} = $output;
    $self->debug_message("VCF conversion output will be at: ".$output);
    $inputs{reference_build_id}   = $refbuild_id;
    $inputs{aligned_reads_sample} = $target_sample;
    $inputs{control_aligned_reads_sample} = $control_sample if $control_sample;

    my $workflow = Workflow::Model->create(
        name => 'Multi-Vcf Merge',
        input_properties  => \@prop_names,        
        output_properties => ['output'],
    );
    $workflow->log_dir($output_directory);

    my $pindel2vcf = $workflow->add_operation(
        name => "Pindel2Vcf",
        operation_type => Workflow::OperationType::Command->get("Genome::Model::Tools::Pindel::RunPindel2Vcf"),
    );

    for my $prop_name (@prop_names) {
        $workflow->add_link(
            left_operation  => $workflow->get_input_connector,
            left_property   => $prop_name,
            right_operation => $pindel2vcf,
            right_property  => $prop_name,
        );
    }

    $workflow->add_link(
        left_operation  => $pindel2vcf,
        left_property   => "output_file",
        right_operation => $workflow->get_output_connector,
        right_property  => "output",
    );

    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }
    $self->debug_message("Now launching the vcf-merge workflow.");
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %inputs);

    unless($result){
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("Workflow did not return correctly.");
    }

    return 1;
}

