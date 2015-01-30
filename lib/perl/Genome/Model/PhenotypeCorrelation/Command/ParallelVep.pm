package Genome::Model::PhenotypeCorrelation::Command::ParallelVep;


use Genome;
use Genome::File::Vep::Reader;
use Genome::File::Vep::Writer;
use Workflow::Simple;
use Carp qw/confess/;
use File::Path qw/mkpath/;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::ParallelVep {
    is => ["Command::V2"],
    doc => "Run VEP in parallel by chromosome.",
    has_input => [
        input_vcf => {
            is => "File",
            doc => "The tabix indexed vcf file of variants to annotate",
        },
        ensembl_annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            doc => 'ID of ImportedAnnotation build with the desired ensembl version.',
        },
        log_dir => {
            is => "File",
            doc => "Workflow log directory",
        },
        work_dir => {
            is => "File",
            doc => "Network accessible directory for intermediate results",
        },
        output_file => {
            is => "File",
            is_output => 1,
            doc => "The final merged output file",
        },
    ],
    has_optional_input => Genome::Db::Ensembl::Command::Run::Base->result_user_inputs(),
};

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
This command annotates vcf variants in parallel with Vep.
EOS
}

sub execute {
    my $self = shift;

    mkpath($self->work_dir);
    mkpath($self->log_dir);

    my $log_dir = $self->log_dir;
    my $vcf = $self->input_vcf;
    $self->status_message("Performing parallel vep on $vcf...");

    my $cmd = Genome::Model::Tools::Tabix::ListChromosomes->create(
        input_file => $vcf,
        suppress_output => 1,
        );
    $cmd->execute;
    my @chromosomes = $cmd->chromosomes;

    print "Got " . @chromosomes . " chromosomes\n";
    my @intermediate_files;
    my %inputs = (
        input_vcf => $vcf,
        ensembl_annotation_build_id => $self->ensembl_annotation_build->id,
        map { ($_, $self->$_) } keys %{ Genome::Db::Ensembl::Command::Run::Base->result_user_inputs() },
    );
    my @fixed_input_keys = keys %inputs;

    my @outputs = map {"job_${_}_result"} 0..$#chromosomes;

    for my $i (0..$#chromosomes) {
        my $output_file = join('/', $self->work_dir, "file$i.vep");
        my $output_file_var = "output_file_$i";
        my $region_var = "region_$i";
        my $region = $chromosomes[$i];

        push(@intermediate_files, $output_file);
        $inputs{$output_file_var} = $output_file;
        $inputs{$region_var} = $region;
    }
    my $workflow = Workflow::Model->create(
        name => "Parallel VEP annotation",
        input_properties => [keys %inputs],
        output_properties => \@outputs,
    );

    for my $i (0..$#chromosomes) {
        my $region = $chromosomes[$i];
        my $output_file_var = "output_file_$i";
        my $region_var = "region_$i";

        my $op = $workflow->add_operation(
            name => "Parallel VEP annotation (chromosome $region)",
            operation_type => Workflow::OperationType::Command->get(
                'Genome::Model::PhenotypeCorrelation::Command::VepWrapper'
            ),
        );
        for my $property (@fixed_input_keys) {
            $workflow->add_link(
                left_operation => $workflow->get_input_connector,
                left_property => $property,
                right_operation => $op,
                right_property => $property,
            );
        }
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $output_file_var,
            right_operation => $op,
            right_property => "output_file",
        );
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $region_var,
            right_operation => $op,
            right_property => "region",
        );
        $workflow->add_link(
            left_operation => $op,
            left_property => "result",
            right_operation => $workflow->get_output_connector,
            right_property => "job_${i}_result",
        );
        $op->log_dir($log_dir);
    }

    $workflow->log_dir($log_dir);
    my @errors = $workflow->validate;
    if (@errors) {
        confess $self->error_message(join("\n", @errors));
    }

    $self->status_message("Executing workflow...");
    my $output = Workflow::Simple::run_workflow_lsf($workflow, %inputs);
    unless (defined $output) {
        $self->error_message("Workflow execution failed!");
        my @error;
        for (@Workflow::Simple::ERROR) {
            push @error, $_->error;
        }
        $self->error_message(join("\n", @error));
        confess $self->error_message;
    }
    $self->status_message("Workflow completed, merging intermediate results...");

    _merge_files($self->output_file, @intermediate_files);

    return 1;
}

sub _merge_files {
    my ($target, @sources) = @_;
    my @readers = map { new Genome::File::Vep::Reader($_) } @sources;
    confess "Can't merge an empty set of files!" unless @readers;
    my $header = $readers[0]->{header};

    my $writer = new Genome::File::Vep::Writer($target, $header);
    for my $r (@readers) {
        while (my $e = $r->next) {
            $writer->write($e);
        }
    }
}

1;
