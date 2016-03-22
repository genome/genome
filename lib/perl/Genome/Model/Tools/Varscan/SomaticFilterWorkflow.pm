package Genome::Model::Tools::Varscan::SomaticFilterWorkflow;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::Varscan::SomaticFilterWorkflow {
    is => 'Genome::Model::Tools::Varscan',
    has_input => [
    normal_bam => {
        is => 'Text',
        doc => "Path to Normal BAM file",
    },
    tumor_bam => {
        is => 'Text',
        doc => "Path to Tumor BAM file",
    },
    prefix => {
        is => 'Text',
        doc => "Input VarScan SNV prefix. The SNV files are found by looking for \$prefix.\$chr.unfiltered",
    },
    reference => {
        is => 'Text',
        doc => "Reference FASTA file for BAMs",
        example_values => [(Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa')],
    },
    outdir => {
        is => 'FilesystemPath',
        doc => 'Directory where output files will be written',
    },
    bamrc_version => {
        is => 'String',
        doc => 'version of bam-readcount to use',
    },
    ],
};

sub help_brief {
    "This filter is identical to GMT Varscan SomaticParallelFilter. Uses a workflow instead of manual 'bsub'. Filters SNVs only."
}

sub help_synopsis {
    return <<EOS
Filter Varscan somatic SNVs
gmt varscan somatic-parallel-filter --normal-bam nbam --tumor-bam  tbam --output test.filtered \
 --prefix test.snps --reference ~/test.fa --bamrc_version 0.7
EOS
}

sub get_variant_files {
    my $self = shift;
    my $variant_files = shift;
    my $prefix = $self->prefix;
    my $reference = $self->reference;
    my $index_file = "$reference.fai";
    unless(-e $index_file) {
        die $self->error_message("Index file for reference ($index_file) not found!\n");
    }
    my $input = new FileHandle ($index_file);
    while (<$input>) {
        chomp;
        my $line = $_;
        my ($chrom) = split(/\t/, $line);
        next if($chrom =~ 'NT_' || $chrom =~ /GL/ || $chrom =~ /MT/ || $chrom =~ /KI/);
        my $variant_f = join "", $prefix, ".$chrom.unfiltered";
        if(-e $variant_f) {
            push @$variant_files, $variant_f;
        }
    }
    unless(@$variant_files) {
        die $self->error_message("Unable to find variant files of prefix $prefix");
    }
}

sub map_workflow_inputs {
    my $self = shift;
    my $input = shift;
    my @variant_files;
    $self->get_variant_files(\@variant_files);
    foreach my $vf(@variant_files) {
        push @{$input->{variants_file}}, $vf;
        my $formatted_output = File::Spec->join(
            $self->outdir,
            basename($vf) . ".formatted",
        );
        my $processed_somatic = File::Spec->join(
            $self->outdir,
            basename($vf) . ".formatted.Somatic.hc",
        );
        my $processed_germline = File::Spec->join(
            $self->outdir,
            basename($vf) . ".formatted.Germline.hc",
        );
        my $processed_loh = File::Spec->join(
            $self->outdir,
            basename($vf) . ".formatted.LOH.hc",
        );
        my $unfiltered_output = File::Spec->join(
            $self->outdir,
            basename($vf) . ".formatted.unfiltered",
        );
        push @{$input->{format_output}}, $formatted_output;
        push @{$input->{unfiltered}}, $unfiltered_output;
        push @{$input->{processed_somatic}}, $processed_somatic;
        push @{$input->{processed_germline}}, $processed_germline;
        push @{$input->{processed_loh}}, $processed_loh;
    }
    $input->{normal_bam} = $self->normal_bam;
    $input->{tumor_bam} = $self->tumor_bam;
}

sub check_result {
    my $self = shift;
    my $result = shift;
    unless($result){
        die $self->error_message("Workflow did not return correctly.");
    }
}

sub execute_workflow_for_operation {
    my $self = shift;
    my $op = shift;
    my $op_shortname = shift;
    my $inputs = shift;

    my $workflow = Genome::WorkflowBuilder::DAG->create(
        name => $op_shortname . ' workflow',
    );
    $workflow->add_operation($op);

    for my $input (keys %$inputs) {
        $workflow->connect_input(
            input_property => $input,
            destination => $op,
            destination_property => $input,
        );
    }

    $workflow->connect_output(
        source => $op,
        source_property => 'result',
        output_property => 'result',
    );

    $workflow->recursively_set_log_dir($self->outdir);

    $self->status_message("Launching %s workflow now.", $op_shortname);
    my $result = $workflow->execute(inputs => $inputs);
    return $result;
}

sub set_lsf_queue {
    my $self = shift;
    my $w = shift;
    my $lsf_queue = Genome::Config::get('lsf_queue_build_worker_alt');
    $w->lsf_queue($lsf_queue);
}

sub run_format_workflow {
    my $self = shift;
    my $input = shift;

    my $w = Genome::WorkflowBuilder::Command->create(
        name => "Format Varscan SNVs",
        command => 'Genome::Model::Tools::Capture::FormatSnvs',
    );
    $w->parallel_by('variants_file');
    $self->set_lsf_queue($w);

    my $inputs = {
        'variants_file' => $input->{variants_file},
        'outdir' => $self->outdir,
        'append_line' => 1,
    };

    my $result = $self->execute_workflow_for_operation($w, 'format', $inputs);

    $self->check_result($result);
}

sub run_process_workflow {
    my $self = shift;
    my $input = shift;

    my $w = Genome::WorkflowBuilder::Command->create(
        name => "Process Varscan SNVs",
        command => 'Genome::Model::Tools::Varscan::ProcessSomatic',
    );
    $w->parallel_by('status_file');
    $self->set_lsf_queue($w);

    my $inputs = {
        'status_file' => $input->{format_output},
    };

    my $result = $self->execute_workflow_for_operation($w, 'process', $inputs);

    $self->check_result($result);
}

sub run_filter_workflow {
    my $self = shift;
    my $processed_input = shift;
    my $bam = shift;

    my $w = Genome::WorkflowBuilder::Command->create(
        name => "Filter Varscan SNVs",
        command => 'Genome::Model::Tools::Somatic::FilterFalsePositives',
    );
    $w->parallel_by('variant_file');
    $self->set_lsf_queue($w);

    my $inputs = {
        'variant_file' => $processed_input,
        'bam_file' => $bam,
        'outdir' => $self->outdir,
        'reference' => $self->reference,
        'bam_readcount_version' => $self->bamrc_version,
    };

    my $result = $self->execute_workflow_for_operation($w, 'filter', $inputs);

    $self->check_result($result);
}

sub execute {
    local $ENV{WF_USE_FLOW} = 1;
    my $self = shift;
    my %input;
    $self->map_workflow_inputs(\%input);
    $self->run_format_workflow(\%input);
    $self->run_process_workflow(\%input);
    $self->run_filter_workflow($input{processed_somatic}, $input{tumor_bam});
    $self->run_filter_workflow($input{processed_germline}, $input{tumor_bam});
    $self->run_filter_workflow($input{processed_loh}, $input{normal_bam});
    return 1;
}

1;

