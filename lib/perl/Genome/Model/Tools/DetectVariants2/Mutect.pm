package Genome::Model::Tools::DetectVariants2::Mutect;

use warnings;
use strict;

use Genome;
use Data::Dumper;

my $DEFAULT_VERSION = 'test';

class Genome::Model::Tools::DetectVariants2::Mutect {
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    doc => "Produces a list of high confidence somatic snps.",
    has_optional => [
        number_of_chunks => {
            is => 'Integer',
            doc => 'number of chunks to split the genome into for mutect to run.',
        },
    ],
    has_param => [
        lsf_resource => {
            default_value => Genome::Config::get('lsf_resource_dv2_mutect'),
        },
        lsf_queue => {
            default_value => Genome::Config::get('lsf_queue_dv2_workflow'),
        },
    ],
};

sub help_synopsis {
    my $self = shift;
    return <<"EOS"

    gmt detect-variant2 mutect something useful...
EOS
}

sub help_detail {
    return <<EOS
    Provide a tumor and normal BAM file and get a list of somatic SNVs.
EOS
}

sub _detect_variants {
    my $self = shift;

    $self->debug_message("beginning execute");


    my $output_dir = $self->_temp_staging_directory;
    my $scratch_dir = $self->_temp_scratch_directory;
    my @basic_params = (
        '--normal-bam' => $self->control_aligned_reads_input,
        '--tumor-bam' => $self->aligned_reads_input,
        '--version' => $self->version,
        '--reference' => $self->reference_sequence_input,
        '--output-file' => $self->_snv_staging_output,
        '--vcf' => $output_dir . "/"  . $self->_snv_staging_output . ".raw.vcf",
    );
    # Update the default parameters with those passed in.
    my %mutect_params = $self->parse_params($self->params, @basic_params);
    if($self->number_of_chunks && $self->number_of_chunks > 1) {
        #die "Unimplemented\n";
        my $reference = Genome::File::Fasta->create(id => $self->reference_sequence_input) or die "Unable to create Genome::File::Fasta object for chunking of jobs\n";
        my @chunks = $reference->divide_into_chunks($self->number_of_chunks);

        my $dag = Genome::WorkflowBuilder::DAG->create(
            name => 'Mutect Parallel Workflow ' . $self->id,
        );
        my $op = Genome::WorkflowBuilder::Command->create(
            name => 'Mutect x' . scalar(@chunks),
            command => 'Genome::Model::Tools::Mutect::ParallelWrapper',
        );
        $op->parallel_by('chunk_num');
        $dag->add_operation($op);

        my $log_dir = $self->output_directory ."/mutect_by_chunk/";
        if(my $parent_dir = Genome::WorkflowBuilder::DAG->parent_log_dir) {
            $log_dir = $parent_dir;
        }
        $dag->recursively_set_log_dir($log_dir);

        $mutect_params{chunk_num} = [1..scalar(@chunks)];
        $mutect_params{total_chunks} = $self->number_of_chunks;
        $mutect_params{fasta_object} = $reference;
        $mutect_params{basename} = $scratch_dir . "/mutect";
        delete $mutect_params{output_file};
        delete $mutect_params{vcf};
        delete $mutect_params{coverage_file};

        for my $param (keys %mutect_params) {
            $dag->connect_input(
                input_property => $param,
                destination => $op,
                destination_property => $param,
            );
        }

        for my $param (qw( vcf output_file result )) {
            $dag->connect_output(
                source => $op,
                source_property => $param,
                output_property => $param,
            );
        }

        my $output = $dag->execute(inputs => \%mutect_params);
        print Dumper $output,"\n";
        my $merger = Genome::Model::Tools::Mutect::MergeOutputFiles->create(mutect_output_files => $output->{output_file}, merged_file => $self->_snv_staging_output);
        unless($merger->execute()) {
            die "Error merging mutect sub-job output files\n";
        }

        my $vcf_merger = Genome::Model::Tools::Joinx::VcfMerge->create(input_files => $output->{vcf}, output_file => $self->_snv_staging_output . ".raw.vcf", merge_samples => 1);
        unless($vcf_merger->execute()) {
            die "Error merging mutect sub-job vcf files\n";
        }
    }
    else {
        my $mutect = Genome::Model::Tools::Mutect->create(
            %mutect_params,
        );
        unless($mutect->execute()) {
            $self->error_message('Failed to execute Mutect: ' . $mutect->error_message);
            return;
        }

        return 1;
    }
    #How is filtering implemented??
    #Scott overrode _run_bed_converter from DetectVariants2::Base below
    #It runs the runs the BED converter twice passing in a parameter that grabs either the passing (hq) or failing (lq) variants from Strelka itself

    $self->debug_message("ending execute");
    return 1;
}


sub _run_bed_converter {
    my $self = shift;
    my $converter = shift;
    my $source = shift;

    my $hq_output = $source . '.bed';

    my $command1 = $converter->create(
        source => $source,
        output => $hq_output,
        reference_build_id => $self->reference_build_id,
        limit_variants_to => 'hq',
    );

    unless($command1->execute) {
        $self->error_message('Failed to convert ' . $source . ' to the standard format.');
        return;
    }

    my $lq_output = $hq_output;
    $lq_output =~ s/\.hq\./\.lq\./ or die "file did not match expected pattern *.hq.*!: $hq_output";

    my $command2 = $converter->create(
        source => $source,
        output => $lq_output,
        reference_build_id => $self->reference_build_id,
        limit_variants_to => 'lq',
    );

    unless($command2->execute) {
        $self->error_message('Failed to convert ' . $source . ' to the standard format.');
        return;
    }
    return 1;
}


sub has_version {
    my $self = shift;
    my $version = shift;
    unless(defined($version)){
        $version = $self->version;
    }
    my %available_version = Genome::Model::Tools::Mutect->mutect_versions;
    if(exists($available_version{$version})){
        return 1;
    }
    return 0;
}

sub parse_params {
    my ($self, $string, @additional_params) = @_;

    my ($dv2_params, $mutect_params) = split ";", $string;
    if(defined $mutect_params) {
        #then we know we have dv2 params
        #only one is allowed so just parse it out
        my ($param_name, $value) = split " ", $dv2_params; #note that this will break if there are flags
        $param_name =~ s/^--//g;
        $param_name =~ s/-/_/g;
        unless($self->can($param_name)) {
            my $package = __PACKAGE__;
            die "$param_name specified in strategy is not available or not properly handled by $package\n";
        }
        else {
            $self->$param_name($value) if defined $value;
        }
    }
    else {
        $mutect_params = $string;
    }
    my @param_list = split(" ", $mutect_params);
    my ($cmd_class, $params) = Genome::Model::Tools::Mutect->resolve_class_and_params_for_argv(@param_list,@additional_params);
    return %$params;
}

sub _sort_detector_output {
    my $self = shift;
    return $self->SUPER::_sort_detector_output(2);
}

sub _create_temp_directories {
    my $self = shift;
    my $network_temp_staging = File::Temp::tempdir(DIR => $self->output_directory, CLEANUP => 1);
    my $scratch_temp_staging = File::Temp::tempdir(DIR => $self->output_directory, CLEANUP => 1);
    $self->_temp_staging_directory($network_temp_staging);
    $self->_temp_scratch_directory($scratch_temp_staging);
    return 1;
}

1;

