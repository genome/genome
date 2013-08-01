package Genome::Model::Tools::DetectVariants2::Mutect;

use warnings;
use strict;

use Genome;
use Workflow;
use Workflow::Simple;
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
        workflow_log_dir => {
            is => 'Text',
            calculate_from => 'output_directory',
            calculate => q{ return File::Basename::dirname($output_directory) . '/mutect_by_chunk_log/'; },
            is_optional => 1,
            doc => 'workflow log directory of each mutect run',
        },
    ],
    has_param => [
        lsf_resource => {
            default_value => 'rusage[mem=4000] select[type==LINUX64 && maxtmp>100000] span[hosts=1]',
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

    $self->status_message("beginning execute");


    my $output_dir = $self->_temp_staging_directory;
    my @basic_params = (   
        '--normal-bam' => $self->control_aligned_reads_input,
        '--tumor-bam' => $self->aligned_reads_input,
        '--version' => $self->version,
        '--reference' => $self->reference_sequence_input,
        '--output-file' => $self->_snv_staging_output,
        '--vcf' => $output_dir . "/"  . $self->_snv_base_name . ".raw.vcf",
        '--coverage-file' => $output_dir . "/coverage.wig",
    );
    # Update the default parameters with those passed in.
    my %mutect_params = $self->parse_params($self->params, @basic_params);
    if($self->number_of_chunks && $self->number_of_chunks > 1) {
        #die "Unimplemented\n";
        my $reference = Genome::File::Fasta->create(id => $self->reference_sequence_input) or die "Unable to create Genome::File::Fasta object for chunking of jobs\n";
        my @chunks = $reference->divide_into_chunks($self->number_of_chunks);

        my $op = Workflow::Operation->create(
            name => 'Mutect x' . scalar(@chunks),
            operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::Mutect::ParallelWrapper'),
        );
        $op->parallel_by('chunk_num');

        if ($self->workflow_log_dir) {
            unless (-d $self->workflow_log_dir) {
                unless (Genome::Sys->create_directory($self->workflow_log_dir)) {
                    die $self->error_message('Failed to create workflow_log_dir: '. $self->workflow_log_dir);
                }
            }
            $op->log_dir($self->workflow_log_dir);
        }

        $mutect_params{chunk_num} = [1..scalar(@chunks)];
        $mutect_params{total_chunks} = $self->number_of_chunks;
        $mutect_params{fasta_object} = $reference;
        $mutect_params{basename} = $output_dir . "/mutect";
        delete $mutect_params{output_file};
        delete $mutect_params{vcf};
        delete $mutect_params{coverage_file};
        my $output = Workflow::Simple::run_workflow_lsf($op, %mutect_params);
        unless (defined $output) {
            my @error;
            for (@Workflow::Simple::ERROR) {
                push @error, $_->error;
            }
            $self->error_message(join("\n", @error));
            die $self->error_message;
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

    $self->status_message("ending execute");
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
    return 1;
}

sub convert_chunks_to_mutect_params {
    my ($self, @chunks) = @_;
    #this thing is a list of array refs containing array refs
    my @mutect_params;
    for my $chunk (@chunks) {
        my @intervals;
        for my $interval (@$chunk) {
            my ($chr, $start, $stop) = @$interval;
            push @intervals, "--intervals $chr:$start-$stop";
        }
        push @mutect_params, \@intervals;
    }
    return @mutect_params;
}
    

1;

