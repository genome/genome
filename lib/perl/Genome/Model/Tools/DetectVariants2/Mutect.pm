package Genome::Model::Tools::DetectVariants2::Mutect;

use warnings;
use strict;

use Genome;
use Workflow;

my $DEFAULT_VERSION = 'test';

class Genome::Model::Tools::DetectVariants2::Mutect {
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    doc => "Produces a list of high confidence somatic snps.",
    has => [
        number_of_chunks => {
            is => 'Integer',
            doc => 'number of chunks to split the genome into for mutect to run.',
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

    # Update the default parameters with those passed in.
    my $mutect_param_string = parse_params($self->params);

    my %basic_params = (   
        normal_bam => $self->control_aligned_reads_input,
        tumor_bam => $self->aligned_reads_input,
        version => $self->version,
        reference => $self->reference_sequence_input,
        output_file => $output_dir . "/" . $self->_snv_staging_output,
        vcf => $output_dir . "/" . $self->_snv_base_name . ".raw.vcf",
        coverage_file => $output_dir . "/coverage.wig",
        params => $mutect_param_string,

    );

    if($self->number_of_chunks) {
        die "Unimplemented\n";

    }
    else {
        my $mutect = Genome::Model::Tools::Mutect->create(
            %basic_params,
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
    my ($self, $string) = @_;

    my ($dv2_params, $mutect_params) = split ":", $string;
    if(defined $mutect_params) {
        #then we know we have dv2 params
        #only one is allowed so just parse it out
        my ($param_name, $value) = split " ", $dv2_params; #note that this will break if there are flags
        unless($self->$param_name) {
            my $package = __PACKAGE__;
            die "$param_name specified in strategy is not available or not properly handled by $package\n";
        }
        else {
            $self->$param_name($value) if defined $value;
        }
        #then set to return mutect params which will be passed directly to the mutect object elsewhere
        return $mutect_params;
    }
    else {
        return $string;
    }
}

1;

