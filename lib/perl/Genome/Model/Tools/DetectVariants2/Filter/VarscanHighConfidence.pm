package Genome::Model::Tools::DetectVariants2::Filter::VarscanHighConfidence;

use warnings;
use strict;

use File::Copy;
use Genome;

class Genome::Model::Tools::DetectVariants2::Filter::VarscanHighConfidence{
    is => 'Genome::Model::Tools::DetectVariants2::Filter',
    has_optional_input => [
        p_value_for_hc  => { is => 'Number', doc => "P-value threshold for high confidence", is_input => 1, default_value => '0.07'},
        max_normal_freq => { is => 'Number', doc => "Maximum normal frequency for HC Somatic", is_input => 1, default_value => '5'},
        min_tumor_freq  => { is => 'Number', doc => "Minimum tumor freq for HC Somatic", is_input => 1, default_value => '10'},
    ],
    has_param => [
         lsf_queue => {
             default_value => $ENV{GENOME_LSF_QUEUE_DV2_WORKFLOW},
         },
         lsf_resource => {
             default_value => "-M 8000000 -R 'select[type==LINUX64 && mem>8000] rusage[mem=8000]'",
         },
     ],
};

sub _variant_type { 'snvs' };

sub _filter_variants {
    my $self = shift;
    my $varscan_status_file = $self->input_directory."/snvs.hq";
    my $base_name = $self->_temp_staging_directory."/snvs";

    my $vshc = Genome::Model::Tools::Varscan::ProcessSomatic->create(
        status_file => $varscan_status_file,
        output_basename => $base_name,
        p_value_for_hc => $self->p_value_for_hc,
        max_normal_freq => $self->max_normal_freq,
        min_tumor_freq => $self->min_tumor_freq,
    );

    unless( $vshc->execute ){
        die $self->error_message("Execution of gmt varscan process-somatic failed.");
    }

    $self->prepare_output;

    my $lq_output = $self->_temp_staging_directory."/snvs.lq";
    my $hq_output = $self->_temp_staging_directory."/snvs.hq";

    ## Possibly move this to a class method on Filter
    my $hq_cnv_cmd = Genome::Model::Tools::Bed::Convert::Snv::VarscanSomaticToBed->create(
        source => $hq_output,
        output => $hq_output.".bed",
        reference_build_id => $self->reference_build_id,
    );

    unless( $hq_cnv_cmd->execute ){
        die $self->error_message(" Failed to execute command to convert hq output to bed format.");
    }

    my $lq_cnv_cmd = Genome::Model::Tools::Bed::Convert::Snv::VarscanSomaticToBed->create(
        source => $lq_output,
        output => $lq_output.".bed",
        reference_build_id => $self->reference_build_id,
    );

    unless( $lq_cnv_cmd->execute ){
        die $self->error_message(" Failed to execute command to convert lq output to bed format.");
    }

    return 1;
}

# Condition output into one hq (sns.Somatic.hc) and one lq file (cat all other outputs)
sub prepare_output {
    my $self = shift;
    my $somatic_hq = $self->_temp_staging_directory."/snvs.Somatic.hc";
    my $somatic_lq = $self->_temp_staging_directory."/snvs.Somatic.lc";
    my $germline = $self->_temp_staging_directory."/snvs.Germline";
    my $loh = $self->_temp_staging_directory."/snvs.LOH";
    my $other = $self->_temp_staging_directory."/snvs.other";

    my $hq_file = $self->_temp_staging_directory."/snvs.hq";
    my $lq_file = $self->_temp_staging_directory."/snvs.lq";
    my $lq_scratch_file = $self->_temp_scratch_directory."/snvs.lq";

    Genome::Sys->copy_file( $somatic_hq, $hq_file );
    
    # FIXME other is possibly not sorted by position
    my @lq_source = ($somatic_lq, $germline, $loh);
    if (-e $other) {
        push @lq_source, $other;
    }
    my $catcmd = Genome::Model::Tools::Cat->create(
        dest => $lq_scratch_file,
        source => \@lq_source,
    );
    unless( $catcmd->execute ){
        die $self->error_message("Failed to run gmt cat on lq files.");
    }
    my $lq_scratch_file_temp = $lq_scratch_file.".tmp";
    my $result = `sort -k2 -n $lq_scratch_file > $lq_scratch_file_temp`;

    unless(Genome::Model::Tools::Bed::ChromSort->execute( input => $lq_scratch_file_temp, output => $lq_file)){
        die $self->error_message("Failed to chrom sort lq output.");
    }
    return 1; 
}

#check that the natively formatted file matches expectation
sub _check_native_file_counts {
    my $self = shift;
    my $total_input = shift;

    my $hq_output_file = $self->output_directory."/".$self->_variant_type.".hq.bed";
    my $detector_style_file = $self->output_directory."/".$self->_variant_type.".hq";

    my $total_output = $self->line_count($hq_output_file);
    my $detector_style_output = $self->line_count($detector_style_file);
    my $detector_fh = Genome::Sys->open_file_for_reading($detector_style_file);
    my $first_line = $detector_fh->getline;
    my @fields = split("\t", $first_line);
    if($fields[1] && $fields[1] !~ /\d+/) { #if second column doesn't hold a position
        $detector_style_output--; #header line present
    }
    unless(($total_output - $detector_style_output) == 0){
        die $self->error_message("Total lines of detector-style output did not match total output lines. Output lines: $total_output \t Detector-style output lines: $detector_style_output");
    }

    return 1;
}


1;
