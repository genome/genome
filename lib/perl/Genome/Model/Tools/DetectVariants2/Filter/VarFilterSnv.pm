package Genome::Model::Tools::DetectVariants2::Filter::VarFilterSnv;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::DetectVariants2::Filter::VarFilterSnv{
    is  => ['Genome::Model::Tools::DetectVariants2::Filter'],
    doc => 'Filters snvs generated from samtools mpileup',
};


sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants2 filter var-filter-snv
EOS
}


sub help_detail {
    return <<EOS 
Filters out snvs 
EOS
}


sub _variant_type { 'snvs' };

sub _filter_variants {
    my $self = shift;

    my $vcf_input_file              = $self->input_directory . '/vars.vcf';
    my $filtered_snv_output_file    = $self->_temp_staging_directory . '/snvs.hq';
    my $fail_filter_snv_output_file = $self->_temp_staging_directory . '/snvs.lq';

    unless (-e $vcf_input_file) {
        $self->error_message("snv input file $vcf_input_file does not exist.");
        die $self->error_message;
    }

    # is input_file .vcf ?
    if (-s $vcf_input_file) {
        my $var_filter = Genome::Model::Tools::Sam::VarFilter->create(
            input_var_file         => $vcf_input_file,
            snv_out_file           => $filtered_snv_output_file,
            filtered_snv_out_file  => $fail_filter_snv_output_file,
            input_var_format       => 'vcf',
            use_version            => $self->detector_version, 
        );
        unless($var_filter->execute) {
            $self->error_message("Running var-filter failed.");
            return;
        }
        my $convert = Genome::Model::Tools::Bed::Convert::Snv::SamtoolsToBed->create( 
            source => $filtered_snv_output_file, 
            output => $self->_temp_staging_directory . "/snvs.hq.bed",
        );

        unless($convert->execute){
            $self->error_message("Failed to convert filter output to bed.");
            die $self->error_message;
        }

        my $convert_lq = Genome::Model::Tools::Bed::Convert::Snv::SamtoolsToBed->create( 
            source => $fail_filter_snv_output_file, 
            output => $self->_temp_staging_directory . "/snvs.lq.bed",
        );

        unless($convert_lq->execute){
            $self->error_message("Failed to convert failed-filter output to bed.");
            die $self->error_message;
        }
    }
    else {
        #FIXME use Genome::Sys... might need to add a method there
        my $hq_bed = $self->_temp_staging_directory . "/snvs.hq.bed";
        my $lq_bed = $self->_temp_staging_directory . "/snvs.lq.bed";
        my $cmd = "touch \"$hq_bed\" \"$lq_bed\" \"$filtered_snv_output_file\" \"$fail_filter_snv_output_file\"";
        Genome::Sys->shellcmd(
            cmd => $cmd
        );
    }
    return 1;
}


sub _check_native_file_counts { #To overwrite the base method
    my $self = shift;

    my $hq_output_file = $self->output_directory."/".$self->_variant_type.".hq.bed";
    my $lq_output_file = $self->output_directory."/".$self->_variant_type.".lq.bed";
    my $hq_detector_style_file = $self->output_directory."/".$self->_variant_type.".hq";
    my $lq_detector_style_file = $self->output_directory."/".$self->_variant_type.".lq";

    my $total_hq_ct = $self->line_count($hq_output_file);
    my $total_lq_ct = $self->line_count($lq_output_file);

    #line count does not include header.
    chomp(my $total_hq_detector_ct = qx(grep -vP '^#' $hq_detector_style_file | wc -l));
    chomp(my $total_lq_detector_ct = qx(grep -vP '^#' $lq_detector_style_file | wc -l));

    unless ($total_hq_ct == $total_hq_detector_ct) {
        die $self->error_message("HQ snv line counts do not match. Bed output: $total_hq_ct \t Detector-style output: $total_hq_detector_ct");
    }
    unless ($total_lq_ct == $total_lq_detector_ct) {
        die $self->error_message("LQ snv line counts do not match. Bed output: $total_lq_ct \t Detector-style output: $total_lq_detector_ct");    
    }

    return 1;
}


1;
