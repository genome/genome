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


sub _check_native_file_counts {
    my $self = shift;
    return $self->_check_native_file_counts_vcf;
}


1;
