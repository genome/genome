package Genome::Model::Tools::DetectVariants2::Filter::VarFilterIndel;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::DetectVariants2::Filter::VarFilterIndel{
    is  => ['Genome::Model::Tools::DetectVariants2::Filter'],
    doc => 'Filters indels generated from samtools mpileup',
};

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants2 filter var-filter-indel
EOS
}

sub help_detail {
    return <<EOS 
Filters out indels
EOS
}

sub _variant_type { 'indels' };

sub _filter_variants {
    my $self = shift;
    my $vcf_input_file = $self->input_directory . '/vars.vcf';

    my $hq_file = $self->_temp_staging_directory . '/indels.hq';   #pass the filter
    my $lq_file = $self->_temp_staging_directory . '/indels.lq';   #fail the filter
    my $hq_bed  = $hq_file . '.bed';
    my $lq_bed  = $lq_file . '.bed';

    unless (-e $vcf_input_file) {
        die $self->error_message("snv input file $vcf_input_file does not exist.");
    }

    my @touch_files;
    # is input_file .vcf ?
    if (-s $vcf_input_file) {
        my $var_filter = Genome::Model::Tools::Sam::VarFilter->create(
            input_var_file          => $vcf_input_file,
            indel_out_file          => $hq_file,
            filtered_indel_out_file => $lq_file,
            input_var_format        => 'vcf',
            use_version             => $self->detector_version, 
        );
        unless ($var_filter->execute) {
            die $self->error_message("Running var-filter failed.");
        }

        if (-z $hq_file) {
            push @touch_files, $hq_bed;
        }
        else {
            my $convert = Genome::Model::Tools::Bed::Convert::Indel::SamtoolsToBed->create( 
                source => $hq_file, 
                output => $hq_bed,
            );
            unless ($convert->execute) {
                die $self->error_message("Failed to convert filter output to bed.");
            }
        }

        if (-z $lq_file) {
            push @touch_files, $lq_bed;
        }
        else {
            my $convert_lq = Genome::Model::Tools::Bed::Convert::Indel::SamtoolsToBed->create( 
                source => $lq_file, 
                output => $lq_bed,
            );
            unless ($convert_lq->execute) {
                die $self->error_message("Failed to convert failed-filter output to bed.");
            }
        }
    }
    else {
        push @touch_files, $hq_file, $hq_bed, $lq_file, $lq_bed;
    }

    if (@touch_files) {
        my $cmd = "touch @touch_files";
        Genome::Sys->shellcmd(cmd => $cmd);
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

    my $hq_detector_fh = Genome::Sys->open_file_for_reading($hq_detector_style_file) or return;
    my $lq_detector_fh = Genome::Sys->open_file_for_reading($lq_detector_style_file) or return;

    my $total_hq_detector_ct = _indel_ct($hq_detector_fh);
    my $total_lq_detector_ct = _indel_ct($lq_detector_fh);
    
    unless ($total_hq_ct == $total_hq_detector_ct) {
        die $self->error_message("HQ indel line counts do not match. Bed output: $total_hq_ct \t Detector-style output: $total_hq_detector_ct");
    }
    unless ($total_lq_ct == $total_lq_detector_ct) {
        die $self->error_message("LQ indel line counts do not match. Bed output: $total_lq_ct \t Detector-style output: $total_lq_detector_ct");    
    }

    return 1;
}


sub _indel_ct {
    my $fh = shift;
    my $ct = 0;

    while (my $line = $fh->getline) {
        next if $line =~ /^#/;
        my @columns = split /\s+/, $line;
        my @indels  = split /,/, $columns[4];
        $ct += scalar @indels;
    }
        
    $fh->close;
    return $ct;
}

#overwrite base _generate_vcf method because that method 
#does not handle vcf indel correctly. This method uses
#hq vcf not hq bed as filter_file. There should not be two vcf
#indel lines with exact same chr and pos. 
sub _generate_vcf {
    my $self = shift;
    my $ori_vcf = $self->previous_result->output_dir.'/'.$self->_variant_type.'.vcf.gz';
    my $out_vcf = $self->output_directory.'/'.$self->_variant_type.'.vcf.gz';
    my $hq_vcf  = $self->output_directory.'/'.$self->_variant_type.'.hq';

    unless(-s $ori_vcf){
        $self->debug_message("Skipping VCF generation, no vcf if the previous result: $ori_vcf");
        return 1;
    }

    my $filter_name = $self->get_module_name_from_class_name(ref $self || $self);
    my $filter_description = $self->filter_description;
    my $vcf_filter_cmd = Genome::Model::Tools::Vcf::VcfFilter->create(
        output_file => $out_vcf,
        vcf_file    => $ori_vcf,
        filter_file => $hq_vcf,
        filter_keep => 1,
        filter_name => $filter_name,
        filter_description => $filter_description,
    );
        
    unless($vcf_filter_cmd->execute){
        die $self->error_message("Could not complete call to gmt vcf vcf-filter!");
    }
    return 1;
}

1;
