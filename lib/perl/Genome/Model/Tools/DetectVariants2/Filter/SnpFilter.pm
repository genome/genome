package Genome::Model::Tools::DetectVariants2::Filter::SnpFilter;

use strict;
use warnings;

use Genome;
use IO::File;

class Genome::Model::Tools::DetectVariants2::Filter::SnpFilter{
    is => ['Genome::Model::Tools::DetectVariants2::Filter'],
    doc => 'Filters out snvs that are around indels',
    has => [
        min_mapping_quality => {
            type => 'String',
            default => '40',
            is_optional => 1,
            is_input => 1,
            doc => 'minimum average mapping quality threshold for a high quality call',
        },
    ],
};

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants2 filter snp-filter
EOS
}

sub help_detail {
    return <<EOS 
Filters out snvs that are around indels
EOS
}

sub _variant_type { 'snvs' };

sub _filter_variants {
    my $self = shift;

    my $snv_input_file = $self->input_directory . "/snvs.hq";
    my $filtered_snv_output_file = $self->_temp_staging_directory . "/snvs.hq";
    my $fail_filter_snv_output_file = $self->_temp_staging_directory . "/snvs.lq";

    unless (-e $snv_input_file) {
        $self->error_message("snv input file $snv_input_file does not exist.");
        die $self->error_message;
    }

    # This is where samtools would have put an indel file if one was generated
    my $filtered_indel_file = $self->detector_directory . "/indels_all_sequences.filtered";
    unless (-e $filtered_indel_file ) {
        $filtered_indel_file = $self->_generate_indels_for_filtering;
    }

    if (-s $snv_input_file) {
        my $snp_filter = Genome::Model::Tools::Sam::SnpFilter->create(
            snp_file   => $snv_input_file,
            out_file   => $filtered_snv_output_file,
            lq_output  => $fail_filter_snv_output_file,
            indel_file => $filtered_indel_file,
            min_mapping_quality => $self->min_mapping_quality,
        );
        unless($snp_filter->execute) {
            $self->error_message("Running sam snp-filter failed.");
            return;
        }
        my $convert = Genome::Model::Tools::Bed::Convert::Snv::SamtoolsToBed->create( 
            source => $filtered_snv_output_file, 
            output => $self->_temp_staging_directory . "/snvs.hq.bed");

        unless($convert->execute){
            $self->error_message("Failed to convert filter output to bed.");
            die $self->error_message;
        }

        my $convert_lq = Genome::Model::Tools::Bed::Convert::Snv::SamtoolsToBed->create( 
            source => $fail_filter_snv_output_file, 
            output => $self->_temp_staging_directory . "/snvs.lq.bed");

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

sub _generate_indels_for_filtering {
    my $self = shift;

    # TODO grab samtools version and parameters by parsing the path of the input directory
    my $version = $self->detector_version;
    my $parameters = $self->detector_params;

    my $sam_pathname = Genome::Model::Tools::Sam->path_for_samtools_version($version);
    my $bam_file = $self->aligned_reads_input;
    my $ref_seq_file = $self->reference_sequence_input;
    my $samtools_cmd = "$sam_pathname pileup -c $parameters -f $ref_seq_file %s $bam_file > %s";

    my $indel_output_file = $self->input_directory . "/indels_all_sequences";
    my $filtered_indel_file = $self->_temp_staging_directory . "/indels_all_sequences.filtered";

    my $indel_cmd = sprintf($samtools_cmd, '-i', $indel_output_file);
    my $rv = Genome::Sys->shellcmd(
        cmd => $indel_cmd,
        input_files => [$bam_file, $ref_seq_file],
        output_files => [$indel_output_file],
        allow_zero_size_output_files => 1,
    );
    unless($rv) {
        $self->error_message("Running samtools indel failed.\nCommand: $indel_cmd");
        return;
    }

    if (-e $indel_output_file and not -s $indel_output_file) {
        $self->warning_message("No indels detected.");
    }

    if (-s $indel_output_file) {
        my %indel_filter_params = ( indel_file => $indel_output_file, out_file => $filtered_indel_file );
        # for capture data we do not know the proper ceiling for depth
        # this was previously set if capture_set_input was. Capture set input is not currently supported, though
        # $indel_filter_params{max_read_depth} = 1000000;
        my $indel_filter = Genome::Model::Tools::Sam::IndelFilter->create(%indel_filter_params);
        unless($indel_filter->execute) {
            $self->error_message("Running sam indel-filter failed.");
            return;
        }
    }
    else {
        Genome::Sys->write_file($filtered_indel_file);
    }


    return $filtered_indel_file;
}

1;
