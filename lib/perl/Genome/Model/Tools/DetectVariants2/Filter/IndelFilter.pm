package Genome::Model::Tools::DetectVariants2::Filter::IndelFilter;

use strict;
use warnings;

use Genome;
use IO::File;

class Genome::Model::Tools::DetectVariants2::Filter::IndelFilter{
    is => ['Genome::Model::Tools::DetectVariants2::Filter'],
    doc => 'Filters out indels that are around indels',
    has_optional_input => [
        max_read_depth => {
            is  => 'Integer',
            doc => 'maximum read depth',
            default => 100,
        },
        min_win_size   => {
            is  => 'Integer',
            doc => 'minimum distance between two adjacent indels',
            default => 10,
        },
        scaling_factor => {
            is  => 'Integer',
            doc => 'scaling factor in score calculation',
            default => 100,
        },
    ]
};

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants2 filter indel-filter
EOS
}

sub help_detail {
    return <<EOS 
Samtools indel filter
EOS
}

sub _variant_type { 'indels' };

sub _filter_variants {
    my $self = shift;

    my $indel_input_file = $self->input_directory . "/indels.hq";
    my $filtered_indel_output_file = $self->_temp_staging_directory . "/indels.hq";
    my $fail_filter_indel_output_file = $self->_temp_staging_directory . "/indels.lq";

    unless (-e $indel_input_file) {
        $self->error_message("indel input file $indel_input_file does not exist.");
        die $self->error_message;
    }

    my @params;

    for my $param (qw(max_read_depth min_win_size scaling_factor)) {
        if($self->$param) {
            push @params,
                $param => $self->$param;
        }
    }

    if (-s $indel_input_file) {
        my $indel_filter = Genome::Model::Tools::Sam::IndelFilter->create(
            indel_file => $indel_input_file,
            out_file => $filtered_indel_output_file,
            lq_file  => $fail_filter_indel_output_file,
            @params,
        );
        unless($indel_filter->execute) {
            $self->error_message("Running sam snp-filter failed.");
            return;
        }
        my $convert = Genome::Model::Tools::Bed::Convert::Indel::SamtoolsToBed->create( 
            source => $filtered_indel_output_file, 
            output => $self->_temp_staging_directory . "/indels.hq.bed");

        unless($convert->execute){
            $self->error_message("Failed to convert filter output to bed.");
            die $self->error_message;
        }

        my $scratch_lq_bed = $self->_temp_scratch_directory . "/indels.lq.bed";
        my $convert_lq = Genome::Model::Tools::Bed::Convert::Indel::SamtoolsToBed->create( 
            source => $fail_filter_indel_output_file, 
            output => $scratch_lq_bed,
        );

        unless($convert_lq->execute){
            $self->error_message("Failed to convert failed-filter output to bed.");
            die $self->error_message;
        }

        my $v = 1; #sort all different versions of lq file
        while(-e (my $scratch_bed = $self->_temp_scratch_directory . '/indels.lq.v' . $v . '.bed')) { #is this retarded?
            my $to_sort = $scratch_bed;
            my $target;
            if($target = readlink($scratch_bed)) {
                $to_sort = $self->_temp_scratch_directory . '/' . $target;
            }

            my $scratch_dir = $self->_temp_scratch_directory;
            my $staging_dir = $self->_temp_staging_directory;

            my $output_bed = $to_sort;
            $output_bed =~ s/$scratch_dir/$staging_dir/;
            unless(-e $output_bed) {
                my $sort_lq_bed = Genome::Model::Tools::Joinx::Sort->create(
                    input_files => [$to_sort],
                    output_file => $output_bed,
                );

                unless($sort_lq_bed->execute){
                    $self->error_message("Failed to sort failed-filter output bed.");
                    die $self->error_message;
                }
            }

            if($target) {
                my $link_output = $scratch_bed;
                $link_output =~ s/$scratch_dir/$staging_dir/;
                symlink($target, $link_output);
            }

            $v++;
        }

        if(my $target = readlink($scratch_lq_bed)) {
            symlink($target, $self->_temp_staging_directory . '/indels.lq.bed');
        } else {
            my $sort_lq_bed = Genome::Model::Tools::Joinx::Sort->create(
                input_files => [$scratch_lq_bed],
                output_file => $self->_temp_staging_directory . '/indels.lq.bed',
            );

            unless($sort_lq_bed->execute){
                $self->error_message("Failed to sort failed-filter output bed.");
                die $self->error_message;
            }
        }
    } else {
        #FIXME use Genome::Sys... might need to add a method there 
        my $hq_bed = $self->_temp_staging_directory . "/indels.hq.bed";
        my $lq_bed = $self->_temp_staging_directory . "/indels.lq.bed";
        my $cmd = "touch \"$hq_bed\" \"$lq_bed\" \"$filtered_indel_output_file\" \"$fail_filter_indel_output_file\"";
        Genome::Sys->shellcmd(
            cmd => $cmd
        );
    }


    return 1;
}

sub _check_native_file_counts {
    my $self = shift;
    my $total_input = shift;

    #TODO To make an accurate check, would need to track how many samtools lines are converted to multiple lines in the BED output
    return 1;
}

1;
