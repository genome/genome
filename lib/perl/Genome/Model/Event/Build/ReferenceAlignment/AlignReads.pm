package Genome::Model::Event::Build::ReferenceAlignment::AlignReads;

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Build::ReferenceAlignment::AlignReads {
    is => ['Genome::Model::Event'],
	has => [
        _calculate_total_read_count => { via => 'instrument_data'},
        #make accessors for common metrics
        (
            map {
                $_ => { via => 'metrics', to => 'value', where => [name => $_], is_mutable => 1 },
            }
            qw/total_read_count/
        ),
    ],
};

sub results_class {
    my $self = shift;
    my $model = $self->model;
    my $processing_profile = $model->processing_profile;
    my $read_aligner_name = $processing_profile->read_aligner_name;
    return 'Genome::InstrumentData::AlignmentResult::' . 
        Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($read_aligner_name);
}

sub lsf_queue {
    my $self = shift;

    if (Genome::Config->can('should_use_alignment_pd') && Genome::Config->should_use_alignment_pd($self->model)) {
        return $ENV{GENOME_LSF_QUEUE_ALIGNMENT_PROD};
    } else {
        return $ENV{GENOME_LSF_QUEUE_ALIGNMENT_DEFAULT};
    }
}

sub bsub_rusage {
    my $self = shift;
    my $delegate = $self->results_class;
    my $rusage = $delegate->required_rusage(
        instrument_data => $self->instrument_data,
        reference_build => $self->model->reference_sequence_build,
        aligner_params  => $self->model->processing_profile->read_aligner_params,
        queue           => $self->lsf_queue,
    );
    return $rusage;
}

sub metrics_for_class {
    my $class = shift;
    return();
    my $delegate = $class->results_class;
    my @metric_names = $delegate->metric_names();
    return @metric_names;
}

sub fwd_reads_passed_quality_filter_count {
    my $self = shift;
    return $self->get_metric_value('fwd_reads_passed_quality_filter_count');
}

sub _calculate_fwd_reads_passed_quality_filter_count {
    my $self = shift;
    my $fwd_reads_passed_quality_filter_count;
    
    unless($self->instrument_data->is_paired_end) {
        $fwd_reads_passed_quality_filter_count = 0;
    } else {
        my ($fwd_fastq) = $self->instrument_data->dump_sanger_fastq_files(paired_end_as_fragment => 1);
        unless(-f $fwd_fastq) {
            $self->error_message("Problem calculating metric (didn't find forward FASTQ)...this doesn't mean the step failed");
            return 0;
        }
        $fwd_reads_passed_quality_filter_count = $self->_count_reads_in_fastq_files($fwd_fastq);
    }
    return $fwd_reads_passed_quality_filter_count;
}

sub rev_reads_passed_quality_filter_count {
    my $self = shift;
    return $self->get_metric_value('rev_reads_passed_quality_filter_count');
}

sub _calculate_rev_reads_passed_quality_filter_count {
    my $self = shift;
    my $rev_reads_passed_quality_filter_count;
    
    unless($self->instrument_data->is_paired_end) {
        $rev_reads_passed_quality_filter_count = 0;
    } else {
        my ($rev_fastq) = $self->instrument_data->dump_sanger_fastq_files(paired_end_as_fragment => 2);
        unless(-f $rev_fastq) {
            $self->error_message("Problem calculating metric (didn't find reverse FASTQ)...this doesn't mean the step failed");
            return 0;
        }
        $rev_reads_passed_quality_filter_count = $self->_count_reads_in_fastq_files($rev_fastq);
    }
    return $rev_reads_passed_quality_filter_count;
}

sub total_reads_passed_quality_filter_count {
    my $self = shift;
    return $self->get_metric_value('total_reads_passed_quality_filter_count');
}

sub _calculate_total_reads_passed_quality_filter_count {
    my $self = shift;
    my $total_reads_passed_quality_filter_count;

    if($self->instrument_data->is_paired_end) {
        $total_reads_passed_quality_filter_count = $self->fwd_reads_passed_quality_filter_count + $self->rev_reads_passed_quality_filter_count;
    } else {
        my @f = grep {-f $_ } $self->instrument_data->dump_sanger_fastq_files;
        unless (@f) {
            $self->error_message("Problem calculating metric (didn't find FASTQ(s))...this doesn't mean the step failed");
            return 0;
        }
        $total_reads_passed_quality_filter_count = $self->_count_reads_in_fastq_files(@f);
    }
    return $total_reads_passed_quality_filter_count;
}

sub total_bases_passed_quality_filter_count {
    my $self = shift;
    return $self->get_metric_value('total_bases_passed_quality_filter_count');
}

sub _calculate_total_bases_passed_quality_filter_count {
    my $self = shift;

    my $total_bases_passed_quality_filter_count;
    if($self->instrument_data->is_paired_end) {
        $total_bases_passed_quality_filter_count =
            $self->fwd_reads_passed_quality_filter_count * $self->instrument_data->fwd_read_length +
            $self->rev_reads_passed_quality_filter_count * $self->instrument_data->rev_read_length;
    } else {
        $total_bases_passed_quality_filter_count = $self->total_reads_passed_quality_filter_count * $self->instrument_data->read_length;
    }
    return $total_bases_passed_quality_filter_count;
}

sub fwd_poorly_aligned_read_count {
    my $self = shift();
    
    return $self->get_metric_value('fwd_poorly_aligned_read_count');
}

sub _calculate_fwd_poorly_aligned_read_count {
    my $self = shift();
    my $fwd_poorly_aligned_read_count;
    
    unless($self->instrument_data->is_paired_end) {
        $fwd_poorly_aligned_read_count = 0;
    } else {
        my $include_flags = 0x40 | 0x4; #is first read & unmapped read (note: also use 0x8 if need to check corresponding rev read is unmapped)
        
        $fwd_poorly_aligned_read_count = $self->_count_lines_in_bam_file(include_flags => $include_flags);
    }
    
    return $fwd_poorly_aligned_read_count;
}

sub rev_poorly_aligned_read_count {
    my $self = shift();
    
    return $self->get_metric_value('rev_poorly_aligned_read_count');
}

sub _calculate_rev_poorly_aligned_read_count {
    my $self = shift();
    my $rev_poorly_aligned_read_count;
    
    unless($self->instrument_data->is_paired_end) {
        $rev_poorly_aligned_read_count = 0;
    } else {
        my $include_flags = 0x80 | 0x4; #is second read & unmapped read (note: also use 0x8 if need to check corresponding fwd read is unmapped)
        
        $rev_poorly_aligned_read_count = $self->_count_lines_in_bam_file(include_flags => $include_flags);
    }
    
    return $rev_poorly_aligned_read_count;
}

sub poorly_aligned_read_count {
    my $self = shift;
    
    return $self->get_metric_value('poorly_aligned_read_count');
}

sub _calculate_poorly_aligned_read_count {
    my $self = shift;
    
    my $poorly_aligned_read_count;
    if($self->instrument_data->is_paired_end) {
        $poorly_aligned_read_count = $self->fwd_poorly_aligned_read_count + $self->rev_poorly_aligned_read_count;
    } else {
        my $include_flags = 0x4; #is unmapped read
        
        $poorly_aligned_read_count = $self->_count_lines_in_bam_file(include_flags => $include_flags);
    }
    return $poorly_aligned_read_count;
}

sub fwd_aligned_read_count {
    my $self = shift;
    return $self->get_metric_value('fwd_aligned_read_count');
}

sub _calculate_fwd_aligned_read_count {
    my $self = shift;
    my $fwd_aligned_read_count;
    unless($self->instrument_data->is_paired_end) {
        $fwd_aligned_read_count = 0;
    } else {
        my $exclude_flags = 0x4; #unmapped read
        my $include_flags = 0x40; #is first read
        
        $fwd_aligned_read_count = $self->_count_lines_in_bam_file(include_flags => $include_flags, exclude_flags => $exclude_flags);
    }
    return $fwd_aligned_read_count;
}

sub rev_aligned_read_count {
    my $self = shift;
    return $self->get_metric_value('rev_aligned_read_count');
}

sub _calculate_rev_aligned_read_count {
    my $self = shift;
    my $rev_aligned_read_count;
    unless($self->instrument_data->is_paired_end) {
        $rev_aligned_read_count = 0;
    } else {
        my $exclude_flags = 0x4; #unmapped read
        my $include_flags = 0x80; #is second read
        
        $rev_aligned_read_count = $self->_count_lines_in_bam_file(include_flags => $include_flags, exclude_flags => $exclude_flags); 
    }
    return $rev_aligned_read_count;
}

sub aligned_read_count {
    my $self = shift;
    return $self->get_metric_value('aligned_read_count');
}

sub _calculate_aligned_read_count {
    my $self = shift;

    my $aligned_read_count;
    if($self->instrument_data->is_paired_end) {
        $aligned_read_count = $self->fwd_aligned_read_count + $self->rev_aligned_read_count;
    } else {
        my $exclude_flags = 0x4; #unmapped read
        
        $aligned_read_count = $self->_count_lines_in_bam_file(exclude_flags => $exclude_flags);
    }
    return $aligned_read_count;
}

sub aligned_base_pair_count {
    my $self = shift;
    return $self->get_metric_value('aligned_base_pair_count');
}

sub _calculate_aligned_base_pair_count {
    my $self = shift;

    my $aligned_base_pair_count;
    
    if($self->instrument_data->is_paired_end) {
        $aligned_base_pair_count =
            $self->fwd_aligned_read_count * $self->instrument_data->fwd_read_length +
            $self->rev_aligned_read_count * $self->instrument_data->rev_read_length;
    } else {
        $aligned_base_pair_count = $self->aligned_read_count * $self->instrument_data->read_length;
    }
    return $aligned_base_pair_count;
}

sub fwd_unaligned_read_count {
    my $self = shift;
    return $self->get_metric_value('fwd_unaligned_read_count');
}

sub _calculate_fwd_unaligned_read_count {
    my $self = shift;
    my $fwd_unaligned_read_count = $self->fwd_poorly_aligned_read_count;
    return $fwd_unaligned_read_count;
}

sub rev_unaligned_read_count {
    my $self = shift;
    return $self->get_metric_value('rev_unaligned_read_count');
}

sub _calculate_rev_unaligned_read_count {
    my $self = shift;
    my $rev_unaligned_read_count = $self->rev_poorly_aligned_read_count;
    return $rev_unaligned_read_count;
}

sub unaligned_read_count {
    my $self = shift;
    return $self->get_metric_value('unaligned_read_count');
}

sub _calculate_unaligned_read_count {
    my $self = shift;
    my $unaligned_read_count;
    if($self->instrument_data->is_paired_end) {
        $unaligned_read_count = $self->fwd_unaligned_read_count + $self->rev_unaligned_read_count;
    } else {
        $unaligned_read_count = $self->poorly_aligned_read_count;
    }
    return $unaligned_read_count;
}

sub unaligned_base_pair_count {
    my $self = shift;
    return $self->get_metric_value('unaligned_base_pair_count');
}

sub _calculate_unaligned_base_pair_count {
    my $self = shift;
    
    my $unaligned_base_pair_count;
    if($self->instrument_data->is_paired_end) {
        $unaligned_base_pair_count =
            $self->fwd_unaligned_read_count * $self->instrument_data->fwd_read_length +
            $self->rev_unaligned_read_count * $self->instrument_data->rev_read_length;
    } else {
        $unaligned_base_pair_count = $self->unaligned_read_count * $self->instrument_data->read_length;
    }
    return $unaligned_base_pair_count;
}

sub total_base_pair_count {
    my $self = shift;
    return $self->get_metric_value('total_base_pair_count');
}

sub _calculate_total_base_pair_count {
    my $self = shift;

    my $total_base_pair_count = $self->instrument_data->total_bases_read;
    return $total_base_pair_count;
}

sub _count_reads_in_fastq_files {
    my $self = shift;
    my @fastq_files = @_;
    
    my @wcs = `wc -l @fastq_files`;
    
    unless(@wcs) {
        $self->error_message('Unable to count reads in FASTQ files (No output from `wc` returned.): ' . join(', ', @fastq_files));
        return;
    }

    my $wc;
    if(scalar @wcs > 1) {
        ($wc) = grep { /total/ } @wcs;
        $wc =~ s/total//;
        $wc =~ s/\s//g;
    } else {
        ($wc) = $wcs[0] =~ /^(\d+)/; #no separate total if only one input, just grab count
    }

    if ($wc % 4) {
        $self->error_message("run has a line count of $wc, which is not divisible by four!");
        return;
    }
    my $total_reads_in_fastq_files = $wc/4; #four lines per read in FASTQ file format
    
    return $total_reads_in_fastq_files;
}

sub _count_lines_in_bam_file {
    my $self = shift;
    my %options = @_;
    
    my $alignment = $self->instrument_data_input->alignment_set;
    
    my $cmd = Genome::Model::Tools::Sam->path_for_samtools_version($alignment->samtools_version);
    $cmd .= ' view ';
    
    if($options{include_flags}) {
        $cmd .= '-f ' . $options{include_flags} . ' ';
    }
    if($options{exclude_flags}) {
        $cmd .= '-F ' . $options{exclude_flags} . ' ';
    }
    
    $cmd .= $alignment->alignment_file;
    $cmd .= ' | wc -l';
    
    my $result = `$cmd`;
    my ($wc) = $result =~ /^(\d+)/;
    
    unless($wc) {
        $self->error_message('Unable to count reads in BAM file (No output from `wc` returned.)');
        return;
    }
    
    return $wc;
}

sub execute {
    my $self = shift;

    return $self->_process_and_link_alignments_to_build('get_or_create');
}


sub shortcut {
    my $self = shift;

    return $self->_process_and_link_alignments_to_build('get');
}

sub _process_and_link_alignments_to_build {
    my $self = shift;
    my $mode = shift;

    $DB::single = $DB::stopper;
    
    unless (-d $self->build_directory) {
        die "Missing build directory???";
    } 

    $self->debug_message(
        "Build directory: " . $self->build_directory
    );

    my $instrument_data_input = $self->instrument_data_input;
   
    my $build = $self->build;
    my $model = $build->model;
    my $processing_profile = $model->processing_profile;

    $self->debug_message("Finding or generating alignments for " . $instrument_data_input->__display_name__);
    my @alignments;
    my @errors;

    my %segment_info;
    if (defined $self->instrument_data_segment_id) {
        $segment_info{instrument_data_segment_id} = $self->instrument_data_segment_id;
        $segment_info{instrument_data_segment_type} = $self->instrument_data_segment_type;
    }
    
    if ($mode eq 'get_or_create') {
        @alignments = $processing_profile->generate_results_for_instrument_data_input($instrument_data_input, %segment_info); 
        unless (@alignments) {
            $self->error_message("Error finding or generating alignments!:\n" .  join("\n",$instrument_data_input->error_message));
            push @errors, $self->error_message;
        }
    } elsif ($mode eq 'get') {
        @alignments = $processing_profile->results_for_instrument_data_input_with_lock($instrument_data_input, %segment_info);
        unless (@alignments) {
            return undef; 
        }
    } else {
        $self->error_message("process/link alignments mode unknown: $mode");
        die $self->error_message;
    }
    
    if (@errors) {
        $self->error_message(join("\n",@errors));
        return 0;
    }

    for my $alignment (@alignments) {
        my $link = $alignment->add_user(user => $build, label => 'uses');
        if ($link) {
            $self->debug_message("Linked alignment " . $alignment->id . " to the build");
        }
        else {
            $self->error_message(
                "Failed to link the build to the alignment " 
                . $alignment->__display_name__ 
                . "!"
            );
            # TODO: die, but not for now
        }
    }

    $self->debug_message("Generating alignments...");
    $self->generate_metric($self->metrics_for_class);

    $self->debug_message("Verifying...");
    unless ($self->verify_successful_completion) {
        $self->error_message("Error verifying completion!");
        return 0;
    }
    
    $self->debug_message("Complete!");
    return 1;
}

sub verify_successful_completion {
    my $self = shift;

    unless (-d $self->build_directory) {
    	$self->error_message("Build directory does not exist: " . $self->build_directory);
        return 0;
    }

    my $instrument_data_input = $self->instrument_data_input;
    my %segment_info;
    if (defined $self->instrument_data_segment_id) {
        $segment_info{instrument_data_segment_id} = $self->instrument_data_segment_id;
        $segment_info{instrument_data_segment_type} = $self->instrument_data_segment_type;
    }
    my @alignments = $self->model->processing_profile->results_for_instrument_data_input($instrument_data_input,%segment_info);
    
    my @errors;
    for my $alignment (@alignments) {
        unless ($alignment->verify_alignment_data) {
            $self->error_message('Failed to verify alignment data: '.  join("\n",$alignment->error_message));
            push @errors, $self->error_message;
        }
    }
    if (@errors) {
        $self->error_message(join("\n",@errors));
        return 0;
    }

    return 1;
}

# these were only ever used for maq, and are probably completely obsolete

sub contaminated_read_count {
    my $self = shift;
    return $self->get_metric_value('contaminated_read_count');
}

sub _calculate_contaminated_read_count {
    my $self = shift;
 
    my $instrument_data_input = $self->instrument_data_input;
    my ($alignment) = $self->build->alignment_results_for_instrument_data($instrument_data_input->value);
    my @f = $alignment->aligner_output_file_paths;
    @f = grep($_ !~ 'sanitized', @f);
    
    my $contaminated_read_count = 0;
    for my $f (@f) {
        my $fh = IO::File->new($f);
        $fh or die "Failed to open $f to read.  Error returning value for contaminated_read_count.\n";
        my $n;
        while (my $row = $fh->getline) {
            if ($row =~ /\[ma_trim_adapter\] (\d+) reads possibly contains adaptor contamination./) {
                $n = $1;
                last;
            }
        }
        unless (defined $n) {
            #$self->warning_message("No adaptor information found in $f!");
            next;
        }
        $contaminated_read_count += $n;
    }
    
    return $contaminated_read_count;
}


  
1;

