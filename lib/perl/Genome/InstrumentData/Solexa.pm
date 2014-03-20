package Genome::InstrumentData::Solexa;

use strict;
use warnings;

use Genome;
use File::Basename;

use List::MoreUtils qw(uniq);

class Genome::InstrumentData::Solexa {
    is => ['Genome::InstrumentData', 'Genome::Searchable'],
    has_constant => [
        sequencing_platform => { value => 'solexa' },
    ],
    has_optional => [
        project_name => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'project_name' ],
            is_mutable => 1,
        },
        target_region_set_name => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'target_region_set_name' ],
            is_mutable => 1,
        },
        flow_cell_id => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'flow_cell_id' ],
            is_mutable => 1,
        },
        # TODO Need to remove, depends on LIMS tables
        flow_cell => {
            is => 'Genome::InstrumentData::FlowCell',
            id_by => 'flow_cell_id',
        },
        lane => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'lane' ],
            is_mutable => 1,
        },
        read_length => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'read_length' ],
            is_mutable => 1,
        },
        old_filt_error_rate_avg => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'filt_error_rate_avg' ],
            is_mutable => 1,
        },
        old_rev_filt_error_rate_avg => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'rev_filt_error_rate_avg' ],
            is_mutable => 1,
        },
        old_fwd_filt_error_rate_avg => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'fwd_filt_error_rate_avg' ],
            is_mutable => 1,
        },
        rev_filt_aligned_clusters_pct => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'rev_filt_aligned_clusters_pct' ],
            is_mutable => 1,
        },
        fwd_filt_aligned_clusters_pct => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'fwd_filt_aligned_clusters_pct' ],
            is_mutable => 1,
        },
        rev_seq_id => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'rev_seq_id' ],
            is_mutable => 1,
        },
        fwd_seq_id => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'fwd_seq_id' ],
            is_mutable => 1,
        },
        rev_read_length => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'rev_read_length' ],
            is_mutable => 1,
        },
        fwd_read_length => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'fwd_read_length' ],
            is_mutable => 1,
        },
        rev_kilobases_read => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'rev_kilobases_read' ],
            is_mutable => 1,
        },
        fwd_kilobases_read => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'fwd_kilobases_read' ],
            is_mutable => 1,
        },
        rev_run_type => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'rev_run_type' ],
            is_mutable => 1,
        },
        fwd_run_type => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'fwd_run_type' ],
            is_mutable => 1,
        },
        run_type => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'run_type' ],
            is_mutable => 1,
        },
        gerald_directory => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'gerald_directory' ],
            is_mutable => 1,
        },
        old_median_insert_size => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'median_insert_size' ],
            is_mutable => 1,
        },
        old_sd_above_insert_size => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'sd_above_insert_size' ],
            is_mutable => 1,
        },
        old_sd_below_insert_size => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'sd_below_insert_size' ],
            is_mutable => 1,
        },
        is_external => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'is_external' ],
            is_mutable => 1,
        },
        archive_path => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'archive_path' ],
            is_mutable => 1,
        },
        bam_path => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'bam_path' ],
            is_mutable => 1,
        },
        adaptor_path => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'adaptor_path' ],
            is_mutable => 1,
        },
        rev_clusters => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'rev_clusters' ],
            is_mutable => 1,
        },
        fwd_clusters => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'fwd_clusters' ],
            is_mutable => 1,
        },
        clusters => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'clusters' ],
            is_mutable => 1,
        },
        read_count => {
            calculate => q| my $reads = $self->clusters; $reads *= 2 if $self->is_paired_end; return $reads; |,
        },
        analysis_software_version => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'analysis_software_version' ],
            is_mutable => 1,
        },
        index_sequence => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'index_sequence' ],
            is_mutable => 1,
        },
        gc_bias_path => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'gc_bias_path' ],
            is_mutable => 1,
        },
        fastqc_path => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'fastqc_path' ],
            is_mutable => 1,
        },

        #TODO These three columns will point to "read_length" or whatever name is decided
        #(see also https://gscweb.gsc.wustl.edu/wiki/Software_Development:Illumina_Indexed_Runs_Warehouse_Schema)
        _sls_read_length => { calculate => q| return $self->read_length + 1| },
        _sls_fwd_read_length => { calculate => q| return $self->fwd_read_length + 1| },
        _sls_rev_read_length => { calculate => q| return $self->rev_read_length + 1| },
        cycles => { calculate => q| return $self->read_length + 1| },

        short_name => {
            is => 'Text',
            calculate_from => ['run_name'],
            calculate => q|($run_name =~ /_([^_]+)$/)[0]|,
            doc => 'The essential portion of the run name which identifies the run.  The rest is redundant information',
        },
        is_paired_end => {
            calculate_from => ['run_type'],
            calculate => q|
                if (defined($run_type) and $run_type =~ m/^Paired$/) {
                    return 1;
                }
                return 0;
            |,
        },
        # TODO Replace with new Genome::Project
        project => {
            is => "Genome::Site::TGI::Project",
            calculate => q|Genome::Site::TGI::Project->get(name => $self->project_name)|
        },

        # TODO: move into the base class (we want to eliminate this subclass and not have the separate Imported subclass
        median_insert_size => {
            calculate_from => ['old_median_insert_size'],
            calculate => q|$old_median_insert_size or $self->get_default_alignment_metrics("median_insert_size")|,
        },
        sd_insert_size => {
            calculate => q|$self->get_default_alignment_metrics("sd_insert_size")|,
        },
        read_1_pct_mismatch => {
            calculate => q|$self->get_default_alignment_metrics("read_1_pct_mismatch")|,
        },
        read_2_pct_mismatch => {
            calculate => q|$self->get_default_alignment_metrics("read_2_pct_mismatch")|,
        },

        # for backward-compatibility
        sd_above_insert_size => {
            calculate_from => ['old_sd_above_insert_size','sd_insert_size'],
            calculate => q|$old_sd_above_insert_size or $sd_insert_size|
        },
        sd_below_insert_size => {
            calculate_from => ['old_sd_below_insert_size','sd_insert_size'],
            calculate => q|$old_sd_below_insert_size or $sd_insert_size|
        },
        fwd_filt_error_rate_avg => {
            calculate_from => ['old_fwd_filt_error_rate_avg','read_1_pct_mismatch'],
            calculate => q|$old_fwd_filt_error_rate_avg or $read_1_pct_mismatch|
        },
        rev_filt_error_rate_avg => {
            calculate_from => ['old_rev_filt_error_rate_avg','read_2_pct_mismatch'],
            calculate => q|$old_rev_filt_error_rate_avg or $read_2_pct_mismatch|
        },
        filt_error_rate_avg => {
            calculate_from => ['old_filt_error_rate_avg','read_1_pct_mismatch','read_2_pct_mismatch'],
            calculate => q|
                if ($old_filt_error_rate_avg) {
                    $old_filt_error_rate_avg;
                }
                elsif ($read_1_pct_mismatch and $read_2_pct_mismatch) {
                    ($read_1_pct_mismatch + $read_2_pct_mismatch)/2;
                }
                elsif (not $read_1_pct_mismatch) {
                    $read_2_pct_mismatch;
                }
                else {
                    $read_1_pct_mismatch;
                }
            |
        },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return if not $self;
    # some properties previously only on imported data are set here
    for my $property_name (qw/import_format/) {
        my $resolver_name = 'resolve_' . $property_name;
        my $value = eval { $self->$resolver_name };
        if ($@) {
            $self->warning_message("Failed to resolve import_format for " . $self->__display_name__ . ": $@");
            next;
        }
        $self->$property_name($value);
    }
    return $self;
}

sub __display_name__ {
    my $self = $_[0];
    my $sample = $self->sample;
    my($flow_cell_id, $subset_name, $id, $sample_display_name)
            = map { defined($_) ? $_ : '' }
                    ( $self->flow_cell_id, $self->subset_name, $self->id, $sample && $sample->__display_name__);
    return $flow_cell_id . '/' . $subset_name . " (" . $id . ") for " . $sample_display_name;
}

sub _calculate_paired_end_kb_usage {
    my $self = shift;
    my $HEADER_LENGTH = shift;
    $HEADER_LENGTH = $HEADER_LENGTH + 5; # Adding 5 accounts for newlines in the FQ file.
    # If data is paired_end fwd_read_length, rev_read_length, fwd_clusters and rev_clusters
    # should all be defined.
    if (!defined($self->fwd_read_length) or $self->fwd_read_length <= 0) {
        $self->error_message("Instrument data fwd_read_length is either undefined or less than 0.");
        die;
    } elsif (!defined($self->rev_read_length) or $self->rev_read_length <= 0) {
        $self->error_message("Instrument data rev_read_length is either undefined or less than 0.");
        die;
    } elsif (!defined($self->fwd_clusters) or $self->fwd_clusters <= 0) {
        $self->error_message("Instrument data fwd_clusters is either undefined or less than 0.");
        die;
    } elsif (!defined($self->rev_clusters) or $self->rev_clusters <= 0) {
        $self->error_message("Instrument data rev_clusters is either undefined or less than 0.");
        die;
    }

    my $fwd = (($self->fwd_read_length + $HEADER_LENGTH) * $self->fwd_clusters)*2;
    my $rev = (($self->rev_read_length + $HEADER_LENGTH) * $self->rev_clusters)*2;
    my $total = ($fwd + $rev) / 1024.0;
    return $total;
}

sub _calculate_non_paired_end_kb_usage {
    my $self = shift;
    my $HEADER_LENGTH = shift;
    $HEADER_LENGTH = $HEADER_LENGTH + 5; # adding 5 accounts for newlines in the FQ file.
    # We will take the max of fwd_read_length or rev_read_length and the max of fwd_clusters and rev_clusters
    # to make sure that even strange data won't cause overly low space allocation
    # Get the max read length.
    my $max_read_length;
    if ( defined($self->read_length) and $self->read_length > 0 ) {
        $max_read_length = $self->read_length;
    } elsif ( defined($self->fwd_read_length) and defined($self->rev_read_length) ) {
        if ( $self->fwd_read_length > $self->rev_read_length ) {
            $max_read_length = $self->fwd_read_length;
        } else {
            $max_read_length = $self->rev_read_length;
        }
    } elsif ( defined($self->fwd_read_length) and $self->fwd_read_length > 0) {
        $max_read_length = $self->fwd_read_length;
    } elsif ( defined($self->rev_read_length) and $self->rev_read_length > 0) {
        $max_read_length = $self->rev_read_length;
    } else {
        $self->error_message("No valid read length value found in instrument data");
        die;
    }
    # Get the max cluster length.
    my $max_clusters;
    if ( defined($self->clusters) and $self->clusters > 0 ) {
        $max_clusters = $self->clusters;
    } elsif ( defined($self->fwd_clusters) and defined($self->rev_clusters) ) {
        if ( $self->fwd_clusters > $self->rev_clusters ) {
            $max_clusters = $self->fwd_clusters;
        } else {
            $max_clusters = $self->rev_clusters;
        }
    } elsif ( defined($self->fwd_clusters) and $self->fwd_clusters > 0) {
        $max_clusters = $self->fwd_clusters;
    } elsif ( defined($self->rev_clusters) and $self->rev_clusters > 0) {
        $max_clusters = $self->rev_clusters;
    } else {
        $self->error_message("No valid number of clusters value found in instrument data");
        die;
    }

    my $total_b = (($max_read_length + $HEADER_LENGTH) * $max_clusters)*2;
    my $total = $total_b / 1024.0;
    return $total;
}

sub calculate_alignment_estimated_kb_usage {
    my $self = shift;
    # Different aligners will require different levels of overhead, so this should return
    # approximately the total size of the instrument data files.
    # In a FQ file, every read means 4 lines, 2 of length $self->read_length,
    # and 2 of read identifier data. Therefore we can expect the total file size to be
    # about (($self->read_length + header_size) * $self->fwd_clusters)*2 / 1024

    # This is length of the id tag in the FQ file, it looks something like: @HWUSI-EAS712_6171L:2:1:1817:12676#TGACCC/1
    my $HEADER_LENGTH = 45;

    if ($self->is_paired_end) {
        return int $self->_calculate_paired_end_kb_usage($HEADER_LENGTH);
    } else {
        return int $self->_calculate_non_paired_end_kb_usage($HEADER_LENGTH);
    }
}

#< Dump to File System >#
sub dump_illumina_fastq_files {
    my $self = shift;

    unless ($self->resolve_quality_converter eq 'sol2phred') {
        $self->error_message("This instrument data is not natively Illumina formatted, cannot dump");
        die $self->error_message;
    }

    return $self->_unprocessed_fastq_filenames(@_);
}

sub dump_solexa_fastq_files {
    my $self = shift;

    unless ($self->resolve_quality_converter eq 'sol2sanger') {
        $self->error_message("This instrument data is not natively Solexa formatted, cannot dump");
        die $self->error_message;
    }

    return $self->_unprocessed_fastq_filenames(@_);
}

sub dump_sanger_fastq_files {
    my $self = shift;

    if (defined $self->bam_path && -s $self->bam_path) {
       $self->debug_message("Now using a bam instead");
       return $self->dump_fastqs_from_bam(@_);
    }

    my @illumina_fastq_pathnames = $self->_unprocessed_fastq_filenames(@_);

    my %params = @_;

    my $requested_directory = delete $params{directory} || Genome::Sys->base_temp_directory;

    my @converted_pathnames;
    my $counter = 0;
    for my $illumina_fastq_pathname (@illumina_fastq_pathnames) {
        my $converted_fastq_pathname;
        if ($self->resolve_quality_converter eq 'sol2sanger') {
            $converted_fastq_pathname = $requested_directory . '/' . $self->id . '-sanger-fastq-'. $counter . ".fastq";
            $self->debug_message("Applying sol2sanger quality conversion.  Converting to $converted_fastq_pathname");
            unless (Genome::Model::Tools::Fastq::Sol2sanger->execute(
                                                                    fastq_file => $illumina_fastq_pathname,
                                                                    sanger_fastq_file => $converted_fastq_pathname)) {
                $self->error_message('Failed to execute sol2sanger quality conversion $illumina_fastq_pathname $converted_fastq_pathname.');
                die($self->error_message);
            }
        } elsif ($self->resolve_quality_converter eq 'sol2phred') {
            $converted_fastq_pathname = $requested_directory . '/' . $self->id . '-sanger-fastq-'. $counter . ".fastq";
            $self->debug_message("Applying sol2phred quality conversion.  Converting to $converted_fastq_pathname");

            unless (Genome::Model::Tools::Fastq::Sol2phred->execute(fastq_file => $illumina_fastq_pathname,
                                                                    phred_fastq_file => $converted_fastq_pathname)) {
                $self->error_message('Failed to execute sol2phred quality conversion.');
                die($self->error_message);
            }
        } elsif ($self->resolve_quality_converter eq 'none') {
            $self->debug_message("No quality conversion required.");
            $converted_fastq_pathname = $illumina_fastq_pathname;
        } else {
            $self->error_message("Undefined quality converter requested, I can't proceed");
            die $self->error_message;
        }
        unless (-e $converted_fastq_pathname && -f $converted_fastq_pathname && -s $converted_fastq_pathname) {
            $self->error_message('Failed to validate the conversion of solexa fastq file '. $illumina_fastq_pathname .' to sanger quality scores');
            die($self->error_message);
        }
        $counter++;

        if (($converted_fastq_pathname ne $illumina_fastq_pathname ) &&
            ($illumina_fastq_pathname =~ m/\/tmp\//)) {

            $self->debug_message("Removing original unconverted file from temp space to save disk space:  $illumina_fastq_pathname");
            unlink $illumina_fastq_pathname;
        }
        push @converted_pathnames, $converted_fastq_pathname;
    }

    return @converted_pathnames;

}

sub dump_trimmed_fastq_files {
    my $self = shift;
    my %params = @_;

    my $data_directory = $params{directory} || Genome::Sys->create_temp_directory;

    my $segment_params = $params{segment_params};
    my $discard_fragments = $params{discard_fragments} || 0;

    my $trimmer_name = $params{trimmer_name};
    my $trimmer_version = $params{trimmer_version};
    my $trimmer_params = $params{trimmer_params};

    unless($trimmer_name) {
        return $self->dump_sanger_fastq_files(%$segment_params, directory => $data_directory, discard_fragments => $discard_fragments);
    }

    $self->debug_message('Trimmer name: '.$trimmer_name);
    $self->debug_message('Trimmer version: '.($trimmer_version // 'NULL'));
    $self->debug_message('Trimmer params: '.($trimmer_params // 'NULL'));
    my @sx_cmd_parts = $self->_convert_trimmer_to_sx_commands(
        trimmer_name => $params{trimmer_name},
        trimmer_version => $params{trimmer_version},
        trimmer_params => $params{trimmer_params},

    );
    if ( @sx_cmd_parts ) {
        $self->debug_message('SX processing detected!');

        # Inputs
        my @fastq_pathnames = $self->dump_sanger_fastq_files(%$segment_params); # dies on error
        $sx_cmd_parts[0] .= ' --input '.join(',', map { 'file='. $_ .':type=sanger' } @fastq_pathnames);
        $self->debug_message('SX inputs: '.join(' ', @fastq_pathnames));

        # Outputs and file names
        my (@trimmed_fastq_pathnames, @output);
        if ($self->is_paired_end) {
            @trimmed_fastq_pathnames =
                map { $data_directory . '/trimmed-sanger-fastq-' . $_ . '.fastq' }
                    ('read1','read2','fragment');
            @output = (
                $trimmed_fastq_pathnames[0].':name=fwd',
                $trimmed_fastq_pathnames[1].':name=rev',
                $trimmed_fastq_pathnames[2].':name=sing',
            );

        }
        else {
            @trimmed_fastq_pathnames =
                map { $data_directory . '/trimmed-sanger-fastq-' . $_ .'.fastq'}
                    ('fragment');
            @output = $trimmed_fastq_pathnames[0];
        }
        $sx_cmd_parts[$#sx_cmd_parts] .= ' --output '.join(',', @output);
        $self->debug_message('SX outputs: '.join(' ', @trimmed_fastq_pathnames));

        # Run SX
        # Assemble command
        my $sx_cmd = join(' | ', map { 'gmt sx '.$_ } @sx_cmd_parts);
        $self->debug_message('SX command: '.$sx_cmd);
        # Run
        my $cmd_ran_ok = eval{ Genome::Sys->shellcmd(cmd => $sx_cmd); };
        if ( not $cmd_ran_ok ) {
            die $self->error_message('Failed to run SX trim command! '.$sx_cmd);
        }
        $self->debug_message('Run SX command...OK');

        # Remove original dumped fastqs
        for my $input_fastq_pathname (@fastq_pathnames) {
            if ($input_fastq_pathname =~ m/^\/tmp/) {
                $self->debug_message("Removing original file from after trimming to save space: $input_fastq_pathname");
                unlink($input_fastq_pathname);
            }
        }

        # in paired end trimming, only return trimmed files with reads, check for errors
        my @paths;
        if (@trimmed_fastq_pathnames == 3){
            my $paired_with_size = grep { -s $_ } @trimmed_fastq_pathnames[0,1];
            if ($paired_with_size == 0){
                $self->debug_message("paired end trimmed files have no size, skipping");
            }elsif($paired_with_size == 1){
                die $self->error_message("only one trimmed pair file with size, trimming produced bad result!");
            }else{
                push @paths, @trimmed_fastq_pathnames[0,1];
            }
            if (-s $trimmed_fastq_pathnames[2]){
                push @paths, $trimmed_fastq_pathnames[2];
            }else{
                $self->debug_message("fragment trimmed file has no size, skipping");
            }
        } else {
            @paths = @trimmed_fastq_pathnames;
        }
        return @paths;
    }
    # if the above did not work, we have a legacy trimmer.

    # DO __NOT__ ADD TO THE CONDITIONAL LOGIC HERE
    # MAKE A TRIMMER IN THE SX API, AND FALL THROUGH TO THE "ELSE" BLOCK
    # EVENTUALLY, ALL OF THIS IF STATMENT NEEDS TO GO AWAY -SSMITH

    $self->debug_message('Preceding w/ non SX processing');
    my @fastq_pathnames = $self->dump_sanger_fastq_files(%$segment_params); # dies on error
    my @trimmed_fastq_pathnames;
    #if the trimmer supports paired end, we just run it once, otherwise we need to loop over the fastqs
    if(@fastq_pathnames == 2 && $trimmer_name eq 'far' && $trimmer_version >= '2.0') {
        my $trimmed_input_fastq_path = $data_directory . '/trimmed-sanger-fastq';
        my $trimmer = Genome::Model::Tools::Far::Trimmer->create(
            params => $trimmer_params,
            use_version => $trimmer_version,
            source => $fastq_pathnames[0],
            source2 => $fastq_pathnames[1],
            target => $trimmed_input_fastq_path,
            far_output => $data_directory .'/far_output_report.txt',
        );
        unless ($trimmer) {
            $self->error_message('Failed to create fastq trim command');
            die($self->error_message);
         }
        unless ($trimmer->execute) {
            $self->error_message('Failed to execute fastq trim command '. $trimmer->command_name);
            die($self->error_message);
         }
        push @trimmed_fastq_pathnames, grep {$_ !~ /single/} glob "$trimmed_input_fastq_path*fastq";
        unless (@trimmed_fastq_pathnames){
            die $self->error_message("Failed to get expected trimmed output files");
        }
    }
    else {
        my $counter = 0;
        for my $input_fastq_pathname (@fastq_pathnames) {
            if($trimmer_name eq 'trimq2_shortfilter') {
                $self->error_message('Trimmer ' . $trimmer_name . ' is not currently supported by this module.');
                die $self->error_message;
             }

            my $trimmed_input_fastq_pathname = $data_directory . '/trimmed-sanger-fastq-' . $counter;
            my $trimmer;
            if ($trimmer_name eq 'fastx_clipper') {
                $trimmer = Genome::Model::Tools::Fastx::Clipper->create(
                    params => $trimmer_params,
                    use_version => $trimmer_version,
                    input_file => $input_fastq_pathname,
                    output_file => $trimmed_input_fastq_pathname,
                );
            } elsif ($trimmer_name eq 'far') {
                $trimmer = Genome::Model::Tools::Far::Trimmer->create(
                    params => $trimmer_params,
                    use_version => $trimmer_version,
                    source => $input_fastq_pathname,
                    target => $trimmed_input_fastq_pathname,
                    far_output => $data_directory .'/far_output_report.txt',
                )
            }
            elsif ($trimmer_name eq 'trim5') {
                $trimmer = Genome::Model::Tools::Fastq::Trim5->create(
                    length => $trimmer_params,
                    input => $input_fastq_pathname,
                    output => $trimmed_input_fastq_pathname,
                );
            }
            elsif ($trimmer_name eq 'bwa_style') {
                my ($trim_qual) = $trimmer_params =~ /--trim-qual-level\s*=?\s*(\S+)/;
                $trimmer = Genome::Model::Tools::Fastq::TrimBwaStyle->create(
                    trim_qual_level => $trim_qual,
                    fastq_file      => $input_fastq_pathname,
                    out_file        => $trimmed_input_fastq_pathname,
                    qual_type       => 'sanger',  #hardcoded for now
                    report_file     => $data_directory.'/trim_bwa_style.report.'.$counter,
                );
            }
            elsif ($trimmer_name =~ /trimq2_(\S+)/) {
                #This is for trimq2 no_filter style
                #move trimq2.report to alignment directory

                my %trimq2_params = (
                    fastq_file  => $input_fastq_pathname,
                    out_file    => $trimmed_input_fastq_pathname,
                    report_file => $data_directory.'/trimq2.report.'.$counter,
                    trim_style  => $1,
                );
                my ($qual_level, $string) = $self->_get_trimq2_params($trimmer_params);

                my ($primer_sequence) = $trimmer_params =~ /--primer-sequence\s*=?\s*(\S+)/;
                $trimq2_params{trim_qual_level} = $qual_level if $qual_level;
                $trimq2_params{trim_string}     = $string if $string;
                $trimq2_params{primer_sequence} = $primer_sequence if $primer_sequence;
                $trimq2_params{primer_report_file} = $data_directory.'/trim_primer.report.'.$counter if $primer_sequence;

                $trimmer = Genome::Model::Tools::Fastq::Trimq2::Simple->create(%trimq2_params);
            }
            elsif ($trimmer_name eq 'random_subset') {
                my $seed_phrase = $self->run_name .'_'. $self->id;
                $trimmer = Genome::Model::Tools::Fastq::RandomSubset->create(
                    input_read_1_fastq_files => [$input_fastq_pathname],
                    output_read_1_fastq_file => $trimmed_input_fastq_pathname,
                    limit_type => 'reads',
                    limit_value => $trimmer_params,
                    seed_phrase => $seed_phrase,
                );
            }
            elsif ($trimmer_name eq 'normalize') {
                my ($read_length,$reads) = split(':',$trimmer_params);
                my $trim = Genome::Model::Tools::Fastq::Trim->execute(
                    read_length => $read_length,
                    orientation => 3,
                    input => $input_fastq_pathname,
                    output => $trimmed_input_fastq_pathname,
                );
                unless ($trim) {
                    die('Failed to trim reads using test_trim_and_random_subset');
                }
                my $random_input_fastq_pathname = $data_directory . '/random-sanger-fastq-' . $counter;
                $trimmer = Genome::Model::Tools::Fastq::RandomSubset->create(
                    input_read_1_fastq_files => [$trimmed_input_fastq_pathname],
                    output_read_1_fastq_file => $random_input_fastq_pathname,
                    limit_type  => 'reads',
                    limit_value => $reads,
                    seed_phrase => $self->run_name .'_'. $self->id,
                );
                $trimmed_input_fastq_pathname = $random_input_fastq_pathname;
            }
            else {
                $self->error_message(
                    sprintf(
                        "Unknown read trimmer_name %s.",
                        $trimmer_name,
                    )
                );
                die($self->error_message);
            }

            # error check and runn the legacy trimmer

            unless ($trimmer) {
                $self->error_message('Failed to create fastq trim command');
                die($self->error_message);
            }

            unless ($trimmer->execute) {
                $self->error_message('Failed to execute fastq trim command '. $trimmer->command_name);
                die($self->error_message);
            }

            # this legacy trimmer has post-execute work (which should have been in the module)
            if ($trimmer_name eq 'normalize') {
                my @empty = ();
                $trimmer->_index(\@empty);
            }

            if ($input_fastq_pathname =~ m/^\/tmp/) {
                $self->debug_message("Removing original file from before trimming to save space: $input_fastq_pathname");
                unlink($input_fastq_pathname);
            }

            push @trimmed_fastq_pathnames, $trimmed_input_fastq_pathname;
            $counter++;
        }
    }
    return @trimmed_fastq_pathnames;
}

sub _convert_trimmer_to_sx_commands {
    # return array of sx cmds
    # return undef if not sx
    # confess on error
    my ($self, %trimmer) = @_;

    my $trimmer_name = delete $trimmer{trimmer_name};
    Carp::confess('No trimmer name provided to convert trimmer to SX command parts!') if not $trimmer_name;
    my $trimmer_version = delete $trimmer{trimmer_version};
    my $trimmer_params = delete $trimmer{trimmer_params};
    Carp::confess('Unknown params sent to convert trimmer to sx commands! '.Data::Dumper::Dumper(\%trimmer)) if %trimmer;

    # New SX! Support running multiple SX functons as a command.
    my $is_sx_cmd = $trimmer_name =~ s/^gmt sx //;
    if ( $is_sx_cmd ) {
        Carp::confess('Cannot have trimmer version w/ SX processing. Include the version in the sx command inside the trimmer name.') if $trimmer_version;
        Carp::confess('Cannot have trimmer params w/ SX processing. Include the params in the sx command inside the trimmer name.') if $trimmer_params;
        return split(/\s*\|\s*/, $trimmer_name)
    }

    # Old SX way. Only supports one trimmer.
    my $sx_class_meta = eval{
        return if not $trimmer_name;
        my $class_name = 'Genome::Model::Tools::Sx::Trim';
        my @words = split(' ', $trimmer_name);
        for my $word (@words) {
            my @parts = map { ucfirst($_) } split('-',$word);
            $class_name .= "::" . join('',@parts);
        }
        $class_name->__meta__;
    };

    # Not SX - OK
    return if not $sx_class_meta;

    # Convert the trimmer params to strings params
    my @params;
    if ( $trimmer_params ) {
        if($trimmer_params =~ '=>') {
            @params = eval("no strict; no warnings; $trimmer_params");
            if ( not @params ) {
                $self->error_message("Invalid params ($trimmer_params) to convert to command line params! $@");
                return;
            }
        } elsif ( $trimmer_params =~ /--\S+?[= ]\S+/ ) {
            while($trimmer_params =~ /--(\S+)?[= ](\S+)/g) {
                my ($key, $value) = ($1, $2);
                #for compatibility with pre-SX far processing profiles
                if($trimmer_name eq 'far') {
                    next if $key eq 'format';
                }

                push @params, $key => $value;
            }
        } else {
            Carp::confess("Unknown params ($trimmer_params) to convert trimmer to sx command!");
        }
    }

    # Add version to the params
    push @params, 'version', $trimmer_version if defined $trimmer_version;

    # Add -- to the odd, and quotes to the even, handle logicals
    my $sx_cmd = 'trim '.$trimmer_name;
    for ( my $i = 0; $i < @params; $i += 2 ) {
        my $property_name = $params[$i];
        $property_name =~ s/\-/_/g;
        my $property_meta = $sx_class_meta->property_meta_for_name($property_name);
        if ( not $property_meta ) {
            my $starts_with_no_property_name = $property_name;
            $starts_with_no_property_name =~ s/^no//;
            $property_meta = $sx_class_meta->property_meta_for_name($starts_with_no_property_name);
        }
        Carp::confess("Failed to get property '$property_name' from sx class ".$sx_class_meta->class_name) if not $property_meta;
        my $param_name = $params[$i];
        $param_name =~ s/_/-/g;
        my $value = $params[$i + 1];
        if ( $property_meta->data_type eq 'Boolean' ) {
            $sx_cmd .= ' --'.$param_name if $value;
        }
        else {
            $sx_cmd .= ' --'.$param_name." '".$value."'";
        }
    }

    # Return the sx cmd
    return $sx_cmd;
}

sub _get_trimq2_params {
    my $self = shift;
    my $trimmer_params = shift;

    #for trimq2_shortfilter, input something like "32:#" (length:string) in processing profile as trimmer_params, for trimq2_smart1, input "20:#" (quality_level:string)
    if ($trimmer_params !~ /:/){
        $trimmer_params = '::';
    }
    my ($first_param, $string) = split /\:/, $trimmer_params;

    return ($first_param, $string);
}

sub _unprocessed_fastq_filenames {
    my $self = shift;
    my @fastqs;
    if ($self->is_external) {
        @fastqs = $self->resolve_external_fastq_filenames(@_);
    } else {
        @fastqs = @{$self->resolve_fastq_filenames(@_)};
    }
    return @fastqs;
}

sub desc {
    my $self = shift;
    return $self->full_name .'('. $self->id .')';
}

sub read1_fastq_name {
    my $self = shift;
    my $lane = $self->lane;

    return "s_${lane}_1_sequence.txt";
}

sub read2_fastq_name {
    my $self = shift;
    my $lane = $self->lane;

    return "s_${lane}_2_sequence.txt";
}

sub fragment_fastq_name {
    my $self = shift;
    my $lane = $self->lane;

    return "s_${lane}_sequence.txt";
}

sub resolve_fastq_filenames {
    my $self = shift;
    my $lane = $self->lane;
    my $desc = $self->desc;

    my %params = @_;
    my $paired_end_as_fragment = delete $params{'paired_end_as_fragment'};
    my $requested_directory = delete $params{'directory'} || Genome::Sys->base_temp_directory;

    my @illumina_output_paths;
    my @errors;

    # First check the archive directory and second get the gerald directory
    for my $dir_type (qw(archive_path gerald_directory)) {
        $self->debug_message("Now trying to get fastq from $dir_type for $desc");

        my $directory = $self->$dir_type;
        $directory = $self->validate_fastq_directory($directory, $dir_type);
        next unless $directory;

        if ($dir_type eq 'archive_path') {
            $directory = $self->dump_illumina_fastq_archive($requested_directory); #need eval{} here ?
            $directory = $self->validate_fastq_directory($directory, 'dump_fastq_dir');
            next unless $directory;
        }

        eval {
            #handle fragment or paired-end data
            if ($self->is_paired_end) {
                if (!$paired_end_as_fragment || $paired_end_as_fragment == 1) {
                    if (-e "$directory/" . $self->read1_fastq_name) {
                        push @illumina_output_paths, "$directory/" . $self->read1_fastq_name;
                    }
                    elsif (-e "$directory/Temp/" . $self->read1_fastq_name) {
                        push @illumina_output_paths, "$directory/Temp/" . $self->read1_fastq_name;
                    }
                    else {
                        die "No illumina forward data in directory for lane $lane! $directory";
                    }
                }
                if (!$paired_end_as_fragment || $paired_end_as_fragment == 2) {
                    if (-e "$directory/" . $self->read2_fastq_name) {
                        push @illumina_output_paths, "$directory/" . $self->read2_fastq_name;
                    }
                    elsif (-e "$directory/Temp/" . $self->read2_fastq_name) {
                        push @illumina_output_paths, "$directory/Temp/" . $self->read2_fastq_name;
                    }
                    else {
                        die "No illumina reverse data in directory for lane $lane! $directory";
                    }
                }
            }
            else {
                if (-e "$directory/" . $self->fragment_fastq_name) {
                    push @illumina_output_paths, "$directory/" . $self->fragment_fastq_name;
                }
                elsif (-e "$directory/Temp/" . $self->fragment_fastq_name) {
                    push @illumina_output_paths, "$directory/Temp/" . $self->fragment_fastq_name;
                }
                else {
                    die "No fragment illumina data in directory for lane $lane! $directory";
                }
            }
        };

        push @errors, $@ if $@;
        last if @illumina_output_paths;
    }
    unless (@illumina_output_paths) {
        $self->error_message("No fastq files were found for $desc");
        $self->error_message(join("\n",@errors)) if @errors;
        die $self->error_message;
    }
    return \@illumina_output_paths;
}


sub dump_illumina_fastq_archive {
    my ($self, $dir) = @_;

    my $archive = $self->archive_path;
    $dir = Genome::Sys->base_temp_directory unless $dir;

    #Prevent unarchiving multiple times during execution
    #Hopefully nobody passes in a $dir expecting to overwrite another set of FASTQs coincidentally from the same lane number
    my $already_dumped = 0;

    if($self->is_paired_end) {
        if (-s $dir . '/' . $self->read1_fastq_name and -s $dir . '/' . $self->read2_fastq_name) {
            $already_dumped = 1;
        }
    }
    else {
        if (-s $dir . '/' . $self->fragment_fastq_name) {
            $already_dumped = 1;
        }
    }

    unless($already_dumped) {
        my $cmd = "tar -xzf $archive --directory=$dir";
        unless (Genome::Sys->shellcmd(
            cmd => $cmd,
            input_files => [$archive],
        )) {
            $self->error_message('Failed to run tar command '. $cmd);
            return;
            #die($self->error_message); Should try to get fastq from gerald_directory instead of dying
        }
    }
    return $dir;
}

sub validate_fastq_directory {
    my ($self, $dir, $dir_type) = @_;

    my $msg_base = "$dir_type : $dir";

    unless ($dir) {
        $self->error_message("$msg_base is null");
        return;
    }

    unless (-e $dir) {
        $self->error_message("$msg_base not existing in file system");
        return;
    }

    unless ($dir_type eq 'archive_path') {
        my @files = glob("$dir/*"); #In scalar context, a glob functions as an iterator--we instead want to check the number of files
        unless (scalar @files) {
            $self->error_message("$msg_base is empty");
            return;
        }
    }

    return $dir;
}


sub resolve_external_fastq_filenames {
    my $self = shift;

    my @fastq_pathnames;
    my $fastq_pathname = $self->create_temp_file_path('fastq');
    unless ($fastq_pathname) {
        die "Failed to create temp file for fastq!";
    }
    return ($fastq_pathname);
}

sub _calculate_total_read_count {
    my $self = shift;

    if($self->is_external) {
        my $data_path_object = Genome::MiscAttribute->get(entity_id => $self->id, property_name=>'full_path');
        my $data_path = $data_path_object->value;
        my $lines = `wc -l $data_path`;
        return $lines/4;
    }
    if ($self->clusters <= 0) {
        die('Impossible value '. $self->clusters .' for clusters field for solexa lane '. $self->id);
    }

    return $self->clusters;
}

sub native_qual_type {

    # old stuff needed sol2sanger, new stuff all uses sol2phred, but
    # we dont care what the version is anymore

    my $self = shift;

    my %analysis_software_versions = (
        'GAPipeline-0.3.0'       => 'illumina',
        'GAPipeline-0.3.0b1'     => 'illumina',
        'GAPipeline-0.3.0b2'     => 'illumina',
        'GAPipeline-0.3.0b3'     => 'illumina',
        'GAPipeline-1.0'         => 'illumina',
        'GAPipeline-1.0-64'      => 'illumina',
        'GAPipeline-1.0rc4'      => 'illumina',
        'GAPipeline-1.1rc1p4'    => 'illumina',
        'SolexaPipeline-0.2.2.5' => 'illumina',
        'SolexaPipeline-0.2.2.6' => 'illumina',
    );

    my $analysis_software_version = $self->analysis_software_version;
    unless ($analysis_software_version) {
        die('No analysis_software_version found for instrument data '. $self->id);
    }
    return $analysis_software_versions{$analysis_software_version} || 'sanger';
}

sub resolve_quality_converter {
    my $self = shift;
    if ($self->native_qual_type eq 'illumina') {
        return "sol2sanger";
    }
    return "sol2phred";
}

sub resolve_import_format {
    my $self = shift;
    if ($self->bam_path) {
        return 'bam'
    }
    elsif ($self->native_qual_type eq 'solexa') {
        return 'solexa fastq'
    }
    elsif ($self->native_qual_type eq 'illumina') {
        return 'illumina fastq'
    }
    elsif ($self->native_qual_type eq 'sanger') {
        return 'sanger fastq'
    }
    else {
        die "cannot resolve import format: native qual type is " . $self->native_qual_type;
    }
}

sub NEW_resolve_quality_converter {
    # TODO: this is the correct logic but it breaks tests,
    # implying our tests are wrong.  Investigate and fix.
    # This mirrors the ::Imported module which does the work correctly.
    # We almost never use this method in production because we use BAMs
    # for reads which were always sanger format.
    my $self = shift;
    if ($self->native_qual_type eq 'solexa') {
        return "sol2sanger";
    }
    elsif ($self->native_qual_type eq 'illumina') {
        return "sol2phred";
    }
    elsif($self->native_qual_type eq 'sanger') {
        return 'none';
    }
    else {
        Carp::confess("unresolved errors with handling non-sanger fastq files: FIXME");
    }
}

sub resolve_adaptor_file {
    my $self = shift;

    #these are constants and should probably be defined in class properties...TODO
    my $dna_primer_file = '/gscmnt/sata114/info/medseq/adaptor_sequences/solexa_adaptor_pcr_primer';
    my $rna_primer_file = '/gscmnt/sata114/info/medseq/adaptor_sequences/solexa_adaptor_pcr_primer_SMART';

    my $adaptor_file;
    if ( $self->sample_type eq 'rna' ) {
        $adaptor_file = $rna_primer_file;
    }
    else {
        $adaptor_file = $dna_primer_file;
    }
    unless (-f $adaptor_file) {
        $self->error_message('Specified adaptor file'. $adaptor_file .' does not exist.');
        die($self->error_message);
    }
    return $adaptor_file;
}

sub create_mock {
    my $class = shift;
    my $self = $class->SUPER::create_mock(@_);
    return unless $self;

    for my $method (qw/
        dump_sanger_fastq_files
        resolve_fastq_filenames
        _calculate_total_read_count
        resolve_adaptor_file
        run_identifier
        dump_fastqs_from_bam
    /) {
        my $ref = $class->can($method);
        die "Unknown method $method on " . $class . ".  Cannot make a pass-through for mock object!" unless $ref;
        $self->mock($method,$ref);
    }

    if (not defined $self->subset_name) {
        if ($self->index_sequence) {
            $self->subset_name($self->lane . '-' . $self->index_sequence);
        }
        else {
            $self->subset_name($self->lane);
        }
    }

    return $self;
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    my $expected;
    if (defined($self->lane) and defined($self->index_sequence)) {
        $expected = $self->lane . '-' . $self->index_sequence;
    }
    elsif (defined($self->lane)) {
        $expected = $self->lane;
    }
    else {
        $expected = '';
    }
    my $subset_name = $self->subset_name;
    $subset_name    = '' unless (defined $subset_name);
    if ($subset_name ne $expected) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => ['subset_name'],
            desc => 'Subset name for Illumina NGS instrument data should be "lane-indexsequence" or just "lane" for unindexed data. '
                    . 'Expected "' . $expected . '" got "' . $subset_name . '".'
        );
    }
    return @errors;
}

sub run_start_date_formatted {
    my $self = shift;

    my ($y, $m, $d) = $self->run_name =~ m/^(\d{2})(\d{2})(\d{2})_.*$/;

    my $dt_format = UR::Time->config('datetime');
    #UR::Time->config(datetime => '%a %b %d %T %Z %Y');
    UR::Time->config(datetime => '%Y-%m-%d');
    my $dt = UR::Time->numbers_to_datetime(0, 0, 0, $d, $m, "20$y");
    UR::Time->config(datetime => $dt_format);

    return $dt;
}

sub total_bases_read {
    my $self = shift;
    my $filter = shift; # optional?
    if(!defined($filter))
    {
        $filter = 'both';
    }
    my $total_bases; # unused?


    my $count;
    if ($self->is_paired_end) {
        # this changed in case we only want the fwd or rev counts...
        $count += ($self->fwd_read_length * $self->fwd_clusters)  unless $filter eq 'reverse-only';
        $count += ($self->rev_read_length * $self->rev_clusters) unless $filter eq 'forward-only';
    } else {
        $count += ($self->read_length * $self->clusters);
    }

    return $count;
}

sub run_identifier {
    my $self = shift;
    return $self->flow_cell_id;
}

sub resolve_median_insert_size {
    my $self = shift;

    #try eland metrics first
    if ($self->median_insert_size) {
        return $self->median_insert_size;
    }
    else {
        #try bwa metrics second
        return $self->get_default_alignment_metrics('median_insert_size');
        #Need try library_insert_size last ?
    }
}

sub resolve_sd_insert_size {
    my $self = shift;

    if ($self->sd_above_insert_size) {
        return $self->sd_above_insert_size;
    }
    else {
        return $self->get_default_alignment_metrics('sd_insert_size');
    }
}

sub get_default_alignment_metrics { #means BWA
    my ($self, $metric_name) = @_;
    my $metrics = $self->get_default_alignment_metrics_hash;

    unless ($metrics) {
        $self->error_message('Failed to get default alignment metrics hash for instrument data: '.$self->id);
        return;
    }

    return $metrics->{$metric_name} if exists $metrics->{$metric_name};
    return;
}

sub get_default_alignment_results {
    my $self = shift;

    my @alignment_results = $self->alignment_results_from_analysis_projects;
    return @alignment_results if @alignment_results;

    # Get alignment results for this inst data and the default aligner name, newest first
    my $pp = Genome::ProcessingProfile::ReferenceAlignment->default_profile;
    @alignment_results = sort {$b->best_guess_date cmp $a->best_guess_date} Genome::InstrumentData::AlignmentResult->get(
        instrument_data_id => $self->id,
        aligner_name       => $pp->read_aligner_name,
    );
    return unless @alignment_results;

    #Prefer alignment results created by apipe-builder models (via AQID autocron)
    my @ars_by_apipe_builder;
    my @ars_by_others;

    for my $alignment_result (@alignment_results) {
        #filter out the bad results and get results only with the new picard bwa alignment metrics
        next if $alignment_result->test_name;
        next unless grep{$_->metric_name eq 'read_1_pct_mismatch'}$alignment_result->metrics;

        my @builds = grep {$_->isa('Genome::Model::Build')} map{$_->user} $alignment_result->users;
        unless (@builds) {
            my @sr = grep {$_->isa('Genome::SoftwareResult')} map{$_->user} $alignment_result->users;
            if (@sr) {
                #look for build with one degree of indirection
                @builds = grep {$_->isa('Genome::Model::Build')} map{$_->user} map{$_->users} @sr;
            }
        }

        next unless @builds;
        my @models_by_apipe_builder = grep {$_->run_as =~ /^apipe-builder/} map { $_->model } @builds;

        if (@models_by_apipe_builder) {
            push @ars_by_apipe_builder, $alignment_result;
        }
        else {
            push @ars_by_others, $alignment_result;
        }
    }

    # Prefer alignment results run by apipe-builder
    if (@ars_by_apipe_builder) {
        my @ars_run_by_apipe_builder = grep {$_->output_dir =~ /\-apipe\-builder\-/} @ars_by_apipe_builder;
        if (@ars_run_by_apipe_builder) {
            return @ars_run_by_apipe_builder;
        }
        else {
            return @ars_by_apipe_builder;
        }
    }
    # lowest priority, in some rare cases, no LIMS eland metrics, no apipe-builder metrics.
    # Maybe we should not allow this to be used as default result.
    return @ars_by_others;
}

sub alignment_results_from_analysis_projects {
    my $self = shift;

    my @aps = $self->analysis_projects;
    return unless scalar(@aps);

    my @inputs = Genome::Model::Input->get(value_id => $self->id, 'model.analysis_projects.id' => [map $_->id, @aps]);
    return unless scalar(@inputs);

    my @alignment_results;
    for my $input (@inputs) {
        my $build = $input->model->last_complete_build;
        next unless $build;
        next unless Genome::Model::Build::Input->get(value_id => $self->id, build_id => $build->id);

        push @alignment_results,
            grep { $_->instrument_data_id eq $self->id }
            grep { $_->isa('Genome::InstrumentData::AlignmentResult') }
            $build->all_results;
    }

    return uniq(@alignment_results);
}

#This method is used in GSC::IndexIllumina to get bwa alignment metrics to retire eland
sub get_default_alignment_metrics_hash {
    my $self = shift;
    my @ar   = $self->get_default_alignment_results;

    my %metrics;

    for my $ar (@ar) {
        map{$metrics{$_->metric_name} = $_->metric_value}$ar->metrics;
        last;
    }
    return \%metrics;
}

sub is_capture {
    my $self = shift;

    return defined $self->target_region_set_name;
}


1;


