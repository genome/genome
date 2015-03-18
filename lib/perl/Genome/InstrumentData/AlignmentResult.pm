package Genome::InstrumentData::AlignmentResult;

use Genome;
use Sys::Hostname;
use IO::File;
use File::Path;
use Path::Class;
use YAML;
use Time::HiRes;
use POSIX qw(ceil);
use File::Copy;
use File::stat;
use Carp qw(confess);
use File::Basename;

use Genome::Utility::Instrumentation;

use warnings;
use strict;

class Genome::InstrumentData::AlignmentResult {
    is_abstract => 1,
    is => ['Genome::SoftwareResult::Stageable', 'Genome::SoftwareResult::WithNestedResults'],
    sub_classification_method_name => '_resolve_subclass_name',
    has => [
        instrument_data         => {
                                    is => 'Genome::InstrumentData',
                                    id_by => 'instrument_data_id'
                                },
        reference_build         => {
                                    is => 'Genome::Model::Build::ImportedReferenceSequence',
                                    id_by => 'reference_build_id',
                                },
        annotation_build        => {
                                    is => 'Genome::Model::Build::ImportedAnnotation',
                                    id_by => 'annotation_build_id',
                                    is_optional => 1,
                                },
        reference_name          => { via => 'reference_build', to => 'name', is_mutable => 0, is_optional => 1 },
        annotation_name         => { via => 'annotation_build', to => 'name', is_mutable => 0, is_optional => 1 },

        aligner                 => {
                                    calculate_from => [qw/aligner_name aligner_version aligner_params/],
                                    calculate => q|no warnings; "$aligner_name $aligner_version $aligner_params"|
                                },

        trimmer                 => {
                                    calculate_from => [qw/trimmer_name trimmer_version trimmer_params/],
                                    calculate => q|no warnings; "$trimmer_name $trimmer_version $trimmer_params"|
                                },

        filter                 => {
                                    calculate_from => [qw/filter_name filter_params force_fragment/],
                                    calculate => q|no warnings; "$filter_name $filter_params $force_fragment"|
                                },

        _disk_allocation        => { is => 'Genome::Disk::Allocation', is_optional => 1, is_many => 1, reverse_as => 'owner' },

    ],
    has_input => [
        instrument_data_id      => {
                                    is => 'Number',
                                    doc => 'the local database id of the instrument data (reads) to align',
                                },
        instrument_data_segment_type => {
                                    is => 'String',
                                    doc => 'Type of instrument data segment to limit within the instrument data being aligned (e.g. "read_group")',
                                    is_optional => 1,
        },
        instrument_data_segment_id => {
                                    is => 'String',
                                    doc => 'Identifier for instrument data segment to limit within the instrument data being aligned (e.g. read group ID)',
                                    is_optional => 1,
        },
        reference_build_id      => {
                                    is => 'Number',
                                    doc => 'the reference to use by id',
                                },
        annotation_build_id     => {
                                    is => 'Number',
                                    doc => 'the annotation to use by id',
                                    is_optional => 1,
                                },
    ],
    has_param => [
        aligner_name            => {
                                    is => 'Text',
                                    doc => 'the name of the aligner to use, maq, blat, newbler etc.',
                                },
        aligner_version         => {
                                    is => 'Text',
                                    doc => 'the version of the aligner to use, i.e. 0.6.8, 0.7.1, etc.',
                                    is_optional=>1,
                                },
        aligner_params          => {
                                    is => 'Text',
                                    is_optional=>1,
                                    doc => 'any additional params for the aligner in a single string',
                                },
        force_fragment          => {
                                    is => 'Boolean',
                                    is_optional=>1,
                                    doc => 'Force this run to be treated as a fragment run, do not do pairing',
                                },
        filter_name             => {
                                    is => 'Text',
                                    doc => 'Filter strategy to use',
                                    is_optional=>1,
                                },
        filter_params           => {
                                    is => 'Text',
                                    doc => 'Filter params to use',
                                    is_optional=>1,
                                },
        trimmer_name            => {
                                    is => 'Text',
                                    doc => 'Trimmer strategy to use',
                                    is_optional=>1,
                                },
        trimmer_version         => {
                                    is => 'Text',
                                    doc => 'Trimmer version to use',
                                    is_optional=>1,
                                },
        trimmer_params          => {
                                    is => 'Text',
                                    is_optional=>1,
                                    doc => 'Trimmer parameters',
                                },
        samtools_version        => {
                                    is=>'Text',
                                    is_optional=>1,
                                    doc=>'Version of samtools to use when creating BAM files',
                                },
        bedtools_version        => {
                                    is => 'Text',
                                    is_optional => 1,
                                    doc => "Version of bedtools to use for BAM > fastq conversion",
                                },
        picard_version          => {
                                    is=>'Text',
                                    is_optional=>1,
                                    doc=>'Version of picard to use when creating bam files',
                                },
        n_remove_threshold      => {
                                    is => 'Number',
                                    is_optional=>1,
                                    doc=>'If set, strips reads containing runs of this many Ns'
                                }
    ],
    has_metric => [
        cigar_md_error_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                        doc=>'The number of alignments with CIGAR / MD strings that failed to be parsed completely.'
                                },
        total_read_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        total_base_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        total_aligned_read_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        total_aligned_base_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        total_unaligned_read_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        total_unaligned_base_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        total_duplicate_read_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        total_duplicate_base_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        total_inserted_base_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        total_deleted_base_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        total_hard_clipped_read_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        total_hard_clipped_base_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        total_soft_clipped_read_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        total_soft_clipped_base_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        paired_end_read_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        paired_end_base_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        read_1_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        read_1_base_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        read_2_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        read_2_base_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        mapped_paired_end_read_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        mapped_paired_end_base_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        proper_paired_end_read_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        proper_paired_end_base_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        singleton_read_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        singleton_base_count => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
        bam_size  => {
                                        is=>'Number',
                                        is_optional=>1,
                                },
    ],
    has_transient => [
        _revivified_bam_file_path => { is => 'Text', is_optional=>1},
        temp_scratch_directory  => {
                                    is=>'Text',
                                    doc=>'Temp scratch directory',
                                    is_optional=>1,
                                },
        _input_fastq_pathnames => { is => 'ARRAY', is_optional => 1 },
        _input_bfq_pathnames   => { is => 'ARRAY', is_optional => 1 },
        _fastq_read_count      => { is => 'Number',is_optional => 1 },
        _sam_output_fh         => { is => 'IO::File',is_optional => 1 },
        _is_inferred_paired_end => { is => 'Boolean', is_optional=>1},
        _extracted_bam_path    => { is => 'String', is_optional=>1},
        _flagstat_file         => { is => 'Text', is_optional=>1},
        _temporary_input_files => { is => 'ARRAY', is_optional => 1, default_value => []},
    ],
};


sub __display_name__ {
    my $self = shift;

    my $test_name = $self->test_name;
    my $instrument_data = $self->instrument_data;
    my $instrument_data_segment_id = $self->instrument_data_segment_id;
    my $reference_build = $self->reference_build;
    my $annotation_build = ($self->can('annotation_build') ? $self->annotation_build : ()); # TODO: pull up and normalize

    my $name;
    if ($test_name) {
        $name = 'TEST: >>' . $test_name . '<< ';
    }
    else {
        $name = '';
    }

    $name .= $self->aligner_name
        . ' ' . $self->aligner_version
        . ' [' . $self->aligner_params . ']'
        . ' on ' . $instrument_data->__display_name__
        . (defined $instrument_data_segment_id ? " ($instrument_data_segment_id)" : '')
        . ($reference_build ? ' against ' . $reference_build->__display_name__ : '')
        . ($annotation_build ? ' annotated by ' . $annotation_build->__display_name__ : '')
        . " (" . $self->id . ")";

    return $name;
}

sub required_arch_os {
    # override in subclasses if 64-bit is not required
    'x86_64'
}

sub required_rusage {
    # override in subclasses
    # e.x.: "-R 'span[hosts=1] rusage[tmp=50000:mem=12000]' -M 1610612736";
    ''
}

sub lsf_queue {
    my $self = shift;

    # if class overrode then use that
    if ($self->__lsf_queue) {
        return $self->__lsf_queue;
    }

    if (Genome::Config->can('should_use_alignment_pd') && Genome::Config->should_use_alignment_pd($self->model)) {
        return $ENV{GENOME_LSF_QUEUE_ALIGNMENT_PROD};
    }

    return $ENV{GENOME_LSF_QUEUE_ALIGNMENT_DEFAULT};
}

sub required_rusage_for_building_index {
    # override if necessary in subclasses.
    my $class = shift;
    my %p = @_;
    my $reference_build = $p{reference_build};

    my $select = "select[mem>=10000 && tmp>=15000]";
    my $rusage = "rusage[mem=10000, tmp=15000]";
    my $options = "-M 10000000";

    return sprintf("-R '%s %s' %s", $select, $rusage, $options);
}


sub _working_dir_prefix {
    "alignment";
}


sub extra_metrics {
    # this will probably go away: override in subclasses if the aligner has custom metrics
    ()
}

sub _resolve_subclass_name {
    my $class = shift;

    if (ref($_[0]) and $_[0]->isa(__PACKAGE__)) {
        my $aligner_name = $_[0]->aligner_name;
        return join('::', 'Genome::InstrumentData::AlignmentResult', $class->_resolve_subclass_name_for_aligner_name($aligner_name));
    }
    elsif (my $aligner_name = $class->define_boolexpr(@_)->value_for('aligner_name')) {
        return join('::', 'Genome::InstrumentData::AlignmentResult', $class->_resolve_subclass_name_for_aligner_name($aligner_name));
    }
    return;
}

sub _resolve_subclass_name_for_aligner_name {
    my ($class,$aligner_name) = @_;
    my @type_parts = split(/[ -]/,$aligner_name);

    my @sub_parts = map { ucfirst } @type_parts;
    my $subclass = join('',@sub_parts);

    return $subclass;
}

sub create {
    my $class = shift;

    if ($class eq __PACKAGE__ or $class->__meta__->is_abstract) {
        # this class is abstract, and the super-class re-calls the constructor from the correct subclass
        return $class->SUPER::create(@_);
    }

    # STEP 1: verify the architecture on which we're running
    my $actual_os = Genome::Config->arch_os();
    $class->debug_message("OS is $actual_os");
    my $required_os = $class->required_arch_os;
    $class->debug_message("Required OS is $required_os");
    unless ($required_os eq $actual_os) {
        die $class->error_message("This logic can only be run on a $required_os machine!  (running on $actual_os)");
    }

    # STEP 2: the base class handles all locking, etc., so it may hang while waiting for a lock
    my $self = $class->SUPER::create(@_);
    return unless $self;

    if (my $output_dir = $self->output_dir) {
        if (-d $output_dir) {
            $self->debug_message("BACKFILL DIRECTORY: $output_dir!");
            return $self;
        }
    }

    # STEP 3: ENSURE WE WILL PROBABLY HAVE DISK SPACE WHEN ALIGNMENT COMPLETES
    # TODO: move disk_group, estimated_size, allocation and promotion up into the software result logic
    my $estimated_kb_usage = $self->estimated_kb_usage;
    $self->debug_message("Estimated disk for this data set: " . $estimated_kb_usage . " kb");
    $self->debug_message("Check for available disk...");
    my @available_volumes = Genome::Disk::Volume->get(disk_group_names => $ENV{GENOME_DISK_GROUP_ALIGNMENTS});
    $self->debug_message("Found " . scalar(@available_volumes) . " disk volumes");
    my $unallocated_kb = 0;
    for my $volume (@available_volumes) {
        $unallocated_kb += $volume->unallocated_kb;
    }
    $self->debug_message("Available disk: " . $unallocated_kb . " kb");
    my $factor = 20;
    unless ($unallocated_kb > ($factor * $estimated_kb_usage)) {
        $self->error_message("NOT ENOUGH DISK SPACE!  This step requires $factor x as much disk as the job will use to be available before starting.");
        die $self->error_message();
    }

    # STEP 4: PREPARE THE STAGING DIRECTORY
    $self->debug_message("Prepare working directories...");
    $self->_prepare_working_and_staging_directories;
    $self->debug_message("Staging path is " . $self->temp_staging_directory);
    $self->debug_message("Working path is " . $self->temp_scratch_directory);

    # STEP 5: PREPARE REFERENCE SEQUENCES
    $self->debug_message("Preparing the reference sequences...");
    unless($self->_prepare_reference_sequences) {
        $self->error_message("Reference sequences are invalid.  We can't proceed:  " . $self->error_message);
        die $self->error_message();
    }

    eval {

        # STEP 6: EXTRACT/COLLECT THE INPUTS TO THE ALIGNER
        my @inputs = $self->collect_inputs;
        unless (@inputs) {
            $self->error_message("Failed to gather input files: " . $self->error_message);
            die $self->error_message;
        }

        # STEP 7: PREPARE THE ALIGNMENT FILE (groups file, sequence dictionary)
        # this also prepares the bam output pipe and crams the alignment headers through it.
        $self->debug_message("Preparing the all_sequences.sam in scratch");
        unless ($self->prepare_scratch_sam_file) {
            $self->error_message("Failed to prepare the scratch sam file with groups and sequence dictionary");
            die $self->error_message;
        }

        # STEP 7: RUN THE ALIGNER
        $self->debug_message("Running aligner...");
        unless ($self->run_aligner(@inputs)) {
            $self->error_message("Failed to collect inputs and/or run the aligner!");
            die $self->error_message;
        }

        # STEP 8: CREATE BAM IN STAGING DIRECTORY
        if ($self->supports_streaming_to_bam) {
            $self->close_out_streamed_bam_file;
        } else {
            $self->debug_message("Constructing a BAM file (if necessary)...");
            unless( $self->create_BAM_in_staging_directory()) {
                $self->error_message("Call to create_BAM_in_staging_directory failed.\n");
                die $self->error_message;
            }
        }
    };

    if ($@) {
        my $error = $@;
        $self->error_message("Oh no!  Caught an exception while in the critical point where the BAM pipe was open: $@");
        if (defined $self->_sam_output_fh) {
            eval {
                $self->_sam_output_fh->close;
            };
            if ($@) {
                $error .= " ... and the input filehandle failed to close due to $@";
            }
        }

        die $error;
    }

    # STEP 9-10, validate BAM file (if necessary)
    $self->debug_message("Postprocessing & Sanity Checking BAM file (if necessary)...");
    unless ($self->postprocess_bam_file()) {
        $self->error_message("Postprocess BAM file failed");
        die $self->error_message;
    }

    # STEP 11: COMPUTE ALIGNMENT METRICS
    $self->debug_message("Computing alignment metrics...");
    $self->_compute_alignment_metrics();

    # STEP 12: PREPARE THE ALIGNMENT DIRECTORY ON NETWORK DISK
    $self->debug_message("Preparing the output directory...");
    $self->debug_message("Staging disk usage is " . $self->_staging_disk_usage . " KB");
    my $output_dir = $self->output_dir || $self->_prepare_output_directory;
    $self->debug_message("Alignment output path is $output_dir");

    # STEP 13: PROMOTE THE DATA INTO ALIGNMENT DIRECTORY
    $self->debug_message("Moving results to network disk...");
    my $product_path;
    unless($product_path= $self->_promote_data) {
        $self->error_message("Failed to de-stage data into alignment directory " . $self->error_message);
        die $self->error_message;
    }

    # STEP 14: RESIZE THE DISK
    $self->_reallocate_disk_allocation;

    $self->status_message("Alignment complete.");
    return $self;
}

sub delete {
    my $self = shift;

    my $name = $self->__display_name__;
    my $class_name = $self->class;

    # Find all the SoftwareResultUser objects that the one being deleted uses
    my @uses = Genome::SoftwareResult::User->get(user => $self);
    my @child_objects = map { $_->software_result } @uses;
    map { $_->delete } @uses;

    # find child objects for which there are no more users.
    for my $child (@child_objects) {
        my @users = Genome::SoftwareResult::User->get(software_result => $child, active => 1);
        $child->delete if !@users;
    }

    # find qc result and delete it
    my @qc_results = Genome::InstrumentData::AlignmentResult::Merged::BamQc->get(alignment_result_id => $self->id);
    for my $qc_result ( @qc_results ) {
        $qc_result->delete;
    }

    return $self->SUPER::delete(@_);
}

sub scratch_sam_file_path {
    my $self = shift;
    return File::Spec->catfile($self->temp_scratch_directory, "all_sequences.sam");
}

sub final_staged_bam_path {
    my $self = shift;
    return File::Spec->catfile($self->temp_staging_directory, "all_sequences.bam");
}

sub prepare_scratch_sam_file {
    my $self = shift;

    my $scratch_sam_file = $self->scratch_sam_file_path;

    unless($self->construct_groups_file) {
        $self->error_message("failed to create groups file");
        die $self->error_message;
    }

    my $groups_input_file = $self->temp_scratch_directory . "/groups.sam";

    my $seq_dict = $self->get_or_create_sequence_dictionary();
    unless (-s $seq_dict) {
        $self->error_message("Failed to get sequence dictionary");
        die $self->error_message;
    }

    my @input_files = ($seq_dict, $groups_input_file);

    $self->debug_message("Cat-ing together: ".join("\n",@input_files). "\n to output file ".$scratch_sam_file);
    my $cat_rv = Genome::Sys->cat(input_files=>\@input_files,output_file=>$scratch_sam_file);
    if ($cat_rv ne 1) {
        $self->error_message("Error during cat of alignment sam files! Return value $cat_rv");
        die $self->error_message;
    }
    else {
        $self->debug_message("Cat of sam files successful.");
    }

    if ($self->supports_streaming_to_bam) {
        my $ref_list  = $self->reference_build->full_consensus_sam_index_path($self->samtools_version);
        my $sam_cmd = sprintf("| %s view -S -b -o %s - ", Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version), $self->temp_scratch_directory . "/raw_all_sequences.bam");
        $self->debug_message("Opening $sam_cmd");

        $self->_sam_output_fh(IO::File->new($sam_cmd));
        unless ($self->_sam_output_fh()) {
            $self->error_message("We support streaming for this alignment module, but can't open a pipe to $sam_cmd");
            die $self->error_message;
        }

        my $temp_fh = IO::File->new($scratch_sam_file);
        unless ($temp_fh) {
            $self->error_message("Can't open temp sam header for reading.");
            die $self->error_message;
        }

        binmode $temp_fh;
        while (my $line = <$temp_fh>) {
            $self->_sam_output_fh->print($line);
        }
    }
    return 1;
}

# Override and return 1 if the aligner module can handle the read group
# extraction more efficiently than copying the whole bam (e.g., by piping).
sub can_extract_read_groups { 0 }

sub requires_fastqs_to_align {
    my $self = shift;
    # disqualify if the aligner can't take a bam
    return 1 unless ($self->accepts_bam_input);

    # read groups are ok here.
    return 1 if (defined $self->instrument_data_segment_id && !$self->instrument_data_segment_type eq 'read_group');

    # n-remove, complex filters, trimmers, instrument data disqualify bam processing
    return 1 if ($self->n_remove_threshold);
    return 1 if ($self->filter_name && ($self->filter_name ne 'forward-only' && $self->filter_name ne 'reverse-only'));
    return 1 if ($self->trimmer_name);

    # obviously we need fastq if we don't have a bam
    return 1 unless (defined $self->instrument_data->bam_path && -e $self->instrument_data->bam_path);

    return 0;
}

sub _extract_input_read_group_bam {
    my $self = shift;
    my $file = sprintf("%s/%s.rg_extracted_%s.bam", $self->temp_scratch_directory,  $self->instrument_data_id, $self->instrument_data_segment_id);

    my $cmd = Genome::Model::Tools::Sam::ExtractReadGroup->create(
        input         => $self->instrument_data->bam_path,
        output        => $file,
        name_sort     => 1,
        use_version   => $self->samtools_version,
        read_group_id => $self->instrument_data_segment_id,
    );

    unless ($cmd->execute) {
        $self->error_message($cmd->error_message);
        return;
    }

    $self->_extracted_bam_path($file);
    return $file;
}

sub collect_inputs {
    my $self = shift;

    $self->debug_message("Unpacking reads...");

    if ($self->requires_fastqs_to_align) {
        $self->debug_message('Requires fastqs to align');
        return $self->_extract_input_fastq_filenames;
    }

    # snag the bam from the instrument data
    my $instr_data = $self->instrument_data;
    my $bam_file = $instr_data->bam_path;
    unless ($bam_file) {
        $self->error_message('Instrument data has no bam_path!');
        return;
    }
    unless (-e $bam_file) {
        $self->error_message("BAM not found: $bam_file");
        return;
    }

    # maybe we want to extract a read group from that bam and deal with that instead...
    if (defined $self->instrument_data_segment_id && !$self->can_extract_read_groups) {
        $self->debug_message('Extract input read group bam');
        $bam_file = $self->_extract_input_read_group_bam;
        unless ($bam_file) {
            $self->error_message(sprintf('Failed to extract read group (%s) into temporary BAM.', $self->instrument_data_segment_id));
            return;
        }
    }

    $self->debug_message("Checking if this read group is paired end...");
    my $paired = $instr_data->is_paired_end;

    #One extra check to see whether the instrument data is really not
    #paired_end or just set wrong in production
    unless ($paired) {
        my $output_file = $bam_file . '.flagstat';
        unless (-s $output_file) {
            $output_file = $self->temp_scratch_directory . '/import_bam.flagstat';
            die unless $self->_create_bam_flagstat($bam_file, $output_file);
        }
        my $stats = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($output_file);
        die $self->error_message('Failed to get flagstat data on input bam: '. $bam_file) unless $stats;

        my $percent_paired = $stats->{reads_paired_in_sequencing} / $stats->{total_reads};

        # Boundaries were arbitrarily chosen, feel free to adjust as a matter
        # of policy.

        if ($percent_paired >= 0.9) {
            die $self->error_message('flagstat on input bam: '. $bam_file.' infers paired_end on instrument_data: '.$instr_data->id);
        }
    }


    return ("$bam_file:1", "$bam_file:2") if $paired;
    return ("$bam_file:0");
}

sub run_aligner {
    my ($self, @inputs) = @_;

    $self->debug_message("Got " . scalar(@inputs) . " input files");
    if (@inputs > 3) {
        $self->error_message("We don't support aligning with more than 3 inputs (the first 2 are treated as PE and last 1 is treated as SE)");
        die $self->error_message;
    }

    # Perform N-removal if requested

    if ($self->n_remove_threshold) {
        $self->debug_message("Running N-remove.  Threshold is " . $self->n_remove_threshold);

        my @n_removed_fastqs;

        for my $input_pathname (@inputs) {
            my $n_removed_file = $input_pathname . ".n-removed.fastq";
            my $n_remove_cmd = Genome::Model::Tools::Fastq::RemoveN->create(n_removed_file=>$n_removed_file, n_removal_threshold=>$self->n_remove_threshold, fastq_file=>$input_pathname);
            unless ($n_remove_cmd->execute) {
                $self->error_message("Error running RemoveN: " . $n_remove_cmd->error_message);
                die $self->error_message;
            }

            my $passed = $n_remove_cmd->passed_read_count();
            my $failed = $n_remove_cmd->failed_read_count();
            $self->debug_message("N removal complete: Passed $passed reads & Failed $failed reads");
            if ($passed > 0) {
                push @n_removed_fastqs, $n_removed_file;
            }

            if ($input_pathname =~ m/^\/tmp/) {
                $self->debug_message("Removing original file before N removal to save space: $input_pathname");
                unlink($input_pathname);
            }
        }
        if (@inputs == 1 && @n_removed_fastqs == 2) {
            $self->debug_message("NOTE: An entire side of the read pairs was filtered away after n-removal.  We'll be running in SE mode from here on out.");
        }

        if (@inputs == 0) {
            $self->error_message("All reads were filtered away after n-removal.  Nothing to do here, bailing out.");
            die $self->error_message;
        }

        @inputs = @n_removed_fastqs;
    }

    # STEP 7: DETERMINE HOW MANY PASSES OF ALIGNMENT ARE REQUIRED
    my @passes;
    if (
        defined($self->filter_name)
        and (
            $self->filter_name eq 'forward-only'
            or $self->filter_name eq 'reverse-only'
        )
    ) {
        my $filter_name = $self->filter_name;
        if (@inputs == 3) {
            die "cannot handle PE and SE data together with $filter_name only data"
        }
        elsif ($filter_name eq 'forward-only') {
            @passes = ( [ shift @inputs ] );
        }
        elsif ($filter_name eq 'reverse-only') {
            @passes = ( [ pop @inputs ] );
        }
        $self->debug_message("Running the aligner with the $filter_name filter.");
    }
    elsif ($self->force_fragment) {
        $self->debug_message("Running the aligner in force-fragment mode.");
        @passes = map { [ $_ ] } @inputs;
    }
    elsif (@inputs == 3) {
        $self->debug_message("Running aligner twice: once for PE & once for SE");
        @passes = ( [ $inputs[0], $inputs[1] ], [ $inputs[2] ] );
    }
    elsif (@inputs == 2) {
        $self->debug_message("Running aligner in PE mode");
        @passes = ( \@inputs );
    }
    elsif (@inputs == 1) {
        $self->debug_message("Running aligner in SE mode");
        @passes = ( \@inputs );
    }

    # STEP 8: RUN THE ALIGNER, APPEND TO all_sequences.sam IN SCRATCH DIRECTORY
    my $fastq_rd_ct = 0;

    if ($self->requires_fastqs_to_align) {
        for my $pass (@passes) {
            for my $file (@$pass) {
                my $line = `wc -l $file`;
                my ($wc_ct) = $line =~ /^(\d+)\s/;
                unless ($wc_ct) {
                    $self->error_message("Fail to count reads in FASTQ file: $file");
                    return;
                }
                if ($wc_ct % 4) {
                    $self->error_message("run has a line count of $wc_ct, which is not divisible by four!");
                    return;
                }
                $fastq_rd_ct += $wc_ct/4;
            }
        }
    }
    else {
        $fastq_rd_ct = $self->determine_input_read_count_from_bam;
    }

    unless ($fastq_rd_ct) {
        $self->error_message("Failed to get a read count in input files before aligning.");
        return;
    }
    $self->_fastq_read_count($fastq_rd_ct);

    for my $pass (@passes) {
        $self->debug_message("Aligning @$pass...");
        # note that the _run_aligner method must _append_ to any existing all_sequences.sam file
        # in case it is not being run on the first pass
        unless ($self->_run_aligner(@$pass)) {
            if (@$pass == 2) {
                $self->error_message("Failed to run aligner on first PE pass");
                die $self->error_message;
            }
            elsif (@$pass == 1) {
                $self->error_message("Failed to run aligner on final SE pass");
                die $self->error_message;
            }
            else {
                $self->error_message("Failed to run aligner on odd number of passes??");
                die $self->error_message;
            }
        }
    }

    for (@inputs) {
       if ($_ =~ m/^\/tmp\/.*\.fastq$/) {
           $self->debug_message("Unlinking fastq file to save space now that we've aligned: $_");
           unlink($_);
       }
    }

    return 1;
}

sub determine_input_read_count_from_bam {
    my $self = shift;

    my $bam_file    = $self->_extracted_bam_path || $self->instrument_data->bam_path;
    my $output_file = $self->temp_scratch_directory . "/input_bam.flagstat";

    die unless $self->_create_bam_flagstat($bam_file, $output_file);

    $self->_flagstat_file($output_file);
    my $stats = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($output_file);
    unless ($stats) {
        $self->status_message('Failed to get flagstat data  on input sequences from '.$output_file);
        return;
    }

    my $total_reads = 0;

    if ($self->filter_name) {
        my $filter_name = $self->filter_name;

        if ($filter_name eq 'forward-only') {
            $total_reads += $stats->{reads_marked_as_read1};
        } elsif ($filter_name eq 'reverse-only') {
            $total_reads += $stats->{reads_marked_as_read2};
        } else {
            $self->error_message("don't know how to handle $filter_name when counting reads in the bam.");
        }
    } else {
        $total_reads += $stats->{total_reads};
    }

    return $total_reads;
}


sub close_out_streamed_bam_file {
    my $self = shift;
    $self->debug_message("Closing bam file...");
    $self->_sam_output_fh->flush;
    $self->_sam_output_fh->close;
    $self->_sam_output_fh(undef);

    if ($self->requires_fixmate) {
        $self->debug_message("Sorting by name to do fixmate...");
        my $bam_file = $self->temp_scratch_directory . "/raw_all_sequences.bam";
        my $final_bam_file = $self->final_staged_bam_path;

        my $samtools = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);

        my $tmp_file = $bam_file.'.sort';
        my $rv = system "$samtools sort -n $bam_file $tmp_file";
        $self->error_message("Sort by name failed") and return if $rv or !-s $tmp_file.'.bam';
        $self->debug_message("unlinking original bam file $bam_file.");
        unlink $bam_file;

        # TODO: run htseq here
        # We need a way to have down-stream steps run before their predecessor cleans-up.

        $self->debug_message("Now running fixmate");
        $rv = system "$samtools fixmate $tmp_file.bam $tmp_file.fixmate";
        $self->error_message("fixmate failed") and return if $rv or !-s $tmp_file.'.fixmate';
        unlink "$tmp_file.bam";

        $self->debug_message("Now putting things back in chr/pos order");
        $rv = system "$samtools sort $tmp_file.fixmate $tmp_file.fix";
        $self->error_message("Sort by position failed") and return if $rv or !-s $tmp_file.'.fix.bam';

        unlink "$tmp_file.fixmate";
        unlink $bam_file;

        move "$tmp_file.fix.bam", $final_bam_file;
    } else {
        $self->debug_message("Skipping fixmate...");
        my $bam_file = $self->temp_scratch_directory . "/raw_all_sequences.bam";

        # TODO: run htseq here
        # We need a way to have down-stream steps run before their predecessor cleans-up.

        my $final_bam_file = $self->final_staged_bam_path;
        move $bam_file, $final_bam_file;
    }
    return 1;
}

sub create_BAM_in_staging_directory {
    my $self = shift;
    # STEP 9: CONVERT THE ALL_SEQUENCES.SAM into ALL_SEQUENCES.BAM
    unless($self->_process_sam_files) {
        $self->error_message("Failed to process sam files into bam files. " . $self->error_message);
        die $self->error_message;
    }

    return 1;
}

sub postprocess_bam_file {
    my $self = shift;

    my $bam_file    = $self->final_staged_bam_path;
    my $output_file = $bam_file . '.flagstat';

    #STEPS 8:  CREATE BAM.FLAGSTAT
    $self->debug_message("Creating all_sequences.bam.flagstat ...");
    die unless $self->_create_bam_flagstat($bam_file, $output_file);

    #STEPS 9: VERIFY BAM IS NOT TRUNCATED BY FLAGSTAT
    $self->debug_message("Verifying the bam...");
    unless ($self->_verify_bam) {
        $self->error_message('Fail to verify the bam');
        die $self->error_message;
    }

    #request by RT#62311 for submission and data integrity
    $self->debug_message('Creating all_sequences.bam.md5 ...');
    unless ($self->_create_bam_md5) {
        $self->error_message('Fail to create bam md5');
        die $self->error_message;
    }

    $self->debug_message("Indexing BAM file ...");
    unless($self->_create_bam_index) {
        $self->error_message('Fail to create bam md5');
        die $self->error_message;
    }
    return 1;
}

sub _use_alignment_summary_cpp { return 1; };

sub _compute_alignment_metrics {
    my $self = shift;
    my $bam = $self->final_staged_bam_path;

    if ($self->_use_alignment_summary_cpp){
        my $out = `bash -c "LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/gsc/scripts/opt/genome_legacy_code/lib:/gsc/pkg/boost/boost_1_42_0/lib /gsc/scripts/opt/genome_legacy_code/bin/alignment-summary-v1.2.6 --bam=\"$bam\" --ignore-cigar-md-errors"`;
        unless ($? == 0) {
            $self->error_message("Failed to compute alignment metrics.");
            die $self->error_message;
        }
        my $res = YAML::Load($out);
        unless (ref($res) eq "HASH") {
            $self->error_message("Failed to parse YAML hash from alignment_summary_cpp output.");
            die $self->error_message;
        }
        # ehvatum TODO: stop using samtools flagstat to detect truncation
        #$self->alignment_file_truncated     ($res->{truncated});
        $self->cigar_md_error_count         ($res->{cigar_md_error});
        $self->total_read_count             ($res->{total});
        $self->total_base_count             ($res->{total_bp});
        $self->total_aligned_read_count     ($res->{total_aligned});
        $self->total_aligned_base_count     ($res->{total_aligned_bp});
        $self->total_unaligned_read_count   ($res->{total_unaligned});
        $self->total_unaligned_base_count   ($res->{total_unaligned_bp});
        $self->total_duplicate_read_count   ($res->{total_duplicate});
        $self->total_duplicate_base_count   ($res->{total_duplicate_bp});
        $self->total_inserted_base_count    ($res->{total_inserted_bp});
        $self->total_deleted_base_count     ($res->{total_deleted_bp});
        $self->total_hard_clipped_read_count($res->{total_hard_clipped});
        $self->total_hard_clipped_base_count($res->{total_hard_clipped_bp});
        $self->total_soft_clipped_read_count($res->{total_soft_clipped});
        $self->total_soft_clipped_base_count($res->{total_soft_clipped_bp});
        $self->paired_end_read_count        ($res->{paired_end});
        $self->paired_end_base_count        ($res->{paired_end_bp});
        $self->read_1_count                 ($res->{read_1});
        $self->read_1_base_count            ($res->{read_1_bp});
        $self->read_2_count                 ($res->{read_2});
        $self->read_2_base_count            ($res->{read_2_bp});
        $self->mapped_paired_end_read_count ($res->{mapped_paired_end});
        $self->mapped_paired_end_base_count ($res->{mapped_paired_end_bp});
        $self->proper_paired_end_read_count ($res->{proper_paired_end});
        $self->proper_paired_end_base_count ($res->{proper_paired_end_bp});
        $self->singleton_read_count         ($res->{singleton});
        $self->singleton_base_count         ($res->{singleton_bp});

        $self->bam_size(stat($bam)->size); #store this for per lane bam recreation
        Genome::Utility::Instrumentation::inc('alignment_result.read_count', $self->total_read_count);
    }

    return 1;
}

sub _create_bam_index {
    my $self = shift;
    my $bam_file    = $self->final_staged_bam_path;

    unless (-s $bam_file) {
        $self->error_message('BAM file ' . $bam_file . ' does not exist or is empty');
        return;
    }

    my $cmd = Genome::Model::Tools::Sam::IndexBam->create(
        bam_file    => $bam_file,
        use_version => $self->samtools_version,
    );

    unless ($cmd->execute) {
        $self->error_message("Failed to index bam file $bam_file !");
        return;
    }

    return 1;
}

sub _create_bam_flagstat {
    my ($self, $bam_file, $output_file) = @_;

    unless (-s $bam_file) {
        $self->error_message('BAM file ' . $bam_file . ' does not exist or is empty');
        return;
    }

    if (-e $output_file) {
        $self->warning_message('Flagstat file '.$output_file.' exists. Now overwrite');
        unlink $output_file;
    }

    my $cmd = Genome::Model::Tools::Sam::Flagstat->create(
        bam_file       => $bam_file,
        output_file    => $output_file,
        use_version    => $self->samtools_version,
        include_stderr => 1,
    );

    unless ($cmd and $cmd->execute) {
        $self->error_message("Failed to create or execute flagstat command on bam: $bam_file");
        return;
    }
    return 1;
}


sub _verify_bam {
    my $self = shift;

    my $bam_file  = $self->final_staged_bam_path;
    my $flag_file = $bam_file . '.flagstat';
    my $flag_stat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flag_file);

    unless($flag_stat) {
        $self->status_message('Fail to get flagstat data from '.$flag_file);
        return;
    }

    if (exists $flag_stat->{errors}) {
        my @errors = @{$flag_stat->{errors}};

        for my $error (@errors) {
            if ($error =~ 'Truncated file') {
                $self->error_message('Bam file: ' . $bam_file . ' appears to be truncated');
                return;
            }
            else {
                $self->status_message('Continuing despite error messages from flagstat: ' . $error);
            }
        }
    }

    if (!exists $flag_stat->{total_reads} || $flag_stat->{total_reads} == 0) {
        $self->error_message("Bam file $bam_file has no reads reported (neither aligned nor unaligned).");
        return;
    }

    unless ($self->_check_read_count($flag_stat->{total_reads})) {
        $self->error_message("Bam file $bam_file failed read_count checking");
        return;
    }

    return 1;
}

sub _check_read_count {
    my ($self, $bam_rd_ct) = @_;
    my $fq_rd_ct = $self->_fastq_read_count;
    my $check = "Read count from bam: $bam_rd_ct and fastq: $fq_rd_ct";

    unless ($fq_rd_ct == $bam_rd_ct) {
        $self->error_message("$check does not match.");
        return;
    }
    $self->debug_message("$check matches.");
    return 1;
}


sub _create_bam_md5 {
    my $self = shift;

    my $bam_file = $self->final_staged_bam_path;
    my $md5_file = $bam_file . '.md5';
    my $cmd      = "md5sum $bam_file > $md5_file";

    my $rv  = Genome::Sys->shellcmd(
        cmd                        => $cmd,
        input_files                => [$bam_file],
        output_files               => [$md5_file],
        skip_if_output_is_present  => 0,
    );
    $self->error_message("Fail to run: $cmd") and return unless $rv == 1;
    return 1;
}


sub _promote_validated_data {
    my $self = shift;

    #my $container_dir = File::Basename::dirname($self->output_dir);
    my $staging_dir = $self->temp_staging_directory;
    my $output_dir  = $self->output_dir;

    $self->debug_message("Now de-staging data from $staging_dir into $output_dir");

    my $copy_cmd = sprintf("cp -rL %s/* %s/", $staging_dir, $output_dir);
    $self->debug_message("Running cp: $copy_cmd");
    my $copy_exit_code = system($copy_cmd);

    if ($copy_exit_code != 0) {

        $self->debug_message("Copy failed, attempting rsync");

        my $rsync_cmd = sprintf("rsync -avzL %s/* %s/", $staging_dir, $output_dir);

        $self->debug_message("Running Rsync: $rsync_cmd");
        my $rsync_exit_code = system($rsync_cmd);

        unless ($rsync_exit_code == 0) {
            $self->error_message("Did not get a valid return from rsync, exit code was $rsync_exit_code for call $rsync_cmd.  Cleaning up and bailing out.");
            rmtree($output_dir);
            die $self->error_message;
        }
    }

    $self->_disk_allocation->set_files_read_only;

    $self->debug_message("Files in $output_dir: \n" . join "\n", glob($output_dir . "/*"));

    return $output_dir;
}

sub _process_sam_files {
    my $self = shift;
    my $groups_input_file;

    # if a bam file is already staged at the end of _run_aligner, trust it to be correct.
    if (-e $self->final_staged_bam_path) {
        return 1;
    }

    my $sam_input_file = $self->scratch_sam_file_path;

    unless (-e $sam_input_file) {
        $self->error_message("$sam_input_file is nonexistent.  Can't convert!");
        die $self->error_message;
    }

    # things which don't produce sam natively must provide an unaligned reads file.
    my $unaligned_input_file = $self->temp_scratch_directory . "/all_sequences_unaligned.sam";

    if (-s $unaligned_input_file) {
        $self->debug_message("Looks like there are unaligned reads not in the main input file.  ");
        my @input_files = ($sam_input_file, $unaligned_input_file);
        $self->debug_message("Cat-ing the unaligned list $unaligned_input_file to the sam file $sam_input_file");
        my $cat_rv = Genome::Sys->cat(input_files=>[$unaligned_input_file],output_file=>$sam_input_file,append_mode=>1);
        if ($cat_rv ne 1) {
            $self->error_message("Error during cat of alignment sam files! Return value $cat_rv");
            die $self->error_message;
        } else {
            $self->debug_message("Cat of sam files successful.");
        }

        unlink($unaligned_input_file);
    }

    my $per_lane_sam_file_rg = $sam_input_file;

    if ($self->requires_read_group_addition) {
        $per_lane_sam_file_rg = $self->temp_scratch_directory . "/all_sequences_rg.sam";
        my $add_rg_cmd = Genome::Model::Tools::Sam::AddReadGroupTag->create(
            input_file     => $sam_input_file,
            output_file    => $per_lane_sam_file_rg,
            read_group_tag => $self->read_and_platform_group_tag_id,
        );

        unless ($add_rg_cmd->execute) {
            $self->error_message("Adding read group to sam file failed!");
            die $self->error_message;
        }
        $self->debug_message("Read group add completed, new file is $per_lane_sam_file_rg");

        $self->debug_message("Removing non-read-group combined sam file: " . $sam_input_file);
        unlink($sam_input_file);
    }

    #For the sake of new bam flagstat that need MD tags added. Some
    #aligner like maq doesn't output MD tag in sam file, now add it
    my $final_sam_file;

    if ($self->fillmd_for_sam) {
        my $sam_path = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);
        my $ref_seq  = $self->reference_build->full_consensus_path('fa');
        $final_sam_file = $self->temp_scratch_directory . '/all_sequences.fillmd.sam';

        my $cmd = "$sam_path fillmd -S $per_lane_sam_file_rg $ref_seq 1> $final_sam_file 2>/dev/null";

        my $rv  = Genome::Sys->shellcmd(
            cmd                          => $cmd,
            input_files                  => [$per_lane_sam_file_rg, $ref_seq],
            output_files                 => [$final_sam_file],
            skip_if_output_is_present    => 0,
        );
        $self->error_message("Fail to run: $cmd") and return unless $rv == 1;
        unlink $per_lane_sam_file_rg;
    }
    else {
        $final_sam_file = $per_lane_sam_file_rg;
    }

    my $ref_list  = $self->reference_build->full_consensus_sam_index_path($self->samtools_version);
    unless ($ref_list) {
        $self->error_message("Failed to get MapToBam ref list: $ref_list");
        return;
    }

    my $per_lane_bam_file = $self->final_staged_bam_path;
    my %params = (
        bam_file => $per_lane_bam_file,
        sam_file => $final_sam_file,
        keep_sam => 0,
        fix_mate => $self->requires_fixmate,
        index_bam => 1,
        ref_list => $ref_list,
        use_version => $self->samtools_version,
    );

    my $to_bam = Genome::Model::Tools::Sam::SamToBam->create(
        %params,
    );

    unless($to_bam->execute) {
        $self->error_message("There was an error converting the Sam file $final_sam_file to $per_lane_bam_file.");
        die $self->error_message;
    }

    $self->debug_message("Conversion successful.  File is: $per_lane_bam_file");
    return 1;
}


sub _gather_params_for_get_or_create {
    my $class = shift;
    my $bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);

    my $aligner_name = $bx->value_for('aligner_name');
    my $subclass = join '::', 'Genome::InstrumentData::AlignmentResult', $class->_resolve_subclass_name_for_aligner_name($aligner_name);

    my %params = $bx->params_list;
    my %is_input;
    my %is_param;
    my $class_object = $subclass->__meta__;
    for my $key ($subclass->property_names) {
        my $meta = $class_object->property_meta_for_name($key);
        if ($meta->{is_input} && exists $params{$key}) {
            $is_input{$key} = $params{$key};
        } elsif ($meta->{is_param} && exists $params{$key}) {
            $is_param{$key} = $params{$key};
        }

    }

    my %software_result_params = ( subclass_name => $subclass );

    return {
        software_result_params => \%software_result_params,
        subclass => $subclass,
        inputs=>\%is_input,
        params=>\%is_param,
    };

}


sub _prepare_working_and_staging_directories {
    my $self = shift;

    unless ($self->_prepare_staging_directory) {
        $self->error_message("Failed to prepare staging directory");
        return;
    }
    my $hostname = hostname;
    my $user = $ENV{'USER'};

    my $scratch_basedir = sprintf("scratch-%s-%s-%s-%s", $hostname, $user, $$, $self->id);
    my $scratch_tempdir =  Genome::Sys->create_temp_directory($scratch_basedir);
    $self->temp_scratch_directory($scratch_tempdir);
    unless($scratch_tempdir) {
        die "failed to create a temp scrach directory for working files";
    }

    return 1;
}


sub _staging_disk_usage {

    my $self = shift;
    my $usage;
    unless ($usage = Genome::Sys->disk_usage_for_path($self->temp_staging_directory)) {
        $self->error_message("Failed to get disk usage for staging: " . Genome::Sys->error_message);
        die $self->error_message;
    }

    return $usage;
}

sub estimated_kb_usage {
    30000000;
    #die "unimplemented method: please define estimated_kb_usage in your alignment subclass.";
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $instrument_data = $self->instrument_data;
    my $staged_basename = File::Basename::basename($self->temp_staging_directory);
    # TODO: the first subdir is actually specified by the disk management system.
    my $directory = join('/', 'alignment_data',$instrument_data->id,$staged_basename);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_ALIGNMENTS};
}


sub _extract_input_fastq_filenames {
    my $self = shift;
    my $instrument_data = $self->instrument_data;

    my %segment_params;

    if (defined $self->instrument_data_segment_type) {
        # sanity check this can be segmented
        if (! $instrument_data->can('get_segments') && $instrument_data->get_segments > 0) {
            $self->error_message("requested to align a given segment, but this instrument data either can't be segmented or has no segments.");
            die $self->error_message;
        }

        # only read groups for now
        if ($self->instrument_data_segment_type ne 'read_group') {
            $self->error_message("specified a segment type we don't support, " . $self->instrument_data_segment_type . ". we only support read group at present.");
            die $self->error_message;
        }

        if (defined $self->filter_name) {
            $self->error_message("filtering reads is currently not supported with segmented inputs, FIXME.");
            die $self->error_message;
        }

        $segment_params{read_group_id} = $self->instrument_data_segment_id;
    }

    my @input_fastq_pathnames;
    if ($self->_input_fastq_pathnames) {
        @input_fastq_pathnames = @{$self->_input_fastq_pathnames};
        my $errors;
        for my $input_fastq (@input_fastq_pathnames) {
            unless (-e $input_fastq && -f $input_fastq && -s $input_fastq) {
                $self->error_message('Missing or zero size sanger fastq file: '. $input_fastq);
                die($self->error_message);
            }
        }
    }
    else {
        # FIXME - getting a warning about undefined string with 'eq'
        if (! defined($self->filter_name)) {
            $self->debug_message('No special filter for this input');
        }
        elsif ($self->filter_name eq 'forward-only') {
            # forward reads only
            $self->debug_message('forward-only filter applied for this input');
        }
        elsif ($self->filter_name eq 'reverse-only') {
            # reverse reads only
            $self->debug_message('reverse-only filter applied for this input');
        }
        else {
            die 'Unsupported filter: "' . $self->filter_name . '"!';
        }

        my $temp_directory = Genome::Sys->create_temp_directory;

        my $trimmer_name = $self->trimmer_name;

        @input_fastq_pathnames = $instrument_data->dump_trimmed_fastq_files(
            segment_params => \%segment_params,
            trimmer_name => $self->trimmer_name,
            trimmer_version => $self->trimmer_version,
            trimmer_params => $self->trimmer_params,
            directory => $temp_directory,
        );
        if (!@input_fastq_pathnames){
            $self->error_message("no input_fastq_pathnames returned from dump_trimmed_fastq_files");
        }

        my @report_files = glob($temp_directory ."/*report*");
        my @length_distributions = glob($temp_directory ."/*.lengthdist");
        for my $report_file (@report_files, @length_distributions) {
            my $staged_path = $report_file;
            my $staging_directory = $self->temp_staging_directory;
            $staged_path =~ s/$temp_directory/$staging_directory/;
            Genome::Sys->copy_file($report_file, $staged_path);
            unlink($report_file);
        }
        $self->_input_fastq_pathnames(\@input_fastq_pathnames);
        $self->add_to_temporary_input_files_queue(@input_fastq_pathnames);
    }
    return @input_fastq_pathnames;
}

sub _prepare_reference_sequences {
    my $self = shift;
    my $reference_build = $self->reference_build;

    my $ref_basename = File::Basename::fileparse($reference_build->full_consensus_path('fa'));
    my $reference_fasta_path = sprintf("%s/%s", $reference_build->data_directory, $ref_basename);

    unless(-e $reference_fasta_path) {
        $self->error_message("Alignment reference path $reference_fasta_path does not exist");
        die $self->error_message;
    }

    return 1;
}

sub get_or_create_sequence_dictionary {
    my $self = shift;

    my $species = "unknown";
    $self->debug_message("Sample id: ".$self->instrument_data->sample_id);
    my $sample = Genome::Sample->get($self->instrument_data->sample_id);
    if ( defined($sample) ) {
        $species =  $sample->species_name;
        if (!$species || $species eq ""  ) {
            $species = "unknown";
        }
    }

    $self->debug_message("Species from alignment: ".$species);

    my $ref_build = $self->reference_build;
    my $seq_dict = $ref_build->get_sequence_dictionary("sam",$species,$self->picard_version);
    return $seq_dict;
}

sub _sam_header_extra {
    my $self = shift;

    my $instr_data = $self->instrument_data;

    my $seems_paired = $self->_is_inferred_paired_end;
    my $paired = defined $seems_paired ? $seems_paired : $instr_data->is_paired_end;
    my $description_for_header = $paired ? "paired end" : "fragment";
    my $aligner_command_line = $self->aligner_params_for_sam_header;

    my $id_tag = $self->read_and_platform_group_tag_id;
    my $pu_tag = sprintf("%s.%s", $instr_data->run_identifier, $instr_data->subset_name);
    my $lib_tag = $instr_data->library_name;
    my $date_run_tag = $instr_data->run_start_date_formatted;
    my $sample_tag = $instr_data->sample_name;
    my $aligner_version_tag = $self->aligner_version;
    my $aligner_cmd  =  $aligner_command_line;
    my $platform = $instr_data->sequencing_platform;
    $platform = ($platform eq 'solexa' ? 'illumina' : $platform);

    my @rg_data = (
        "\@RG",
        "ID:$id_tag",
        "PL:$platform",
        "PU:$pu_tag",
        "LB:$lib_tag",
        "DS:$description_for_header",
        "DT:$date_run_tag",
        "SM:$sample_tag",
        "CN:WUGSC",
        );

    my @pg_data = (
        "\@PG",
        "ID:$id_tag",
        "VN:$aligner_version_tag",
        "CL:$aligner_cmd",
        );

    return {
        RG => join("\t", @rg_data),
        PG => join("\t", @pg_data),
        };
}

sub construct_groups_file {
    my $self = shift;
    my $output_file = shift || $self->temp_scratch_directory . "/groups.sam";

    my $extra = $self->_sam_header_extra;

    my $rg_tag = $extra->{RG};
    my $pg_tag = $extra->{PG};
    die "Failed to generate \@RG line for sam header" unless $rg_tag;
    die "Failed to generate \@PG line for sam header" unless $pg_tag;

    $self->debug_message("RG: $rg_tag");
    $self->debug_message("PG: $pg_tag");

    my $fh = IO::File->new($output_file, "a")
        || die "failed opening groups file for writing";

    $fh->printf("%s\n", $rg_tag);
    $fh->printf("%s\n", $pg_tag);
    $fh->close;

    unless (-s $output_file) {
        $self->error_message("Failed to create groups file");
        die $self->error_message;
    }

    return 1;
}

sub read_and_platform_group_tag_id {
    my $self = shift;
    my $id = $self->instrument_data->id;

    if (defined $self->instrument_data_segment_id) {
        $id .= "-". $self->instrument_data_segment_id;
    }

    return $id;
}

sub aligner_params_for_sam_header {
    die "You must implement aligner_params_for_sam_header in your AlignmentResult subclass. This specifies the parameters used to align the reads";
}

sub fillmd_for_sam {
    #Maybe this can be set to return 0 as default.
    die 'Must implement fillmd_for_sam in AlignmentResult subclass. return either 1 or 0';
}

sub verify_alignment_data {
    return 1;
}

# This will generally only return something until the first merge has deleted the original per-lane bam file
sub alignment_bam_file_paths {
    return glob(File::Spec->join(shift->output_dir, "*.bam"));
}

# This method will recreate the per-lane bam file and return the path.
# This must be provided an allocation into which the bam will go (so that it can be cleaned up in the parent process once it is done).
# The calling process is responsible for cleaning up the allocation after we are done with it.
sub revivified_alignment_bam_file_paths {
    my $self = shift;
    my %p = Params::Validate::validate(@_, {disk_allocation => { isa => 'Genome::Disk::Allocation'}});

    # If we have a merged alignment result, the per-lane bam can be regenerated and we will do so now
    # This is less efficient than using the in-place per-lane bam. However, it aids us in terms of
    # contention, and in our transition period (allowing us to delete all per-lane bams now).
    if ($self->get_merged_alignment_results) {
        return $self->_revivified_bam_file_path if defined $self->_revivified_bam_file_path;
    } elsif (my @bams = $self->alignment_bam_file_paths) {
        return @bams;
    }

    my $revivified_bam = File::Spec->join($p{disk_allocation}->absolute_path, 'all_sequences.bam');
    my $merged_bam    = $self->get_merged_bam_to_revivify_per_lane_bam;

    unless ($merged_bam and -s $merged_bam) {
        die $self->error_message('Failed to get valid merged bam to recreate per lane bam '.$self->id);
    }

    my $cmd = Genome::InstrumentData::AlignmentResult::Command::RecreatePerLaneBam->create(
        merged_bam          => $merged_bam,
        per_lane_bam        => $revivified_bam,
        instrument_data_id  => $self->read_and_platform_group_tag_id,
        samtools_version    => $self->samtools_version,
        picard_version      => $self->picard_version,
        bam_header          => $self->bam_header_path,
        comparison_flagstat => $self->flagstat_path,
    );

    unless ($cmd->execute) {
        die $self->error_message('Failed to execute RecreatePerLaneBam for '.$self->id);
    }

    if (-s $revivified_bam) {
        # Cache the path of this revivified bam for future access
        $self->_revivified_bam_file_path($revivified_bam);
        return ($revivified_bam);
    }
    else {
        die $self->error_message("After running RecreatePerLaneBam, no per-lane bam (%s) exists still!", $revivified_bam);
    }
}


sub bam_header_path {
    return File::Spec->join(shift->output_dir, 'all_sequences.bam.header');
}


sub flagstat_path {
    return File::Spec->join(shift->output_dir, 'all_sequences.bam.flagstat');
}


sub get_merged_bam_to_revivify_per_lane_bam {
    my $self = shift;
    my $merged_result = $self->get_smallest_merged_alignment_result($self->get_unarchived_merged_alignment_results);

    unless ($merged_result) {
        $merged_result = $self->get_smallest_merged_alignment_result($self->get_merged_alignment_results);
        unless ($merged_result) {
            die $self->error_message('Failed to get archived merged result for per lane alignment '.$self->id);
        }
        $merged_result->_auto_unarchive;
    }

    my $merged_bam = $merged_result->merged_alignment_bam_path;
    unless (-s $merged_bam) {
        die $self->error_message("Merged bam (%s) does not exist for merged result id (%s)", $merged_bam, $merged_result->id);
    }
    return $merged_bam;
}

sub get_merged_alignment_results {
    my $self = shift;
    # Always load from the database, since other merged results may have committed since we updated the UR cache
    my @results = Genome::InstrumentData::AlignmentResult::Merged->load(
        'inputs.value_id' => $self->instrument_data_id,
        test_name => $self->test_name,
    );
    my @filtered_results = $self->filter_non_database_objects(@results);
    return $self->filter_non_matching_results(@filtered_results);
}

# This was refactored out to override in the test - mock objects break this logic
sub filter_non_database_objects {
    my ($self, @results) = @_;
    my @db_results;
    for my $result (@results) {
        if (UR::Context->current->object_exists_in_underlying_context($result)) {
            push @db_results, $result;
        }
    }
    return @db_results;
}

# This uses a sort of 'backwards lookup' from merged alignment results going back to per-lane alignment results.
# This is nice because it allows us to keep the logic in one place rather than mirroring it here.
sub filter_non_matching_results {
    my ($self, @merged_results) = @_;

    my @matching_results;
    for my $merged_result (@merged_results) {
        my @individual_results = $merged_result->collect_individual_alignments($self->_user_data_for_nested_results);
        if (@individual_results) {
            push @matching_results, $merged_result if grep{$_->id eq $self->id}@individual_results;
        }
    }

    return @matching_results;
}

sub get_unarchived_merged_alignment_results {
    my $self = shift;
    my @merged = $self->get_merged_alignment_results;
    my @unarchived;
    for my $merged (@merged) {
        unless ( grep { $_->is_archived } $merged->disk_allocations ) {
            push @unarchived, $merged;
        }
    }
    return @unarchived;
}

sub get_smallest_merged_alignment_result {
    my ($self, @merged_alignment_results) = @_;
    return unless @merged_alignment_results;
    my $smallest_result = $merged_alignment_results[0];

    for my $result (@merged_alignment_results) {
        my @current_instrument_data = $result->instrument_data;
        my @smallest_instrument_data = $smallest_result->instrument_data;
        if ( scalar(@current_instrument_data) < scalar(@smallest_instrument_data) ) {
            $smallest_result = $result;
        }
    }

    return $smallest_result;
}

sub requires_read_group_addition {
    return 1;
}

sub supports_streaming_to_bam {
    0;
}

sub requires_fixmate {
    1;
}

sub accepts_bam_input {
    0;
}

sub aligner_params_required_for_index {
    0;
}

sub aligner_name_for_aligner_index {
    my $self = shift;
    return $self->aligner_name;
}

sub get_reference_sequence_index {
    my $self = shift;
    my $build = shift || $self->reference_build;
    my @overrides = @_;

    my $index = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_with_lock(
        aligner_name => $self->aligner_name_for_aligner_index,
        aligner_version => $self->aligner_version,
        aligner_params => $self->aligner_params,
        reference_build => $build,
        users => $self->_user_data_for_nested_results,
        @overrides
        );

    if (!$index) {
        die $self->error_message(sprintf("No reference index prepared for %s with params %s and reference build %s", $self->aligner_name, $self->aligner_params, $self->reference_build->id));
    }

    return $index;
}

# this behavior was in the alignment class earlier and was set as changeable in SoftwareResult::Stageable.
sub _needs_symlinks_followed_when_syncing {
    1;
}


# note: this may be completely wrong. fix later!
sub _derive_insert_size_bounds {
    my ($self, $up_default, $low_default) = @_;

    my $median = $self->instrument_data->resolve_median_insert_size;
    my $stddev = $self->instrument_data->resolve_sd_insert_size;

    my ($upper, $lower);

    if (defined $median && defined $stddev) {
        $upper = $median + $stddev*5;
        $lower = $median - $stddev*5;
    }

    if (!defined $upper || $upper <= 0) {
        $self->debug_message("Calculated upper bound on insert size is undef or less than 0, defaulting to $up_default");
        $upper = $up_default;
    }
    if (!defined $lower || not $median || $lower < 0 || $lower > $upper) {
        $self->debug_message("Calculated lower bound on insert size is undef or invalid, defaulting to $low_default");
        $lower = $low_default;
    }
    return ($lower, $upper);
}

sub show_temporary_input_files_queue {
    my $self = shift;

    my @paths = $self->temporary_input_files_queue;

    unless(@paths) {
        $self->debug_message("Paths in Temporary Storage Queue: None");
        return 1;
    }

    $self->debug_message("Paths in Temporary Storage Queue:");
    for (my $i = 0; $i < @paths; $i++) {
        if ($paths[$i]->is_dir) {
            $self->debug_message('---> [%d] (Directory) %s', $i, $paths[$i]->stringify);
        } else {
            $self->debug_message('---> [%d] (File) %s', $i, $paths[$i]->stringify);
        }
    }

    return 1;
}

sub clear_temporary_input_files_queue {
    my $self = shift;

    while (my $path = shift @{$self->_temporary_input_files}) {
        if ($path->is_dir) {
            $self->debug_message(
                "[delete] Temprorary Storage Queue:  "
                . "(Directory) - '$path'"
            );
            $path->rmtree;
        }
        else {
            $self->debug_message(
                "[delete] Temprorary Storage Queue:  "
                . "(File) - '$path'"
            );
            $path->remove;
        }
    }

    return 1;
}

sub temporary_input_files_queue {
    my ($self, @args) = @_;
    return @{$self->_temporary_input_files}
}

sub add_to_temporary_input_files_queue {
    my ($self, @input_paths) = @_;

    for my $path (@input_paths) {
        my $p;
        if (-d $path) {
            $self->debug_message(
                "[push] Temprorary Storage Queue:  "
                . "(Directory) - '$path'"
            );
            $p = Path::Class::Dir->new("$path");
        }
        else {
            $self->debug_message(
                "[push] Temprorary Storage Queue:  "
                . "(File) - '$path'"
            );
            $p = Path::Class::File->new("$path");
        }
        push(@{$self->_temporary_input_files}, $p);
    }

    return 1;
}

sub lock_bam_file_access {
    my $self = shift;
    unless ($self->get_merged_alignment_results) {
        my @bams = $self->alignment_bam_file_paths;
        unless (@bams) {
            die $self->error_message("Alignment result with class (%s) and id (%s) has neither ".
                "merged results nor valid bam paths. This likely means that this alignment result ".
                "needs to be removed and realigned because data has been lost. ".
                "Please create an apipe-support ticket for this.", $self->class, $self->id);
            #There is no way to recreate per lane bam if merged bam
            #does not exist and per lane bam is removed. Software
            #result of this per lane alignment needs to be removed and
            #this per lane instrument data needs to be realigned
            #with that aligner.
        }

        my $lock_var = File::Spec->join('genome', __PACKAGE__, 'lock-per-lane-alignment-'.$self->id);
        my $lock = Genome::Sys->lock_resource(
            resource_lock => $lock_var,
            scope         => 'site',
            max_try       => 288, # Try for 48 hours every 10 minutes
            block_sleep   => 600,
        );
        die $self->error_message("Unable to acquire the lock for per lane alignment result id (%s) !", $self->id) unless $lock;

        # If the build before us successfully created a merged alignment result, we no longer need a lock
        # If it failed, we will add an observer just as the first build did.
        if ($self->get_merged_alignment_results) {
            Genome::Sys->unlock_resource(resource_lock => $lock);
        } else {
            # The problem here is if we commit BEFORE merge is done, we unlock too early.
            # However, if we unlock any other way we may fail to unlock more often and leave old locks.
            UR::Context->process->add_observer(
                aspect   => 'commit',
                once => 1,
                callback => sub {
                    Genome::Sys->unlock_resource(resource_lock => $lock);
                }
            );
        }
    }
}

1;

