package Genome::InstrumentData::AlignmentResult;

use Genome;
use Sys::Hostname;
use IO::File;
use File::Path;
use YAML;
use Time::HiRes;
use POSIX qw(ceil);
use File::Copy;
use Carp qw(confess);

use Genome::Utility::Instrumentation;

use warnings;
use strict;

class Genome::InstrumentData::AlignmentResult {
    is_abstract => 1,
    is => 'Genome::SoftwareResult::Stageable',
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
                                    is => 'Text', default_value => 'maq',
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
    ],
    has_transient => [
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
        . ($instrument_data_segment_id ? " ($instrument_data_segment_id)" : '')
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
    # e.x.: "-R 'select[model!=Opteron250 && type==LINUX64] span[hosts=1] rusage[tmp=50000:mem=12000]' -M 1610612736";
    ''
}

sub required_rusage_for_building_index {
    # override if necessary in subclasses.
    my $class = shift;
    my %p = @_;
    my $reference_build = $p{reference_build};

    my $select = "select[mem>=10000 && tmp>=15000]";
    my $rusage = "rusage[mem=10000, tmp=15000]";
    my $options = "-M 10000000 -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER}";

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
    $class->status_message("OS is $actual_os");
    my $required_os = $class->required_arch_os;
    $class->status_message("Required OS is $required_os");
    unless ($required_os eq $actual_os) {
        die $class->error_message("This logic can only be run on a $required_os machine!  (running on $actual_os)");
    }

    # STEP 2: the base class handles all locking, etc., so it may hang while waiting for a lock
    my $self = $class->SUPER::create(@_);
    return unless $self;

    if (my $output_dir = $self->output_dir) {
        if (-d $output_dir) {
            $self->status_message("BACKFILL DIRECTORY: $output_dir!");
            return $self;
        }
    }

    # STEP 3: ENSURE WE WILL PROBABLY HAVE DISK SPACE WHEN ALIGNMENT COMPLETES
    # TODO: move disk_group, estimated_size, allocation and promotion up into the software result logic
    my $estimated_kb_usage = $self->estimated_kb_usage;
    $self->status_message("Estimated disk for this data set: " . $estimated_kb_usage . " kb");
    $self->status_message("Check for available disk...");
    my @available_volumes = Genome::Disk::Volume->get(disk_group_names => "info_alignments");
    $self->status_message("Found " . scalar(@available_volumes) . " disk volumes");
    my $unallocated_kb = 0;
    for my $volume (@available_volumes) {
        $unallocated_kb += $volume->unallocated_kb;
    }
    $self->status_message("Available disk: " . $unallocated_kb . " kb");
    my $factor = 20;
    unless ($unallocated_kb > ($factor * $estimated_kb_usage)) {
        $self->error_message("NOT ENOUGH DISK SPACE!  This step requires $factor x as much disk as the job will use to be available before starting.");
        die $self->error_message();
    }

    # STEP 4: PREPARE THE STAGING DIRECTORY
    $self->status_message("Prepare working directories...");
    $self->_prepare_working_and_staging_directories;
    $self->status_message("Staging path is " . $self->temp_staging_directory);
    $self->status_message("Working path is " . $self->temp_scratch_directory);

    # STEP 5: PREPARE REFERENCE SEQUENCES
    $self->status_message("Preparing the reference sequences...");
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
        $self->status_message("Preparing the all_sequences.sam in scratch");
        unless ($self->prepare_scratch_sam_file) {
            $self->error_message("Failed to prepare the scratch sam file with groups and sequence dictionary");
            die $self->error_message;
        }

        # STEP 7: RUN THE ALIGNER
        $self->status_message("Running aligner...");
        unless ($self->run_aligner(@inputs)) {
            $self->error_message("Failed to collect inputs and/or run the aligner!");
            die $self->error_message;
        }

        # STEP 8: CREATE BAM IN STAGING DIRECTORY
        if ($self->supports_streaming_to_bam) {
            $self->close_out_streamed_bam_file;
        } else {
            $self->status_message("Constructing a BAM file (if necessary)...");
            unless( $self->create_BAM_in_staging_directory()) {
                $self->error_message("Call to create_BAM_in_staging_directory failed.\n");
                die $self->error_message;
            }
        }
    };

    if ($@) {
        my $error = $@;
        $self->status_message("Oh no!  Caught an exception while in the critical point where the BAM pipe was open: $@");
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
    $self->status_message("Postprocessing & Sanity Checking BAM file (if necessary)...");
    unless ($self->postprocess_bam_file()) {
        $self->error_message("Postprocess BAM file failed");
        die $self->error_message;
    }

    # STEP 11: COMPUTE ALIGNMENT METRICS
    $self->status_message("Computing alignment metrics...");
    $self->_compute_alignment_metrics();

    # STEP 12: PREPARE THE ALIGNMENT DIRECTORY ON NETWORK DISK
    $self->status_message("Preparing the output directory...");
    $self->status_message("Staging disk usage is " . $self->_staging_disk_usage . " KB");
    my $output_dir = $self->output_dir || $self->_prepare_output_directory;
    $self->status_message("Alignment output path is $output_dir");

    # STEP 13: PROMOTE THE DATA INTO ALIGNMENT DIRECTORY
    $self->status_message("Moving results to network disk...");
    my $product_path;
    unless($product_path= $self->_promote_data) {
        $self->error_message("Failed to de-stage data into alignment directory " . $self->error_message);
        die $self->error_message;
    }

    # STEP 14: RESIZE THE DISK
    # TODO: move this into the actual original allocation so we don't need to do this
    $self->status_message("Resizing the disk allocation...");
    if ($self->_disk_allocation) {
        my %params;
        $params{allow_reallocate_with_move} = 0;
        $params{allow_reallocate_with_move} = 1 if $self->_disk_allocation->kilobytes_requested < 10_485_760; # 10GB
        unless (eval { $self->_disk_allocation->reallocate(%params) }) {
            $self->warning_message("Failed to reallocate my disk allocation: " . $self->_disk_allocation->id);
        }
        $self->output_dir($self->_disk_allocation->absolute_path); #update if was moved
    }

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
        my @users = Genome::SoftwareResult::User->get(software_result => $child);
        $child->delete if !@users;
    }

    # find qc result and delete it
    my @qc_results = Genome::InstrumentData::AlignmentResult::Merged::BamQc->get(alignment_result_id => $self->id);
    for my $qc_result ( @qc_results ) {
        $qc_result->delete;
    }

    return $self->SUPER::delete(@_);
}

sub prepare_scratch_sam_file {
    my $self = shift;

    my $scratch_sam_file = $self->temp_scratch_directory . "/all_sequences.sam";

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

    $self->status_message("Cat-ing together: ".join("\n",@input_files). "\n to output file ".$scratch_sam_file);
    my $cat_rv = Genome::Sys->cat(input_files=>\@input_files,output_file=>$scratch_sam_file);
    if ($cat_rv ne 1) {
        $self->error_message("Error during cat of alignment sam files! Return value $cat_rv");
        die $self->error_message;
    }
    else {
        $self->status_message("Cat of sam files successful.");
    }

    if ($self->supports_streaming_to_bam) {
        my $ref_list  = $self->reference_build->full_consensus_sam_index_path($self->samtools_version);
        my $sam_cmd = sprintf("| %s view -S -b -o %s - ", Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version), $self->temp_scratch_directory . "/raw_all_sequences.bam");
        $self->status_message("Opening $sam_cmd");

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
        input  => $self->instrument_data->bam_path,
        output => $file,
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

    $self->status_message("Unpacking reads...");

    if ($self->requires_fastqs_to_align) {
        $self->status_message('Requires fastqs to align');
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
    if (defined $self->instrument_data_segment_id) {
        $self->status_message('Extract input read group bam');
        $bam_file = $self->_extract_input_read_group_bam;
        unless ($bam_file) {
            $self->error_message(sprintf('Failed to extract read group (%s) into temporary BAM.', $self->instrument_data_segment_id));
            return;
        }
    }

    # Some old imported bam does not have is_paired_end set, patch for now
    $self->status_message("Checking if this read group is paired end...");
    my $paired = $instr_data->is_paired_end;

    # Should !$paired be !defined($paired)?
    if ((!$paired || defined $self->instrument_data_segment_id) && $instr_data->can('import_format') && $instr_data->import_format eq 'bam') {
        if (defined $self->instrument_data_segment_id) {
            $self->status_message(sprintf('Inferring paired end status for "%s" segment...', $self->instrument_data_segment_id));
        } else {
            $self->status_message('Inferring paired end status...');
        }

        my $output_file = $bam_file . '.flagstat';
        unless (-s $output_file) {
            $output_file = $self->temp_scratch_directory . '/import_bam.flagstat';
            die unless $self->_create_bam_flagstat($bam_file, $output_file);
        }
        my $stats = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($output_file);
        die $self->error_message('Failed to get flagstat data on input bam: '. $bam_file) unless $stats;

        # I guess ideally we would align both paired end and fragment if this
        # were not 100% but we can't do that yet. I don't know at what point we
        # should flip $paired to true but I don't think it should be based on a
        # raw count nor based on non-zero.
        my $percent_paired = $stats->{reads_paired_in_sequencing} / $stats->{total_reads};

        # Boundaries were arbitrarily chosen, feel free to adjust as a matter
        # of policy.
        if ($percent_paired > 0.1 && $percent_paired < 0.9) {
            die $self->error_message(sprintf(
                'Trying to infer paired end status on mixed data (%.2f%% paired), not sure which to choose!',
                100 * $percent_paired));
        } elsif ($percent_paired >= 0.9) {
            $self->status_message('Setting paired end status to true.');
            $self->_is_inferred_paired_end(1);
            $paired = 1;
        } else {
            $self->status_message('Setting paired end status to false.');
            $self->_is_inferred_paired_end(0);
            $paired = 0;
        }
    }

    return ("$bam_file:1", "$bam_file:2") if $paired;
    return ("$bam_file:0");
}

sub run_aligner {
    my ($self, @inputs) = @_;

    $self->status_message("Got " . scalar(@inputs) . " input files");
    if (@inputs > 3) {
        $self->error_message("We don't support aligning with more than 3 inputs (the first 2 are treated as PE and last 1 is treated as SE)");
        die $self->error_message;
    }

    # Perform N-removal if requested
    
    if ($self->n_remove_threshold) {
        $self->status_message("Running N-remove.  Threshold is " . $self->n_remove_threshold);

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
            $self->status_message("N removal complete: Passed $passed reads & Failed $failed reads");
            if ($passed > 0) {
                push @n_removed_fastqs, $n_removed_file;
            }

            if ($input_pathname =~ m/^\/tmp/) {
                $self->status_message("Removing original file before N removal to save space: $input_pathname");
                unlink($input_pathname);
            }
        }
        if (@inputs == 1 && @n_removed_fastqs == 2) {
            $self->status_message("NOTE: An entire side of the read pairs was filtered away after n-removal.  We'll be running in SE mode from here on out.");
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
        $self->status_message("Running the aligner with the $filter_name filter.");
    }
    elsif ($self->force_fragment) {
        $self->status_message("Running the aligner in force-fragment mode.");
        @passes = map { [ $_ ] } @inputs;
    }
    elsif (@inputs == 3) {
        $self->status_message("Running aligner twice: once for PE & once for SE");
        @passes = ( [ $inputs[0], $inputs[1] ], [ $inputs[2] ] );
    }
    elsif (@inputs == 2) {
        $self->status_message("Running aligner in PE mode");
        @passes = ( \@inputs );
    }
    elsif (@inputs == 1) {
        $self->status_message("Running aligner in SE mode");
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
        $self->status_message("Aligning @$pass...");
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
           $self->status_message("Unlinking fastq file to save space now that we've aligned: $_");
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
    $self->status_message("Closing bam file...");
    $self->_sam_output_fh->flush;
    $self->_sam_output_fh->close;
    $self->_sam_output_fh(undef);

    $self->status_message("Sorting by name to do fixmate...");
    my $bam_file = $self->temp_scratch_directory . "/raw_all_sequences.bam";
    my $final_bam_file = $self->temp_staging_directory . "/all_sequences.bam";
    my $samtools = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);

    my $tmp_file = $bam_file.'.sort';
    #402653184 bytes = 3 Gb
    my $rv = system "$samtools sort -n -m 402653184 $bam_file $tmp_file";
    $self->error_message("Sort by name failed") and return if $rv or !-s $tmp_file.'.bam';
    $self->status_message("unlinking original bam file $bam_file.");
    unlink $bam_file;

    # TODO: run htseq here
    # We need a way to have down-stream steps run before their predecessor cleans-up.

    $self->status_message("Now running fixmate");
    $rv = system "$samtools fixmate $tmp_file.bam $tmp_file.fixmate";
    $self->error_message("fixmate failed") and return if $rv or !-s $tmp_file.'.fixmate';
    unlink "$tmp_file.bam";

    $self->status_message("Now putting things back in chr/pos order");
    $rv = system "$samtools sort -m 402653184 $tmp_file.fixmate $tmp_file.fix";
    $self->error_message("Sort by position failed") and return if $rv or !-s $tmp_file.'.fix.bam';

    unlink "$tmp_file.fixmate";
    unlink $bam_file;

    move "$tmp_file.fix.bam", $final_bam_file;
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

    my $bam_file    = $self->temp_staging_directory . '/all_sequences.bam';
    my $output_file = $bam_file . '.flagstat';

    #STEPS 8:  CREATE BAM.FLAGSTAT
    $self->status_message("Creating all_sequences.bam.flagstat ...");
    die unless $self->_create_bam_flagstat($bam_file, $output_file);

    #STEPS 9: VERIFY BAM IS NOT TRUNCATED BY FLAGSTAT
    $self->status_message("Verifying the bam...");
    unless ($self->_verify_bam) {
        $self->error_message('Fail to verify the bam');
        die $self->error_message;
    }

    #request by RT#62311 for submission and data integrity
    $self->status_message('Creating all_sequences.bam.md5 ...');
    unless ($self->_create_bam_md5) {
        $self->error_message('Fail to create bam md5');
        die $self->error_message;
    }

    $self->status_message("Indexing BAM file ...");
    unless($self->_create_bam_index) {
        $self->error_message('Fail to create bam md5');
        die $self->error_message;
    }
    return 1;
}

sub _use_alignment_summary_cpp { return 1; };

sub _compute_alignment_metrics {
    my $self = shift;
    my $bam = $self->temp_staging_directory . "/all_sequences.bam";

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

        Genome::Utility::Instrumentation::inc('alignment_result.read_count', $self->total_read_count);
    }

    return 1;
}

sub _create_bam_index {
    my $self = shift;
    my $bam_file    = $self->temp_staging_directory . '/all_sequences.bam';

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

    my $bam_file  = $self->temp_staging_directory . '/all_sequences.bam';
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
    $self->status_message("$check matches.");
    return 1;
}


sub _create_bam_md5 {
    my $self = shift;

    my $bam_file = $self->temp_staging_directory . '/all_sequences.bam';
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

    $self->status_message("Now de-staging data from $staging_dir into $output_dir");

    my $copy_cmd = sprintf("cp -rL %s/* %s/", $staging_dir, $output_dir);
    $self->status_message("Running cp: $copy_cmd");
    my $copy_exit_code = system($copy_cmd);

    if ($copy_exit_code != 0) {

        $self->status_message("Copy failed, attempting rsync");

        my $rsync_cmd = sprintf("rsync -avzL %s/* %s/", $staging_dir, $output_dir);

        $self->status_message("Running Rsync: $rsync_cmd");
        my $rsync_exit_code = system($rsync_cmd);

        unless ($rsync_exit_code == 0) {
            $self->error_message("Did not get a valid return from rsync, exit code was $rsync_exit_code for call $rsync_cmd.  Cleaning up and bailing out.");
            rmtree($output_dir);
            die $self->error_message;
        }
    }

    chmod 02775, $output_dir;
    for my $subdir (grep { -d $_  } glob("$output_dir/*")) {
        chmod 02775, $subdir;
    }

    # Make everything in here read-only
    for my $file (grep { -f $_  } glob("$output_dir/*")) {
        chmod 0444, $file;
    }

    $self->status_message("Files in $output_dir: \n" . join "\n", glob($output_dir . "/*"));

    return $output_dir;
}

sub _process_sam_files {
    my $self = shift;
    my $groups_input_file;

    # if a bam file is already staged at the end of _run_aligner, trust it to be correct.
    if (-e $self->temp_staging_directory . "/all_sequences.bam") {
        return 1;
    }

    my $sam_input_file = $self->temp_scratch_directory . "/all_sequences.sam";

    unless (-e $sam_input_file) {
        $self->error_message("$sam_input_file is nonexistent.  Can't convert!");
        die $self->error_message;
    }

    # things which don't produce sam natively must provide an unaligned reads file.
    my $unaligned_input_file = $self->temp_scratch_directory . "/all_sequences_unaligned.sam";

    if (-s $unaligned_input_file) {
        $self->status_message("Looks like there are unaligned reads not in the main input file.  ");
        my @input_files = ($sam_input_file, $unaligned_input_file);
        $self->status_message("Cat-ing the unaligned list $unaligned_input_file to the sam file $sam_input_file");
        my $cat_rv = Genome::Sys->cat(input_files=>[$unaligned_input_file],output_file=>$sam_input_file,append_mode=>1);
        if ($cat_rv ne 1) {
            $self->error_message("Error during cat of alignment sam files! Return value $cat_rv");
            die $self->error_message;
        } else {
            $self->status_message("Cat of sam files successful.");
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
        $self->status_message("Read group add completed, new file is $per_lane_sam_file_rg");

        $self->status_message("Removing non-read-group combined sam file: " . $sam_input_file);
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

    my $per_lane_bam_file = $self->temp_staging_directory . "/all_sequences.bam";
    my %params = (
        bam_file => $per_lane_bam_file,
        sam_file => $final_sam_file,
        keep_sam => 0,
        fix_mate => 1,
        index_bam => 1,
        ref_list => $ref_list,
        use_version => $self->samtools_version,
    );
    
    if ($self->aligner_name =~ /rtg/){
        #adukes - fix_mate screws up bitflags with rtg alignment, this is probably not the ideal spot for this...
        $params{fix_mate} = 0;
    }

    my $to_bam = Genome::Model::Tools::Sam::SamToBam->create(
        %params,
    );

    unless($to_bam->execute) {
        $self->error_message("There was an error converting the Sam file $final_sam_file to $per_lane_bam_file.");
        die $self->error_message;
    }

    $self->status_message("Conversion successful.  File is: $per_lane_bam_file");
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

    #my $inputs_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($subclass, %is_input);
    #my $params_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($subclass, %is_param);

    my %software_result_params = (#software_version=>$params_bx->value_for('aligner_version'),
        #params_id=>$params_bx->id,
        #inputs_id=>$inputs_bx->id,
        subclass_name=>$subclass);

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

    my $scratch_basedir = sprintf("scratch-%s-%s-%s", $hostname, $user, $$);
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
    return "info_alignments";
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
            $self->status_message('No special filter for this input');
        }
        elsif ($self->filter_name eq 'forward-only') {
            # forward reads only
            $self->status_message('forward-only filter applied for this input');
        }
        elsif ($self->filter_name eq 'reverse-only') {
            # reverse reads only
            $self->status_message('reverse-only filter applied for this input');
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
    }
    return @input_fastq_pathnames;
}


sub input_bfq_filenames {
    my $self = shift;
    my @input_fastq_pathnames = @_;

    my @input_bfq_pathnames;
    if ($self->_input_bfq_pathnames) {
        @input_bfq_pathnames = @{$self->_input_bfq_pathnames};
        for my $input_bfq (@input_bfq_pathnames) {
            unless (-s $input_bfq) {
                $self->error_message('Missing or zero size sanger bfq file: '. $input_bfq);
                die $self->error_message;
            }
        }
    }
    else {
        my $counter = 0;
        for my $input_fastq_pathname (@input_fastq_pathnames) {
            my $input_bfq_pathname = Genome::Sys->create_temp_file_path('sanger-bfq-'. $counter++);
            #Do we need remove sanger fastq here ?
            unless (Genome::Model::Tools::Maq::Fastq2bfq->execute(
                    fastq_file => $input_fastq_pathname,
                    bfq_file   => $input_bfq_pathname,
                )) {
                $self->error_message('Failed to execute fastq2bfq quality conversion.');
                die $self->error_message;
            }
            unless (-s $input_bfq_pathname) {
                $self->error_message('Failed to validate the conversion of sanger fastq file '. $input_fastq_pathname .' to sanger bfq.');
                die $self->error_message;
            }
            push @input_bfq_pathnames, $input_bfq_pathname;
        }
        $self->_input_bfq_pathnames(\@input_bfq_pathnames);
    }
    return @input_bfq_pathnames;
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
    $self->status_message("Sample id: ".$self->instrument_data->sample_id);
    my $sample = Genome::Sample->get($self->instrument_data->sample_id);
    if ( defined($sample) ) {
        $species =  $sample->species_name;
        if (!$species || $species eq ""  ) {
            $species = "unknown";
        }
    }

    $self->status_message("Species from alignment: ".$species);

    my $ref_build = $self->reference_build;
    my $seq_dict = $ref_build->get_sequence_dictionary("sam",$species,$self->picard_version);
    return $seq_dict;
}

sub construct_groups_file {
    my $self = shift;
    my $output_file = shift || $self->temp_scratch_directory . "/groups.sam";

    my $aligner_command_line = $self->aligner_params_for_sam_header;
    my $instr_data = $self->instrument_data;

    my $insert_size_for_header;
    if ($instr_data->can('resolve_median_insert_size') && $instr_data->resolve_median_insert_size) {
        $insert_size_for_header= $instr_data->resolve_median_insert_size;
    }
    else {
        $insert_size_for_header = 0;
    }

    my $paired = defined $self->_is_inferred_paired_end ? $self->_is_inferred_paired_end : $instr_data->is_paired_end;
    my $description_for_header = $paired ? "paired end" : "fragment";

    # build the header
    my $id_tag       = $self->read_and_platform_group_tag_id;
    my $pu_tag       = sprintf("%s.%s", $instr_data->run_identifier, $instr_data->subset_name);
    my $lib_tag      = $instr_data->library_name;
    my $date_run_tag = $instr_data->run_start_date_formatted;
    my $sample_tag   = $instr_data->sample_name;
    my $aligner_version_tag = $self->aligner_version;
    my $aligner_cmd  =  $aligner_command_line;

    my $platform = $instr_data->sequencing_platform;
    $platform = ($platform eq 'solexa' ? 'illumina' : $platform);

    #@RG     ID:2723755796   PL:illumina     PU:30945.1      LB:H_GP-0124n-lib1      PI:0    DS:paired end   DT:2008-10-03   SM:H_GP-0124n   CN:WUGSC
    #@PG     ID:0    VN:0.4.9        CL:bwa aln -t4
    my $rg_tag = "\@RG\tID:$id_tag\tPL:$platform\tPU:$pu_tag\tLB:$lib_tag\tPI:$insert_size_for_header\tDS:$description_for_header\tDT:$date_run_tag\tSM:$sample_tag\tCN:WUGSC\n";
    my $pg_tag = "\@PG\tID:$id_tag\tVN:$aligner_version_tag\tCL:$aligner_cmd\n";

    $self->status_message("RG: $rg_tag");
    $self->status_message("PG: $pg_tag");

    my $header_groups_fh = IO::File->new(">>".$output_file) || die "failed opening groups file for writing";
    print $header_groups_fh $rg_tag;
    print $header_groups_fh $pg_tag;
    $header_groups_fh->close;

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

sub alignment_bam_file_paths {
    my $self = shift;

    return glob($self->output_dir . "/*.bam");
}

sub requires_read_group_addition {
    return 1;
}

sub supports_streaming_to_bam {
    0;
}

sub accepts_bam_input {
    0;
}

sub aligner_params_required_for_index {
    0;
}

sub get_reference_sequence_index {
    my $self = shift;
    my $build = shift || $self->reference_build;
    my @overrides = @_;
    my $index = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_with_lock(aligner_name=>$self->aligner_name, aligner_version=>$self->aligner_version, aligner_params=>$self->aligner_params, reference_build=>$build, @overrides);

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
        $self->status_message("Calculated upper bound on insert size is undef or less than 0, defaulting to $up_default");
        $upper = $up_default;
    }
    if (!defined $lower || not $median || $lower < 0 || $lower > $upper) {
        $self->status_message("Calculated lower bound on insert size is undef or invalid, defaulting to $low_default");
        $lower = $low_default;
    }
    return ($lower, $upper);
}

1;

