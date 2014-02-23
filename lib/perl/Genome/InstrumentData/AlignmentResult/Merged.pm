package Genome::InstrumentData::AlignmentResult::Merged;

use strict;
use warnings;

use Sys::Hostname;
use File::stat;
use File::Path 'rmtree';
use List::MoreUtils qw{ uniq };

use Genome;

use Genome::Utility::Text; #quiet warning about deprecated use of autoload

class Genome::InstrumentData::AlignmentResult::Merged {
    is => 'Genome::InstrumentData::AlignedBamResult',
    has => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            calculate => q{
                return Genome::InstrumentData->get([$self->instrument_data_id]);
            }
        },
        reference_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            id_by => 'reference_build_id',
        },
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            id_by => 'annotation_build_id',
            is_optional => 1,
        },
        reference_name => {
            via => 'reference_build',
            to => 'name',
            is_mutable => 0,
            is_optional => 1
        },
        _disk_allocation => {
            is => 'Genome::Disk::Allocation',
            is_optional => 1,
            is_many => 1,
            reverse_as => 'owner'
        },
    ],

    has_param => [
        #the parameters from individual alignments that are used to produce this result)
        #Note: samtools_version is used in this module, so if it is removed from AlignmentResult will need to be explicitly added below
        (map
            {$_->property_name => { is => $_->data_type, doc => $_->doc, is_optional => $_->is_optional} }
            Genome::InstrumentData::AlignmentResult->__meta__->_legacy_properties(via => 'params')
        ),
        merger_name => {
            is => 'Text',
            doc => 'The name of the merge program to use (e.g. "samtools")',
        },
        merger_params => {
            is => 'Text',
            is_optional => 1,
            doc => 'Additional parameters to pass to the merge program',
        },
        merger_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'Version of the merge program to use',
        },
        duplication_handler_name => {
            is => 'Text',
            is_optional => 1,
            doc => 'The name of the program to use for marking or removing duplicate reads',
        },
        duplication_handler_params => {
            is => 'Text',
            is_optional => 1,
            doc => 'Additional parameters to pass to the dpulication handler',
        },
        duplication_handler_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'Version of the duplication handler to use',
        },
        filter_name => {
            is => 'Text',
            is_many => 1,
            doc => 'Filters for any of the individual alignments (if applicable)',
        },
    ],

    has_input => [
        reference_build_id => {
            is => 'Number',
            doc => 'the reference to use by id',
        },
        annotation_build_id => {
            is => 'Number',
            doc => 'the annotation to use by id',
            is_optional => 1,
        },
        instrument_data_id => {
            is => 'Number',
            doc => 'the local database ids of the instrument data (reads) for this merged alignment',
            is_many => 1,
        },
        instrument_data_segment => {
            is => 'Text',
            is_many => 1,
            doc => 'Segments for individual alignments (if applicable)',
        },
    ],
    has_transient_optional => [
        temp_staging_directory  => {
            is => 'Text',
            doc => 'A directory to use for staging the alignment data while working.',
        },
        temp_scratch_directory  => {
            is => 'Text',
            doc => 'A directory for working files not intended to be kept',
        },
        bam_file => {
          is => 'Text',
          is_output => 1,
          via => '__self__',
          to => 'merged_alignment_bam_path',
        },
    ],

    has_calculated => [
        _final_bam_file => {
            is => 'Text', calculate_from => ['temp_staging_directory', 'id',],
            calculate => q{ return join('/', $temp_staging_directory, $id . '.bam'); },
        },
        merged_alignment_bam_path => {
            is => 'Text', calculate_from => ['output_dir', 'id'],
            calculate => q{ return join('/', $output_dir, $id . '.bam'); }
        },
        merged_alignment_bam_flagstat => {
            is => 'Text', calculate_from => ['merged_alignment_bam_path'],
            calculate => q{ return $merged_alignment_bam_path . '.flagstat' }
        },
    ],

    doc => 'Represents merged (and possibly deduplicated) instrument data',
};

sub create {
    my $class = shift;

    #This will do some locking and the like for us.
    my $self = $class->SUPER::create(@_);
    return unless ($self);

    my $rv = eval {

        #TODO In a future version collect relevant alignments from other merged alignment results when available
        $self->debug_message('Collecting alignments for merger...');
        my @alignments = $self->collect_individual_alignments;

        $self->debug_message('Preparing directories...');
        $self->_prepare_output_directory; #This gets a disk allocation
        my @tmp_dirs = $self->_prepare_working_directories(\@alignments); #need to keep these in scope while in use

        my $bams_per_library = {};
        my $libraries = {};

        my @bams_for_final_merge;
        if(defined $self->duplication_handler_name) {
            #handle duplicates on a per-library basis
            for my $alignment (@alignments) {
                my $library = $alignment->instrument_data->library;

                push @{ $bams_per_library->{$library->id} }, $alignment->alignment_bam_file_paths;
                $libraries->{$library->id} = $library;
            }

            for my $library_id (keys %$bams_per_library) {
                my $library = $libraries->{$library_id};
                my $sanitized_library_name = Genome::Utility::Text::sanitize_string_for_filesystem(join('-',$library->name, $library->id));
                my $library_merged_bam = join('/', $self->temp_scratch_directory, $sanitized_library_name . '.bam');
                my $per_library_post_duplication_bam = join('/', $self->temp_scratch_directory, $sanitized_library_name . '-post_dup.bam');

                $self->debug_message('Merging alignments for library ' . $library->__display_name__ . '...');
                $self->merge_alignments($bams_per_library->{$library_id}, $library_merged_bam);
                $self->debug_message('Handling duplicates for library' . $library->__display_name__ . '...');
                $self->handle_duplicates($library_merged_bam, $per_library_post_duplication_bam);
                push @bams_for_final_merge, $per_library_post_duplication_bam;
            }
        } else {
            #just collect the BAMs for a merge
            for my $alignment (@alignments) {
                push @bams_for_final_merge, $alignment->alignment_bam_file_paths;
            }
        }

        $self->debug_message('Merging per-library bams...');
        my $final_bam = $self->_final_bam_file;
        $self->merge_alignments(\@bams_for_final_merge, $final_bam);

        $self->debug_message("Indexing the final BAM file...");
        my $index_cmd = Genome::Model::Tools::Sam::IndexBam->create(
            bam_file    => $final_bam,
            use_version => $self->samtools_version,
        );
        my $index_cmd_rv = $index_cmd->execute;

        if($index_cmd_rv ne 1) {
            #not failing here because this is not a critical error.  this can be regenerated manually if needed.
            $self->warning_message('Failed to create bam index for ' . $final_bam);
        }

        $self->create_bam_md5($final_bam);

        $self->_promote_validated_data;
        @tmp_dirs = (); #clear out temp directories

        for my $alignment (@alignments) {
            $alignment->add_user(user => $self, label => 'uses');
        }

        return 1;
    };
    if(my $error = $@) {
        $self->_cleanup;
        die $error;
    } elsif ($rv ne 1) {
        $self->error_message('Unexpected return value: ' . $rv);
        $self->_cleanup;
        die $self->error_message;
    }

    $self->debug_message("Resizing the disk allocation...");
    if ($self->_disk_allocation) {
        unless ($self->_disk_allocation->reallocate) {
            $self->warning_message("Failed to reallocate disk allocation: " . $self->_disk_allocation->id);
        }
    }

    

    $self->debug_message('All processes completed.');

    return $self;
}

sub collect_individual_alignments {
    my $self = shift;

    my @instrument_data = $self->instrument_data();
    my %params;
    for my $property (Genome::InstrumentData::AlignmentResult->__meta__->_legacy_properties(via => 'params')) {
        my $property_name = $property->property_name;
        next if grep($_ eq $property_name, ('filter_name', 'test_name')); #these are handled below
        $params{$property_name} = $self->$property_name;
    }

    my @filters = $self->filter_name;
    my $filters = {};
    for my $filter_string (@filters) {
        my ($id, $filter) = split(':', $filter_string);
        $filters->{$id} = $filter;
    }

    my @segments = $self->instrument_data_segment;

    my $segments = {};

    for my $segment_string (@segments) {
        my ($id, $segment_id, $segment_type) = split(':', $segment_string);
        $segments->{$id}{$segment_type} ||= [];
        push @{$segments->{$id}{$segment_type} }, $segment_id;
    }

    my @alignments;
    my @not_found;
    for my $i (@instrument_data) {
        my @segment_params;
        if($segments->{$i->id}) {
            for my $type (keys %{ $segments->{$i->id} }) {
                my $segment_ids = $segments->{$i->id}{$type};
                for my $segment_id (@$segment_ids) {
                    push @segment_params, {
                        'instrument_data_segment_type' => $type,
                        'instrument_data_segment_id' => $segment_id,
                    };
                }
            }
        } else {
            push @segment_params, {
                'instrument_data_segment_type' => undef,
                'instrument_data_segment_id' => undef
            };
        }

        for my $segment_param (@segment_params) {
            my $alignment = Genome::InstrumentData::AlignmentResult->get_with_lock(
                %params,
                reference_build_id => $self->reference_build_id,
                annotation_build_id => ($self->annotation_build_id || undef),
                instrument_data_id => $i->id,
                filter_name => ($filters->{$i->id} || undef),
                test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
                %$segment_param,
            );

            if($alignment) {
                push @alignments, $alignment;
            } else {
                push @not_found, $i;
            }
        }
    }

    if(scalar @not_found) {
        $self->error_message(
            'Failed to find individual alignments for all instrument_data. Missing: ' .
            join(', ', map($_->__display_name__, @not_found) )
        );
        die $self->error_message;
    }

    return @alignments;
}

sub required_rusage {
    return ''; #FIXME This needs to be filled in
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

sub resolve_allocation_kilobytes_requested {
    my $self = shift;
    my @alignments = $self->collect_individual_alignments;
    return $self->estimated_kb_usage(\@alignments);
}

sub resolve_allocation_subdirectory {
    return $_[0]->resolve_alignment_subdirectory;
}

sub estimated_kb_usage {
    my $self = shift;
    my $alignments = shift;

    my @bams;
    for my $alignment (@$alignments) {
        my @aln_bams = $alignment->alignment_bam_file_paths;
        unless (@aln_bams) {
            $self->debug_message("alignment $alignment has no bams at " . $alignment->output_dir);
        }
        push @bams, @aln_bams;
    }
    my $total_size;
    
    unless (@bams) {
        die "No bams?";
    }

    for (@bams) {
        my $size = stat($_)->size;
        $self->debug_message("BAM has size: " . $size);
        $total_size += $size;
    }

    #take the total size plus a 10% safety margin
    # 2x total size; full build merged bam, full build deduped bam
    $total_size = sprintf("%.0f", ($total_size/1024)*1.1); 
    $total_size = ($total_size * 2);

    return $total_size;
}

sub resolve_alignment_subdirectory {
    my $self = shift;

    my $hostname = hostname;
    my $user = $ENV{'USER'};
    my $base_dir = sprintf("merged-alignment-%s-%s-%s-%s", $hostname, $user, $$, $self->id);
    # TODO: the first subdir is actually specified by the disk management system.
    my $directory = join('/', 'build_merged_alignments', $base_dir);
    return $directory;
}

sub _prepare_working_directories {
    my $self = shift;

    return $self->temp_staging_directory if $self->temp_staging_directory;

    my $output_dir = $self->output_dir;

    #file sizes are so large that /tmp/ would be exhausted--stage files to the allocation itself instead
    my $staging_tempdir = File::Temp->newdir( 
        "staging-XXXXX",
        DIR     => $output_dir, 
        CLEANUP => 1,
    );

    my $scratch_tempdir = File::Temp->newdir( 
        "scratch-XXXXX",
        DIR     => $output_dir, 
        CLEANUP => 1,
    );

    # fix permissions on this temp dir so others can clean it up later if need be
    chmod(0775,$staging_tempdir);
    chmod(0775,$scratch_tempdir);

    $self->temp_staging_directory($staging_tempdir->dirname);
    $self->temp_scratch_directory($scratch_tempdir->dirname);

    return ($staging_tempdir, $scratch_tempdir);
}

sub _promote_validated_data {
    my $self = shift;

    my $staging_dir = $self->temp_staging_directory;
    my $output_dir  = $self->output_dir;

    $self->debug_message("Now de-staging data from $staging_dir into $output_dir");

    for my $staged_file (glob("$staging_dir/*")) {
        my $destination = $staged_file;
        $destination =~ s/$staging_dir/$output_dir/;
        rename($staged_file, $destination);
    }

    chmod 02775, $output_dir;
    for my $subdir (grep { -d $_  } glob("$output_dir/*")) {
        chmod 02775, $subdir;
    }

    # Make everything in here read-only
    for my $file (grep { -f $_  } glob("$output_dir/*")) {
        chmod 0444, $file;
    }

    $self->debug_message("Files in $output_dir: \n" . join "\n", glob($output_dir . "/*"));

    return $output_dir;
}

sub merge_alignments {
    my $self = shift;

    my $input_bams = shift;
    my $output_path = shift;

    my $rv = $self->_run_merger($input_bams, $output_path);
    unless($rv) {
        die $self->error_message('Failed to merge.');
    }

    $self->verify_result($input_bams, [$output_path]);
}

sub handle_duplicates {
    my $self = shift;

    my $input_bam = shift;
    my $output_path = shift;

    my $rv = $self->_run_duplication_handler($input_bam, $output_path);
    unless($rv) {
        die $self->error_message('Failed to handle duplicates.');
    }

    $self->verify_result([$input_bam], [$output_path]);    
}

sub verify_result {
    my $self = shift;

    my $input_bams = shift;
    my $output_bams = shift;

    my $input_total = 0;
    my $output_total = 0;

    for my $bam_file (@$input_bams) {
        my $bam_total = $self->_bam_flagstat_total($bam_file);
        unless($bam_total) {
            $self->error_message('Could not verify.  Error in ' . $bam_file);
            die $self->error_message;
        }
        $input_total += $bam_total;
    }

    for my $bam_file (@$output_bams) {
        my $bam_total = $self->_bam_flagstat_total($bam_file);
        unless($bam_total) {
            $self->error_message('Could not verify.  Error in ' . $bam_file);
            die $self->error_message;
        }
        $output_total += $bam_total;
    }

    $self->debug_message('input count: ' . $input_total);
    $self->debug_message('output count: ' . $output_total);

    if($input_total eq $output_total) {
        $self->debug_message('Counts match.  Verification OK.');
        return 1;
    } else {
        $self->error_message('Counts do not match. Verification failed.');
        die $self->error_message;
    }
}

sub _bam_flagstat_total {
    my $self      = shift;
    my $bam_file  = shift;
    my $flag_file = $bam_file . '.flagstat';

    unless(Genome::Sys->check_for_path_existence($bam_file)) {
        $self->error_message('BAM file not found: ' . $bam_file);
        die $self->error_message;
    }

    unless (-s $flag_file) {
        my $cmd = Genome::Model::Tools::Sam::Flagstat->create(
            use_version    => $self->samtools_version,
            bam_file       => $bam_file,
            output_file    => $flag_file,
            include_stderr => 1,
        );

        unless($cmd and $cmd->execute) {
            $self->error_message("Fail to create or execute flagstat command on bam file: $bam_file");
            return;
        }
    }
    my $flagstat_data = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flag_file);

    unless($flagstat_data) {
        $self->error_message('No output from samtools flagstat');
        return;
    }

    if(exists $flagstat_data->{errors}) {
        for my $error (@{ $flagstat_data->{errors} }) {
            if($error =~ m/Truncated file/) {
                $self->error_message('Flagstat output for ' . $bam_file . ' indicates possible truncation.');
                return;
            }
        }
    }
    my $total = $flagstat_data->{total_reads};

    $self->debug_message('flagstat for ' . $bam_file . ' reports ' . $total . ' in total');    
    return $total;
}

sub _run_merger {
    my $self = shift;

    my $input_bams = shift;
    my $output_path = shift;

    my $comment = $self->_resolve_comment;

    my $merger_module = $self->_resolve_merger_module;
    my $merger = $merger_module->create(
        input_bams => $input_bams,
        output_path => $output_path,
        parameters => $self->merger_params,
        version => $self->merger_version,
        scratch_directory => $self->temp_scratch_directory,
        samtools_version => $self->samtools_version,
        ($comment? (include_comment => $comment) : ()),
    );

    unless($merger->execute()) {
        $self->error_message('Failed to execute merger.');
        return;
    }

    return 1;
}

sub _resolve_comment {
    my $self = shift;

    my @i = $self->instrument_data;
    my $includes_imported = grep { $_->isa('Genome::InstrumentData::Imported') } @i;

    return if $includes_imported;

    return 'This sequence data was produced in Saint Louis, Missouri, USA.';
}

sub _resolve_merger_module {
    my $self = shift;

    return join('::', 'Genome::InstrumentData::Command::AlignmentResult::Merged::Merger', Genome::Utility::Text::string_to_camel_case($self->merger_name));
}

sub _resolve_duplication_handler_module {
    my $self = shift;

    return join('::', 'Genome::InstrumentData::Command::AlignmentResult::Merged::DuplicationHandler', Genome::Utility::Text::string_to_camel_case($self->duplication_handler_name));
}

sub _run_duplication_handler {
    my $self = shift;

    my $input_bam = shift;
    my $output_path = shift;

    my $comment = $self->_resolve_comment;

    my $duplication_handler_module = $self->_resolve_duplication_handler_module;
    my $duplication_handler = $duplication_handler_module->create(
        input_bam => $input_bam,
        output_path => $output_path,
        parameters => $self->duplication_handler_params,
        version => $self->duplication_handler_version,
        scratch_directory => $self->temp_scratch_directory,
        log_file => $self->_resolve_duplication_log_name($input_bam),
        metrics_file => $self->_resolve_duplication_metrics_name($input_bam),
        ($comment? (include_comment => $comment) : ()),
    );

    unless($duplication_handler->execute()) {
        $self->error_message('Failed to execute duplication handler.');
        return;
    }

    return 1;
}

sub _resolve_duplication_log_name {
    my $self = shift;
    my $bam_path = shift;

    my $scratch_dir = $self->temp_scratch_directory;
    my $staging_dir = $self->temp_staging_directory;

    $bam_path =~ s/$scratch_dir/$staging_dir/;
    $bam_path =~ s/.bam$//;

    return $bam_path . '.log';
}

sub _resolve_duplication_metrics_name {
    my $self = shift;
    my $bam_path = shift;

    my $scratch_dir = $self->temp_scratch_directory;
    my $staging_dir = $self->temp_staging_directory;

    $bam_path =~ s/$scratch_dir/$staging_dir/;
    $bam_path =~ s/.bam$//;

    return $bam_path . '.metrics';    
}

sub create_bam_md5 {
    my $self = shift;

    my $bam_file = shift;
    my $md5_file = $bam_file.'.md5';
    my $cmd = "md5sum $bam_file > $md5_file";

    $self->debug_message("Creating md5 file for the BAM file...");

    Genome::Sys->shellcmd(
        cmd                        => $cmd, 
        input_files                => [$bam_file],
        output_files               => [$md5_file],
        skip_if_output_is_present  => 0,
    ); 

    return 1;
}

sub _cleanup {
    my $self = shift;

    return unless $self->_disk_allocation;

    $self->debug_message('Now deleting allocation with owner_id = ' . $self->id);
    my $allocation = $self->_disk_allocation;
    if ($allocation) {
        my $path = $allocation->absolute_path;
        unless (rmtree($path)) {
            $self->error_message("could not rmtree $path");
            return;
       }
       $allocation->deallocate; 
    }
}

sub bowtie_version { return shift->scalar_property_from_underlying_alignment_results('bowtie_version'); }

sub scalar_property_from_underlying_alignment_results {
    my ($self, $property) = @_;
    my @properties;
    my @bad_alignments;
    for ($self->collect_individual_alignments) {
        if ($_->can($property)) {
            push @properties, $_->$property;
        } else {
          push @bad_alignments, $_->id;
        }
    }

    die("The following alignment results do not have the property $property: " . join(',',  @bad_alignments)) if @bad_alignments;

    @properties = uniq @properties;
    if ( scalar(@properties) > 1 ) {
        my @sorted_properties = sort @properties;
        warn("Merged alignment results contains mismatched $property values, picking " . $sorted_properties[0]);
        return $sorted_properties[0];
    } elsif (scalar(@properties) ) {
        return $properties[0];
    } else {
        die("Could not determine $property value for your alignments!");
    }
}


1;
