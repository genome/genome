package Genome::FeatureList;

use strict;
use warnings;

use Genome;
use Genome::Utility::File::Mode qw(mode);
use Memoize qw();
use Text::ParseWords;

class Genome::FeatureList {
    table_name => 'model.feature_list',
    id_by => [
        id => { is => 'Text', len => 64 },
    ],
    has => [
        name => { is => 'Text', len => 1500 },
        format => {
            is => 'Text',
            len => 64,
            valid_values => [ "1-based", "true-BED", "multi-tracked", "multi-tracked 1-based", "unknown" ],
            doc => 'Indicates whether the file follows the BED spec.',
        },
        file_content_hash => {
            is => 'Text',
            len => 32,
            doc => 'MD5 of the BED file (to ensure integrity)',
        },
        is_multitracked => {
            is => 'Boolean',
            calculate_from => 'format',
            calculate => q( return scalar ($format =~ /multi-tracked/); ),
        },
        is_1_based => {
            is => 'Boolean',
            calculate_from => 'format',
            calculate => q( return scalar ($format =~ /1-based/); ),
        },
    ],
    has_optional => [
        source => {
            is => 'Text',
            len => 64,
            doc => 'Provenance of this feature list. (e.g. Agilent)',
        },
        reference_id => {
            is => 'Text',
            len => 32,
            doc => 'ID of the reference sequence build for which the features apply',
        },
        reference => {
            is => 'Genome::Model::Build::ReferenceSequence',
            id_by => 'reference_id',
        },
        reference_name => { via => 'reference', to => 'name' },
        subject_id => {
            is => 'Text',
            len => 32,
            doc => 'ID of the subject to which the features are relevant',
        },
        subject => { is => 'Genome::Model::Build', id_by => 'subject_id' },
        disk_allocation => {
            is => 'Genome::Disk::Allocation',
            reverse_as => 'owner',
            is_many => 1,
        },
        output_dir => {
            is => 'Text',
            via => 'disk_allocation',
            to => 'absolute_path',
        },
        file_path => {
            is => 'Text',
            calculate_from => 'disk_allocation',
            calculate => q(
                  if($disk_allocation) {
                    my $directory = $disk_allocation->absolute_path;
                       return join('/', $directory, $self->id . '.bed');
                  }
               ),
        },
        content_type => {
            is => 'Text',
            len => 255,
            valid_values => [ "exome", "targeted", "validation", "roi", undef ],
            doc => 'The type of list (used for determining automated processing)',
        },
        description => {
            is => 'Text',
            len => 255,
            doc => 'General description of the BED file',
        },
    ],
    has_optional_transient => [
        _merged_bed_file_path => {
            is => 'Text',
            doc => 'The path to the temporary dumped copy of the merged post-processed BED file',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
    doc => 'A feature-list is, generically, a set of coördinates for some reference',
};

sub __display_name__ {
    my $self = shift;

    return $self->name . ' (' . $self->id . ')';
}

sub create {
    my $class = shift;
    my %params = @_;
    my $file = delete $params{file_path};

    my $self = $class->SUPER::create(%params);

    if($file and Genome::Sys->check_for_path_existence($file)) {
        my $allocation = Genome::Disk::Allocation->allocate(
            disk_group_name => $ENV{GENOME_DISK_GROUP_REFERENCES},
            allocation_path => 'feature_list/' . $self->id,
            kilobytes_requested => ( int((-s $file) / 1024) + 1),
            owner_class_name => $self->class,
            owner_id => $self->id
        );

        unless ($allocation) {
            $self->delete;
            return;
        }

        my $retval = eval {
            Genome::Sys->copy_file($file, $self->file_path);
        };
        if($@ or not $retval) {
            $self->error_message('Copy failed: ' . ($@ || 'returned' . $retval) );
            $self->delete;
            return;
        }
    }

    #TODO If this is in __errors__, gets called too soon--still would be nice to have
    unless($self->verify_file_md5) {
        $self->error_message('MD5 of copy does not match supplied value!');
        $self->delete;
        return;
    }

    my $result = eval{
        mode($self->file_path)->rm_all_writable;
    };
    if($@ or !$result){
        $self->error_message("Could not modify file permissions for: ".$self->file_path);
    }

    #This wouldn't be necessary, but the differing block sizes between disks sometimes make the estimate
    #above off by a small amount.
    $self->disk_allocation->reallocate;

    return $self;
}

sub delete {
    my $self = shift;

    #If we commit the delete, need to get rid of the allocation.
    my $upon_delete_callback = $self->_cleanup_allocation_sub;
    $self->ghost_class->add_observer(aspect=>'commit', callback=>$upon_delete_callback);

    return $self->SUPER::delete(@_);
}

sub _cleanup_allocation_sub {
    my $self = shift;

    my $id = $self->id;
    my $class_name = $self->class;
    return sub {
        print "Now deleting allocation with owner_id = $id\n";
        my $allocation = Genome::Disk::Allocation->get(owner_id => $id, owner_class_name => $class_name);
        if ($allocation) {
            $allocation->deallocate;
        }
    };
}

sub verify_file_md5 {
    my $self = shift;

    my $bed_file = $self->file_path;

    my $md5_sum = Genome::Sys->md5sum($bed_file);
    if($md5_sum eq $self->file_content_hash) {
        return $md5_sum;
    } else {
        return;
    }
}

sub get_one_based_file {
    my $self = shift;

    if ($self->format eq 'unknown') {
        $self->error_message("Cannot convert format of BED file with unknown format");
        die $self->error_message;
    }
    
    my $new_format;
    my $new_description = $self->description;
    if ($self->is_1_based) {
        return $self->file_path;
    }
    else {
        return $self->transform_zero_to_one_based($self->file_path, $self->is_multitracked);
    }
}
sub transform_zero_to_one_based {
    my $class = shift;
    my $file = shift;
    my $is_multitracked = shift;
    my $bed_file_content;

    my $fh = Genome::Sys->open_file_for_reading($file);
    my $line_no = 0;
    while(my $line = <$fh>) {
        chomp($line);
        $line_no++;

        if($is_multitracked) {
            if ($line =~ /^track/) {
                $bed_file_content .= "$line\n";
                next;
            }
        }

        my @entry = $class->_parse_entry($line, $line_no);
        $entry[1]++;
        $bed_file_content .= join("\t",@entry) ."\n";
    }
    my $temp_file = Genome::Sys->create_temp_file_path;
    Genome::Sys->write_file($temp_file, $bed_file_content);
    return $temp_file;
}

#The raw "BED" file we import will be in one many BED-like formats.
#The output of this method is the standardized "true-BED" representation
sub processed_bed_file_content {
    my $self = shift;
    my %args = @_;

    my $track_name = delete($args{track_name});
    unless (defined($track_name)) {
        $track_name = 'target_region';
    }
    my $short_name = delete($args{short_name});
    unless (defined($short_name)) {
        $short_name = 1;
    }
    if($self->format eq 'unknown'){
        $self->error_message('Cannot process BED file with unknown format');
        die $self->error_message;
    }

    my $file = $self->file_path;
    unless($self->verify_file_md5) {
        $self->error_message('MD5 mismatch! BED file modified or corrupted?');
        die $self->error_message;
    }

    my $fh = Genome::Sys->open_file_for_reading($file);

    my $print = 1;
    my $bed_file_content;
    my $name_counter = 0;
    my $line_no = 0;
    while(my $line = <$fh>) {
        chomp($line);
        $line_no++;

        if(!$bed_file_content && index($line, 'browser') == 0) {
            next;
        }

        if($self->is_multitracked) {
            if($line =~ /^track /) {
                my %track_attrs = $self->parse_track_definition($line);
                my $name = $self->_standardize_track_name($track_attrs{name});
                if ($name eq 'probes') {
                    if ($track_name eq 'tiled_region') {
                        $print = 1;
                    } else {
                        $print = 0;
                    }
                    next;
                } elsif ($name eq 'targets') {
                    if ($track_name eq 'target_region') {
                        $print = 1;
                    } else {
                        $print = 0;
                    }
                    next;
                } else {
                    $self->warning_message("Unknown track name (line $line_no '$line'). Including regions.");
                    $print = 1;
                    next;
                }
            }
        }

        if ($print) {
            my @entry = $self->_parse_entry($line, $line_no);
            if (!defined($entry[3])) {
                $entry[3] = $entry[0] .':'. $entry[1] .'-'. $entry[2];
            }
            $entry[0] =~ s/chr//g;
            if ($entry[0] =~ /random/) { next; }

            # Correct for 1-based start positions in imported BED files,
            # unless at zero already(which means we shouldn't be correcting the position anyway...)
            if ($self->is_1_based) {
                if($entry[1] == 0) {
                    $self->error_message('BED file was imported as 1-based but contains a 0 in the start position!');
                    die($self->error_message);
                }
                $entry[1]--;
            }
            #Bio::DB::Sam slows down dramatically when large names are used, so just number the regions sequentially
            if ($short_name) {
                $entry[3] = 'r' . $name_counter++;
            }
            $bed_file_content .= join("\t",@entry[0..3]) ."\n";
        }
    }
    return $bed_file_content;
}

sub _parse_entry {
    my ($self, $line, $line_no) = @_;

    my @entry = split("\t", $line);
    unless (scalar(@entry) >= 3) {
        die $self->error_message("At least three fields are required in BED format files.  Error with line $line_no:\n$line\n\n");
    }

    return @entry;
}

sub processed_bed_file {
    my $self = shift;
    my %args = @_;
    if($self->format eq 'unknown'){
        $self->error_message('Cannot process BED file with unknown format');
        die $self->error_message;
    }

    my $content = $self->processed_bed_file_content(%args);
    my $temp_file = Genome::Sys->create_temp_file_path();
    Genome::Sys->write_file($temp_file, $content);
    return $temp_file;
}
Memoize::memoize('processed_bed_file');

sub generate_merged_bed_file {
    my $self = shift;
    my %args = @_;

    my $processed_bed_file = $self->processed_bed_file(%args);
    my $output_file = Genome::Sys->create_temp_file_path( $self->id . '.merged.bed' );

    my %merge_params = (
        input_file => $processed_bed_file,
        output_file => $output_file,
        report_names => 1,
        #All files should have zero-based start postitions at this point
        maximum_distance => 0,
    );

    my $merge_command = Genome::Model::Tools::BedTools::Merge->create(%merge_params);
    unless($merge_command) {
        $self->error_message('Failed to create merge command.');
        die $self->error_message;
    }
    unless ($merge_command->execute) {
        $self->error_message('Failed to merge BED file with params '. Data::Dumper::Dumper(%merge_params) . ' ' . $merge_command->error_message);
        die $self->error_message;
    }

    return $output_file;
}

sub merged_bed_file {
    my $self = shift;
    my %args = @_;
    if ($self->format eq 'unknown'){
        $self->error_message('Cannot merge BED file with unknown format');
        die $self->error_message;
    }

    unless($self->_merged_bed_file_path) {
        my $temp_file = $self->generate_merged_bed_file(%args);

        $self->_merged_bed_file_path($temp_file);
    }

    return $self->_merged_bed_file_path;
}

sub generate_converted_bed_file {
    my $self = shift;
    my %args = @_;

    my $merge = delete($args{merge});
    my $reference = delete($args{reference});

    my $original_file_path;
    if ($merge) {
        $original_file_path = $self->merged_bed_file(%args);
    } else {
        $original_file_path = $self->processed_bed_file(%args);
    }
    my $original_md5 = Genome::Sys->md5sum($original_file_path);

    my $result_users = delete($args{result_users});
    $result_users->{'feature-list'} = $self;
    $result_users->{requestor}    ||= $reference;
    $result_users->{sponsor}      ||= Genome::Sys->current_user;

    my $sr = Genome::Model::Build::ReferenceSequence::ConvertedBedResult->get_or_create(
        source_reference => $self->reference,
        target_reference => $reference,
        source_bed => $original_file_path,
        source_md5 => $original_md5,
        users => $result_users,
    );

    my $converted_file_path = delete($args{file_path});
    if (defined $converted_file_path) {
        Genome::Sys->create_symlink($sr->target_bed, $converted_file_path);
        return $converted_file_path;
    } else {
        return $sr->target_bed;
    }
}

sub converted_bed_file {
    my $self = shift;
    my %args = @_;
    if ($self->format eq 'unknown'){
        $self->error_message('Cannot convert BED file with unknown format');
        die $self->error_message;
    }
    unless ($self->reference) {
        $self->error_message('Cannot convert BED file without an associated reference');
        die $self->error_message;
    }

    return $self->generate_converted_bed_file(%args);
}

sub _resolve_param_value_from_text_by_name_or_id {
    my $class = shift;
    my $param_arg = shift;

    #First try default behaviour of looking up by name or id
    my @results = Command::V2->_resolve_param_value_from_text_by_name_or_id($class, $param_arg);

    #If that didn't work, and the argument is a filename, see if we have a feature list matching the provided file.
    if(!@results and -f $param_arg) {
        my $md5 = Genome::Sys->md5sum($param_arg);
        @results = Genome::FeatureList->get(file_content_hash => $md5);

        @results = grep( !Genome::Sys->diff_file_vs_file($param_arg, $_->file_path), @results);
    }

    return @results;
}

# Early detection for the common problem that the bed reference_name is not set correctly
sub _check_bed_list_is_on_correct_reference {
    my $self = shift;
    my $bed_file = $self->file_path;
    # Accept 0 or 1 return codes since grep returns 1 if it does not find anything
    my @chr_lines = Genome::Sys->capture([0,1], "grep --max-count=1 chr $bed_file");

    if (@chr_lines and not ($self->reference_name =~ m/nimblegen/) ) {
        die $self->error_message("It looks like your bed has 'chr' chromosomes but does not have a 'nimblegen' reference name (It is currently %s).\n".
            "This will result in your variant sets being filtered down to nothing. An example of a fix to this situation: \n".
            "genome feature-list update '%s' --reference nimblegen-human-buildhg19 (if your reference is hg19)", $self->reference_name, $self->name);
    }

    return 1;
}

sub resolve_bed_for_reference {
    my ($self, $reference_sequence_build) = @_;

    $self->_check_bed_list_is_on_correct_reference;

    my $bed_file;
    if($self->reference->id eq $reference_sequence_build->id) {
        $bed_file = $self->file_path;
    } else {
        $bed_file = $self->converted_bed_file(
            reference => $reference_sequence_build,
        );
    }
    return $bed_file;
}

sub get_tabix_and_gzipped_bed_file {
    my $self = shift;

    if ($self->format eq 'unknown') {
        die $self->error_message("Cannot convert format of BED file with unknown format");
    }

    my $processed_bed = $self->processed_bed_file(short_name => 0);
    my $sorted_processed_bed = Genome::Sys->create_temp_file_path;
    Genome::Model::Tools::Joinx::Sort->execute(
        input_files => [$processed_bed],
        output_file => $sorted_processed_bed,
    );
    return $self->gzip_and_tabix_bed($sorted_processed_bed);
}

sub gzip_and_tabix_bed {
    my $class = shift;
    my $file = shift;

    Genome::Sys->validate_file_for_reading($file);
    my $gzipped_file = Genome::Sys->create_temp_file_path;

    unless ( Genome::Sys->gzip_file($file, $gzipped_file) ) {
        die $class->error_message("Failed to gzip file ($file) to ($gzipped_file)");
    }

    my $tabix_cmd = Genome::Model::Tools::Tabix::Index->create(
        input_file => $gzipped_file,
        preset => 'bed',
    );

    unless ($tabix_cmd->execute) {
        die $class->error_message("Failed to tabix index file ($gzipped_file)");
    }

    return $gzipped_file;
}

sub chromosome_list {
    my $self = shift;
    my $indexed_bed = $self->get_tabix_and_gzipped_bed_file;
    my $list_command = Genome::Model::Tools::Tabix::ListChromosomes->execute(
        input_file => $indexed_bed,
    );

    my @chromosomes = $list_command? ($list_command->chromosomes) : ();
    unless (@chromosomes) {
        die $self->error_message("Failed to tabix list-chromosomes on indexed bed (%s) of bed (%s)", $indexed_bed, $self->processed_bed_file(short_name => 0));
    }

    return @chromosomes;
}

sub get_target_track_only {
    my ($self, $target_track) = @_;

    unless ($self->is_multitracked) {
        die $self->error_message("get_target_track_only can only be called on multitracked feature lists");
    }

    my ($ofh, $output_path)  = Genome::Sys->create_temp_file;
    my $ifh = Genome::Sys->open_file_for_reading($self->file_path);

    my $correct_track = 0;
    my $target_track_is_empty = 1;
    while (my $line = $ifh->getline) {
        if ( my ($name) = $line =~ m/^track name=(\w+) description=/) {
            if ($name eq $target_track) {
                $correct_track = 1;
            } else {
                $correct_track = 0;
            }
        } else {
            if ($correct_track) {
                $ofh->print($line);
                $target_track_is_empty = 0;
            }
        }
    }

    if ($target_track_is_empty) {
        die $self->error_message("No data found for target track ($target_track)! Does it exist in the file?");
    }

    return $output_path;
}

use constant STANDARDIZED_TRACK_NAMES => {
    # nimblegen lingo:
    'tiled_region'  => 'probes',
    'target_region' => 'targets',

    # agilent lingo:
    'Target Regions' => 'targets',
    'Probes'         => 'probes',
};

sub _standardize_track_name {
    my $class = shift;
    my $track_name = shift;


    if(exists STANDARDIZED_TRACK_NAMES->{$track_name}) {
        return STANDARDIZED_TRACK_NAMES->{$track_name};
    }

    return $track_name;
}

# input a track definition such as:
# track name=foo description="this is a track" some_arbitrary_attribute=yes
# and parse it into a hash of attribute/value pairs:
# %attrs = (
#    name                     => 'foo',
#    description              => 'this is a track',
#    some_arbitrary_attribute => 'yes',
# );
sub parse_track_definition {
    my $class = shift;
    my $track_def = shift;
    $track_def =~ s/^track\s*//;

    my %track_attrs;
    my @attr_pairs = parse_line('\s+', 1, $track_def); # break on spaces, ignoring spaces inside quotes
    foreach (grep {defined} @attr_pairs) {
        my ($attr, $value) = parse_line('=', 0, $_);
        die "badly formed track definition in BED file, could not parse track definition: $track_def" unless ($attr && $value);
        $track_attrs{$attr} = $value;
    }
    return %track_attrs;
}

sub _derive_format {
    my $class = shift;
    my ($is_1_based, $is_multitracked) = @_;

    return 'unknown' unless defined($is_1_based) && defined($is_multitracked);

    return 'multi-tracked 1-based' if $is_multitracked && $is_1_based;
    return '1-based' if !$is_multitracked && $is_1_based;
    return 'multi-tracked' if $is_multitracked && !$is_1_based;
    return 'true-BED' if !$is_multitracked && !$is_1_based;

    die 'Logic error in _derive_format.';
}


1;
