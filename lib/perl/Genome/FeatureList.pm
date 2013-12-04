package Genome::FeatureList;

use strict;
use warnings;

use Genome;

class Genome::FeatureList {
    table_name => 'model.feature_list',
    id_by => [
        id => {
            is => 'Text',
            len => 64,
        },
    ],
    has => [
        name => {
            is => 'Text',
            len => 200,
        },
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
        reference_name => {
            via => 'reference',
            to => 'name',
        },
        subject_id => {
            is => 'NUMBER',
            len => 10,
            doc => 'ID of the subject to which the features are relevant',
        },
        subject => {
            is => 'Genome::Model::Build',
            id_by => 'subject_id',
        },
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
            is => 'VARCHAR2',
            len => 255,
            valid_values => [ "exome", "targeted", "validation", "roi", undef ],
            doc => 'The type of list (used for determining automated processing)',
        },
        description => {
            is => 'VARCHAR2',
            len => 255,
            doc => 'General description of the BED file',
        },
    ],
    has_optional_transient => [
        #TODO These could be pre-computed and stored in the allocation rather than re-generated every time
        _processed_bed_file_path => {
            is => 'Text',
            doc => 'The path to the temporary dumped copy of the post-processed BED file',
        },
        _merged_bed_file_path => {
            is => 'Text',
            doc => 'The path to the temporary dumped copy of the merged post-processed BED file',
        },
        _converted_bed_file_paths => {
            is => 'HASH',
            doc => 'The paths to the temporary dumped copies of merged post-processed BED files converted to other references',
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
    #set the newly copied file to read only
    my $result = eval{
        chmod 0444, $self->file_path;
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
    while(my $line = <$fh>) {
        chomp($line);
        if($is_multitracked) {
            if ($line =~ /^track/) {
                $bed_file_content .= "$line\n";
                next;
            }
        }
        my @entry = split("\t",$line);
        unless (scalar(@entry) >= 3) {
            my $error_message = 'At least three fields are required in BED format files.  Error with line: '. $line;
            die($error_message);
        }
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
    while(my $line = <$fh>) {
        chomp($line);
        if($self->is_multitracked) {
            if ($line =~ /^track name=tiled_region/ or $line =~ /^track name=probes/) {
                if ($track_name eq 'tiled_region') {
                    $print = 1;
                } else {
                    $print = 0;
                }
                next;
            } elsif ($line =~ /^track name=target_region/ or $line =~ /^track name=targets/) {
                if ($track_name eq 'target_region') {
                    $print = 1;
                } else {
                    $print = 0;
                }
                next;
            } elsif ($line =~ /^track\s.*?name=/) {
                $self->warning_message('Unknown track name. Including regions.');
                $print = 1;
                next;
            }
        }
        if ($print) {
            my @entry = split("\t",$line);
            unless (scalar(@entry) >= 3) {
                $self->error_message('At least three fields are required in BED format files.  Error with line: '. $line);
                die($self->error_message);
            }
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

sub processed_bed_file {
    my $self = shift;
    my %args = @_;
    if($self->format eq 'unknown'){
        $self->error_message('Cannot process BED file with unknown format');
        die $self->error_message;
    }

    unless($self->_processed_bed_file_path) {
        my $content = $self->processed_bed_file_content(%args);
        my $temp_file = Genome::Sys->create_temp_file_path( $self->id . '.processed.bed' );
        Genome::Sys->write_file($temp_file, $content);
        $self->_processed_bed_file_path($temp_file);
    }

    return $self->_processed_bed_file_path;
}

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

    my $converted_file_path = delete($args{file_path}) || Genome::Sys->create_temp_file_path( $self->id . '.merged.' . $reference->id . '.converted.bed' );
    my $original_file_path;
    if ($merge) {
        $original_file_path = $self->merged_bed_file(%args);
    } else {
        $original_file_path = $self->processed_bed_file(%args);
    }
    my $converted_bed_file = Genome::Model::Build::ReferenceSequence::Converter->convert_bed($original_file_path, $self->reference, $converted_file_path, $reference);

    unless(-s $converted_bed_file) {
        $self->error_message('Could not convert to requested reference!');
        die $self->error_message;
    }

    return $converted_bed_file;
}

sub converted_bed_file {
    my $self = shift;
    my %args = @_;
    my $reference = $args{reference};
    #TODO: add short name and merge options
    if ($self->format eq 'unknown'){
        $self->error_message('Cannot convert BED file with unknown format');
        die $self->error_message;
    }

    unless ($self->reference) {
        $self->error_message('Cannot convert BED file without an associated reference');
        die $self->error_message;
    }

    my $converted_bed_files = $self->_converted_bed_file_paths;

    my $result;
    if(exists $converted_bed_files->{$reference->id}) {
        $result = $converted_bed_files->{$reference->id};
        unless(-e $result) { # cached result may have been deleted
            $result = $self->generate_converted_bed_file(%args);
        }
    } else {
        $result = $self->generate_converted_bed_file(%args);
    }
    $converted_bed_files->{$reference->id} = $result;
    $self->_converted_bed_file_paths($converted_bed_files);

    return $result;
}

sub _resolve_param_value_from_text_by_name_or_id {
    my $class = shift;
    my $param_arg = shift;

    #First try default behaviour of looking up by name or id
    my @results = Genome::Command::Base->_resolve_param_value_from_text_by_name_or_id($class, $param_arg);

    #If that didn't work, and the argument is a filename, see if we have a feature list matching the provided file.
    if(!@results and -f $param_arg) {
        my $md5 = Genome::Sys->md5sum($param_arg);
        @results = Genome::FeatureList->get(file_content_hash => $md5);

        @results = grep( !Genome::Sys->diff_file_vs_file($param_arg, $_->file_path), @results);
    }

    return @results;
}

1;
