package Genome::FeatureList::Command::Merge;

use strict;
use warnings;

use feature qw(say);

use List::MoreUtils qw(any uniq);
use Genome;

class Genome::FeatureList::Command::Merge {
    is => 'Command::V2',
    has_input => [
        source_lists => {
            is => 'Genome::FeatureList',
            is_many => 1,
            doc => 'feature-lists to combine',
            shell_args_position => 1,
        },
        name => {
            is => 'Text',
            doc => 'name for the new feature-list',
        },
        reference => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'the (combined) reference for the merged feature-list, if the lists do not share a reference',
            is_optional => 1,
        },
    ],
    doc => 'command to combine several feature-lists into one',
};

sub help_brief {
    return 'Merge feature-lists together';
}

sub help_synopsis {
    return <<EOS
genome feature-list merge --name=ComboOf1And2 MyList1 MyList2
EOS
;
}

sub help_detail {
    return <<EOHELP
This command is used to combine several feature-lists into a new feature-list.
This process is carried out by merging only the tracks shared by all of the BED files
into a new BED file.

If the references do not match, supply a new reference sequence that is compatible
with all of the source references.
EOHELP
}

sub execute {
    my $self = shift;

    my @inputs = $self->source_lists;
    unless(@inputs > 1) {
        die $self->error_message('Must provide at least two FeatureLists to combine. Got %s.', scalar(@inputs));
    }

    unless($self->reference) {
        $self->reference($inputs[0]->reference);
    }
    $self->status_message('Using reference: %s', $self->reference->name);
    $self->validate_list_references_are_compatible($self->reference, @inputs);

    my $is_1_based = $self->validate_consistent_starting_base(@inputs);
    my $is_multitracked = $self->validate_consistent_multitrackedness(@inputs);

    my $new_bed_content = $self->combine_bed_content(@inputs);
    my $new_bed_md5 = Genome::Sys->md5sum_data($new_bed_content);

    my $merged_list = $self->find_existing_list($new_bed_md5, $self->name);
    unless($merged_list) {
        my $temp_file = Genome::Sys->create_temp_file_path;
        Genome::Sys->write_file($temp_file, $new_bed_content);

        my $cmd = Genome::FeatureList::Command::Create->create(
            name => $self->name,
            reference => $self->reference,
            file_path => $temp_file,
            content_type => 'targeted',
            format => Genome::FeatureList->_derive_format($is_1_based, $is_multitracked),
            description => 'generated with `genome feature-list merge`',
        );
        unless($merged_list = $cmd->execute) {
            die $self->error_message('Failed to create FeatureList');
        }
    }

    say $merged_list->id;
    return 1;
}

sub validate_properties_required_for_merge {
    my $class = shift;
    my @feature_lists = @_;

    my $error_count = 0;
    for my $list (@feature_lists) {
        unless($list->reference) {
            $class->error_message(
                'FeatureList %s has no reference set.',
                $list->__display_name__
            );
            $error_count++;
        }

        if($list->format eq 'unknown') {
            $class->error_message(
                'FeatureList %s has an unknown format.',
                $list->__display_name__
            );
            $error_count++;
        }
    }

    if($error_count) {
        die 'Cannot proceed without required properties set.';
    }

    return 1;
}

sub validate_list_references_are_compatible {
    my $class = shift;
    my $reference = shift;
    my @feature_lists = @_;

    my $error_count = 0;
    for my $list (@feature_lists) {
        unless($list->reference->is_compatible_with($reference)) {
            $class->error_message(
                'FeatureList %s has an incompatible reference, %s.',
                $list->__display_name__,
                $list->reference->name
            );
            $error_count++;
        }
    }

    if($error_count) {
        die 'To combine FeatureLists with different references, first define a new reference compatible with the input references and supply it as the --reference parameter.';
    }

    return 1;
}

sub validate_consistent_starting_base {
    my $class = shift;
    my @feature_lists = @_;
    my $is_1_based = $feature_lists[0]->is_1_based;

    my $error_count = 0;
    for my $list (@feature_lists) {
        unless($list->is_1_based eq $is_1_based) {
            $class->error_message(
                'FeatureList %s and %s have different starting bases.',
                $feature_lists[0]->__display_name__,
                $list->__display_name__
            );
            $error_count++;
        }
    }

    if($error_count) {
        die 'Combining FeatureLists with different starting bases is not currently supported.';
    }

    return $is_1_based;
}

sub validate_consistent_multitrackedness {
    my $class = shift;
    my @feature_lists = @_;
    my $is_multitracked = $feature_lists[0]->is_multitracked;

    my $error_count = 0;
    for my $list (@feature_lists) {
        unless($list->is_multitracked eq $is_multitracked) {
            $class->error_message(
                'FeatureList %s and %s are not both multi-tracked.',
                $feature_lists[0]->__display_name__,
                $list->__display_name__
            );
            $error_count++;
        }
    }

    if($error_count) {
        die $class->error_message('FeatureLists to be combined must either all be single tracked or all be multi-tracked.');
    }

    unless($is_multitracked) {
        die $class->error_message('Non-tracked BED files are not currently supported by this command.');
    }

    return $is_multitracked;
}

sub find_existing_list {
    my $class = shift;
    my ($bed_md5, $name) = @_;


    my @existing_lists = Genome::FeatureList->get(file_content_hash => $bed_md5);
    return unless @existing_lists;

    my $list_to_use;

    if(@existing_lists == 1) {
        $list_to_use = $existing_lists[0];
    } else {
        my @matches = grep { $_->name eq $name } @existing_lists;
        if(@matches == 1) {
            $list_to_use = $matches[0];
        }
    }

    unless($list_to_use) {
        $class->error_message(
            'Multiple matching existing FeatureLists found: %s',
            join(' ', map $_->id, @existing_lists)
        );
        die $class->error_message('To proceed make sure the --name parameter matches one of the existing lists.');
    }

    $class->status_message('Found existing list: %s', $list_to_use->__display_name__);
    return $list_to_use;
}


#### Based on GSC::BEDFile ####

sub combine_bed_content {
    my $class = shift;
    my @feature_lists = @_;

    # Example:
    # %track_content = ($bed_file_name => $content_hash)
    # where $content_hash = { $track_name => $track_data_array }

    my %track_content;
    foreach my $list (@feature_lists) {
        $class->debug_message('Attempting to parse tracks from %s', $list->name);
        my $file_path = $list->file_path;
        my $bed_data = Genome::Sys->read_file($file_path);
        my $content_hash = $class->hash_bed_content_by_tracks($bed_data);
        $track_content{$file_path} = $content_hash;
    }

    return $class->_combine_bed_content(\%track_content);
}

sub _combine_bed_content {
    my $class = shift;
    my $track_content = shift;

    foreach my $name (keys %$track_content) {
        $class->standardize_track_names($track_content->{$name});
    }

    my @approved_track_names = qw/targets probes/;

    my %completeness_count; # account for whether a track is provided in all input files or not
    foreach my $name (keys %$track_content) {
        my $c = $track_content->{$name};
        foreach my $track (@approved_track_names) {
            if (exists $c->{$track}) {
                $completeness_count{$track}++;
            }
        }
    }

    # in order to use information from a track, the track must be included in all input files. if not included, disregard the whole track
    my $combined_content;
    my $num_files = scalar keys %$track_content;
    foreach my $track (sort @approved_track_names) {
        if ($completeness_count{$track} == $num_files) {
            my @track_content;
            foreach my $file (keys %$track_content) {
                push @track_content, @{$track_content->{$file}->{$track}};
            }
            $combined_content .= "track name=$track\n" . join ("\n", uniq @track_content) . "\n";
        }
        else {
            $class->warning_message("track '$track' is not provided in all input files (" . join (', ', sort keys %$track_content) . "). Cannot combine its content");
        }
    }

    unless ($combined_content) {
        die $class->error_message("Could not find any complete track information that was provided in all input files! Must have at least one complete track");
    }

    return $combined_content;
}

sub standardize_track_names {
    my $class = shift;
    my $track_content = shift;
    my $accept_weird_tracks = shift; # typically no during normal processing / pooling. set to yes for some exploratory / information gathering purposes such as get_comparable_bed_file_content

    my %standardized_track_names = %{Genome::FeatureList->STANDARDIZED_TRACK_NAMES};

    my %approved_track_names = map {$_ => 1} uniq values %standardized_track_names;
    $approved_track_names{no_track} = 1 if ($accept_weird_tracks);

    # consolidate several naming conventions into one standardized track name set: (probes, targets) 
    foreach my $non_standard_track (keys %standardized_track_names) {
        my $standard_track = $standardized_track_names{$non_standard_track};
        if (exists $track_content->{$non_standard_track}) {
            push @{$track_content->{$standard_track}}, @{$track_content->{$non_standard_track}};
            delete $track_content->{$non_standard_track};
        }
    }

    my @track_names_to_ignore = ( 'Covered' );
    delete $track_content->{$_} foreach @track_names_to_ignore;  # just delete ignored tracks completely from the hash

    # make sure all tracks were named in a way known to us
    unless ($accept_weird_tracks) {
        if (any {not exists $approved_track_names{$_}} keys %$track_content) {
            die $class->error_message("Found unexpected extra tracks in BED file! " . join (', ', keys %$track_content));
        }
    }

    return $track_content;
}

# takes in raw file content and returns a hash keyed by track names
sub hash_bed_content_by_tracks {
    my $class = shift;
    my $content = shift;
    return unless $content;
    my $accept_weird_tracks = shift; # typically no during normal processing / pooling. set to yes for some exploratory / information gathering purposes such as get_comparable_bed_file_content

    my %ret;
    my $current_track;
    foreach my $line (split /\n/, $content) {
        next if $line =~ /^browser/;
        if ($line =~ /^track/) {
            my %track_attrs = Genome::FeatureList->parse_track_definition($line);
            $current_track = $track_attrs{name};
            unless ($current_track) {
                die $class->error_message("badly formed BED file track definition, track with no 'name' information: " . Data::Dumper::Dumper(\%track_attrs));
            }
            # don't pass along info about empty tracks after all
            #$ret{$current_track} = [];
        }
        else {
            unless ($current_track) {
                if ($accept_weird_tracks) {
                    $current_track = 'no_track';
                }
                else {
                    die $class->error_message("badly formed BED file, missing track definition");
                }
            }
            push @{$ret{$current_track}}, $line;
        }
    }
    return \%ret;
}

sub is_suitable_for_pooling {
    my $class = shift;
    my $content = shift;
    return unless $content;

    $content = eval {$class->hash_bed_content_by_tracks($content) };
    return unless $content;

    $content = eval {$class->standardize_track_names($content) };
    return 1 if $content;
    return;
}

# for purposes of comaparing BED files, get a sorted, standardized version of just the "important" 
# bits (the first three columns)
sub get_comparable_bed_file_content {
    my $class = shift;
    my $content = shift;
    $content = $class->hash_bed_content_by_tracks($content, 'accept_weird_tracks');
    $content = $class->standardize_track_names($content, 'accept_weird_tracks');

    my $functional_content;
    foreach my $track (keys %$content) {
        foreach my $row ( @{$content->{$track}} ) {
            $row =~ s/^chr//; # remove leading 'chr' if present'
            my @cols = split "\t", $row;
            $row = join "\t", @cols[0..2]; # ignore optional columns
        }
        $functional_content .= "track name=" . '"' . $track . '"' . "\n" . join ("\n", uniq sort @{$content->{$track}}) . "\n";
    }
    return $functional_content;
}

# Compare only the important bits. If we ignore the optional "description" fields, etc, is this BED
# file *functionally* the same as another?
sub loosely_verify_same_content {
    my $class = shift;
    my ($content1, $content2) = @_;
    return $class->get_comparable_bed_file_content($content1) eq $class->get_comparable_bed_file_content($content2);
}

sub get_bed_track_names {
    my $class = shift;
    my $content = shift;
    if (ref $content) {
        $content = $content->generate_or_retrieve_bed_content;
    }
    my @track_info = grep {/^track/} split /\n/, $content;

    my @track_names;
    foreach my $line (@track_info) {
        my %track_attrs = Genome::FeatureList->parse_track_definition($line);
        $track_attrs{name} ||= 'BAD';
        push @track_names, $track_attrs{name};
    }
    return sort @track_names;
}

sub combine_capture_set_names {
    my $class = shift;
    my $roi = shift;
    my @bed = @_;
    die "ROI passed in with multiple BED files, this cannot happen" if ($roi && @bed > 1);

    my $cs_name = join " + " , map {$_->generate_or_retrieve_capture_set_name($roi)} sort {$a->id <=> $b->id} @bed;
    $cs_name =~ s/\s+/ /; # make sure there aren't any double spaces
    return $cs_name;
}

sub combine_bed_filenames {
    my $class = shift;
    my $roi = shift;
    my @bed = @_;
    die "ROI passed in with multiple BED files, this cannot happen" if ($roi && @bed > 1);

    my $bed_name = join ("__", 
                         map {$_ =~ s/\.bed$//i; $_} 
                         map {$_->generate_or_retrieve_bed_filename($roi)} 
                         sort {$a->id <=> $b->id} @bed);

    $bed_name .= '.bed';
    return $bed_name;
}

####/ GSC::BEDFile ####


1;
