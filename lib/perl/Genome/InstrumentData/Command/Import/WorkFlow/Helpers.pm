package Genome::InstrumentData::Command::Import::WorkFlow::Helpers;

use strict;
use warnings;

use Genome;

require Carp;
require File::Basename;
require File::Copy;
require Filesys::Df;
require List::MoreUtils;
require LWP::Simple;
use Regexp::Common;

class Genome::InstrumentData::Command::Import::WorkFlow::Helpers { 
    is => 'UR::Singleton',
};

#<MOVE>#
sub move_path {
    my ($self, $from, $to) = @_;
    $self->status_message('Move path...');

    my $from_sz = -s $from;
    $self->status_message("From: $from");
    $self->status_message("To: $to");
    my $move_ok = File::Copy::move($from, $to);
    if ( not $move_ok ) {
        $self->error_message('Move failed!');
        return;
    }
    my $to_sz = -s $to;
    if ( not $to_sz or $to_sz != $from_sz ) {
        $self->error_message("Move succeeded, but destination size is different from original! $to_sz vs $from_sz");
        return;
    }

    $self->status_message('Move path...done');
    return 1;
}
#<>#

#<FILES>#
sub source_file_retrieval_method {
    my ($self, $file) = @_;

    Carp::confess('No source file to get retrieval method!') if not $file;

    if ( $file =~ m#^https://cghub.ucsc.edu# ) {
        return 'cg hub';
    }
    elsif ( $file =~ m#^(http|https|file)\://# ) {
        return 'remote url';
    }
    else {
        return 'local disk';
    }
}
sub source_files_retrieval_method {
    my ($self, @files) = @_;

    Carp::confess('No source files to get retrieval method!') if not @files;

    my @source_file_retrieval_methods;
    for my $file ( @files ) {
        push @source_file_retrieval_methods, $self->source_file_retrieval_method($file);
    }

    @source_file_retrieval_methods = List::MoreUtils::uniq(@source_file_retrieval_methods);
    if ( @source_file_retrieval_methods > 1 ) {
        Carp::confess('Mixed file retrieval methods for source files! '. join(' ', sort @source_file_retrieval_methods));
    }

    return $source_file_retrieval_methods[0];
}

sub kilobytes_required_for_processing_of_source_files {
    my ($self, @source_files) = @_;

    Carp::confess('No source files to get kb needed!') if not @source_files;

    my %formats_and_multipliers = (
        bam => 3,
        fastq => 2,
        sra => 4,
    );

    my $kb_required = 0;
    for my $source_file ( @source_files ) {
        my $size = $self->file_size($source_file);
        if ( not $size ) {
            $self->error_message('Source file does not have any size! '.$source_file);
            return;
        }

        my $format = $self->source_file_format($source_file);
        return if not $format;

        my $kb_required_for_source_file = int($size / 1024) + 1; #convert to kb, +1 for ceiling
        $kb_required_for_source_file *= 3 if $source_file =~ /\.gz$/; # assume ~30% compression rate for gzipped fasta/q

        my $multiplier = $formats_and_multipliers{$format};
        $kb_required_for_source_file *= $multiplier; # for convresion to bam
        $kb_required += $kb_required_for_source_file;
    }

    return $kb_required;
}

sub file_size {
    my ($self, $source_file) = @_;

    Carp::confess('No file to get size!') if not $source_file;

    my $size;
    if ( $source_file =~ /^http/ ) {
        $size = $self->remote_file_size($source_file);
    }
    else {
        $size = -s $source_file;
    }

    return $size;
}

sub remote_file_size {
    my ($self, $remote_file) = @_;

    Carp::confess('No remote file to get size!') if not $remote_file;

    my $agent = LWP::UserAgent->new;
    my $response = $agent->head($remote_file);
    if ( not $response->is_success ) {
        $self->error_message($response->message) if $response->message;
        $self->error_message('HEAD failed for remote file!');
        return;
    }

    return $response->headers->content_length;
}

sub source_file_format {
    my ($self, $source_file) = @_;

    Carp::confess('No source file to get format!') if not $source_file;

    my %suffixes_to_original_format = (
        fastq => 'fastq',
        fq => 'fastq',
        fasta => 'fasta',
        fa => 'fasta',
        fna => 'fasta',
        bam => 'bam',
        sra => 'sra',
    );
    $source_file =~ s/\Q.$_\E$// for Genome::InstrumentData::Command::Import::WorkFlow::ArchiveToFastqs->types;
    my ($source_file_base_name, $path, $suffix) = File::Basename::fileparse(
        $source_file, keys %suffixes_to_original_format
    );
    if ( not $suffix or not $suffixes_to_original_format{$suffix} ) {
        $self->error_message('Unrecognized source file format! '.$source_file_base_name);
        return;
    }

    return $suffixes_to_original_format{$suffix};
}

sub is_source_file_archived {
    my ($self, $source_file) = @_;

    Carp::confess('No source file to determined if archived!') if not $source_file;

    return grep { 
        $source_file =~ /\Q$_\E$/
    } (qw/ .tgz .tar .tar.gz /);
}

sub verify_adequate_disk_space_is_available_for_source_files {
    my ($self, %params) = @_;
    $self->debug_message('Verify adequate disk space is available...');

    my $working_directory = delete $params{working_directory};
    Carp::confess('No tmp dir to verify adequate temp space is avaliable!') if not $working_directory;
    $self->debug_message("Tmp dir: $working_directory");
    my $source_files = delete $params{source_files};
    $self->debug_message("Source files: ".join(' ', @$source_files));
    Carp::confess('No source files to verify temp space!') if not $source_files;

    my $df = eval{ Filesys::Df::df($working_directory); };
    if( not $df ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to get "df" for temp dir! '.$working_directory);
        return;
    }
    $self->debug_message("Available Kb: ".$df->{bavail});

    my $kb_required = $self->kilobytes_required_for_processing_of_source_files(@$source_files);
    return if not $kb_required;
    $self->debug_message("Required Kb: ".$kb_required);

    my $remaining_space = $df->{bavail} - $kb_required;
    if ( $remaining_space < 1024 ) { # 1 Mb
        $self->error_message("There is not enough space in $working_directory to process source files!");
        return;
    }

    $self->debug_message('Verify adequate disk space is available...done');
    return 1;
}
#<>#

#<SAMTOOLS>#
sub load_or_run_flagstat {
    my ($self, $bam_path) = @_;
    $self->debug_message('Load or run flagstat...');

    Carp::confess('No bam path given to run flagstat!') if not $bam_path;
    Carp::confess('Bam path given to run flagstat does not exist!') if not -s $bam_path;

    my $flagstat_path = $bam_path.'.flagstat';
    if ( -s $flagstat_path ) {
        return $self->load_flagstat_for_bam_path($bam_path);
    }
    else {
        return $self->run_flagstat($bam_path);
    }
}

sub run_flagstat {
    my ($self, $bam_path) = @_;
    $self->debug_message('Run flagstat...');

    Carp::confess('No bam path given to run flagstat!') if not $bam_path;
    Carp::confess('Bam path given to run flagstat does not exist!') if not -s $bam_path;

    my $flagstat_path = $bam_path.'.flagstat';
    $self->debug_message("Bam path: $bam_path");
    $self->debug_message("Flagstat path: $flagstat_path");
    my $cmd = "samtools flagstat $bam_path > $flagstat_path";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $flagstat_path ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run flagstat!');
        return;
    }

    my $flagstat = $self->load_flagstat_for_bam_path($bam_path);
    return if not $flagstat;

    $self->debug_message('Run flagstat...done');
    return $flagstat;
}

sub load_flagstat_for_bam_path {
    my ($self, $bam_path) = @_;
    $self->debug_message('Load flagstat for bam path...');

    Carp::confess('No bam path to derive flagstat path to load!') if not $bam_path;

    my $flagstat_path = $bam_path.'.flagstat';
    $self->debug_message('Flagstat path: '.$flagstat_path);
    Carp::confess('Flagstat file is empty!') if not -s $flagstat_path;
    my $flagstat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_path);
    if ( not $flagstat ) {
        $self->error_message('Failed to load flagstat file!');
        return;
    }
    $flagstat->{path} = $flagstat_path;
    $flagstat->{is_paired_end} = $self->is_bam_paired_end($bam_path);

    $self->debug_message('Flagstat output:');
    $self->debug_message( join("\n", map { ' '.$_.': '.$flagstat->{$_} } sort keys %$flagstat) );

    $self->debug_message('Load flagstat for bam path...done');
    return $flagstat;
}

sub is_bam_paired_end {
    # Assumes bam is sorted with secondary alignments and duplicate reads removed!
    # Only checks first 10K reads
    my ($self, $bam_path, $flagstat) = @_;

    Carp::confess('No bam path given to is_bam_paired_end!') if not $bam_path;
    Carp::confess('Bam path given to is_bam_paired_end does not exist!') if not -s $bam_path;

    if ( $flagstat ) {
        return 0 if $flagstat->{reads_paired_in_sequencing} == 0; # no reads marked as read 1/2
        return 0 if $flagstat->{reads_marked_as_read1} != $flagstat->{reads_marked_as_read2}; # uneven read 1/2
    }

    my $bam_fh = IO::File->new("samtools view $bam_path | head -10000 |");
    if ( not $bam_fh ) {
        $self->error_message('Failed to open file handle to samtools command!');
        return;
    }

    my (%seen, $line1, $line2, $name1, $name2);
    my $is_paired_end = 1;
    while ( $line1 = $bam_fh->getline ) {
        $line2 = $bam_fh->getline;
        if ( not $line2 ) {
            # missing next line
            $is_paired_end = 0;
            last;
        }
        my ($name1) = split(/\t/, $line1);
        my ($name2) = split(/\t/, $line2);
        if ( $name1 ne $name2 ) {
            # different templates - maybe not sorted, but we shouldn't be imported these
            $is_paired_end = 0;
            last;
        }
        if ( exists $seen{$name1} ) {
            # template has multiple entries
            $is_paired_end = 0;
            last;
        }
        $seen{$name1} = 1;
        # check flags?
    }
    $bam_fh->close;

    return $is_paired_end;
}

sub validate_bam {
    my ($self, $bam_path, $flagstat_path) = @_;
    $self->debug_message('Validate bam...');

    my $flagstat = $self->load_or_run_flagstat($bam_path, $flagstat_path);
    return if not $flagstat;

    if ( not $flagstat->{total_reads} > 0 ) {
        $self->error_message('Flagstat determined that there are no reads in bam! '.$bam_path);
        return;
    }

    if ( $flagstat->{reads_marked_as_read1} > 0 and $flagstat->{reads_marked_as_read2} > 0 ) {
        if ( $flagstat->{reads_marked_as_read1} != $flagstat->{reads_marked_as_read2} ) {
            $self->error_message('Flagstat indicates that there are not equal pairs in bam! '.$bam_path);
            return;
        }
    }

    return $flagstat;
}

sub load_headers_from_bam {
    my ($self, $bam_path) = @_;
    $self->debug_message('Load headers...');

    Carp::confess('No bam path given to load headers!') if not $bam_path;
    Carp::confess('Bam path given to load headers does not exist!') if not -s $bam_path;

    $self->debug_message("Bam path: $bam_path");
    my $headers_fh = IO::File->new("samtools view -H $bam_path |");
    if ( not $headers_fh ) {
        $self->error_message('Failed to open file handle to samtools command!');
        return;
    }

    my $headers = {};
    while ( my $line = $headers_fh->getline ) {
        chomp $line;
        my ($type, $tags) = split(/\t/, $line, 2);
        push @{$headers->{$type}}, $tags;
    }
    $headers_fh->close;

    $self->debug_message('Load headers...done');
    return $headers;
}
 
sub read_groups_from_headers {
    my ($self, $rg_headers) = @_;
    $self->debug_message('Read groups from headers...');

    Carp::confess('No read group headers given to read groups from headers!') if not $rg_headers;
    Carp::confess('Invalid read group headers given to read groups from headers! '.Data::Dumper::Dumper($rg_headers)) unless ref($rg_headers) eq 'ARRAY';

    my %read_groups_from_headers;
    return \%read_groups_from_headers if not @$rg_headers;

    for my $rg_header ( @$rg_headers ) {
        my %tags = map { split(':', $_, 2) } split(/\t/, $rg_header);
        my $rg_id = delete $tags{ID};
        if ( not defined $rg_id ) {
            $self->error_message("No ID tag in read group header! \@RG\t$rg_header");
            return;
        }
        $read_groups_from_headers{ $rg_id } = join("\t", map { $_.':'.$tags{$_} } sort keys %tags);
    }

    $self->debug_message('Read groups from headers...done');
    return \%read_groups_from_headers;
}

sub load_read_groups_from_bam {
    my ($self, $bam_path) = @_;
    $self->debug_message('Load read groups from bam...');

    my $headers = $self->load_headers_from_bam($bam_path);
    return if not $headers;

    my $read_groups_from_headers = $self->read_groups_from_headers($headers->{'@RG'} || []);
    return if not $read_groups_from_headers;

    $self->debug_message('Load read groups from bam...done');
    return [ sort keys %$read_groups_from_headers ];
}

sub headers_to_string {
    my ($self, $orig_headers) = @_;
    $self->debug_message('Header string from headers...');

    Carp::confess('No headers given to headers to string!') if not $orig_headers;
    Carp::confess('Invalid headers given to headers to string! '.Data::Dumper::Dumper($orig_headers)) if ref($orig_headers) ne 'HASH';

    my %headers = %$orig_headers;
    my $string;
    for my $type (qw/ @HD @SQ @RG @PG @CO /) {
        my $tags = delete $headers{$type};
        next if not $tags;
        $string .= join("\n", map { $type."\t".$_ } @$tags)."\n";
    }

    for my $type ( keys %headers ) {
        my $tags = delete $headers{$type};
        $string .= join("\n", map { $type."\t".$_ } @$tags)."\n";
    }

    $self->debug_message('Header string from headers...done');
    return $string;
}
#<>#

#<READ COUNT>#
sub load_read_count_for_fastq_paths {
    my ($self, @fastq_paths) = @_;

    Carp::confess('No fastq paths to load read counts!') if not @fastq_paths;

    my %fastq_paths_and_read_counts;
    for my $fastq_path ( @fastq_paths ) {
        my $line_count_path = $fastq_path.'.count';
        if ( not -s $line_count_path ) {
            $self->error_message('No line count path for fastq path! '.$fastq_path);
            return;
        }
        my $read_count = $self->load_read_count_from_line_count_path($line_count_path);
        return if not defined $read_count;
        $fastq_paths_and_read_counts{$fastq_path} = $read_count;
    }

    return \%fastq_paths_and_read_counts;
}

sub load_read_count_from_line_count_path {
    my ($self, $line_count_path) = @_;

    my $line_count = eval{ Genome::Sys->read_file($line_count_path); };
    if ( not defined $line_count ) {
        $self->error_message('Failed to open line count file! '.$@);
        return;
    }

    $line_count =~ s/\s+//g;
    if ( $line_count !~ /^\d+$/ ) {
        $self->error_message('Invalid line count! '.$line_count);
        return;
    }

    if ( $line_count == 0 ) {
        $self->error_message('Read count is 0!');
        return;
    }

    if ( $line_count % 4 != 0 ) {
        $self->error_message('Line count is not divisible by 4! '.$line_count);
        return;
    }

    return $line_count / 4;
}
#<>#

#<KEY VALUE PAIRS TO HASH>#
sub key_value_pairs_to_hash {
    my ($self, @key_value_pairs) = @_;

    my %properties;
    for my $key_value_pair ( @key_value_pairs ) {
        my ($label, $value) = split('=', $key_value_pair);
        if ( not defined $value or $value eq '' ) {
            $self->error_message('Failed to parse with instrument data property label/value! '.$key_value_pair);
            return;
        }
        if ( exists $properties{$label} and $value ne $properties{$label} ) {
            $self->error_message(
                "Multiple values for instrument data property! $label => ".join(', ', sort $value, $properties{$label})
            );
            return;
        }
        $properties{$label} = $value;
    }

    return \%properties;
}
#<>#

#<MD5>#
sub md5_path_for {
    Carp::confess('No path given to get md5 path!') if not $_[1];
    return $_[1].'.md5';
}

sub original_md5_path_for {
    Carp::confess('No path given to get original md5 path!') if not $_[1];
    return $_[1].'.md5-orig';
}

sub load_or_run_md5 {
    my ($self, $path, $md5_path) = @_;
    $self->debug_message('Load or run MD5...');

    Carp::confess('No path given to run MD5!') if not $path;
    Carp::confess('Path given to run MD5 on does not exist!') if not -s $path;

    $md5_path ||= $self->md5_path_for($path);
    my $md5;
    if ( -s $md5_path ) {
        $md5 = $self->load_md5($md5_path);
    }
    else {
        $md5 = $self->run_md5($path, $md5_path);
    }

    return $md5;
}

sub run_md5 {
    my ($self, $path, $md5_path) = @_;
    $self->debug_message('Run MD5...');

    Carp::confess('No path given to run md5!') if not $path;
    Carp::confess('Path given to run md5 does not exist!') if not -s $path;

    $md5_path ||= $self->md5_path_for($path);
    $self->debug_message("Path: $path");
    $self->debug_message("MD5 path: $md5_path");
    my $cmd = "md5sum $path > $md5_path";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $md5_path ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run md5!');
        return;
    }

    my $md5 = $self->load_md5($md5_path);
    return if not $md5;

    $self->debug_message('Run MD5...done');
    return $md5;
}

sub load_md5 {
    my ($self, $md5_path) = @_;
    $self->debug_message('Load MD5...');

    Carp::confess('No MD5 path to load!') if not $md5_path;

    $self->debug_message('MD5 path: '.$md5_path);
    my $fh = Genome::Sys->open_file_for_reading($md5_path);
    if ( not $fh ) {
        $self->error_message('Failed to open MD5 path!');
        return;
    }

    my $line = $fh->getline;
    $fh->close;
    if ( not $line ) {
        $self->error_message('Failed to read line from MD5 path!');
        return;
    }
    chomp $line;

    my ($md5) = split(/\s+/, $line);
    if ( not $md5 ) {
        $self->error_message('Failed to get MD5 from line! '.$line);
        return;
    }

    $self->debug_message('MD5: '.$md5);

    $self->debug_message('Load MD5...done');
    return $md5;
}

sub were_original_path_md5s_previously_imported {
    my ($self, %params) = @_;

    my $md5s = delete $params{md5s};
    Carp::confess('No md5s given to check if previously imported!') if not $md5s or not @$md5s;
    my $downsample_ratio = delete $params{downsample_ratio};
    Carp::confess('Unknown params sent to were_original_path_md5s_previously_imported: '.Data::Dumper::Dumper(\%params)) if %params;

    # Gather original_data_path_md5 attrs
    my @original_data_path_md5_attrs = Genome::InstrumentDataAttribute->get(
        attribute_label => 'original_data_path_md5',
        'attribute_value in' => $md5s,
    );

    # No existing original_data_path_md5s == not previously imported
    return if not @original_data_path_md5_attrs;

    # Get unique inst data from original_data_path_md5 attrs
    my @instrument_data = List::MoreUtils::uniq( map { $_->instrument_data } @original_data_path_md5_attrs);

    # Check if md5 was previously imported
    my @previously_imported_instrument_data;
    my $was_previously_imported = sub{
        my ($downsample_ratio_attr) = @_;
        # if no downsample_ratio given and no downsample_ratio for inst data
        return 1 if not defined $downsample_ratio and not $downsample_ratio_attr;
        # if given a downsample_ratio and there is a downsample_ratio_attr
        #  and the attribute_value matches the given downsample_ratio
        return 1 if defined $downsample_ratio and $downsample_ratio_attr
            and $downsample_ratio == $downsample_ratio_attr->attribute_value;
        # not previously imported
        return;
    };
    for my $instrument_data ( @instrument_data ) {
        my $downsample_ratio_attr = Genome::InstrumentDataAttribute->get(
            instrument_data_id => $instrument_data->id,
            attribute_label => 'downsample_ratio',
        );
        push @previously_imported_instrument_data, $instrument_data if $was_previously_imported->($downsample_ratio_attr);
    }
    return if not @previously_imported_instrument_data;

    $self->error_message(
        'Instrument data was previously'.
            ( $downsample_ratio ? " downsampled by a ratio of $downsample_ratio and" : '' ).
            ' imported! Found existing instrument data: '.
            join(', ', sort map { $_->id } @previously_imported_instrument_data)
    );
    return 1;
}
#<>#

sub remove_paths_and_auxiliary_files {
    my ($self, @paths) = @_;

    Carp::confess('No source paths to remove!') if not @paths;

    for my $path ( @paths ) {
        for my $path_to_rm ( glob($path.'*') ) {
            unlink $path_to_rm;
        }
    }

    return 1;
}

sub work_flow_operation_class_from_name {
    my ($self, $name) = @_;
    die 'No name given to get work flow operation class!' if not $name;
    return 'Genome::InstrumentData::Command::Import::WorkFlow::'
        . join('', map { ucfirst } split(' ', $name));
}

sub add_operation_to_workflow_by_class {
    my ($self, $wf, $command_class_name) = @_;

    die 'No work flow given to add opertion!' if not $wf;
    die 'No class given to add operation to work flow!' if not $command_class_name;

    my $operation_type = Workflow::OperationType::Command->create(command_class_name => $command_class_name);
    if ( not $operation_type ) {
        $self->error_message("Failed to create work flow operation for ".$command_class_name);
        return;
    }

    my $name = $command_class_name->command_name_brief;
    $name =~ s/\-/ /g;
    my $operation = $wf->add_operation(
        name => $name,
        operation_type => $operation_type,
    );

    return $operation;
}


sub add_operation_to_workflow_by_name {
    my ($self, $wf, $name) = @_;

    die 'No work flow given to add opertion!' if not $wf;
    die 'No name given to add operation to work flow!' if not $name;

    return $self->add_operation_to_workflow_by_class( $wf, $self->work_flow_operation_class_from_name($name) );
}

sub insert_extension_into_bam_path {
    my ($self, $bam_path, $ext) = @_;

    die 'No bam path given to insert extension to bam path!' if not $bam_path;
    die 'No extension given to insert extension to bam path!' if not $ext;

    if ( not $bam_path =~ s/\.bam$// ) {
        die 'Failed to insert extension into bam path! Bam path does not end in .bam! '.$bam_path;
    }

    return join('.', $bam_path, $ext, 'bam');
}

#<VALIDATORS>#
sub is_downsample_ratio_invalid {
    my ($self, $downsample_ratio) = @_;

    Carp::confess('No downsample ratio to check!') if not defined $downsample_ratio;

    if ( $downsample_ratio !~ /$RE{num}{real}/ ) {
        return (
            UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ downsample_ratio /],
                desc => 'Invalid number! '.$downsample_ratio,
            )
        );
    }

    if ( $downsample_ratio <= 0 or $downsample_ratio >= 1 ) {
        return (
            UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ downsample_ratio /],
                desc => 'Must be greater than 0 and less than 1! '.$downsample_ratio,
            )
        );
    }

    return;
}
#<>#

#< InstData Metrics >#
sub update_bam_metrics_for_instrument_data {
    my ($self, $instrument_data) = @_;

    Carp::confess('No instrument data given to update bam for instrument data!') if not $instrument_data;

    my $bam_path = $instrument_data->bam_path;
    Carp::confess('No bam path set to update bam for instrument data!') if not $bam_path;
    Carp::confess('Bam path to update bam for instrument data does not exist!') if not -s $bam_path;

    my $flagstat_path = $bam_path.'.flagstat';
    my $flagstat = $self->load_flagstat_for_bam_path($bam_path);
    return if not $flagstat;

    my $read_length = $self->determine_read_length_in_bam($bam_path);
    return if not defined $read_length;

    my %metrics = (
        bam_path => $bam_path,
        is_paired_end => $flagstat->{is_paired_end},
        read_count => $flagstat->{total_reads},
        read_length => $read_length,
        # TODO add these?
        #base_count => $read_length * $read_count, # might be inaccurate if reads are not the same length
        #fragment_count => $read_count * 2, # needed? add other fragment info?
    );

    for my $name ( keys %metrics ) {
        eval{ $instrument_data->attributes(attribute_label => $name)->delete; }; # remove existing
        $instrument_data->add_attribute(
            attribute_label => $name,
            attribute_value => $metrics{$name},
            nomenclature => 'WUGC',
        );
    }

    return 1;
}

sub determine_read_length_in_bam {
    my ($self, $bam_path) = @_;

    my $read_length = `samtools view $bam_path | head -n 1000 | awk '{print \$10}'  | xargs -I{} expr length {} | perl -mstrict -e 'my \$c = 0; my \$t = 0; while (<>) { chomp; \$c++; \$t += \$_; } printf("%d\\n", \$t/\$c);'`;
    chomp $read_length;

    if ( not $read_length ) {
        $self->error_message('Failed to get read length from bam! '.$bam_path);
        return;
    }

    if ( $read_length !~ /^\d+$/ ) {
        $self->error_message("Non numeric read length ($read_length) from bam! ".$bam_path);
        return;
    }

    if ( $read_length == 0 ) {
        $self->error_message('Read length for bam is 0! '.$bam_path);
        return;
    }

    return $read_length;
}
#<>#

1;

