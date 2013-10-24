package Genome::InstrumentData::Command::Import::WorkFlow::Helpers;

use strict;
use warnings;

use Genome;

require Carp;
require File::Basename;
require File::Copy;
require Filesys::Df;
require LWP::Simple;

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
        $self->error_message("Move succeeded, but destination size is diffeerent from original! $to_sz vs $from_sz");
        return;
    }

    $self->status_message('Move path...done');
    return 1;
}
#<>#

#<FILE SIZE>#
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
            $self->error_message('Source file does have any size! '.$source_file);
            return;
        }

        my $format = $self->source_file_format($source_file);
        return if not $format;

        my $kb_required_for_source_file = int($size / 1024); #convert to kb
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

    my $source_file_base_name = File::Basename::basename($source_file);
    my @parts = split(/\./, $source_file_base_name);
    my $suffix;
    do {
        $suffix = pop @parts;
    } until not defined $suffix or ( $suffix !~ /t?gz/ and $suffix ne 'tar' );

    if ( not $suffix ) {
        $self->error_message("Failed to get suffix from source file! $source_file");
        return;
    }

    my %suffixes_to_original_format = (
        fastq => 'fastq',
        fq => 'fastq',
        fasta => 'fasta',
        fa => 'fasta',
        fna => 'fasta',
        bam => 'bam',
        sra => 'sra',
    );
    my $format = $suffixes_to_original_format{$suffix};
    if ( not $format ) {
        $self->error_message('Unrecognized source file format! '.$source_file_base_name);
        return;
    }

    return $format;
}

sub verify_adequate_disk_space_is_available_for_source_files {
    my ($self, %params) = @_;
    $self->status_message('Verify adequate disk space is available...');

    my $working_directory = delete $params{working_directory};
    Carp::confess('No tmp dir to verify adequate temp space is avaliable!') if not $working_directory;
    $self->status_message("Tmp dir: $working_directory");
    my $source_files = delete $params{source_files};
    $self->status_message("Source files: ".join(' ', @$source_files));
    Carp::confess('No source files to verify temp space!') if not $source_files;

    my $df = eval{ Filesys::Df::df($working_directory); };
    if( not $df ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to get "df" for temp dir! '.$working_directory);
        return;
    }
    $self->status_message("Available Kb: ".$df->{bavail});

    my $kb_required = $self->kilobytes_required_for_processing_of_source_files(@$source_files);
    return if not $kb_required;
    $self->status_message("Required Kb: ".$kb_required);

    my $remaining_space = $df->{bavail} - $kb_required;
    if ( $remaining_space < 1024 ) { # 1 Mb
        $self->error_message("There is not enough space in $working_directory to process source files!");
        return;
    }

    $self->status_message('Verify adequate disk space is available...done');
    return 1;
}
#<>#

#<SAMTOOLS>#
sub load_or_run_flagstat {
    my ($self, $bam_path, $flagstat_path) = @_;
    $self->status_message('Load or run flagstat...');

    Carp::confess('No bam path given to run flagstat!') if not $bam_path;
    Carp::confess('Bam path given to run flagstat does not exist!') if not -s $bam_path;

    $flagstat_path ||= $bam_path.'.flagstat';
    my $flagstat;
    if ( -s $flagstat_path ) {
        $flagstat = $self->load_flagstat($flagstat_path);
    }
    else {
        $flagstat = $self->run_flagstat($bam_path, $flagstat_path);
    }

    return $flagstat;
}

sub run_flagstat {
    my ($self, $bam_path, $flagstat_path) = @_;
    $self->status_message('Run flagstat...');

    Carp::confess('No bam path given to run flagstat!') if not $bam_path;
    Carp::confess('Bam path given to run flagstat does not exist!') if not -s $bam_path;

    $flagstat_path ||= $bam_path.'.flagstat';
    $self->status_message("Bam path: $bam_path");
    $self->status_message("Flagstat path: $flagstat_path");
    my $cmd = "samtools flagstat $bam_path > $flagstat_path";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $flagstat_path ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run flagstat!');
        return;
    }

    my $flagstat = $self->load_flagstat($flagstat_path);
    return if not $flagstat;

    $self->status_message('Run flagstat...done');
    return $flagstat;
}

sub load_flagstat {
    my ($self, $flagstat_path) = @_;
    $self->status_message('Load flagstat...');

    Carp::confess('No flagstat path to load!') if not $flagstat_path;

    $self->status_message('Flagstat path: '.$flagstat_path);
    my $flagstat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_path);
    if ( not $flagstat ) {
        $self->error_message('Failed to load flagstat file!');
        return;
    }

    # FIXME What is paired end?
    if ( $flagstat->{reads_paired_in_sequencing} > 0 and $flagstat->{reads_marked_as_read1} == $flagstat->{reads_marked_as_read2} ) {
        # Only set paired end if read1 and read2 are equal
        $flagstat->{is_paired_end} = 1;
    }
    else {
        $flagstat->{is_paired_end} = 0;
    }

    $self->status_message('Flagstat output:');
    $self->status_message( join("\n", map { ' '.$_.': '.$flagstat->{$_} } sort keys %$flagstat) );

    $self->status_message('Load flagstat...done');
    return $flagstat;
}

sub validate_bam {
    my ($self, $bam_path, $flagstat_path) = @_;
    $self->status_message('Validate bam...');

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
    $self->status_message('Load headers...');

    Carp::confess('No bam path given to load headers!') if not $bam_path;
    Carp::confess('Bam path given to load headers does not exist!') if not -s $bam_path;

    $self->status_message("Bam path: $bam_path");
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

    $self->status_message('Load headers...done');
    return $headers;
}
 
sub read_groups_from_headers {
    my ($self, $rg_headers) = @_;
    $self->status_message('Read groups from headers...');

    Carp::confess('No read group headers given to read groups from headers!') if not $rg_headers;
    Carp::confess('Invalid read group headers given to read groups from headers! '.Data::Dumper::Dumper($rg_headers)) unless ref($rg_headers) eq 'ARRAY';

    my %read_groups_from_headers;
    return \%read_groups_from_headers if not @$rg_headers;

    for my $rg_header ( @$rg_headers ) {
        my %tags = map { split(':', $_, 2) } split(/\t/, $rg_header);
        my $rg_id = delete $tags{ID};
        if ( not $rg_id ) {
            $self->error_message("No ID tag in read group header! \@RG\t$rg_header");
            return;
        }
        $read_groups_from_headers{ $rg_id } = join("\t", map { $_.':'.$tags{$_} } sort keys %tags);
    }

    $self->status_message('Read groups from headers...done');
    return \%read_groups_from_headers;
}

sub load_read_groups_from_bam {
    my ($self, $bam_path) = @_;
    $self->status_message('Load read groups from bam...');

    my $headers = $self->load_headers_from_bam($bam_path);
    return if not $headers;

    my $read_groups_from_headers = $self->read_groups_from_headers($headers->{'@RG'} || []);
    return if not $read_groups_from_headers;

    $self->status_message('Load read groups from bam...done');
    return [ sort keys %$read_groups_from_headers ];
}

sub headers_to_string {
    my ($self, $orig_headers) = @_;
    $self->status_message('Header string from headers...');

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

    $self->status_message('Header string from headers...done');
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
sub load_or_run_md5 {
    my ($self, $path, $md5_path) = @_;
    $self->status_message('Load or run MD5...');

    Carp::confess('No path given to run MD5!') if not $path;
    Carp::confess('Path given to run MD5 on does not exist!') if not -s $path;

    $md5_path ||= $path.'.md5';
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
    $self->status_message('Run MD5...');

    Carp::confess('No path given to run md5!') if not $path;
    Carp::confess('Path given to run md5 does not exist!') if not -s $path;

    $md5_path ||= $path.'.md5';
    $self->status_message("Path: $path");
    $self->status_message("MD5 path: $md5_path");
    my $cmd = "md5sum $path > $md5_path";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $md5_path ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run md5!');
        return;
    }

    my $md5 = $self->load_md5($md5_path);
    return if not $md5;

    $self->status_message('Run MD5...done');
    return $md5;
}

sub load_md5 {
    my ($self, $md5_path) = @_;
    $self->status_message('Load MD5...');

    Carp::confess('No MD5 path to load!') if not $md5_path;

    $self->status_message('MD5 path: '.$md5_path);
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

    $self->status_message('MD5: '.$md5);

    $self->status_message('Load MD5...done');
    return $md5;
}

sub ensure_original_data_path_md5s_were_not_previously_imported {
    my ($self, @md5s) = @_;

    Carp::confess('No md5s given to ensure instrument data were not imported!') if not @md5s;

    my @instrument_data_attr = Genome::InstrumentDataAttribute->get(
        attribute_label => 'original_data_path_md5',
        'attribute_value in' => \@md5s,
    );
    if ( @instrument_data_attr ) {
        $self->error_message(
            "Instrument data was previously imported! Found existing instrument data with MD5s: ".
            join(', ', map { $_->instrument_data_id.' => '.$_->attribute_value } @instrument_data_attr),
        );
        return;
    }

    return 1;
}

#<>#

sub remove_paths_and_auxiliary_files {
    my ($self, @paths) = @_;

    Carp::confess('No source paths to remove!') if not @paths;

    for my $path ( @paths ) {
        for my $ext ( '', '.md5', '.md5-orig', '.flagstat' ) {
            my $path = $path.$ext;
            unlink $path if -e $path;
        }
    }

    return 1;
}

1;

