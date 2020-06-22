package Genome::InstrumentData::Imported;

use strict;
use warnings;

use Genome;

require File::Basename;
use File::stat;
use File::Path;
use Genome::InstrumentData::Microarray;
use Set::Scalar;

class Genome::InstrumentData::Imported {
    is => 'Genome::InstrumentData',
    roles => 'Genome::Role::Searchable',
    has_optional => [
        source => { is => 'Genome::Subject', via => 'sample', to => 'source', },
        source_id => { is=> 'Text', via => 'source', to => 'id', },
        source_name => { is=> 'Text', via => 'source', to => 'name', },
        reference_sequence_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            id_by => 'reference_sequence_build_id',
        },
        full_name => {
            calculate_from => [ 'run_name', 'subset_name' ],
            calculate => q( $subset_name ? "$run_name/$subset_name" : $run_name ),
        },
    ],
    has_optional_attribute => [
        import_date => {
            is => 'DateTime',
        },
        user_name => {
            is => 'Text',
        },
        original_data_path => {
            is => 'DirectoryPath',
        },
        read_count => {
            is => 'Number',
        },
        base_count => {
            is => 'Number',
        },
        fragment_count => {
            is => 'Number',
        },
        fwd_read_length => {
            is => 'Number',
        },
        is_paired_end => {
            is => 'Boolean',
        },
        median_insert_size => {
            is => 'Number',
        },
        read_length => {
            is => 'Number',
        },
        rev_read_length => {
            is => 'Number',
        },
        sd_above_insert_size => {
            is => 'Number',
        },
        target_region_set_name => {
            is => 'Text',
        },
        sra_accession => {
            is => 'Text',
        },
        sra_sample_id => {
            is => 'Text',
        },
        barcode => {
            is => 'Text',
        },
        reference_sequence_build_id => {
            is => 'Number',
        },
        blacklisted_segments => {
            is_many => 1,
            is => 'Text',
        },
        downsample_ratio => {
            is => 'Float',
        },
    ],
};

sub __display_name__ {
    my $self = $_[0];
    return (
        join(' ', map { $self->$_ } qw/sequencing_platform import_format id/)
        . ($self->desc ? ' (' . $self->desc . ')' : '')
    );
}

sub run_name {
    my $self = shift;
    return $self->{run_name} if defined $self->{run_name};
    return $self->id;
}

sub data_directory {
    my $self = shift;

    my $disk_allocation = $self->disk_allocation;
    if ( not $disk_allocation ) {
        die $self->error_message("Could not find disk allocation for ".$self->__display_name__);
    }

    return $disk_allocation->absolute_path;
}

sub calculate_alignment_estimated_kb_usage {
    my $self = shift;
    my $answer;

    # Check for an existing allocation for this instrument data, which would've been created by the importer
    my $allocation = $self->disk_allocation;
    if ($allocation) {
        return int(($allocation->kilobytes_requested/1000) + 100);
    }

    if($self->original_data_path !~ /\,/ ) {
        if (-d $self->original_data_path) {
            my $source_size = Genome::Sys->directory_size_recursive($self->original_data_path);
            $answer = ($source_size/1000)+ 100;
        } else {
            unless ( -e $self->original_data_path) {
                $self->error_message("Could not locate directory or file to import.");
                die $self->error_message;
            }
            my $stat = stat($self->original_data_path);
            $answer = ($stat->size/1000) + 100;
        }
    }
    else {
        my @files = split /\,/  , $self->original_data_path;
        my $stat;
        my $size;
        foreach (@files) {
            if (-s $_) {
                $stat = stat($_);
                $size += $stat->size;
            } else {
                die "file not found - $_\n";
            }
        }
        $answer = ($size/1000) + 100;
    }
    return int($answer);
}


sub create {
    my $class = shift;

    my %params = @_;

    unless (exists $params{import_date}) {
        my $date = UR::Context->current->now;
        $params{import_date} = $date;
    }
    unless (exists $params{user_name}) {
        my $user = Genome::Sys->username;
        $params{user_name} = $user;
    }

    my $self = $class->SUPER::create(%params);
    return $self;
}

sub native_qual_format {
    my $self = shift;
    if ($self->import_format =~ /illumina/) {
        return "illumina";
    }
    elsif ($self->import_format =~ /solexa/) {
        return "solexa";
    }
    elsif ($self->import_format =~ /sanger/) {
        return "sanger";
    }
    return;
}

################## Solexa Only ###################
# aliasing these methods before loading Genome::InstrumentData::Solexa causes it to
# believe Genome::InstrumentData::Solexa is already loaded.  So we load it first...
##################################################
BEGIN: {
    Genome::InstrumentData::Solexa->class;
    no warnings 'once';
    *solexa_dump_sanger_fastq_files= \&Genome::InstrumentData::Solexa::dump_sanger_fastq_files;
    *dump_illumina_fastq_files= \&Genome::InstrumentData::Solexa::dump_illumina_fastq_files;
    *dump_solexa_fastq_files= \&Genome::InstrumentData::Solexa::dump_solexa_fastq_files;
    *dump_illumina_fastq_archive = \&Genome::InstrumentData::Solexa::dump_illumina_fastq_archive;
    *_unprocessed_fastq_filenames= \&Genome::InstrumentData::Solexa::_unprocessed_fastq_filenames;
    *validate_fastq_directory = \&Genome::InstrumentData::Solexa::validate_fastq_directory;
    *resolve_fastq_filenames = \&Genome::InstrumentData::Solexa::resolve_fastq_filenames;
    *fragment_fastq_name = \&Genome::InstrumentData::Solexa::fragment_fastq_name;
    *read1_fastq_name = \&Genome::InstrumentData::Solexa::read1_fastq_name;
    *read2_fastq_name = \&Genome::InstrumentData::Solexa::read2_fastq_name;
    *dump_trimmed_fastq_files = \&Genome::InstrumentData::Solexa::dump_trimmed_fastq_files;
    *_get_trimq2_params = \&Genome::InstrumentData::Solexa::_get_trimq2_params;
    *resolve_median_insert_size = \&Genome::InstrumentData::Solexa::resolve_median_insert_size;
    *resolve_sd_insert_size = \&Genome::InstrumentData::Solexa::resolve_sd_insert_size;
    *get_default_alignment_metrics = \&Genome::InstrumentData::Solexa::get_default_alignment_metrics;
    *get_default_alignment_results = \&Genome::InstrumentData::Solexa::get_default_alignment_results;
    *get_default_alignment_metrics_hash = \&Genome::InstrumentData::Solexa::get_default_alignment_metrics_hash;
    *_convert_trimmer_to_sx_commands = \&Genome::InstrumentData::Solexa::_convert_trimmer_to_sx_commands;
    *is_capture = \&Genome::InstrumentData::Solexa::is_capture;
    *target_region_set = \&Genome::InstrumentData::Solexa::target_region_set;
}

sub dump_sanger_fastq_files {
    my $self = shift;

    if ($self->import_format eq 'bam') {
        return $self->dump_fastqs_from_bam(@_);
    } else {
        return $self->solexa_dump_sanger_fastq_files(@_);
    }
}

sub total_bases_read {
    my $self = shift;
    #guess return zero rather than fail a build for this
    my $read_length = $self->read_length || 0;
    my $read_count = $self->read_count || 0;

    return $read_length * $read_count;
}

# leave as-is for first test,
# ultimately find out what uses this and make sure it really wants clusters
sub _calculate_total_read_count {
    my $self = shift;
    return $self->fragment_count;
}

sub short_run_name {
    my $self = shift;
    unless($self->run_name eq $self->id){
        my (@names) = split('-',$self->run_name);
        return $names[-1];
    }
    return $self->run_name;
}

sub flow_cell_id {
    my $self = shift;
    return $self->short_run_name;
}

sub lane {
    my $self = shift;

    my $lane_attribute = $self->attributes(attribute_label => 'lane');
    if ( $lane_attribute ) {
        return $lane_attribute->attribute_value;
    }

    my $subset_name = $self->subset_name;
    if (($subset_name =~ m/DACC/ || $self->sequencing_platform eq 'solexa') && $subset_name =~/[-\.]/){
        my ($lane) = $subset_name =~ /(\d)[-\.]/;
        return $lane;
    }else{
        return $subset_name;
    }
}

sub run_start_date_formatted {
    return Genome::Model::Tools::Sam->time();
}

sub seq_id {
    my $self = shift;
    return $self->id;
}

sub instrument_data_id {
    my $self = shift;
    return $self->id;
}

sub resolve_quality_converter {
    my $self = shift;

    if ($self->import_format eq "solexa fastq") {
        return 'sol2sanger';
    } elsif ($self->import_format eq "illumina fastq") {
        return 'sol2phred';
    } elsif ($self->import_format eq 'sanger fastq') {
        return 'none';
    } else {
        $self->error_message("cannot resolve quality convertor for import format of type " . $self->import_format);
        die $self->error_message;
    }
}

sub gerald_directory {
    undef;
}

sub desc {
    my $self = shift;
    return $self->description || "[unknown]";
}

sub is_external {
    0;
}

sub resolve_adaptor_file {
 return '/gscmnt/sata114/info/medseq/adaptor_sequences/solexa_adaptor_pcr_primer';
}

sub run_identifier {
 my $self = shift;
 return $self->id;
}

sub _archive_file_name { # private for now...can be public
    my $self = shift;

    my $format = $self->import_format;
    if ( $format =~ /fastq/ ){
        return 'archive.tgz';
    }
    elsif ( $format eq 'bam' ){
        return 'all_sequences.bam';
    }
    elsif ( $format eq 'sff' ){
        return 'all_sequences.sff';
    }
    elsif ( $format eq 'raw sra download' ) {
        my $alloc = $self->disk_allocation;
        die $self->error_message("Couldn't get allocations")
            if not $alloc;
        my $base_dir = $alloc->absolute_path;

        opendir(my $base_dh, $base_dir)
            || die $self->error_message("Couldn't open $base_dir: $!");

        my @base_files = grep {not ($_ eq '.' or $_ eq '..')} readdir ($base_dh);
        closedir($base_dh);

        if ( (scalar(@base_files) != 1) and (not -d $base_files[0]) ) {
            die $self->error_message(
                "Could not understand directory structure of raw sra download in $base_dir.\n"
                . "Expected one directory, but found something else."
            );
        }

        my $sub_dir = $base_files[0];

        opendir(my $sub_dh, "$base_dir/$sub_dir")
            || die $self->error_message("Couldn't open $base_dir/$sub_dir: $!");

        my @sub_files = grep {not ($_ eq '.' or $_ eq '..')} readdir ($sub_dh);
        closedir($sub_dh);

        if (scalar(@sub_files) == 1) {
            my $sra_file = "$sub_dir/$sub_files[0]";

            die $self->error_message("$sra_file had zero size!")
                if (not -s "$base_dir/$sra_file");

            die $self->error_message("$sra_file was not a .sra file!")
                if ($sra_file !~ /\.sra$/);

            return "$sub_dir/$sub_files[0]";
        } elsif (scalar(@sub_files) == 5) {

            my @expected_items = qw(col idx lock md md5);

            for my $item (@expected_items) {
                die $self->error_message("Missing $item in raw sra directory")
                    unless grep {$item eq $_} @sub_files;
            }

            warn $self->status_message(
                "Warning: Found a directory instead of a single .sra file:\n"
                . "$base_dir/$sub_dir."
            );

            return "$sub_dir";

        } else {
            die $self->error_message(
                "Did not find expected directory structure of raw sra download:\n"
                . "$base_dir/$sub_dir."
            );
        }
    }
    else {
        Carp::confess("Unknown import format: $format");
    }
}

sub archive_path {
    my $self = shift;

    my $alloc = $self->disk_allocation;
    return if not $alloc;

    my $file_name = $self->_archive_file_name;
    return $alloc->absolute_path.'/'.$file_name;
}

sub bam_path {
    my ($self) = @_;

    my ($allocation) = $self->disk_allocation;
    unless ($allocation) {
        $self->fatal_message("Found no disk allocation for imported instrument data %s, so cannot find bam!", $self->__display_name__);
    }

    my $bam_file = $allocation->absolute_path . "/all_sequences.bam";

    return $bam_file if -e $bam_file;

    my @possible_bams = glob($allocation->absolute_path . "/*.bam");
    if (@possible_bams == 1) {
        return $possible_bams[0];
    } elsif (@possible_bams > 1) {
        $self->fatal_message('Multiple BAMs found in allocation for imported instrument data: %s.  Name one "all_sequences.bam" to have it used here.', $self->__display_name__);
    }

    return; #some other format besides BAM (e.g., FASTQs)
}

sub get_read_groups_set {
    my ($self) = @_;

    my $bam_file = $self->bam_path;
    unless ($bam_file and -e $bam_file) {
        $self->error_message("Bam file $bam_file doesn't exist");
        die $self->error_message;
    }
    my $cmd = Genome::Model::Tools::Sam::ListReadGroups->create(input=>$bam_file,
            silence_output=>1);
    unless ($cmd->execute) {
        $self->error_message("Failed to run list read groups command for $bam_file");
        die $self->error_message;
    }

    my $read_groups = Set::Scalar->new($cmd->read_groups);
    return $read_groups;
}

sub get_segments {
    my $self = shift;

    my %options = @_;
    my $allow_blacklisted_segments = delete $options{allow_blacklisted_segments};

    my @unknown_options = keys %options;
    if (@unknown_options) {
        die $self->error_message('Unknown option(s): ' . join(', ', @unknown_options));
    }

    unless ($self->import_format eq "bam") {
        return ();
    }

    my $read_groups = $self->get_read_groups_set();
    unless ($allow_blacklisted_segments) {
        my $blacklisted_segments = Set::Scalar->new($self->blacklisted_segments);
        $read_groups = $read_groups - $blacklisted_segments;
    }

    return map {{segment_type=>'read_group', segment_id=>$_}} $read_groups->elements;
}

1;

