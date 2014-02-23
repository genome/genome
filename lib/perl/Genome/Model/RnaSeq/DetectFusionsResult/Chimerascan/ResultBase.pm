package Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::ResultBase;

use strict;
use warnings;

use above 'Genome';
use Genome::Utility::List 'in';
use File::Path qw();
use File::Spec;

class Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::ResultBase {
    is => "Genome::Model::RnaSeq::DetectFusionsResult",
    has_input => [
        original_bam_paths => {
            is => "Text",

            is_many => 1,
            doc => "The path(s) to the original instrument_data BAM files."
        },
        reuse_bam => {
            is => 'Boolean',
        },
        bowtie_version => {
            is => 'Text',
        },
    ],
    has => [
        bedpe_file => {
            is => 'Path',
            calculate_from => 'disk_allocations',
            calculate => q| File::Spec->join($disk_allocations->absolute_path, 'chimeras.bedpe'); |,

        },
    ],
};

sub _run_chimerascan {
    my ($bowtie_version, $c_pargs, $c_opts) = @_;
    die("Must be defined in subclass");
}

sub _chimerascan_index_cmd {
    die("Must be defined in subclass");
}

sub _chimerascan_result_class {
    die("Must be defined in subclass");
}


sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    $self->_prepare_output_directory();
    $self->_prepare_staging_directory();

    my ($fastq1, $fastq2, $qname_sorted_bam) = $self->_resolve_original_files($self->reuse_bam);
    my $index_dir = $self->_resolve_index_dir($self->bowtie_version);

    my $c_pargs = $self->_resolve_chimerascan_positional_arguments($index_dir, $fastq1, $fastq2);

    $self->_prepare_to_run_chimerascan($self->bowtie_version, $self->reuse_bam, $fastq1, $fastq2, $qname_sorted_bam);

    $self->_run_chimerascan($self->bowtie_version, $c_pargs, $self->detector_params);

    # cleanup $fastq1, $fastq2, and $qname_sorted_bam?
    if (-f $fastq1) {
        unless (unlink $fastq1) {
            $self->warning_message("Failed to unlink fastq1: $fastq1: $!");
        }
    }
    if (-f $fastq2) {
        unless (unlink $fastq2) {
            $self->warning_message("Failed to unlink fastq2: $fastq2: $!");
        }
    }

    # This is the queryname sorted BAM only used when --reuse-bam=1.  It's a copy of the original alignment result(ex. TopHat)
    if (-f $qname_sorted_bam) {
        unless (unlink $qname_sorted_bam) {
            $self->warning_message("Failed to unlink qname_sorted_bam: $qname_sorted_bam: $!");
        }
    }

    my $unsorted_aligned_reads_bam = File::Spec->join($self->temp_staging_directory, 'aligned_reads.bam');
    # This is an unsorted Chimerascan BAM file that is no longer needed since a sorted and indexed version exists
    if (-f $unsorted_aligned_reads_bam) {
        unless (unlink $unsorted_aligned_reads_bam) {
            $self->warning_message("Failed to unlink unsorted aligned reads BAM: $unsorted_aligned_reads_bam: $!");
        }
    }

    $self->_promote_data();
    $self->_remove_staging_directory();
    $self->_reallocate_disk_allocation();

    return $self;
}

sub _resolve_original_files {
    my ($self, $reuse_bam) = @_;

    my @fastq_files;
    if ($reuse_bam) {
        return $self->_resolve_original_files_reusing_bam();
    } else {
        return ($self->_resolve_original_fastq_files, undef);
    }
}

sub _resolve_original_fastq_files {
    my ($self) = @_;

    # When do we actually have mulitple original BAMs?  Shouldn't the original BAM file be the merged alignment result?
    my @original_bam_paths = $self->original_bam_paths;
    unless (@original_bam_paths) {
        die("Couldn't find 'original_bam_paths' to make fastq files!");
    }

    # queryname sort each BAM and get fastq1 and fastq2 from /tmp queryname sorted BAMs
    my (@fastq1_files, @fastq2_files);
    my $counter = 0; #make sure the original_fastq names are unique
    for my $bam_path (@original_bam_paths) {
        my $queryname_sorted_bam = File::Spec->join($self->temp_staging_directory,
                'original_queryname_sorted.bam');
        $self->_qname_sort_bam($bam_path, $queryname_sorted_bam);

        my $fastq1 = File::Spec->join($self->temp_staging_directory, "original_fastq1.$counter");
        my $fastq2 = File::Spec->join($self->temp_staging_directory, "original_fastq2.$counter");
        $counter++;

        # make fastqs from the qname sorted bam
        $self->_convert_bam_to_fastqs($queryname_sorted_bam, $fastq1, $fastq2);

        unless(unlink($queryname_sorted_bam)){
            $self->debug_message("Failed to remove temporary sorted bam file: $queryname_sorted_bam.  Reason: $!");
        }

        push @fastq1_files, $fastq1;
        push @fastq2_files, $fastq2;

    }

    # concatinate forward fastqs together
    my $fastq1 = File::Spec->join($self->temp_staging_directory, 'fastq1');
    my $cmd = sprintf('cat %s > %s', join(" ", @fastq1_files), $fastq1);
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [@fastq1_files],
        output_files => [$fastq1],
    );
    
    # concatinate reverse fastqs together
    my $fastq2 = File::Spec->join($self->temp_staging_directory, 'fastq2');
    $cmd = sprintf('cat %s > %s', join(" ", @fastq2_files), $fastq2);
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [@fastq2_files],
        output_files => [$fastq2],
    );

    #cleanup temporary fastq files
    foreach my $fastq (@fastq1_files, @fastq2_files){
        unless(unlink($fastq)){
            $self->debug_message("Failed to remove temporary fastq file: $fastq.  Reason: $!");
        }
    }

    return ($fastq1, $fastq2);
}

sub _prepare_to_run_chimerascan {
    my ($self, $bowtie_version, $reuse_bam, $fastq1, $fastq2, $qname_sorted_bam) = @_;

    if ($reuse_bam) {
        my $tmp_dir = $self->_symlink_in_fastqs($fastq1, $fastq2);
        $self->debug_message("Attempting to reuse BAMs from pipeline");

        my $reusable_dir = File::Spec->join($self->temp_staging_directory,
                'reusable');
        Genome::Sys->create_directory($reusable_dir);

        $self->debug_message("Creating the reheadered bam.");
        my $reheadered_bam = $self->_create_reheadered_bam($reusable_dir,
                $bowtie_version, $qname_sorted_bam);

        $self->debug_message("Creating the 'both mates mapped' bam.");
        my ($both_mates_bam, $sorted_both_mates_bam, $sorted_both_mates_index) =
                $self->_create_both_mates_bam($reusable_dir, $reheadered_bam);

        $self->debug_message("Creating the unmapped fastqs.");
        my ($unmapped_1, $unmapped_2) =
                $self->_create_unmapped_fastqs($reusable_dir, $reheadered_bam);
        unlink($reheadered_bam);

        $self->debug_message("Symlinking reusable parts into Chimerascan working dir.");
        $self->_symlink_in_files($tmp_dir, $both_mates_bam, $unmapped_1,
                $unmapped_2, $sorted_both_mates_bam, $sorted_both_mates_index);
    }
}

sub _resolve_original_files_reusing_bam {
    my ($self) = @_;

    # query name sort the bam
    my $alignment_result = $self->alignment_result;

    my $bam_file = $alignment_result->bam_file || die (
            "Couldn't get BAM file from alignment result (" . $alignment_result->id . ")");

    # No reason to write to staging directory, but rather /tmp that will be cleaned up
    my $tmp_dir = Genome::Sys->create_temp_directory();
    my $queryname_sorted_bam = File::Spec->join($tmp_dir,
            'alignment_result_queryname_sorted.bam');
    $self->_qname_sort_bam($bam_file, $queryname_sorted_bam);

    # make fastqs from the qname sorted bam
    my $fastq1 = $self->temp_staging_directory . "/fastq1";
    my $fastq2 = $self->temp_staging_directory . "/fastq2";

    $self->_convert_bam_to_fastqs($queryname_sorted_bam, $fastq1, $fastq2);
    $alignment_result->add_user(label => 'uses' , user => $self);
    return ($fastq1, $fastq2, $queryname_sorted_bam);
}

sub _convert_bam_to_fastqs {
    my ($self, $bam_path, $fastq1_path, $fastq2_path) = @_;

    my $cmd = Genome::Model::Tools::Picard::StandardSamToFastq->create(
        input                          => $bam_path,
        fastq                          => $fastq1_path,
        second_end_fastq               => $fastq2_path,
        re_reverse                     => 1,
        include_non_pf_reads           => 1,
        include_non_primary_alignments => 0,
        use_version                    => $self->picard_version,
        maximum_memory                 => 16,
        maximum_permgen_memory         => 256,
        max_records_in_ram             => 2_500_000,
        additional_jvm_options         => '-XX:-UseGCOverheadLimit',
    );
    unless ($cmd->execute){
        die ("Failed to convert BAM file to fastq: $bam_path");
    }
    unless (-s $fastq1_path and -s $fastq2_path) {
        die ("Fastq files were not found after conversion from BAM file: $bam_path");
    }
}

sub _qname_sort_bam {
    my ($self, $bam_path, $output_bam_path) = @_;

    my $rv = Genome::Model::Tools::Picard::SortSam->execute(
        input_file             => $bam_path,
        output_file            => $output_bam_path,
        sort_order             => 'queryname',
        max_records_in_ram     => 2_500_000,
        maximum_permgen_memory => 256,
        maximum_memory         => 16,
        use_version            => $self->picard_version,
    );
    unless ($rv) {
        die ("Failed to queryname sort BAM file: $bam_path");
    }

    return $output_bam_path;
}

sub _symlink_in_fastqs {
    my ($self, $fastq1, $fastq2) = @_;

    my $tmp_dir = File::Spec->join($self->temp_staging_directory, 'tmp');
    Genome::Sys->create_directory($tmp_dir);
    Genome::Sys->create_symlink($fastq1, File::Spec->join($tmp_dir, 'reads_1.fq'));
    Genome::Sys->create_symlink($fastq2, File::Spec->join($tmp_dir, 'reads_2.fq'));

    return $tmp_dir;
}

sub _create_both_mates_bam {
    my ($self, $reusable_dir, $reheadered_bam) = @_;

    $self->debug_message("Filtering to only reads where both mates map (for aligned_reads.bam)");
    my $both_mates_bam = File::Spec->join($reusable_dir, 'both_mates.bam');
    my $view_cmd = "samtools view -F 12 -b -h $reheadered_bam > $both_mates_bam";
    Genome::Sys->shellcmd(
        cmd => $view_cmd,
        input_files => [$reheadered_bam],
        output_files => [$both_mates_bam],
    );

    $self->debug_message("Sorting both_mates_bam by position");
    my $sorted_both_mates_bam = File::Spec->join($reusable_dir, 'sorted_both_mates.bam');
    my $rv = Genome::Model::Tools::Picard::SortSam->execute(
        input_file => $both_mates_bam,
        output_file => $sorted_both_mates_bam,
        sort_order => 'coordinate',
        use_version => $self->picard_version,
    );
    unless ($rv) {
        die('Failed to sort by position!');
    }

    # create index
    $self->debug_message("Indexing BAM: $sorted_both_mates_bam");
    my $sorted_both_mates_index = $sorted_both_mates_bam . '.bai';
    my $index_cmd = Genome::Model::Tools::Picard::BuildBamIndex->create(
        input_file => $sorted_both_mates_bam,
        output_file => $sorted_both_mates_index,
        use_version => $self->picard_version,
    );


    $index_cmd->execute();
    unless (-e $sorted_both_mates_index) {
        die "Failed to create BAM index of \"$sorted_both_mates_bam\".";
    }
    return ($both_mates_bam, $sorted_both_mates_bam, $sorted_both_mates_index);
}

sub _create_reheadered_bam {
    my ($self, $reusable_dir, $bowtie_version, $qname_sorted_bam) = @_;

    $self->debug_message("Creating bam re-headered with transcripts from the chimerascan-index");
    my $index = $self->_get_index($bowtie_version);
    my $seqdict_file = $index->get_sequence_dictionary;
    my $new_bam_header = $self->_get_new_bam_header($reusable_dir, $qname_sorted_bam, $seqdict_file);
    my $reheadered_bam = File::Spec->join($reusable_dir, 'reheadered.bam');
    my $rv = Genome::Model::Tools::Picard::ReplaceSamHeader->execute(
        input_file => $qname_sorted_bam,
        output_file => $reheadered_bam,
        header_file => $new_bam_header,
    );
    unless ($rv) {
        die('Failed to reheader alignment result BAM file!');
    }
    unlink($new_bam_header);

    return $reheadered_bam;
}


# Add all the @SQ lines from <suppliment_header_file> to the header of the <source_bam>
# Returns: the filename of the newly created header
sub _get_new_bam_header {
    my ($self, $reusable_dir, $source_bam, $suppliment_header_file) = @_;

    die "Source bam doesn't exist at ($source_bam)" unless -e $source_bam;
    die "suppliment_header not supplied" unless -e $suppliment_header_file;

    my $new_bam_header = File::Spec->join($reusable_dir, 'new_bam_header.seqdict');
    # awk '!x[$1 $2]++' prints uniq lines (based on first two columns) in order they originally appeared
    my $sq_lines_only = q($1 == "@SQ");
    my $cmd = "cat <(samtools view -H $source_bam) " .
              "<(awk '$sq_lines_only' $suppliment_header_file) " .
              "| awk '!x[\$1 \$2]++' > $new_bam_header";
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$source_bam, $suppliment_header_file],
        output_files => [$new_bam_header],
    );

    return $new_bam_header
}


sub _create_unmapped_fastqs {
    my ($self, $reusable_dir, $reheadered_bam) = @_;

    # Run samtools to filter, include only reads where one or both mates didn't map (primary and non-primary).
    $self->debug_message("Filtering to get all the other reads.");
    my $not_both_mates_bam = File::Spec->join($reusable_dir, 'not_both_mates.bam');
    # sam format specifies bitmask 4 = read itself unmapped and 8 = read's mate is unmapped (we want either of them)
    my $gawk_cmd = q{substr($1, 0, 1) == "@" || and($2, 4) || and($2, 8)};
    my $view_cmd = "samtools view -h $reheadered_bam | gawk '$gawk_cmd' | " .
            "samtools view -S -b - > $not_both_mates_bam";
    Genome::Sys->shellcmd(
        cmd => $view_cmd,
        input_files => [$reheadered_bam],
        output_files => [$not_both_mates_bam],
    );

    $self->debug_message("Generating the unaligned fastqs from $not_both_mates_bam");
    my $unmapped_1 = File::Spec->join($reusable_dir, 'unmapped_1');
    my $unmapped_2 = File::Spec->join($reusable_dir, 'unmapped_2');

    my $cmd = Genome::Model::Tools::Picard::StandardSamToFastq->create(
        input                          => $not_both_mates_bam,
        fastq                          => $unmapped_1,
        second_end_fastq               => $unmapped_2,
        re_reverse                     => 1,
        include_non_pf_reads           => 1,
        include_non_primary_alignments => 0,
        use_version                    => $self->picard_version,
        maximum_memory                 => 16,
        maximum_permgen_memory         => 256,
    );
    unless ($cmd->execute()) { die ('Failed to convert unaligned BAM to FASTQ!'); }

    return ($unmapped_1, $unmapped_2);
}

sub _symlink_in_files {
    my ($self, $tmp_dir, $both_mates_bam, $unmapped_1, $unmapped_2,
            $sorted_both_mates_bam, $sorted_both_mates_index) = @_;

    Genome::Sys->create_symlink($both_mates_bam,
            File::Spec->join($self->temp_staging_directory, 'aligned_reads.bam'));

    Genome::Sys->create_symlink($unmapped_1, File::Spec->join($tmp_dir, 'unaligned_1.fq'));
    Genome::Sys->create_symlink($unmapped_2, File::Spec->join($tmp_dir, 'unaligned_2.fq'));

    Genome::Sys->create_symlink($sorted_both_mates_bam,
            File::Spec->join($self->temp_staging_directory, 'sorted_aligned_reads.bam'));
    Genome::Sys->create_symlink($sorted_both_mates_index,
            File::Spec->join($self->temp_staging_directory, 'sorted_aligned_reads.bam.bai'));
}

sub _resolve_chimerascan_positional_arguments {
    my ($self, $index_dir, $fastq1, $fastq2) = @_;

    my $output_directory = $self->temp_staging_directory;

    my @c_pargs = ($index_dir, $fastq1, $fastq2, $output_directory);
    return \@c_pargs;
}

sub _staging_disk_usage {
    #return enough to cover our temp dir (which will be under staging)
    return 60 * 1024 * 1024;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    return 'build_merged_alignments/chimerascan/' . $self->id;
}

sub _resolve_index_dir {
    my ($self, $bowtie_version) = @_;

    my $index = $self->_get_index($bowtie_version);

    if ($index) {
        $self->debug_message(sprintf('Registering software result %s as a user ' .
                'of the generated index (%s)', $self->id, $index->id));
        $index->add_user(user => $self, label => 'uses');
    } else {
        die("Unable to get a Chimerascan/ChimerascanVrl result");
    }

    return $index->output_dir;
}

sub _get_index {
    my ($self, $bowtie_version, $query_underlying_context) = @_;

    (my $index_class = $self->_chimerascan_result_class) =~ s/::Result$/::Index/;
    my %params = (
        test_name => $ENV{GENOME_ALIGNER_INDEX_TEST_NAME} || undef,
        version => $self->version,
        bowtie_version => $bowtie_version,
        reference_build => $self->alignment_result->reference_build,
        annotation_build => $self->annotation_build,
        picard_version => $self->picard_version,
    );

    my $index;
    if ($query_underlying_context) {
        # Force UR to query the datasource instead of using its cache for this lookup.
        my $previous_value = UR::Context->query_underlying_context;
        UR::Context->query_underlying_context(1);
        $index = $index_class->get_with_lock(%params);
        UR::Context->query_underlying_context($previous_value);
    } else {
        $index = $index_class->get_with_lock(%params);
    }

    return $index;
}

1;

