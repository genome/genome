package Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult;

use strict;
use warnings;

use above 'Genome';
use Genome::Utility::List 'in';

# Notes from Chris Miller:
# bsub -oo err.log -q long
# -M 16000000 -R 'select[type==LINUX64 && mem>16000] span[hosts=1]
#   rusage[mem=16000]' -n 8 -J chimera -oo outputdir/chimera.err
# "python /gsc/bin/chimerascan_run.py -v -p 8
#   /gscmnt/sata921/info/medseq/cmiller/annotations/chimeraScanIndex/
#   $fastq1 $fastq2 $outputdir"

class Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult {
    is => "Genome::Model::RnaSeq::DetectFusionsResult",
};

our %INDIRECT_PARAMETER_VALIDATORS = (
        '--bowtie-version' => '_validate_bowtie_version',
        '--reuse-bam' => '_validate_reuse_bam',
);

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    # TODO If process crashes before calling _reallocate_disk_allocation then
    # the output directory will remain at initially allocated size
    # indefinitely.
    $self->_prepare_output_directory();
    $self->_prepare_staging_directory();

    my $options = $self->_resolve_chimerascan_options();
    my $arguments = $self->_resolve_chimerascan_arguments($options->{indirect});
    $self->_prepare_to_run_chimerascan($options, $arguments);
    $self->_run_chimerascan($options, $arguments);

    $self->_promote_data();
    $self->_remove_staging_directory();
    $self->_reallocate_disk_allocation();

    return $self;
}

sub _resolve_chimerascan_options {
    my ($self) = @_;

    my $detector_params = $self->detector_params;
    my %options;
    my ($direct_params, $indirect_params) =
            $self->_preprocess_detector_params($detector_params);
    $options{direct} = $direct_params;
    $options{indirect} = $indirect_params;

    return \%options;
}

sub _resolve_chimerascan_arguments {
    my ($self, $indirect_params) = @_;

    my $bowtie_version = $indirect_params->{'--bowtie-version'};
    my $index_dir = $self->_resolve_index_dir($bowtie_version);
    my $bowtie_path = Genome::Model::Tools::Bowtie->base_path($bowtie_version);

    my ($fastq1, $fastq2) = @{$self->_fastq_files};
    my $output_directory = $self->temp_staging_directory;

    my @arguments = ($index_dir, $fastq1, $fastq2, $output_directory);
    return \@arguments;
}

sub _prepare_to_run_chimerascan {
    my ($self, $options, $arguments) = @_;

    my @can_reuse_bam_in_versions = qw(0.4.3 0.4.5);
    my $bowtie_version = $options->{indirect}->{'--bowtie-version'};
    my $should_reuse_bam = $options->{indirect}->{'--reuse-bam'};

    if ($should_reuse_bam) {
        $self->status_message("Attempting to reuse BAMs from pipeline");
        if (in($self->version, @can_reuse_bam_in_versions)) {
            # These are intermediate files required for downstream conversion and filtering
            $self->_create_reheadered_queryname_sorted_primary_bam_file($bowtie_version);
            # simply create symlinks for the original FASTQ files
            $self->_link_in_converted_fastqs($arguments);
            # includes sorted_aligned_reads and .bai
            $self->_link_in_aligned_reads_bam();
            # convert unaligned read pairs, one end unaligned, from BAM to FASTQ
            $self->_create_unaligned_fastqs();
        } else {
            die(sprintf("This version of chimerascan (%s) doesn't " .
                        "support reusing the alignment BAM.", $self->version));
        }
    }
}

# Add all the @SQ lines from <suppliment_header_file> to the header of the <source_bam>
# Returns: the filename of the newly created header
sub _get_new_bam_header {
    my ($self, $source_bam, $suppliment_header_file) = @_;

    die "Source bam doesn't exist at ($source_bam)" unless -e $source_bam;
    die "suppliment_header not supplied" unless -e $suppliment_header_file;

    my $new_bam_header = File::Spec->join($self->temp_staging_directory, 'new_bam_header.seqdict');
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

sub _create_reheadered_queryname_sorted_primary_bam_file {
    my ($self,$bowtie_version) = @_;

    $self->status_message("Generating queryname sorted BAM with primary alignments only");
    my $alignment_result = $self->alignment_result;

    # Use the alignment_result BAM as the input, the chimerascan index seqdict .sam as the header, and output as the aligned_reads
    my $index = $self->_get_index($bowtie_version);
    my $seqdict_file = $index->get_sequence_dictionary;
    my $input_bam_file = $alignment_result->bam_file;
    my $new_bam_header = $self->_get_new_bam_header($input_bam_file, $seqdict_file);
    my $reheadered_bam_file = $self->temp_staging_directory .'/reheadered_aligned_reads.bam';
    unless (Genome::Model::Tools::Picard::ReplaceSamHeader->execute(
        input_file => $input_bam_file,
        output_file => $reheadered_bam_file,
        header_file => $new_bam_header,
    )) {
        die('Failed to reheader alignment result BAM file!');
    }

    my $primary_bam_file = join("/", $self->temp_staging_directory,'reheadered_primary_aligned_reads.bam');
    # Run samtools to filter, include primary alignments only.
    my $view_cmd = 'samtools view -F 256 -b '. $reheadered_bam_file .' > '. $primary_bam_file;
    Genome::Sys->shellcmd(
        cmd => $view_cmd,
        input_files => [$reheadered_bam_file],
        output_files => [$primary_bam_file],
    );
    unlink($reheadered_bam_file);

    my $queryname_sorted_file = join("/", $self->temp_staging_directory,'reheadered_primary_aligned_reads_queryname_sorted.bam');
    # Queryname sort the BAM for later FilterSamReads to get aligned/unaligned
    unless (Genome::Model::Tools::Picard::SortSam->execute(
        input_file => $primary_bam_file,
        output_file => $queryname_sorted_file,
        sort_order => 'queryname',
        use_version => $self->picard_version,
    )) {
        die('Failed to sort by queryname!');
    }
}

sub _link_in_aligned_reads_bam {
    my ($self) = @_;

    $self->status_message("Generating aligned-reads BAM from primary queryname sorted BAM");
    my $queryname_sorted_file = join("/", $self->temp_staging_directory,'reheadered_primary_aligned_reads_queryname_sorted.bam');
    my $output_file = join("/", $self->temp_staging_directory,'aligned_reads.bam');

    my $cmd = Genome::Model::Tools::Picard::FilterSamReads->create(
        input_file  => $queryname_sorted_file,
        output_file => $output_file,
        filter      => 'includeAligned',
        sort_order  => 'coordinate',
        use_version => $self->picard_version,
    );

    $cmd->execute();
    unless (-e $output_file) {
        die "Error generating BAM of aligned reads from alignment_result.";
    }

    my $sorted_bam = join("/", $self->temp_staging_directory,
            "sorted_aligned_reads.bam");
    $self->status_message("Creating symlink to alignment result BAM (pretending " .
            "also to be the sorted-aligned-reads BAM.");
    Genome::Sys->create_symlink($output_file, $sorted_bam);
    my $index_cmd = Genome::Model::Tools::Picard::BuildBamIndex->create(
        input_file => $sorted_bam,
        output_file => $sorted_bam .'.bai',
        use_version => $self->picard_version,
    );

    $index_cmd->execute();
    unless (-e "$sorted_bam.bai") {
        die "Failed to create BAM index of \"$sorted_bam\".";
    }
}

sub _link_in_converted_fastqs {
    my ($self) = @_;

    my ($fastq1, $fastq2) = @{$self->_fastq_files};
    my $chimerascan_tmp_dir = $self->temp_staging_directory .'/tmp';
    unless (-d $chimerascan_tmp_dir) {
        Genome::Sys->create_directory($chimerascan_tmp_dir);
    }
    my $converted_fastq_prefix = join('/', $chimerascan_tmp_dir, 'reads_');
    $self->status_message("Creating symlink to fastq files (pretending " .
            "to be the converted-fastq files");
    Genome::Sys->create_symlink($fastq1,$converted_fastq_prefix .'1.fq');
    Genome::Sys->create_symlink($fastq2,$converted_fastq_prefix .'2.fq');
}

sub _create_unaligned_fastqs {
    my ($self) = @_;

    $self->status_message("Generating the unaligned fastqs from the " .
            "alignment results BAM");
    my $stage            = $self->temp_staging_directory;

    my $queryname_sorted_file = join("/", $stage,'reheadered_primary_aligned_reads_queryname_sorted.bam');
    unless (-s $queryname_sorted_file) {
        die('Failed to find the primary queryname sorted BAM!');
    }

    my $output_file = $stage .'/unaligned_reads.bam';

    my $filter_cmd = Genome::Model::Tools::Picard::FilterSamReads->create(
        input_file  => $queryname_sorted_file,
        output_file => $output_file,
        filter      => 'excludeAligned',
        use_version => $self->picard_version,
    );

    $filter_cmd->execute();
    unless (-e $output_file) {
        die "Error generating BAM of unaligned reads from the primary queryname sorted BAM.";
    }
    unlink($queryname_sorted_file);
    my $read1_fastq      = join("/", $stage, 'tmp', 'unaligned_1.fq');
    my $read2_fastq      = join("/", $stage, 'tmp', 'unaligned_2.fq');

    my $cmd = Genome::Model::Tools::Picard::StandardSamToFastq->create(
        input            => $output_file,
        fastq            => $read1_fastq,
        second_end_fastq => $read2_fastq,
        re_reverse       => 1,
        include_non_pf_reads => 1,
        include_non_primary_alignments => 0,
        use_version      => $self->picard_version,
        maximum_memory => 16,
        maximum_permgen_memory => 256,
    );

    unless ($cmd->execute()) { die ('Failed to convert unaligned BAM to FASTQ!'); }
}

sub _run_chimerascan {
    my ($self, $options, $arguments) = @_;

    local $ENV{PYTHONPATH} = ($ENV{PYTHONPATH} ? $ENV{PYTHONPATH} . ":" : "") .
            $self->_python_path_for_version($self->version);

    warn Data::Dumper::Dumper(($options, $arguments));
    my $cmd_path = $self->_path_for_version($self->version);
    unless ($cmd_path) {
        die $self->error_message("Failed to find a path for chimerascan for " .
                "version " . $self->version . "!");
    }
    my $direct = $options->{direct};
    my $indirect = $options->{indirect};
    my $bowtie_version = $indirect->{'--bowtie-version'};
    my $bowtie_path = Genome::Model::Tools::Bowtie->base_path($bowtie_version);

    my $arguments_str = join(" ", @{$arguments});
    my $output_directory = $arguments->[-1];
    my $cmd = "python $cmd_path/chimerascan_run.py -v $direct " .
              "--bowtie-path=" . $bowtie_path .
              " $arguments_str > $output_directory/chimera_result.out";

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => $arguments,
    );
}

# return the detector params (sent directly to detector) and a hash of
# indirect parameters
sub _preprocess_detector_params {
    my ($self, $params) = @_;

    my %indirect_parameters;
    for my $name (keys %INDIRECT_PARAMETER_VALIDATORS) {
        # \Q$foo\E ensures that regex symbols are 'quoted'
        if($params and $params =~ m/(\s*\Q$name\E[=\s]([^\s]*)\s*)/) {
            my $str = $1;
            my $val = $2;
            $params =~ s/\Q$str\E/ /;

            my $validation_method_name = $INDIRECT_PARAMETER_VALIDATORS{$name};
            $self->$validation_method_name($val, $params); # dies if invalid

            $indirect_parameters{$name} = $val;
        } else {
            die(sprintf("Couldn't find parameter named \"%s\" in \"%s\"",
                    $name, $params || ''));
        }
    }
    return $params, \%indirect_parameters;
}

sub _validate_bowtie_version {
    my ($self, $val, $params) = @_;

    my $bowtie_version = $val ||
            die("You must supply a bowtie version in the detector parameters " .
                "in the form of \"--bowtie-version=<version>\" Got detector " .
                "parameters: [\"$params\"]");
    my ($major_version) = split(/\./, $bowtie_version);
    if ($major_version ne 0) {
        die("Chimerascan currently only supports bowtie major version 0, " .
            "not $major_version");
    }
}

sub _validate_reuse_bam {
    my ($self, $val, $params) = @_;

    unless ($val eq 1 or $val eq 0) {
        die "You must specify either 1 (true) or 0 (false) for parameter " .
                \"--reuse-bam\", you specified \"$val\"";
    }
}

sub _staging_disk_usage {
    #return enough to cover our temp dir (which will be under staging)
    return 60 * 1024 * 1024;
}

sub _path_for_version {
    my ($class, $version) = @_;

    my ($base_path, $bin_sub, $python_sub) =
            _get_chimerascan_path_for_version($version);
    my $result = join('/', $base_path, $bin_sub);

    if (-e $result) {
        return $result;
    } else {
        die("Binary path ($bin_sub) does not exist under chimerascan path " .
                "($base_path)!");
    }
}

sub _python_path_for_version {
    my ($class, $version) = @_;

    my ($base_path, $bin_sub, $python_sub) =
            _get_chimerascan_path_for_version($version);
    my $result = join('/', $base_path, $python_sub);

    if (-e $result) {
        return $result;
    } else {
        die("Python path ($python_sub) does not exist under chimerascan " .
            "path ($base_path)!");
    }
}

sub _get_chimerascan_path_for_version {
    my ($version) = @_;

    my $path = sprintf("/usr/lib/chimerascan%s", $version);
    if (-e $path) {
        # installed via deb-package method
        my $bin_sub = 'bin';
        my $python_sub = join('/', 'lib', 'python2.6', 'site-packages');
        return $path, $bin_sub, $python_sub;
    } else {
        # fall back to old installation method
        $path = $ENV{GENOME_SW} . "/chimerascan/chimerascan-$version";
        if (-e $path) {
            my $bin_sub = 'chimerascan';
            my $python_sub = join('/', 'build', 'lib.linux-x86_64-2.6');
            return $path, $bin_sub, $python_sub;
        } else {
            die("You requested an unavailable version of Chimerascan. " .
                "Requested: $version");
        }
    }
}

sub _resolve_index_dir {
    my ($self, $bowtie_version) = @_;

    my $index = $self->_get_index($bowtie_version);

    unless ($index) {
        #We want to shell out and create the chimerscan index in a different UR context.
        #That way it is committed, even if we fail and other builds don't need to wait on
        #the overall chimerscan run just to get their index results.
        my $cmd = 'genome model rna-seq detect-fusions chimerascan-index';
        $cmd .= ' --version=' . $self->version;
        $cmd .= ' --bowtie-version=' . $bowtie_version;
        $cmd .= ' --reference-build=' . $self->alignment_result->reference_build->id;
        $cmd .= ' --annotation-build=' . $self->annotation_build->id;

        Genome::Sys->shellcmd(cmd => $cmd);

        # Force UR to query the datasource instead of using its cache for this lookup.
        my $previous_value = UR::Context->query_underlying_context;
        UR::Context->query_underlying_context(1);
        $index = $self->_get_index($bowtie_version);
        UR::Context->query_underlying_context($previous_value);
    }

    if ($index) {
        $self->status_message('Registering software result ' . $self->id . ' (self) as a user " .
                "of the generated index');
        $index->add_user(user => $self, label => 'uses');
    } else {
        die("Unable to get a chimerascan index result");
    }

    return $index->output_dir;
}

sub _get_index {
    my ($self, $bowtie_version) = @_;

    my $index_class = 'Genome::Model::RnaSeq::DetectFusionsResult' .
                      '::ChimerascanResult::Index';
    my $index = $index_class->get_with_lock(
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
        version => $self->version,
        bowtie_version => $bowtie_version,
        reference_build => $self->alignment_result->reference_build,
        annotation_build => $self->annotation_build,
    );

    return $index;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    return 'build_merged_alignments/chimerascan/' . $self->id;
}

1;
