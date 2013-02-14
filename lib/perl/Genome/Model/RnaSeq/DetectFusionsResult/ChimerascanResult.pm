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

    my $tmp_dir = $self->_symlink_in_fastqs();
    if ($should_reuse_bam) {
        $self->status_message("Attempting to reuse BAMs from pipeline");

        if (in($self->version, @can_reuse_bam_in_versions)) {
            my $reusable_dir = File::Spec->join($self->temp_staging_directory, 'reusable');
            Genome::Sys->create_directory($reusable_dir);

            $self->status_message("Creating the 'both mates mapped' bam.");
            my ($both_mates_bam, $both_mates_index, $sorted_reheadered_bam) =
                    $self->_create_both_mates_bam($reusable_dir, $bowtie_version);

            $self->status_message("Creating the unmapped fastqs.");
            my ($unmapped_1, $unmapped_2) =
                    $self->_create_unmapped_fastqs($reusable_dir, $sorted_reheadered_bam);
            unlink($sorted_reheadered_bam);

            $self->status_message("Symlinking reusable parts into Chimerascan working dir.");
            $self->_symlink_in_files($tmp_dir, $both_mates_bam, $unmapped_1, $unmapped_2,
                        $both_mates_index);
        } else {
            die(sprintf("This version of chimerascan (%s) doesn't " .
                        "support reusing the alignment BAM.", $self->version));
        }
    }
}

sub _symlink_in_fastqs {
    my ($self) = @_;

    my ($fq1, $fq2) = @{$self->_fastq_files};
    my $tmp_dir = File::Spec->join($self->temp_staging_directory, 'tmp');
    Genome::Sys->create_directory($tmp_dir);
    Genome::Sys->create_symlink($fq1, File::Spec->join($tmp_dir, 'reads_1.fq'));
    Genome::Sys->create_symlink($fq2, File::Spec->join($tmp_dir, 'reads_2.fq'));

    return $tmp_dir;
}

sub _create_both_mates_bam {
    my ($self, $reusable_dir, $bowtie_version) = @_;

    $self->status_message("Creating bam re-headered with transcripts from the chimerascan-index");
    my $alignment_result = $self->alignment_result;

    my $index = $self->_get_index($bowtie_version);
    my $seqdict_file = $index->get_sequence_dictionary;
    my $input_bam_file = $alignment_result->bam_file;
    my $new_bam_header = $self->_get_new_bam_header($reusable_dir, $input_bam_file, $seqdict_file);
    my $reheadered_bam = File::Spec->join($reusable_dir, 'reheadered.bam');
    my $rv = Genome::Model::Tools::Picard::ReplaceSamHeader->execute(
        input_file => $input_bam_file,
        output_file => $reheadered_bam,
        header_file => $new_bam_header,
    );
    unless ($rv) {
        die('Failed to reheader alignment result BAM file!');
    }
    unlink($new_bam_header);


    $self->status_message("Sorting reheadered BAM by position");
    my $sorted_reheadered_bam = File::Spec->join($reusable_dir, 'sorted_reheadered.bam');
    $rv = Genome::Model::Tools::Picard::SortSam->execute(
        input_file => $reheadered_bam,
        output_file => $sorted_reheadered_bam,
        sort_order => 'coordinate',
        use_version => $self->picard_version,
    );
    unless ($rv) {
        die('Failed to sort by position!');
    }
    unlink($reheadered_bam);


    $self->status_message("Filtering to only reads where both mates map (for aligned_reads.bam)");
    my $both_mates_bam = File::Spec->join($reusable_dir, 'both_mates_mapped_sorted_aligned_reads.bam');
    my $view_cmd = "samtools view -F 12 -b -h $sorted_reheadered_bam > $both_mates_bam";
    Genome::Sys->shellcmd(
        cmd => $view_cmd,
        input_files => [$sorted_reheadered_bam],
        output_files => [$both_mates_bam],
    );

    # create index
    $self->status_message("Indexing BAM: $both_mates_bam");
    my $both_mates_index = $both_mates_bam . '.bai';
    my $index_cmd = Genome::Model::Tools::Picard::BuildBamIndex->create(
        input_file => $both_mates_bam,
        output_file => $both_mates_index,
        use_version => $self->picard_version,
    );

    $index_cmd->execute();
    unless (-e $both_mates_index) {
        die "Failed to create BAM index of \"$both_mates_bam\".";
    }
    return ($both_mates_bam, $both_mates_index, $sorted_reheadered_bam);
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
    my ($self, $reusable_dir, $sorted_reheadered_bam) = @_;

    # Run samtools to filter, include only alignments where both mates map (primary and non-primary).
    $self->status_message("Filtering to get all the other reads.");
    my $not_both_mates_bam = File::Spec->join($reusable_dir, 'unmapped_sorted_aligned_reads.bam');
    # sam format specifies bitmask 4 = read itself unmapped and 8 = read's mate is unmapped (we want either of them)
    my $gawk_cmd = q{and($2, 4) || and($2, 8) || substr($1, 0, 1) == "@"};
    my $view_cmd = "samtools view -h $sorted_reheadered_bam | gawk '$gawk_cmd' | " .
            "samtools view -S -b - > $not_both_mates_bam";
    Genome::Sys->shellcmd(
        cmd => $view_cmd,
        input_files => [$sorted_reheadered_bam],
        output_files => [$not_both_mates_bam],
    );

    $self->status_message("Generating the unaligned fastqs from $not_both_mates_bam");
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
    my ($self, $tmp_dir, $both_mates_bam, $unmapped_1, $unmapped_2, $both_mates_index) = @_;

    Genome::Sys->create_symlink($both_mates_bam,
            File::Spec->join($self->temp_staging_directory, 'aligned_reads.bam'));

    Genome::Sys->create_symlink($unmapped_1, File::Spec->join($tmp_dir, 'unaligned_1.fq'));
    Genome::Sys->create_symlink($unmapped_2, File::Spec->join($tmp_dir, 'unaligned_2.fq'));

    Genome::Sys->create_symlink($both_mates_bam,
            File::Spec->join($self->temp_staging_directory, 'sorted_aligned_reads.bam'));
    Genome::Sys->create_symlink($both_mates_index,
            File::Spec->join($self->temp_staging_directory, 'sorted_aligned_reads.bam.bai'));
}

sub _run_chimerascan {
    my ($self, $options, $arguments) = @_;

    local $ENV{PYTHONPATH} = ($ENV{PYTHONPATH} ? $ENV{PYTHONPATH} . ":" : "") .
            $self->_python_path_for_version($self->version);

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
