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
    my $should_reuse_bam = $options->{indirect}->{'--reuse-bam'};

    if ($should_reuse_bam) {
        $self->status_message("Attempting to reuse BAMs from pipeline");
        if (in($self->version, @can_reuse_bam_in_versions)) {
            # includes sorted_aligned_reads and .bai
            $self->_link_in_aligned_reads_bam();

            $self->_link_in_converted_fastqs($arguments);
            $self->_create_unaligned_fastqs();
        } else {
            die(sprintf("This version of chimerascan (%s) doesn't " .
                        "support reusing the alignment BAM.", $self->version));
        }
    }
}

sub _link_in_aligned_reads_bam {
    my ($self) = @_;

    $self->status_message("Generating aligned-reads BAM from merged BAM");
    my $alignment_result = $self->alignment_result;
    my $input_file = $alignment_result->bam_file;
    my $output_file = join("/", $self->temp_staging_directory,
            "aligned_reads.bam");

    my $cmd = Genome::Model::Tools::Picard::FilterSamReads->create(
        input_file  => $input_file,
        output_file => $output_file,
        filter      => 'includeAligned',
    );

    $cmd->execute();
    unless (-e $output_file) {
        die "Error generating BAM of unaligned reads from alignment_result.";
    }

    my $sorted_bam = join("/", $self->temp_staging_directory,
            "sorted_aligned_reads.bam");
    $self->status_message("Creating symlink to aligned-reads BAM (pretending " .
            "also to be the sorted-aligned-reads BAM.");
    Genome::Sys->create_symlink($output_file, $sorted_bam);

    $cmd = Genome::Model::Tools::Sam::IndexBam->create(
        bam_file => $sorted_bam,
    );

    $cmd->execute();
    unless (-e "$sorted_bam.bai") {
        die "Failed to create BAM index of \"$sorted_bam\".";
    }
}

sub _link_in_converted_fastqs {
    my ($self) = @_;

    my ($fastq1, $fastq2) = @{$self->_fastq_files};
    my $converted_fastq_prefix = join("/", $self->temp_staging_directory,
            "tmp", "reads_");
    $self->status_message("Creating symlink to fastq files (pretending " .
            "to be the converted-fastq files");
    Genome::Sys->create_symlink($fastq1, "$converted_fastq_prefix" . "1");
    Genome::Sys->create_symlink($fastq2, "$converted_fastq_prefix" . "2");
}

sub _create_unaligned_fastqs {
    my ($self) = @_;

    $self->status_message("Generating the unaligned fastqs from the " .
            "alignment results BAM");
    my $stage            = $self->temp_staging_directory;
    my $alignment_result = $self->alignment_result;
    my $bam_file         = $alignment_result->bam_file;
    my $output_directory = join("/", $stage, "unaligned_fastqs");

    my $cmd = Genome::Model::Tools::Sam::BamToUnalignedFastq->create(
        bam_file         => $bam_file,
        output_directory => $output_directory,
    );
    $cmd->execute();

    for my $i (1..2) {
        my @pair_fqs = $self->_get_pair_fqs($output_directory, $i);
        unless (scalar(@pair_fqs)) {
            die "Couldn't find any pair $i fastqs in $output_directory!";
        }
        my $destination = join("/", $stage, "tmp", "unaligned_$i.fq");
        $self->status_message(sprintf("Found %d lanes of data for pair %d",
                scalar(@pair_fqs), $i));
        if (scalar(@pair_fqs) == 1) {
            # just move the file
            my $source = join("/", $output_directory, $pair_fqs[0]);
            Genome::Sys->shellcmd("mv $source $destination",
                    input_files  => [$source],
                    output_files => [$destination],
            );
        } else {
            # cat the files (different lanes) together
            my $cmd = "cat";
            my @sources;
            for my $pair_fq (@pair_fqs) {
                my $source = join("/", $output_directory, $pair_fq);
                push(@sources, $source);
                $cmd = $cmd . " $source";
            }
            $cmd = $cmd . " > $destination";
            Genome::Sys->shellcmd($cmd,
                    input_files  => \@sources,
                    output_files => [$destination],
            );
        }
        unless (-s $destination) {
            die "Failed to create unaligned_fastq at \"$destination\"";
        }
        Genome::Sys->remove_directory_tree($output_directory);
    }
}

sub _get_pair_fqs {
    my ($self, $dir, $num) = @_;

    my @filenames = split(/\n/, `ls $dir`);
    my @results;
    for my $filename (@filenames) {
        if ($filename =~ m/s_._([12])_sequence\.txt/) {
            my $pair_num = $1;
            push(@results, $filename) if $pair_num eq $num;
        } else {
            die "Couldn't determine the pair number based " .
                    "on filename: $filename";
        }
    }
    return @results;
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
