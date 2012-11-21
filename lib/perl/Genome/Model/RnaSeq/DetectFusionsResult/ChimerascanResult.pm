package Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult;

use strict;
use warnings;

use above 'Genome';

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

    #We want to shell out and create the chimerscan index in a different UR context.
    #That way it is committed, even if we fail and other builds don't need to wait on
    #the overall chimerscan run just to get their index results.
    my $cmd = 'genome model rna-seq detect-fusions chimerascan-index';
    $cmd .= ' --version=' . $self->version;
    $cmd .= ' --bowtie-version=' . $bowtie_version;
    $cmd .= ' --reference-build=' . $self->alignment_result->reference_build->id;
    $cmd .= ' --annotation-build=' . $self->annotation_build->id;

    Genome::Sys->shellcmd(cmd => $cmd);

    my $index_class = 'Genome::Model::RnaSeq::DetectFusionsResult' .
                      '::ChimerascanResult::Index';
    my $index = $index_class->get_with_lock(
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
        version => $self->version,
        bowtie_version => $bowtie_version,
        reference_build => $self->alignment_result->reference_build,
        annotation_build => $self->annotation_build,
    );

    if ($index) {
        $self->status_message( 'Registering ' . $self->build->id . ' as a user of the generated index' );
        $index->add_user( user => $self->build, label => 'uses' );
    } else {
        die("Unable to get a chimerascan index result");
    }

    return $index->output_dir;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    return 'build_merged_alignments/chimerascan/' . $self->id;
}

1;
