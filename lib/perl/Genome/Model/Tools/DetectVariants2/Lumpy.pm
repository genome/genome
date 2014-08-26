#lumpy.pm
#this program is ment to run the bam files specified.  it will run both a paired end and a single read file before depositing the data into a histo file in the lumpy-sv_results dir

package Genome::Model::Tools::DetectVariants2::Lumpy;
use warnings;
use strict;

use Genome;
use File::Basename;
use IPC::System::Simple;

my @FULL_CHR_LIST = (1 .. 22, 'X', 'Y', 'MT');

class Genome::Model::Tools::DetectVariants2::Lumpy {is => 'Genome::Model::Tools::DetectVariants2::Detector',};

sub _detect_variants {
    my $self     = shift;

    my @pe_cmds;
    my @sr_cmds;
    for my $input_bam ($self->aligned_reads_input, $self->control_aligned_reads_input) {
        next unless defined($input_bam);
        for my $bam (split_bam_by_readgroup($input_bam)) {
            if ($self->pe_param) {
                push(@pe_cmds, $self->paired_end_parameters_for_bam($bam));
            }
            if ($self->sr_param) {
                push(@sr_cmds, $self->split_read_parameters_for_bam($bam));
            }
        }
    }
    my $pe_cmd = join(",@pe_cmds");
    my $sr_cmd = join(",@sr_cmds");

    my @cmd = $self->_open_params();
    splice @cmd, 1, 0, "@pe_cmds", "@sr_cmds";
    my $cmmd = "@cmd";

    my $run = Genome::Sys->shellcmd(
        cmd                          => $cmmd,
        output_files                 => [$self->_sv_staging_output],
        allow_zero_size_output_files => 1,
    );
}

sub split_bam_by_readgroup {
    my $bam_file  = shift;

    my $split_dir = Genome::Sys->create_temp_directory();
    my $split_bam_basename = File::Spec->join($split_dir, 'split_bam');
    my $command = "bamtools split -in $bam_file -stub $split_bam_basename -tag RG";
    Genome::Sys->shellcmd(cmd => $command);
    return glob("$split_bam_basename*");
}

sub extract_paired_end_reads {
    my $bam = shift;

    my $filtered_bam  = Genome::Sys->create_temp_file_path();
    my $command = "samtools view -b -F 1294 $bam -o $filtered_bam";
    Genome::Sys->shellcmd(
        cmd                          => $command,
        allow_zero_size_output_files => 1,
    );
    return $filtered_bam;
}

sub extract_split_reads {
    my $self = shift;
    my $bam = shift;
    my $filtered_bam = Genome::Sys->create_temp_file_path();
    my $extract_split_reads_bwamen_script = $self->lumpy_script_for_extract_split_reads_bwamem();
    my $command = join(
        '|',
        "samtools view -h $bam",
        "$extract_split_reads_bwamen_script -i stdin",
        "java -Xmx8g -XX:MaxPermSize=256m -cp /gsc/scripts/lib/java/samtools/picard-tools-1.82/SamFormatConverter.jar net.sf.picard.sam.SamFormatConverter I=/dev/stdin O=$filtered_bam"
    );
    Genome::Sys->shellcmd(
        cmd                          => $command,
        allow_zero_size_output_files => 1,
    );
    return $filtered_bam;
}

sub paired_end_parameters_for_bam {
    my $self = shift;
    my $bam  = shift;

    my $filtered_bam = extract_paired_end_reads($bam);
    my %metrics = $self->calculate_metrics($bam);

    return sprintf(
        ' -pe bam_file:%s,histo_file:%s,mean:%s,stdev:%s,read_length:150,%s',
        $filtered_bam,
        $metrics{histogram},
        $metrics{mean},
        $metrics{standard_deviation},
        $self->pe_param
    );
}

sub split_read_parameters_for_bam {
    my $self = shift;
    my $bam = shift;

    my $filtered_bam = $self->extract_split_reads($bam);
    return sprintf(
        " -sr bam_file:%s,%s",
        $filtered_bam,
        $self->sr_param
    );
}

sub calculate_metrics {
    my $self = shift;
    my $bam  = shift;

    my $histogram = Genome::Sys->create_temp_file_path();
    my $pairend_distro_script = $self->lumpy_script_for_pairend_distro;
    my @commands   = qq(samtools view $bam | tail -n+100 | $pairend_distro_script -r1 100 -X 4 -N 10000 -o $histogram);
    my $output = IPC::System::Simple::capture(@commands);

    if ($output =~ m/mean:([-+]?[0-9]*\.?[0-9]*)\s+stdev:([-+]?[0-9]*\.?[0-9]*)/) {
        my $mean = $1;
        my $standard_deviation = $2;
        my %metrics = (
            mean               => $mean,
            standard_deviation => $standard_deviation,
            histogram          => $histogram,
        );
        return %metrics;
    }
    else {
        die "ERROR couldn't determine mean and standard deviation: $output";
    }
}

sub sr_param {
    my $self   = shift;
    my %params = $self->params_hash();
    return $params{'sr'};
}

sub pe_param {
    my $self   = shift;
    my %params = $self->params_hash();
    return $params{'pe'};
}

sub lumpy_param {
    my $self   = shift;
    my %params = $self->params_hash();
    return $params{'lp'};
}

sub params_hash {
    my $self            = shift;
    my $unparsed_params = $self->params;
    my @params          = split('//', $unparsed_params);
    my %parameters;
    foreach my $place (@params) {
        if ($place =~ m/^\-([a-z]{2}),(.*)$/) {
            $parameters{$1} = $2;
        }
        else {
            die sprintf(
                "You specified the parameters incorrectly. Unparsed parametere were: (%s)  The malformed parameters were: (%s)",
                $unparsed_params, $place);
        }
    }
    return %parameters;
}

sub _open_params {
    my $self      = shift;
    my $lump_text = $self->lumpy_param;
    $lump_text =~ s/,/ /g;
    $lump_text =~ s/:/ /g;
    # print "$lump_text";
    my $executable_path = $self->lumpy_command;
    my $output_files    = $self->_sv_staging_output;
    my @sur_cmd         = ("$executable_path $lump_text ", " > $output_files");
    return @sur_cmd;
}

sub lumpy_directory {
    my $self    = shift;
    my $version = $self->version();

    return _lumpy_directory($version);
}

sub _lumpy_directory {
    my $version = shift;
    return File::Spec->catdir(File::Spec->rootdir, "usr", "lib", "lumpy" . "$version");
}

sub lumpy_command {
    my $self = shift;
    return File::Spec->catfile($self->lumpy_directory(), "bin", "lumpy");
}

sub lumpy_scripts_directory {
    my $self = shift;
    return File::Spec->catfile($self->lumpy_directory(), 'scripts');
}

sub lumpy_script_for {
    my $self        = shift;
    my $script_name = shift;

    die "no script name given" if not $script_name;
    my $script_location = File::Spec->catfile($self->lumpy_scripts_directory(), "$script_name");

    die "script does not exist $script_location" if not -e $script_location;
    return $script_location;
}

sub lumpy_script_for_extract_split_reads_bwamem {
    my $self = shift;
    return $self->lumpy_script_for("extractSplitReads_BwaMem");
}

sub lumpy_script_for_pairend_distro {
    my $self = shift;
    return $self->lumpy_script_for("pairend_distro.py");
}

sub has_version {
    my $self    = shift;
    my $version = shift;
    if (-d _lumpy_directory($version)) {
        return 1;
    }
    else {
        return 0;
    }
}

