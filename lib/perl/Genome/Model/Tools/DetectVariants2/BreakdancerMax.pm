package Genome::Model::Tools::DetectVariants2::BreakdancerMax;

use warnings;
use strict;

use Genome;
use File::Basename;

class Genome::Model::Tools::DetectVariants2::BreakdancerMax{
    is => 'Genome::Model::Tools::DetectVariants2::Detector',
    has => [
        _config_base_name => {
            is => 'Text',
            default_value => 'breakdancer_config',
            is_constant => 1,
        },
        config_file => {
            is_output   => 1,
            is_input    => 1,
            is_optional => 1,
            doc => 'breakdancer config file path if provided not made by this tool',
        },
        _config_staging_output => {
            calculate_from => ['_temp_staging_directory', '_config_base_name'],
            calculate => q{ join("/", $_temp_staging_directory, $_config_base_name); },
            is_optional => 1,
        },
        chromosome => {
            is => 'Text',
            is_input => 1,
            is_optional => 1,
            default_value => 'all',
        },
        version => {
            is => 'Version',
            is_optional => 1,
            is_input => 1,
            default_value =>  Genome::Model::Tools::Breakdancer->default_breakdancer_version,
            valid_values  => [Genome::Model::Tools::Breakdancer->available_breakdancer_versions],
            doc => "Version of breakdancer to use",
        },
        workflow_log_dir => {
            is => 'Text',
            calculate_from => 'output_directory',
            calculate => q{ return File::Basename::dirname($output_directory) . '/breakdancer_by_chromosome_log/'; },
            is_optional => 1,
            doc => 'workflow log directory of per chromosome breakdancer run',
        },
        _bam2cfg_params=> {
            calculate_from => ['params'],
            calculate => q{
                return (split(':', $params))[0];
            },
            is_optional => 1,
            doc => 'This is the property used internally by the tool for bam2cfg parameters. It splits params.',
        },
        _breakdancer_params => {
            calculate_from => ['params'],
            calculate => q{
                return (split(':', $params))[1];
            },
            is_optional => 1,
            doc => 'This is the property used internally by the tool for breakdancer parameters. It splits params.',
        },
        _base_output_directory => {
            is => 'Text',
            is_optional => 1,
            doc => 'Store the base output directory when using per-chromosome output dirs',
        },
    ],
    has_param => [
        lsf_resource => {
            default_value => "-M 25000000 -R 'select[mem>25000] rusage[mem=25000]'",
        },
    ],
};

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants2 breakdancer -aligned-reads-input tumor.bam -control-aligned-reads-input normal.bam --output-dir breakdancer_dir
gmt detect-variants2 breakdancer -aligned-reads-input tumor.bam -control-aligned-reads-input normal.bam --output-dir breakdancer_dir --version 0.0.1r59
EOS
}

sub help_detail {
    return <<EOS
This tool discovers structural variation.  It generates an appropriate configuration based on
the input BAM files and then uses that configuration to run breakdancer.
EOS
}


sub _create_temp_directories {
    my $self = shift;
    local %ENV = %ENV;
    $ENV{TMPDIR} = $self->output_directory;
    return $self->SUPER::_create_temp_directories(@_);
}

sub _promote_staged_data {
    my $self = shift;
    my $output_dir = $self->SUPER::_promote_staged_data(@_);

    #remove the staging directory before the reallocation occurs so that it doesn't get double counted.
    my $staging_dir = $self->_temp_staging_directory;
    Genome::Sys->remove_directory_tree($staging_dir);

    return $output_dir;
}

sub _resolve_output_directory {
    my $self = shift;

    if ($self->_base_output_directory and $self->_base_output_directory ne $self->output_directory) {
        return 1;
    }

    $self->_base_output_directory($self->output_directory);

    if($self->chromosome and $self->chromosome ne 'all') {
        my $output_dir = $self->output_directory;
        unless(-d $output_dir) {
            #This should only happen if a single chromosome was executed directly
            Genome::Sys->create_directory($output_dir);
        }
        #Put per-chromosome outputs in subdirectories to avoid collisions in SoftwareResults
        $self->output_directory($output_dir . '/' . $self->chromosome);
    }

    return 1;
}

sub _detect_variants {
    my $self = shift;

    $self->run_config;
    $self->run_breakdancer;

    return 1;
}

sub run_config {
    my $self = shift;
    my $cfg_file = $self->config_file;

    if ($cfg_file) {
        unless (Genome::Sys->check_for_path_existence($cfg_file)) {
            $self->error_message("Given breakdancer config file $cfg_file is not valid");
            die $self->error_message;
        }
        $self->status_message("Using given breakdancer config file: $cfg_file");
    }
    else {
        $self->status_message("Run bam2cfg to make breakdancer_config file");

        my %params = (
            tumor_bam   => $self->aligned_reads_input,
            output_file => $self->_config_staging_output,
            params      => $self->_bam2cfg_params,
            use_version => $self->version,
        );

        $params{normal_bam} = $self->control_aligned_reads_input
            if $self->control_aligned_reads_input;

        my $bam2cfg = Genome::Model::Tools::Breakdancer::BamToConfig->create(%params);

        unless ($bam2cfg->execute) {
            $self->error_message("Failed to run bam2cfg");
            die;
        }

        $self->config_file($self->_config_staging_output);
        $self->status_message('BreakdancerMax config is created ok');
    }
    return 1;
}


sub run_breakdancer {
    my $self = shift;
    my $bd_params = $self->_breakdancer_params || "";

    #Allow 0 size of config, breakdancer output
    if (-z $self->config_file) {
        $self->warning_message("0 size of breakdancer config file. Probably it is for testing of small bam files");
        my $output_file = $self->_sv_staging_output;
        `touch $output_file`;
        return 1;
    }

    my $is_ctx = 0;
    if ($bd_params =~ /\-o/) {
        my $chr = $self->chromosome;
        if ($chr eq 'all') {
            require Workflow::Simple;

            my $op = Workflow::Operation->create(
                name => 'BreakdancerMax by chromosome',
                operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::DetectVariants2::BreakdancerMax'),
            );

            $op->parallel_by('chromosome');

            if(Workflow::Model->parent_workflow_log_dir) {
                $op->log_dir(Workflow::Model->parent_workflow_log_dir);
            } elsif ($self->workflow_log_dir) {
                unless (-d $self->workflow_log_dir) {
                    unless (Genome::Sys->create_directory($self->workflow_log_dir)) {
                        $self->error_message('Failed to create workflow_log_dir: '. $self->workflow_log_dir);
                        die;
                    }
                }
                $op->log_dir($self->workflow_log_dir);
            }

            my $cfg_file = $self->config_file;

            unless (Genome::Sys->check_for_path_existence($cfg_file)) {
                $self->error_message('prerun breakdancer config file '.$cfg_file.' does not exist');
                die $self->error_message;
            }

            my @chr_list = $self->_get_chr_list;
            if (scalar @chr_list == 0) {
                #FIXME Sometimes samtools idxstats does not get correct
                #stats because of bam's bai file is not created by
                #later samtools version (0.1.9 ?)
                die $self->error_message("chr list from samtools idxstats is empty, using full chr list now");
            }

            $self->status_message('chromosome list is '.join ',', @chr_list);

            my %params = (
                aligned_reads_input         => $self->aligned_reads_input,
                reference_build_id          => $self->reference_build_id,
                output_directory            => $self->_temp_staging_directory,
                config_file => $cfg_file,
                params      => $self->params,
                version     => $self->version,
                chromosome  => \@chr_list,
            );
            $params{control_aligned_reads_input} = $self->control_aligned_reads_input
                if $self->control_aligned_reads_input;

            my $output = Workflow::Simple::run_workflow_lsf($op, %params);

            unless (defined $output) {
                my @error;
                for (@Workflow::Simple::ERROR) {
                    push @error, $_->error;
                }
                $self->error_message(join("\n", @error));
                die $self->error_message;
            }

            my $merge_obj = Genome::Model::Tools::Breakdancer::MergeFiles->create(
                input_files => join(',', map { $self->_temp_staging_directory . '/' . $_ . '/' . $self->_sv_base_name . '.' . $_ } @chr_list),
                output_file => $self->_sv_staging_output,
            );
            my $merge_rv = $merge_obj->execute;
            Carp::confess 'Could not execute breakdancer file merging!' unless defined $merge_rv and $merge_rv == 1;

            return 1;
        }
        else {
            $self->_sv_base_name($self->_sv_base_name . '.' . $chr);
            $bd_params =~ s/\-o/\-o $chr/;
        }
    }
    elsif ($bd_params =~ /\-d/) {
        $is_ctx = 1;
        my $sv_staging_out = $self->_sv_staging_output;
        $bd_params =~ s/\-d/\-d $sv_staging_out/;
    }

    my $breakdancer_path = Genome::Model::Tools::Breakdancer->breakdancer_max_command_for_version($self->version);
    my $cfg_file         = $self->config_file;

    my $cmd = "$breakdancer_path " . $cfg_file . " " . $bd_params . " > "  . $self->_sv_staging_output;

    $self->status_message("EXECUTING BREAKDANCER STEP: $cmd");
    my $return = Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files  => [$cfg_file],
        output_files => [$self->_sv_staging_output],
        allow_zero_size_output_files => 1,
    );

    unless ($return) {
        $self->error_message("Running breakdancer failed using command: $cmd");
        die;
    }

    unless (-s $self->_sv_staging_output) {
        $self->error_message("$cmd output " . $self->_sv_staging_output . " does not exist or has zero size");
        die;
    }

    if ( $is_ctx ) { # validate the output fastqs and write an md5 for them
        $self->status_message('Validate CTX fastqs...');
        my $validate_fastqs = $self->_validate_ctx_fastqs;
        die if not $validate_fastqs;
        $self->status_message('Validate CTX fastqs...OK');
    }

    $self->status_message('breakdancer run finished ok');
    return 1;
}

sub _get_chr_list {
    my $self = shift;

    my $chromosome = $self->chromosome;
    if (not $chromosome) {
        die "Chromosoome not set???!  Expected 'all', or a specific sequence name!";
    }

    my $tmp_idx_dir = File::Temp::tempdir(
        "Normal_bam_idxstats_XXXXX",
        CLEANUP => 1,
        DIR     => '/tmp',  #File::Temp can not remove inside dir for _temp_staging_dir
        #DIR     => $self->_temp_staging_directory,
    );

    my $tmp_idx_file = $tmp_idx_dir . '/normal_bam.idxstats';

    my $idxstats = Genome::Model::Tools::Sam::Idxstats->create(
        bam_file    => $self->aligned_reads_input,
        output_file => $tmp_idx_file,
    );
    unless ($idxstats->execute) {
        $self->error_message("Failed to run samtools idxstats output $tmp_idx_file");
        die;
    }

    my $map_chr_list = $idxstats->map_ref_list($tmp_idx_file);

    if ($chromosome eq 'all') {
        # this is the _real_ <all>, since <all> is a legacy flag meaning only 1-22,X,Y,MT
        # return (grep { $map_chr_list{$_} } @FULL_CHR_LIST);
        return @$map_chr_list;
    }
    else {
        my %map_chr_list = map { $_ => 1 } @$map_chr_list;
        if ($map_chr_list{$chromosome}) {
            return ($chromosome);
        }
        else {
            die "Chromosome $chromosome not found in BAM idxstats!";
        }
    }
}


sub has_version {
    my $self    = shift;
    my $version = shift;

    unless (defined $version) {
        $version = $self->version;
    }
    my @versions = Genome::Model::Tools::Breakdancer->available_breakdancer_versions;
    for my $v (@versions) {
        if ($v eq $version) {
            return 1;
        }
    }
    return 0;
}

sub params_for_detector_result {
    my $self = shift;
    my ($params) = $self->SUPER::params_for_detector_result;

    $params->{chromosome_list} = $self->chromosome;
    return $params;
}

sub _validate_ctx_fastqs {
    my $self = shift;

    # Get the fastqs [only for ctx]
    my @fastqs = grep { -s } glob($self->_sv_staging_output . "*.fastq");
    if ( not @fastqs ) {
        $self->error_message('Expected fastqs in sv staging output! '.$self->_sv_staging_output);
        return;
    }

    my $md5sum;
    for my $fastq ( @fastqs ) {
        # Validate fastq
        my $reader = Genome::Model::Tools::Sx::FastqReader->create(file => $fastq);
        if ( not $reader ) {
            $self->error_message('Failed to open fastq file! '.$fastq);
            return;
        }
        return if not $reader->validate;
        # Add md5sum
        my $md5sum_for_fastq = eval{ Genome::Sys->md5sum($fastq); };
        if ( not $md5sum_for_fastq ) {
            $self->error_message($@) if $@;
            $self->error_message("Failed to get md5sum for fastq $fastq");
            return;
        }
        my $basename = File::Basename::basename($fastq);
        $md5sum .= "$md5sum_for_fastq\t$basename\n";
    }

    my $md5file = $self->_sv_staging_output . '.fastqs.md5';
    my $fh = Genome::Sys->open_file_for_writing($md5file);
    $fh->print($md5sum);
    $fh->close;

    return 1;
}

1;
