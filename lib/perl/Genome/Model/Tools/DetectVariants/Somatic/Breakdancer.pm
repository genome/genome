package Genome::Model::Tools::DetectVariants::Somatic::Breakdancer;

use warnings;
use strict;

use Genome;
use File::Copy;

my $DEFAULT_VERSION = '2010_06_24';
my @FULL_CHR_LIST   = (1..22, 'X', 'Y');

class Genome::Model::Tools::DetectVariants::Somatic::Breakdancer{
    is => 'Genome::Model::Tools::DetectVariants::Somatic',
    has => [
        _config_base_name => {
            is => 'Text',
            default_value => 'breakdancer_config',
            is_input => 1,
        },
        config_file => {
            #calculate_from => ['_config_base_name', 'output_directory'],
            #calculate => q{ join("/", $output_directory, $_config_base_name); },
            is_output   => 1,
            is_input    => 1,
            is_optional => 1,
            doc => 'breakdancer config file path if provided not made by this tool',
        },
        _config_staging_output => {
            calculate_from => ['_temp_staging_directory', '_config_base_name'],
            calculate => q{ join("/", $_temp_staging_directory, $_config_base_name); },
        },
        chromosome => {
            is => 'Text',
            is_input => 1,
            is_optional => 1,
            valid_values => [@FULL_CHR_LIST, 'all'],
            default_value => 'all',
        },
        version => {
            is => 'Version',
            is_optional => 1,
            is_input => 1,
            default_value => $DEFAULT_VERSION,
            doc => "Version of breakdancer to use",
        },
        workflow_log_dir => {
            is => 'Text',
            is_optional => 1,
            doc => 'workflow log directory of per chromosome breakdancer run',
        },
        detect_svs => { value => 1, is_constant => 1, },
        sv_params => {
            is => 'Text',
            is_input => 1,
            is_optional => 1,
            doc => "Parameters to pass to bam2cfg and breakdancer. The two should be separated by a ':'. i.e. 'bam2cfg params:breakdancer params'",
        },
        _bam2cfg_params=> {
            calculate_from => ['sv_params'],
            calculate => q{
                return (split(':', $sv_params))[0];
            },
            doc => 'This is the property used internally by the tool for bam2cfg parameters. It splits sv_params.',
        },
        _breakdancer_params => {
            calculate_from => ['sv_params'],
            calculate => q{
                return (split(':', $sv_params))[1];
            },
            doc => 'This is the property used internally by the tool for breakdancer parameters. It splits sv_params.',
        },
        skip => {
            is => 'Boolean',
            default => '0',
            is_input => 1,
            is_optional => 1,
            doc => "If set to true... this will do nothing! Fairly useless, except this is necessary for workflow.",
        },
        skip_if_output_present => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enable this flag to shortcut this step if the output is already present. Useful for pipelines.',
        },
        ],
    has_param => [ 
        lsf_resource => {
            default_value => "-M 8000000 -R 'select[type==LINUX64 && mem>8000] rusage[mem=8000]'",
        },
        lsf_queue => {
            default_value => 'long'
        }, 
    ],
    # These are params from the superclass' standard API that we do not require for this class (dont show in the help)
    has_constant_optional => [
        snv_params=>{},
        indel_params=>{},
        capture_set_input =>{},
        detect_snvs=>{},
        detect_indels=>{},
    ],
};

my %BREAKDANCER_VERSIONS = (
    '0.0.1r59'   => {
        dir => '/gsc/scripts/pkg/bio/breakdancer/breakdancer-0.0.1r59',
        cfg => 'bam2cfg.pl',
        max => 'BreakDancerMax.pl',
    },
    '2010_02_17' => {
        dir => '/gsc/scripts/pkg/bio/breakdancer/breakdancer-2010_02_17/bin',
        cfg => 'bam2cfg.pl',
        max => 'BreakDancerMax.pl',
    },
    '2010_03_02' => {
        dir => '/gsc/scripts/pkg/bio/breakdancer/breakdancer-2010_03_02/bin',
        cfg => 'bam2cfg.pl',
        max => 'BreakDancerMax.pl',
    },
    '2010_06_24' => {
        dir => $ENV{GENOME_SW} . '/breakdancermax/breakdancer-20100624',
        cfg => 'perl/bam2cfg.pl',
        max => 'cpp/breakdancer_max',
    },
);


sub help_brief {
    "discovers structural variation using breakdancer",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt somatic breakdancer -t tumor.bam -n normal.bam --output-dir breakdancer_dir
gmt somatic breakdancer -t tumor.bam -n normal.bam --output-dir breakdancer_dir --version 0.0.1r59 --skip-if-output-present 
EOS
}

sub help_detail {                           
    return <<EOS 
This tool discovers structural variation.  It generates an appropriate configuration based on
the input BAM files and then uses that configuration to run breakdancer.
EOS
}

sub _should_skip_execution {
    my $self = shift;
    
    if ($self->skip) {
        $self->status_message("Skipping execution: Skip flag set");
        return 1;
    }
    if (($self->skip_if_output_present)&&(-s $self->sv_output)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }
    
    return $self->SUPER::_should_skip_execution;
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
        my $config_path = $self->breakdancer_config_command;
        my $cmd = "$config_path " . $self->_bam2cfg_params .' '.$self->aligned_reads_input . ' ' . $self->control_aligned_reads_input . " > "  . $self->_config_staging_output;

        $self->status_message("EXECUTING CONFIG STEP: $cmd");
        my $return = Genome::Sys->shellcmd(
            cmd => $cmd,
            input_files  => [$self->aligned_reads_input, $self->control_aligned_reads_input],
            output_files => [$self->_config_staging_output],
        );

        unless ($return) {
            $self->error_message("Running breakdancer config failed using command: $cmd");
            die;
        }

        unless (-s $self->_config_staging_output) {
            $self->error_message("$cmd output " . $self->_config_staging_output . " does not exist or has zero size");
            die;
        }
        $self->config_file($self->_config_staging_output);
        $self->status_message('Breakdancer config is created ok');
    }
    return 1;
}


sub run_breakdancer {
    my $self = shift;
    my $bd_params = $self->_breakdancer_params;

    if ($bd_params =~ /\-o/) {
        my $chr = $self->chromosome;
        if ($chr eq 'all') {
            require Workflow::Simple;
        
            my $op = Workflow::Operation->create(
                name => 'Breakdancer by chromosome',
                operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::DetectVariants::Somatic::Breakdancer'),
            );

            $op->parallel_by('chromosome');
            $op->log_dir($self->workflow_log_dir) if $self->workflow_log_dir;

            my $temp_staging_dir = File::Temp::tempdir(
                "breakdancer_by_chromosome_XXXXX",
                DIR     => $self->output_directory,
                CLEANUP => 1,
            ); 
            $self->_temp_staging_directory($temp_staging_dir);
            my $cfg_file = $temp_staging_dir . '/breakdancer_config';

            unless (Genome::Sys->check_for_path_existence($self->config_file)) {
                $self->error_message('prerun breakdancer config file '.$self->config_file.' does not exist');
                die $self->error_message;
            }

            copy $self->config_file, $cfg_file; #Make each chr run sharing the same config file

            unless (Genome::Sys->check_for_path_existence($cfg_file)) {
                $self->error_message("breakdancer config file $cfg_file is not copied ok");
                die $self->error_message;
            }

            my @chr_list = $self->_get_chr_list;
            if (scalar @chr_list == 0) {
                #FIXME Sometimes samtools idxstats does not get correct
                #stats because of bam's bai file is not created by
                #later samtools version (0.1.9 ?)
                $self->warning_message("chr list from samtools idxstats is empty, using full chr list now"); 
                @chr_list = @FULL_CHR_LIST;
            }

            $self->status_message('chromosome list is '.join ',', @chr_list);

            my $output = Workflow::Simple::run_workflow_lsf(
                $op,
                aligned_reads_input         => $self->aligned_reads_input, 
                control_aligned_reads_input => $self->control_aligned_reads_input,
                reference_sequence_input    => $self->reference_sequence_input,
                output_directory            => $self->_temp_staging_directory,
                config_file => $cfg_file,
                sv_params   => $self->sv_params,
                version     => $self->version,
                chromosome  => \@chr_list,
            );
            unless (defined $output) {
                my @error;
                for (@Workflow::Simple::ERROR) {
                    push @error, $_->error;
                }
                $self->error_message(join("\n", @error));
                die $self->error_message;
            }

            my $merge_obj = Genome::Model::Tools::Breakdancer::MergeFiles->create(
                input_files => join(',', map { $self->_temp_staging_directory . '/' . $self->_sv_base_name . '.' . $_ } @chr_list),
                output_file => $self->_temp_staging_directory . '/' . $self->_sv_base_name,
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
        my $sv_staging_out = $self->_sv_staging_output;
        $bd_params =~ s/\-d/\-d $sv_staging_out/;
    }

    my $breakdancer_path = $self->breakdancer_max_command;
    my $cfg_file         = $self->config_file;

    my $cmd = "$breakdancer_path " . $cfg_file . " " . $bd_params . " > "  . $self->_sv_staging_output;

    $self->status_message("EXECUTING BREAKDANCER STEP: $cmd");
    my $return = Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files  => [$self->config_file],
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

    $self->status_message('breakdancer run finished ok');
    return 1;
}

sub _get_chr_list {
    my $self = shift;

    my $tmp_idx_dir = File::Temp::tempdir(
        "Normal_bam_idxstats_XXXXX",
        CLEANUP => 1,
        DIR     => $self->_temp_staging_directory,
    );

    my $tmp_idx_file = $tmp_idx_dir . '/normal_bam.idxstats';

    my $idxstats = Genome::Model::Tools::Sam::Idxstats->create(
        bam_file    => $self->control_aligned_reads_input,
        output_file => $tmp_idx_file,
    );
    unless ($idxstats->execute) {
        $self->error_message("Failed to run samtools idxstats output $tmp_idx_file");
        die;
    }

    my $unmap_chr_list = $idxstats->unmap_ref_list($tmp_idx_file);
    my @chr_list; 

    for my $chr (@FULL_CHR_LIST) {
        push @chr_list, $chr unless grep{$chr eq $_}@$unmap_chr_list;
    }

    return @chr_list;
}


sub breakdancer_path {
    my $self = shift;
    return $self->path_for_breakdancer_version($self->version);
}

sub breakdancer_max_command { 
    my $self = shift;
    return $self->breakdancer_max_command_for_version($self->version);
}

sub breakdancer_config_command { 
    my $self = shift;
    return $self->breakdancer_config_command_for_version($self->version);
}

sub available_breakdancer_versions {
    my $self = shift;
    return keys %BREAKDANCER_VERSIONS;
}

sub path_for_breakdancer_version {
    my ($self, $version) = @_;

    if (defined $BREAKDANCER_VERSIONS{$version}) {
        my $dir = $BREAKDANCER_VERSIONS{$version}->{dir};
        unless (-d $dir) {
            $self->error_message("breakdancer base dir $dir for version $version is not valid");
            die $self->error_message;
        }
        return $dir;
    }
    die 'No path for breakdancer version '. $version;
}

sub breakdancer_max_command_for_version {
    my ($self, $version) = @_;

    if (defined $BREAKDANCER_VERSIONS{$version}->{max}) {
        my $max_cmd = $self->path_for_breakdancer_version($version) . "/" .  $BREAKDANCER_VERSIONS{$version}->{max};
        unless (-s $max_cmd and -x $max_cmd) {
            $self->error_message("breakdancer_max command $max_cmd for version $version is not valid");
            die $self->error_messge;
        }
        return $max_cmd;
    }
    die 'No breakdancer max command for breakdancer version '. $version;
}

sub breakdancer_config_command_for_version {
    my ($self, $version) = @_;

    if (defined $BREAKDANCER_VERSIONS{$version}->{cfg}) {
        my $cfg_cmd = $self->path_for_breakdancer_version($version) . "/" .  $BREAKDANCER_VERSIONS{$version}->{cfg};
        unless (-s $cfg_cmd and -x $cfg_cmd) {
            $self->error_message("breakdancer config command $cfg_cmd for version $version is not valid");
            die $self->error_messge;
        }
        return $cfg_cmd;
    }
    die 'No breakdancer config command for breakdancer version '. $version;
}

sub default_breakdancer_version {
    die "default breakdancer version: $DEFAULT_VERSION is not valid" unless $BREAKDANCER_VERSIONS{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}
 
1;
