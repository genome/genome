package Genome::InstrumentData::AlignmentResult::Bwamem;

use strict;
use warnings;
use Carp qw/confess/;
use Data::Dumper;
use File::Basename;
use Genome;
use Getopt::Long;

class Genome::InstrumentData::AlignmentResult::Bwamem {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => { value => 'bwamem', is_param=>1 },
    ],
    has_transient_optional => [
        _bwa_sam_cmd => { is => 'Text' }
    ]
};

sub required_arch_os { 'x86_64' }

sub required_rusage {
    my $class = shift;
    my %p = @_;
    my $instrument_data = delete $p{instrument_data};
    my $aligner_params  = delete $p{aligner_params};

    my $tmp_mb = $class->tmp_megabytes_estimated($instrument_data);
    my $mem_mb = 1024 * 16; 
    my $cpus = 4;

    if ($aligner_params and $aligner_params =~ /-t\s*([0-9]+)/) {
        $cpus = $1;
    }

    my $mem_kb = $mem_mb*1024;
    my $tmp_gb = $tmp_mb/1024;

    my $user = getpwuid($<);
    my $queue = 'alignment';
    $queue = 'alignment-pd' if (Genome::Config->should_use_alignment_pd);

    my $host_groups;
    $host_groups = qx(bqueues -l $queue | grep ^HOSTS:);
    $host_groups =~ s/\/\s+/\ /;
    $host_groups =~ s/^HOSTS:\s+//;

    my $select  = "select[ncpus >= $cpus && mem >= $mem_mb && gtmp >= $tmp_gb] span[hosts=1]";
    my $rusage  = "rusage[mem=$mem_mb, gtmp=$tmp_gb]";
    my $options = "-M $mem_kb -n $cpus -q $queue";

    my $required_usage = "-R '$select $rusage' $options";

    #check to see if our resource requests are feasible (This uses "maxmem" to check theoretical availability)
    #factor of four is based on current six jobs per host policy this should be revisited later
    my $select_check = "select[ncpus >= $cpus && maxmem >= " . ($mem_mb * 4) . " && maxgtmp >= $tmp_gb] span[hosts=1]";
    my $select_cmd = "bhosts -R '$select_check' $host_groups | grep ^blade";

    my @selected_blades = qx($select_cmd);

    if (@selected_blades) {
        return $required_usage;
    } else {
        die $class->error_message("Failed to find hosts that meet resource requirements ($required_usage). [Looked with `$select_cmd`]");
    }
}

# TODO copied verbatim from normal bwa, but this may be totally off for bwa mem
sub tmp_megabytes_estimated {
    my $class = shift || die;
    my $instrument_data = shift;

    my $default_megabytes = 90000;


    if (not defined $instrument_data) {
        return $default_megabytes;
    } elsif ($instrument_data->bam_path) {
        my $bam_path = $instrument_data->bam_path;

        my $scale_factor = 3.25; # assumption: up to 3x during sort/fixmate/sort and also during fastq extraction (2x) + bam = 3

        my $bam_bytes = -s $bam_path;
        unless ($bam_bytes) {
            die $class->error_message("Instrument Data " . $instrument_data->id  . " has BAM ($bam_path) but has no size!");
        }

        if ($instrument_data->can('get_segments')) {
            my $bam_segments = scalar $instrument_data->get_segments;
            if ($bam_segments > 1) {
                $scale_factor = $scale_factor / $bam_segments;
            }
        }

        return int(($bam_bytes * $scale_factor) / 1024**2);
    } elsif ($instrument_data->can("calculate_alignment_estimated_kb_usage")) {
        my $kb_usage = $instrument_data->calculate_alignment_estimated_kb_usage;
        return int(($kb_usage * 3) / 1024) + 100; # assumption: 2x for the quality conversion, one gets rm'di after; aligner generates 0.5 (1.5 total now); rm orig; sort and merge maybe 2-2.5
    } else {
        return $default_megabytes;
    }

    return;
}

# override this from AlignmentResult.pm to filter reads with secondary alignment flag (0x100)
sub _check_read_count {
    my ($self) = @_;
    my $fq_rd_ct = $self->_fastq_read_count;
    my $sam_path = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);

    my $cmd = "$sam_path view -F 256 -c " . $self->temp_staging_directory . "/all_sequences.bam";
    my $bam_read_count = `$cmd`;
    my $check = "Read count from bam: $bam_read_count and fastq: $fq_rd_ct";

    unless ($fq_rd_ct == $bam_read_count) {
        $self->error_message("$check does not match.");
        return;
    }
    $self->status_message("$check matches.");
    return 1;
}

sub _run_aligner {
    my $self = shift;
    my @input_pathnames = @_;

    # process inputs
    if (@input_pathnames != 1 and @input_pathnames != 2) {
        $self->error_message(
            "Expected 1 or 2 input path names. Got: " . Dumper(\@input_pathnames));
    }
    my $input_filenames = join(' ', @input_pathnames);

    # get temp dir
    my $tmp_dir = $self->temp_scratch_directory;
    my $log_filename = $tmp_dir . '/aligner.log';
    my $raw_sam = $tmp_dir . '/raw_sequences.sam';
    #my $fixed_sam = $tmp_dir . '/fixed_sequences.sam'; # XXX
    my $final_sam = $tmp_dir . '/all_sequences.sam';

    # get refseq info
    my $reference_build = $self->reference_build;
    my $reference_fasta_path = $self->get_reference_sequence_index->full_consensus_path('fa');

    # get params and verify the aligner
    my $params = $self->decomposed_aligner_params;

    my $aligner_version = $self->aligner_version;

    unless (Genome::Model::Tools::Bwa->supports_mem($aligner_version)) {
        die $self->error_message("The pipeline does not support using bwa mem with bwa-$aligner_version.");
    }

    my $command_name = Genome::Model::Tools::Bwa->path_for_bwa_version($aligner_version);

    # run cmd
    my $full_command = sprintf '%s mem %s %s %s 1>> %s 2>> %s',
        $command_name, $params, $reference_fasta_path, $input_filenames, $raw_sam, $log_filename;

    Genome::Sys->shellcmd(
        cmd          => $full_command,
        input_files  => [ @input_pathnames ],
        output_files => [ $raw_sam, $log_filename ],
        skip_if_output_is_present => 0,
    );

    #my $is_paired = @input_pathnames == 2 ? 1 : 0; # XXX
    #my $include_secondary = 1; # XXX
    #my $mark_secondary_as_duplicate = 0; # XXX

    #$self->status_message("Fixing flags and mates in merged sam file."); # XXX

    #$self->_fix_sam($raw_sam, $fixed_sam, $is_paired, $include_secondary, $mark_secondary_as_duplicate); # XXX

    #unlink($raw_sam) || die $self->error_message("Could not unlink $raw_sam."); # XXX

    $self->status_message("Resorting fixed sam file by coordinate.");

    my $picard_sort_cmd = Genome::Model::Tools::Picard::SortSam->create(
        sort_order             => 'coordinate',
        #input_file             => $fixed_sam, # XXX
        input_file             => $raw_sam, # XXX
        output_file            => $final_sam,
        max_records_in_ram     => 2000000,
        maximum_memory         => 8,
        maximum_permgen_memory => 256,
        temp_directory         => $self->temp_scratch_directory,
        use_version            => $self->picard_version,
    );

    # TODO not sure if the following is necessary
    #my $add_rg_cmd = Genome::Model::Tools::Sam::AddReadGroupTag->create(
    #    input_filehandle  => $sorted_sam,
    #    output_filehandle => $final_sam,
    #    read_group_tag    => $self->read_and_platform_group_tag_id,
    #    pass_sam_headers  => 0,
    #);

    unless ($picard_sort_cmd and $picard_sort_cmd->execute) {
        die $self->error_message("Failed to create or execute picard sort command.");
    }

    unlink($raw_sam) || die $self->error_message("Could not unlink $raw_sam."); # XXX
    #unlink($fixed_sam) || die $self->error_message("Could not unlink $fixed_sam."); # XXX

    # verify the bwa mem logfile
    unless ($self->_verify_bwa_bwa_mem_did_happen($log_filename)) {
        die $self->error_message("bwa mem seems to fail based on run log: $log_filename");
    }

    return 1;
}

sub _verify_bwa_bwa_mem_did_happen {
    my ($self, $log_file) = @_;
    # TODO implement this to make sure bwa mem finished; see _verify_bwa_samxe_did_happen in Bwa.pm
    return 1;
}

sub decomposed_aligner_params {
    my $self = shift;
    my $param_string = $self->aligner_params || '';

    my $param_hash = $self->get_aligner_params_hash($param_string);

    my $cpu_count = $self->_available_cpu_count;
    my $processed_param_string = $self->join_aligner_params_hash($param_hash);

    $self->status_message("[decomposed_aligner_params] cpu count is $cpu_count");
    $self->status_message("[decomposed_aligner_params] bwa mem params are: $processed_param_string");

    # Make sure the thread count argument matches the number of CPUs available.
    if ($param_hash->{t} ne $cpu_count) {
        $param_hash->{t} = $cpu_count;
        my $modified_param_string = $self->join_aligner_params_hash($param_hash);
        $self->status_message("[decomposed_aligner_params] autocalculated CPU requirement, bwa mem params modified: $modified_param_string");
    }

    if (not exists $param_hash->{M}) {
        $param_hash->{M} = '';
        my $modified_param_string = $self->join_aligner_params_hash($param_hash);
        $self->status_message("[decomposed_aligner_params] forcing -M, bwa mem params modified: $modified_param_string");
    }

    my $final_param_string = $self->join_aligner_params_hash($param_hash);

    return $final_param_string;
}

sub aligner_params_for_sam_header {
    my $self = shift;

    my $param_string = $self->aligner_params || '';
    my $param_hash = $self->get_aligner_params_hash($param_string);

    delete $param_hash->{t}; # we don't want cpu count to be part of the sam header

    my $modified_param_string = $self->join_aligner_params_hash($param_hash);

    return "bwa mem $modified_param_string";
}

# helper for decomposed_aligner_params and aligner_params_for_sam_header
sub get_aligner_params_hash {
    my $self = shift;
    my $param_string = shift;

    Getopt::Long::Configure("bundling");

    my %param_hash;
    my $rv = Getopt::Long::GetOptionsFromString($param_string,
        't=i' => \$param_hash{t},
        'k=i' => \$param_hash{k},
        'w=i' => \$param_hash{w},
        'd=i' => \$param_hash{d},
        'r=f' => \$param_hash{r},
        'c=i' => \$param_hash{c},
        'P'   => \$param_hash{P},
        'A=i' => \$param_hash{A},
        'B=i' => \$param_hash{B},
        'O=i' => \$param_hash{O},
        'E=i' => \$param_hash{E},
        'L=i' => \$param_hash{L},
        'U=i' => \$param_hash{U},
        'p'   => \$param_hash{p},
        'R=s' => \$param_hash{R},
        'T=i' => \$param_hash{T},
        'a'   => \$param_hash{a},
        'C'   => \$param_hash{C},
        'H'   => \$param_hash{H},
        'M'   => \$param_hash{M},
        'v=i' => \$param_hash{v},
    );

    die $self->error_message("Failed to parse parameter string: $param_string") unless $rv;

    my @switches = qw(a C H M p P);

    for my $key (keys %param_hash) {
        if (not defined $param_hash{$key}) {
            delete $param_hash{$key};
            next;
        }
        if (grep { $key eq $_ } @switches) {
            if ($param_hash{$key} == 1) {
                $param_hash{$key} = '';
            } else {
                delete $param_hash{$key};
            }
        }
    }

    return \%param_hash;
}

# helper for decomposed_aligner_params and aligner_params_for_sam_header
sub join_aligner_params_hash {
    my $self = shift;
    my $param_hash = shift;

    my @param_list;

    for my $key (sort { $a cmp $b } keys %$param_hash) {
        my $val = $param_hash->{$key};
        push @param_list, "-$key";
        push @param_list, $val if $val;
    }

    return join ' ', @param_list;
}

sub fillmd_for_sam {
    return 1;
}

sub requires_read_group_addition {
    return 0;
}

sub supports_streaming_to_bam {
    return 0;
}

sub multiple_reference_mode {
    return 0;
}

sub accepts_bam_input {
    return 0;
}

# Bwa mem should just find a corresponding Bwa index and symlink it. This is the
# best we can do within the existing framework if we don't want to recreate an
# identical index already created by the 'regular' bwa module.
sub prepare_reference_sequence_index {
    my $class = shift;
    my $refindex = shift;

    my $staging_dir = $refindex->temp_staging_directory;

    $class->status_message("Bwa mem version 0.7.2 is looking for a bwa version 0.7.2 index.");

    Genome::Sys->create_symlink($refindex->reference_build->get_sequence_dictionary("sam"), $staging_dir ."/all_sequences.dict" );

    my $bwa_index = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_or_create(
        reference_build_id => $refindex->reference_build_id,
        aligner_name       => 'bwa',
        #aligner_params     => $refindex->aligner_params, # none of the aligner params should affect the index step so I think this okay
        aligner_version    => $refindex->aligner_version,
        test_name          => $ENV{GENOME_ALIGNER_INDEX_TEST_NAME},
    );

    for my $filepath (glob($bwa_index->output_dir . "/*")){
        my $filename = File::Basename::fileparse($filepath);
        next if $filename eq 'all_sequences.fa';
        next if $filename eq 'all_sequences.dict';
        Genome::Sys->create_symlink($filepath, $staging_dir . "/$filename");
    }

    $bwa_index->add_user(
        label => 'uses',
        user  => $refindex
    );

    return 1;
}

1;

