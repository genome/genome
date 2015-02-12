package Genome::InstrumentData::AlignmentResult::BwamemBase;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::AlignmentResult::BwamemBase {
    is => 'Genome::InstrumentData::AlignmentResult',
};

sub required_arch_os { 'x86_64' }

sub aligner_name_for_aligner_index {
    return 'bwa';
}

sub requires_read_group_addition {
    return 0;
}

sub _verify_params_and_inputs {
    my $self = shift;
    my @input_paths = @_;

    if (@input_paths != 1 and @input_paths != 2) {
        $self->error_message(
            "Expected 1 or 2 input path names. Got: " . Dumper(\@input_paths));
    }

    my $reference_fasta_path = $self->_aligner_index_fasta;

    # Verify inputs and outputs.
    for (@input_paths, $reference_fasta_path) {
        die $self->error_message("Missing input '$_'.") unless -e $_;
        die $self->error_message("Input '$_' is empty.") unless -s $_;
    }

    # Verify the aligner and get params.
    my $aligner_version = $self->aligner_version;
    unless (Genome::Model::Tools::Bwa->supports_mem($aligner_version)) {
        die $self->error_message(
            "The pipeline does not support using " .
            "bwa mem with bwa-$aligner_version."
        );
    }
}

# Bwamem already correctly generates pairs, but fixmate is incompatible with
# the supplementary alignment flag (0x800) and will ruin the pairing. So, we
# turn it off.
sub requires_fixmate {
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

# Override _check_read_count() from Genome::InstrumentData::AlignmentResult to
# filter reads with secondary or supplementary alignment flags (0x100 or 0x800)
# when comparing to the fastq.
sub _check_read_count {
    my ($self, $bam_rd_ct) = @_;

    my $param_hash = $self->decomposed_aligner_params;
    my $flag = exists $param_hash->{M} ? 0x100 : 0x800;

    my $sam_path = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);
    my $cmd = "$sam_path view -F $flag -c " . $self->temp_staging_directory . "/all_sequences.bam";
    my $filtered_bam_rd_ct = `$cmd`;

    $self->debug_message("Overriding _check_read_count: filtering flag $flag from bam read count.");
    $self->debug_message("Actual read count: $bam_rd_ct; filtered read count: $filtered_bam_rd_ct");

    return $self->SUPER::_check_read_count($filtered_bam_rd_ct);
}

sub _aligner_index_fasta {
    my $self = shift;
    return $self->get_reference_sequence_index->full_consensus_path('fa');
}

sub _disconnect_from_db {
    my ($self) = @_;

    $self->debug_message("Closing data source db handle...");
    if ($self->__meta__->data_source->has_default_handle) {
        if ($self->__meta__->data_source->disconnect_default_handle) {
            $self->debug_message("Disconnected data source db handle (as expected).");
        } else {
            $self->debug_message("Unable to disconnect data source db handle.");
        }
    } else {
        $self->debug_message("Data source db handle already closed.");
    }
}

sub _check_db_connection {
    my ($self) = @_;

    if ($self->__meta__->data_source->has_default_handle) {
        $self->debug_message("Data source db handle unexpectedly reconnected itself.");
    } else {
        $self->debug_message("Data source db handle still closed (as expected).");
    }
}

sub _verify_bwa_mem_did_happen {
    my ($self, $log_file) = @_;

    unless ($log_file and -e $log_file) {
        $self->error_message("Log file $log_file is does not exist.");
        return;
    }

    unless ($log_file and -s $log_file) {
        $self->error_message("Log file $log_file is empty.");
        return;
    }

    my $line_count = 100;
    my @last_lines = `tail -$line_count $log_file`;

    if (not (
        ($last_lines[-3] =~ /^\[main\] Version:/) and
        ($last_lines[-2] =~ /^\[main\] CMD:/) and
        ($last_lines[-1] =~ /^\[main\] Real time:/) )
    ) {
        $self->error_message("Last lines of $log_file were unexpected. Dumping last $line_count lines.");
        $self->debug_message($_) for @last_lines;
        return;
    }
    return 1;
}

# Generates the param hash from the processing profile aligner_params. Corrects
# the cpu and M flags (if necessary) and the returns the param hash.
sub decomposed_aligner_params {
    my $self = shift;

    my $param_string = $self->aligner_params || '';
    my $param_hash = $self->_param_string_to_hash($param_string);

    # I'm not sure why we need all these debug messages...
    $self->debug_message(
        "[decomposed_aligner_params] unmodified bwa mem params are: "
        . $self->_param_hash_to_string($param_hash));

    $self->_fix_cpu_flag($param_hash);
    $self->_fix_M_flag($param_hash);

    $self->debug_message(
        "[decomposed_aligner_params] final bwa mem params are: "
        . $self->_param_hash_to_string($param_hash));

    return $param_hash;
}

# Gets the param hash using decomposed_aligner_params, strips out the cpu count
# flag, and returns the string using _param_hash_to_string.
sub aligner_params_for_sam_header {
    my $self = shift;

    my $param_hash = $self->decomposed_aligner_params;
    delete $param_hash->{t}; # we don't want cpu count to be in the sam header
    my $param_string = $self->_param_hash_to_string($param_hash);

    return "bwa mem $param_string";
}

# Helper for decomposed_aligner_params. Takes a param string (from the
# processing profile) and generates a param hash.
sub _param_string_to_hash {
    my $class = shift;
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

    die $class->error_message("Failed to parse parameter string: $param_string") unless $rv;

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

# Helper for decomposed_aligner_params. Forces the CPU flag to match the
# available cpu count in a param hash.
sub _fix_cpu_flag {
    my $self = shift;
    my $param_hash = shift;

    my $cpu_count = $self->_available_cpu_count;
    $self->debug_message("[_fix_cpu_flag] cpu count is $cpu_count");

    # Make sure the thread count argument matches the number of CPUs available.
    if ((not exists $param_hash->{t}) or (not defined $param_hash->{t}) or ($param_hash->{t} ne $cpu_count)) {
        $param_hash->{t} = $cpu_count;
        my $modified_param_string = $self->_param_hash_to_string($param_hash);
        $self->debug_message("[_fix_cpu_flag] autocalculated CPU requirement, bwa mem params modified: $modified_param_string");
    }

    return $param_hash;
}

# Helper for decomposed_aligner_params. Forces the M flag if we're on an older
# version of bwamem in a param hash.
sub _fix_M_flag {
    my $self = shift;
    my $param_hash = shift;

    # version check
    my $supports_supplementary_flag = Genome::Model::Tools::Bwa->supports_supplementary_alignment_flag($self->aligner_version);

    if ((not exists $param_hash->{M}) and (not $supports_supplementary_flag)) {
        $param_hash->{M} = ''; # in this instance, '' means -M is added with no argument
        my $modified_param_string = $self->_param_hash_to_string($param_hash);
        $self->debug_message("[_fix_mem_flag] forcing -M, bwa mem params modified: $modified_param_string");
    }

    return $param_hash;
}

# Takes a param hash and turns it into a string; sorts by keys before joining.
sub _param_hash_to_string {
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

# Bwa mem should just find a corresponding Bwa index and symlink it. This is the
# best we can do within the existing framework if we don't want to recreate an
# identical index already created by the 'regular' bwa module.
sub prepare_reference_sequence_index {
    my $class = shift;
    my $refindex = shift;

    my $staging_dir = $refindex->temp_staging_directory;

    my $aligner_version = $refindex->aligner_version;

    $class->debug_message("Bwa mem version $aligner_version is looking for a bwa version $aligner_version index.");

    Genome::Sys->create_symlink($refindex->reference_build->get_sequence_dictionary("sam"), $staging_dir ."/all_sequences.dict" );

    my $bwa_index = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_or_create(
        users              => $refindex->_user_data_for_nested_results,
        reference_build_id => $refindex->reference_build_id,
        aligner_name       => 'bwa',
        #aligner_params    => $refindex->aligner_params, # none of the aligner params should affect the index step so I think this okay
        aligner_version    => $aligner_version,
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
