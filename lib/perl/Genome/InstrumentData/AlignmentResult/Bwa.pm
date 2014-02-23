package Genome::InstrumentData::AlignmentResult::Bwa;

use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use Carp qw/confess/;
use Genome;

class Genome::InstrumentData::AlignmentResult::Bwa {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => { value => 'bwa', is_param=>1 },
    ],
    has_transient_optional => [
         _bwa_sam_cmd => { is=>'Text' }
    ]
};

sub required_arch_os { 'x86_64' }

sub required_rusage {
    my $class = shift;
    my %p = @_;
    my $instrument_data = delete $p{instrument_data};
    my $aligner_params  = delete $p{aligner_params};
    my $reference_build = delete $p{reference_build};

    my $tmp_mb = $class->tmp_megabytes_estimated($instrument_data);
    my $mem_mb = 1024 * 8; 
    my $cpus = 2;

    if ($aligner_params and $aligner_params =~ /-t\s*([0-9]+)/) {
        $cpus = $1;
    }

    if ($reference_build and $reference_build->id eq '107494762') {
        $class->status_message(sprintf 'Doubling memory requirements for alignments against %s.', $reference_build->name);
        $mem_mb *= 2;
    }

    my $mem_kb = $mem_mb*1024;
    my $tmp_gb = $tmp_mb/1024;

    my $user = getpwuid($<);
    my $queue = $ENV{GENOME_LSF_QUEUE_ALIGNMENT_DEFAULT};
    $queue = $ENV{GENOME_LSF_QUEUE_ALIGNMENT_PROD} if (Genome::Config->should_use_alignment_pd);

    my $host_groups;
    $host_groups = qx(bqueues -l $queue | grep ^HOSTS:);
    if ($host_groups =~ /all hosts/) {
        $host_groups = '';
    }
    else {
        $host_groups =~ s/\/\s+/\ /;
        $host_groups =~ s/^HOSTS:\s+//;
    }

    my $select  = "select[ncpus >= $cpus && mem >= $mem_mb && gtmp >= $tmp_gb] span[hosts=1]";
    my $rusage  = "rusage[mem=$mem_mb, gtmp=$tmp_gb]";
    my $options = "-M $mem_kb -n $cpus -q $queue";

    my $required_usage = "-R '$select $rusage' $options";

    Workflow::Dispatcher::Lsf->class;
    no warnings;
    unless ($Workflow::Dispatcher::Lsf::OPENLAVA) {
        #check to see if our resource requests are feasible (This uses "maxmem" to check theoretical availability)
        #factor of four is based on current six jobs per host policy this should be revisited later
        my $select_check = "select[ncpus >= $cpus && maxmem >= " . ($mem_mb * 4) . " && maxgtmp >= $tmp_gb] span[hosts=1]";
        my $select_cmd = "bhosts -R '$select_check' $host_groups";

        my @selected_blades = qx($select_cmd);
        unless (@selected_blades) {
            die $class->error_message("Failed to find hosts that meet resource requirements ($required_usage). [Looked with `$select_cmd`]");
        }
    }

    return $required_usage;
}

sub tmp_megabytes_estimated {
    my $class = shift || die;
    my $instrument_data = shift;

    my $default_megabytes = 90000;


    if (not defined $instrument_data) {
        return $default_megabytes;
    }
    elsif ($instrument_data->bam_path) {
        my $bam_path = $instrument_data->bam_path;

        my $scale_factor = 3.25; # assumption: up to 3x during sort/fixmate/sort and also during fastq extraction (2x) + bam = 3

        my $bam_bytes = -s $bam_path;
        unless ($bam_bytes) {
            #die $class->error_message("Instrument Data " . $instrument_data->id  . " has BAM ($bam_path) but has no size!");
        }

        if ($instrument_data->can('get_segments')) {
            my $bam_segments = scalar $instrument_data->get_segments;
            if ($bam_segments > 1) {
                $scale_factor = $scale_factor / $bam_segments;
            }
        }

        return int(($bam_bytes * $scale_factor) / 1024**2);
    }
    elsif ($instrument_data->can("calculate_alignment_estimated_kb_usage")) {
        my $kb_usage = $instrument_data->calculate_alignment_estimated_kb_usage;
        return int(($kb_usage * 3) / 1024) + 100; # assumption: 2x for the quality conversion, one gets rm'di after; aligner generates 0.5 (1.5 total now); rm orig; sort and merge maybe 2-2.5
    }
    else {
        return $default_megabytes;
    }

    return;
}

sub _all_reference_indices {
    my $self = shift;
    my @overrides = @_;

    my @indices;
    my $b = $self->reference_build;
    if ($self->multiple_reference_mode) {
        do {
            $self->status_message("Getting reference sequence index for build ".$b->__display_name__);
            $self->status_message("...using overrides @overrides\n") if @overrides;
            my $index = $self->get_reference_sequence_index($b,@overrides);
            $self->status_message("Index: " . $index->__display_name__);
            unshift(@indices, $index);
            $b = $b->append_to;
        } while ($b);
    } else {
        # TODO: the above condition works for the single-layer reference too right? -ssmith
        $self->status_message("Getting reference sequence index for build ".$b->__display_name__);
        $self->status_message("...using overrides @overrides\n") if @overrides;
        my $index = $self->get_reference_sequence_index($b,@overrides);
        $self->status_message("Index: " . $index->__display_name__);
        unshift(@indices, $index);
    }
    return @indices;
}

sub _intermediate_result {
    my ($self, $params, $index, @input_files) = @_;

    my @results;
    for my $idx (0..$#input_files) {
        my $path = $input_files[$idx];
        my ($input_pass) = $path =~ m/\.bam:(\d)$/;
        print "INPUT FILE: $path, INPUT PASS: $input_pass\n" . Dumper(\@_);
        if (defined($input_pass)) {
            $path =~ s/\.bam:\d$/\.bam/;
        } else {
            $input_pass = $idx+1;
        }

        my $aligner_version = $self->aligner_version;
        if ($aligner_version =~ /^(.*)-i(.*)/) {
            my $old = $aligner_version;
            $aligner_version = $1;
            $self->warning_message("FOR iBWA (BWA $old), USING (IDENTICAL) $aligner_version FOR INTERMEDIATE RESULTS"); 
        }

        my %intermediate_params = (
            instrument_data_id           => $self->instrument_data->id,
            aligner_name                 => $self->aligner_name,
            aligner_version              => $aligner_version,
            aligner_params               => $params,
            aligner_index_id             => $index->id,
            flagstat_file                => $self->_flagstat_file,
            input_file                   => $path,
            input_pass                   => $input_pass,
            instrument_data_segment_type => $self->instrument_data_segment_type,
            instrument_data_segment_id   => $self->instrument_data_segment_id,
            samtools_version             => $self->samtools_version,
            trimmer_name                 => $self->trimmer_name,
            trimmer_version              => $self->trimmer_version,
            trimmer_params               => $self->trimmer_params,
            test_name                    => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
        );

        my $includes = join(' ', map { '-I ' . $_ } UR::Util::used_libs);
        my $class = 'Genome::InstrumentData::IntermediateAlignmentResult::Command::Bwa';
        my $parameters = join(', ', map($_ . ' => "' . (defined($intermediate_params{$_}) ? $intermediate_params{$_} : '') . '"', sort keys %intermediate_params));

        if(UR::DBI->no_commit()) {
            my $rv = eval "$class->execute($parameters);";
            if(!$rv or $@) {
                my $err = $@;
                die('Failed to generate intermediate result!' . ($err? $err : ' command returned false') );
            }
        } else {
            my $cmd = qq{$^X $includes -e 'use above "Genome"; $class->execute($parameters); UR::Context->commit;' };
            my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd) };
            if(!$rv or $@) {
                my $err = $@;
                die('Failed to generate intermediate result!' . ($err? $err : ' command returned false') );
            }
        }

        my $intermediate_result = Genome::InstrumentData::IntermediateAlignmentResult::Bwa->get_with_lock(%intermediate_params);
        unless ($intermediate_result) {
            confess "Failed to generate IntermediateAlignmentResult for $path, params were: " . Dumper(\%intermediate_params);
        }

        $self->status_message(sprintf("Got/created an intermediate alignment result %s with file path %s", $intermediate_result->id, $intermediate_result->sai_file));
        push(@results, $intermediate_result);
    }

    my @bad_results = grep {!-e $_->sai_file || !-s $_->sai_file} @results;
    if (@bad_results > 0) {
        my $str_bad_ids = join " ", map {$_->id} @bad_results;
        confess sprintf("The following intermediate alignment result(s) have nonexistent or zero-length SAI files -- cannot proceed: %s", $str_bad_ids);
    }

    for my $result (@results) {
        $result->add_user(user => $self, label => 'uses');
    }

    return @results;
}

sub _samxe_cmdline {
    my ($self, $aligner_params, $input_groups, @input_pathnames) = @_;

    my $cmdline;
    if (@input_pathnames == 1) {
        my $params = $aligner_params->{'bwa_samse_params'} || '';
        $self->_bwa_sam_cmd("bwa samse " . $params);

        unless (scalar @$input_groups == 1) {
            die "Multiple reference mode not supported for bwa samse!";
        }

        my $aligner_index = $input_groups->[0]->{aligner_index};
        my $intermediate_results = $input_groups->[0]->{intermediate_results};
        $cmdline = sprintf("samse $params %s %s %s",
            $aligner_index->full_consensus_path('fa'),
            $input_groups->[0]->{intermediate_results}->[0]->sai_file,
            $input_pathnames[0]);

    } elsif (@input_pathnames == 2) {
        my $params = $self->_derive_bwa_sampe_parameters;
        $self->_bwa_sam_cmd("bwa sampe " . $params);

        my @cmdline_inputs;
        for my $group (@$input_groups) {
            my $aligner_index = $group->{aligner_index};
            my $intermediate_results = $group->{intermediate_results};
            push(@cmdline_inputs,
                $aligner_index->full_consensus_path('fa'),
                map { $_->sai_file } @$intermediate_results,
                );
        }

        # fastq/bam input files come after the first set of "ref.fa seq1.sai seq2.sai"
        # insert them where they need to go
        splice(@cmdline_inputs, 3, 0, @input_pathnames);
        $cmdline = "sampe $params -P " . join(' ', @cmdline_inputs);
    } else {
        $self->error_message("Input pathnames should have 2 elements... contents: " . Dumper(\@input_pathnames) );
    }

    return $cmdline;
}

sub _run_aligner {
    my $self = shift;
    my @input_pathnames = @_;

    my $tmp_dir = $self->temp_scratch_directory;

    # get refseq info
    my $reference_build = $self->reference_build;
    my $reference_fasta_path = $self->get_reference_sequence_index->full_consensus_path('fa');

    # decompose aligner params for each stage of bwa alignment
    my %aligner_params = $self->decomposed_aligner_params;

    #### STEP 1: Use "bwa aln" to align each fastq independently to the reference sequence

    my $bwa_aln_params = (defined $aligner_params{'bwa_aln_params'} ? $aligner_params{'bwa_aln_params'} : "");

    my $aligner_version = $self->aligner_version;
    if ($aligner_version =~ /^(.*)-i(.*)/) {
        my $old = $aligner_version;
        $aligner_version = $1;
        $self->warning_message("FOR iBWA (BWA $old), USING (IDENTICAL) $aligner_version FOR THE REFERENCE INDEX");
    }

    my @indices = $self->_all_reference_indices(aligner_version => $aligner_version);
    for (@indices) {
        if ($_->aligner_version =~ /i/) {
            die "got an index for an ibwa aligner???" . Data::Dumper::Dumper($_);
        }
    }
    my @input_groups;
    my @aln_log_files;
    for my $index (@indices) {
        my @results = $self->_intermediate_result($bwa_aln_params, $index, @input_pathnames);
        push(@input_groups, {
                aligner_index => $index,
                intermediate_results => [ @results ],
            }
        );
        push(@aln_log_files, map { $_->log_file } @results);
    }

    #### STEP 2: Use "bwa samse" or "bwa sampe" to perform single-ended or paired alignments, respectively.
    #### Runs once for ALL input files

    map { s/\.bam:\d/.bam/ } @input_pathnames; # strip :[12] suffix from bam files if present

    my $samxe_logfile = $tmp_dir . "/bwa.samxe.log";
    my $samxe_cmdline = $self->_samxe_cmdline(\%aligner_params, \@input_groups, @input_pathnames);
    my $sam_file = $self->temp_scratch_directory . "/all_sequences.sam";

    my $sam_command_line = sprintf("%s $samxe_cmdline 2>> $samxe_logfile",
        Genome::Model::Tools::Bwa->path_for_bwa_version($self->aligner_version));

    unless ($self->_filter_samxe_output($sam_command_line, $sam_file)) {
        die "Failed to process sam command line, error_message is ".$self->error_message;
    }

    unless (-s $sam_file) {
        die "The sam output file $sam_file is zero length; something went wrong.";
    }

    for my $file (@input_pathnames) {
        if ($file =~ m/^\/tmp\//) {
            $self->status_message("removing $file to save space!");
            unlink($file);
        }
    }

    #### STEP 3: verify the samxe_logfile
    unless ($self->_verify_bwa_samxe_did_happen($samxe_logfile)) {
        die $self->error_message("bwa samxe seems to fail based on run log: $samxe_logfile");
    }

    #### STEP 4: Merge log files.

    my $log_input_fileset = join " ",  (@aln_log_files, $samxe_logfile);
    my $log_output_file   = $self->temp_staging_directory . "/aligner.log";
    my $concat_log_cmd = sprintf('cat %s >> %s', $log_input_fileset, $log_output_file);

    Genome::Sys->shellcmd(
        cmd          => $concat_log_cmd,
        input_files  => [ @aln_log_files, $samxe_logfile ],
        output_files => [ $log_output_file ],
        skip_if_output_is_present => 0,
    );

    return 1;
}

sub _filter_samxe_output {
    my ($self, $sam_cmd, $sam_file_name) = @_;

#    my $cmd = "$sam_cmd | grep -v ^@ >> $sam_file_name";
#    print $cmd,"\n\n";
#    $DB::single = 1;
#    Genome::Sys->shellcmd(cmd => $cmd);
#    return 1;

    my $sam_run_output_fh = IO::File->new( $sam_cmd . "|" );
    binmode $sam_run_output_fh;
    $self->status_message("Running $sam_cmd");
    if ( !$sam_run_output_fh ) {
            $self->error_message("Error running $sam_cmd $!");
            return;
    }

    my $sam_out_fh;
    # UGLY HACK: the multi-aligner code redefines this to zero so it can extract sam files.
    if ($self->supports_streaming_to_bam) {
        $sam_out_fh = $self->_sam_output_fh;
        $self->status_message("Streaming output through existing file handle");
    } else {
        $sam_out_fh = IO::File->new(">>" . $self->temp_scratch_directory . "/all_sequences.sam");
        $self->status_message("Opened for output ( " . $self->temp_scratch_directory . "/all_sequences.sam )");
    }
    my $add_rg_cmd = Genome::Model::Tools::Sam::AddReadGroupTag->create(
            input_filehandle     => $sam_run_output_fh,
            output_filehandle    => $sam_out_fh,
            read_group_tag => $self->read_and_platform_group_tag_id,
            pass_sam_headers => 0,
        );

    unless ($add_rg_cmd->execute) {
        $self->error_message("Adding read group to sam file failed!");
        die $self->error_message;
    }

    $sam_run_output_fh->close;

    return 1;
}

#For very few cases, bwa samxe still silently runs and gives bad
#output when ref_seq and sai become not accessible during the middle of
#process that is caused by things like power outage. refer to #81719
sub _verify_bwa_samxe_did_happen {
    my ($self, $log_file) = @_;

    unless ($log_file and -s $log_file) {
        $self->error_message("log file $log_file is not valid");
        return;
    }

    my @last_lines = `tail -50 $log_file`;
    my $rp_ct;

    @last_lines = grep {!/^\[bwtcache_destroy\] \d+ cache waits encountered/} @last_lines;

    # XXX I am trying to avoid switching on $bwa_version in this module, but
    # I'm not entirely sure that storing 'log_format' in GMT::Bwa and defining
    # the logic for interpreting that property here is the correct separation
    # of functionality between GMT::Bwa and this module.
    my $last_line;
    my $log_format = Genome::Model::Tools::Bwa->log_format($self->aligner_version);
    if ($log_format eq 'new') {
        $last_line = $last_lines[-4];
    } else {
        $last_line = $last_lines[-1];
    }

    if ($last_line =~ /(\d+) sequences have been processed/) {
        $rp_ct = $1;
    }
    else {
        $self->error_message("The last line of samxe logfile: $last_line is not expected");
        return;
    }

    #If getting very few fastq (like in unit test), it's very likely to error out.
    my $fq_rd_ct = $self->_fastq_read_count;
    if ($fq_rd_ct and $fq_rd_ct <= 1000) {
        $self->warning_message("fastq read count is $fq_rd_ct. Skip the samxe log check");
        return 1;
    }

    #if ($self->_is_inferred_paired_end) {
        #my $fq_rp_ct = $fq_rd_ct/2;
        #unless ($rp_ct == $fq_rp_ct) {
            #$self->error_message("read pair count conflict: samxe_log shows $rp_ct, input_file shows $fq_rp_ct");
            #return;
        #}
    #}

    chomp (my $fail_ct = qx(grep -c 'fail to infer insert size: too few good pairs' $log_file));

    my $line = qx(grep -m1 'sequences have been processed' $log_file);
    my ($batch_size) = $line =~ /(\d+) sequences/;
    $self->warning_message("The batch size: $batch_size is not 262144 as expected. Check bwa source code")
        unless $batch_size == 262144;
    #hard code fail percentile threshold for now. The percentile calculation is not accurate 
    #because the last batch count is mostly smaller than 262114, but this estimate is close enough 
    #since threshold is arbitrary too.
    my $threshold = 55;
    my $fail_percentile = sprintf("%.2f", $fail_ct * $batch_size * 100 / $rp_ct);

    # NOTE: To get distribution statistics, we're pretending this is a timing.
    Genome::Utility::Instrumentation::timing(
        'alignment_result.bwa.infer_insert_size_fail_percent',
        $fail_percentile);

    if ($fail_percentile > $threshold) {
        $fail_percentile = 100 if $fail_percentile > 100; #a calculation bug
        $self->error_message("samxe failed to infer insert size on $fail_percentile% read pairs. The current threshold is $threshold%");
        return;
    }

    return 1;
}


sub decomposed_aligner_params {
    my $self = shift;
    my $params = $self->aligner_params || ":::";

    my @spar = split /\:/, $params;

    my $bwa_aln_params = $spar[0] || "";

    my $cpu_count = $self->_available_cpu_count;

    $self->status_message("[decomposed_aligner_params] cpu count is $cpu_count");

    $self->status_message("[decomposed_aligner_params] bwa aln params are: $bwa_aln_params");

    # Make sure the thread count argument matches the number of CPUs available.
    if (!$bwa_aln_params || $bwa_aln_params !~ m/-t/) {
        $bwa_aln_params .= "-t$cpu_count";
    } elsif ($bwa_aln_params =~ m/-t/) {
        $bwa_aln_params =~ s/-t ?\d/-t$cpu_count/;
    }

    $self->status_message("[decomposed_aligner_params] autocalculated CPU requirement, bwa aln params modified: $bwa_aln_params");


    return ('bwa_aln_params' => $bwa_aln_params, 'bwa_samse_params' => $spar[1], 'bwa_sampe_params' => $spar[2]);
}

sub aligner_params_for_sam_header {
    my $self = shift;

    my %params = $self->decomposed_aligner_params;
    my $aln_params = $params{bwa_aln_params} || "";

    if ($self->instrument_data->is_paired_end) {
        $self->_derive_bwa_sampe_parameters;
    }
    my $sam_cmd = $self->_bwa_sam_cmd || "";

    return "bwa aln $aln_params; $sam_cmd ";
}

sub _derive_bwa_sampe_parameters {
    my $self = shift;
    my %aligner_params = $self->decomposed_aligner_params;
    my $bwa_sampe_params = (defined $aligner_params{'bwa_sampe_params'} ? $aligner_params{'bwa_sampe_params'} : "");

    # Ignore where we have a -a already specified
    if ($bwa_sampe_params =~ m/\-a\s*(\d+)/) {
        $self->status_message("Aligner params specify a -a parameter ($1) as upper bound on insert size.");
    } else {
        # come up with an upper bound on insert size.
        my $instrument_data = $self->instrument_data;
        my $sd_above        = $instrument_data->resolve_sd_insert_size || 0;
        my $median_insert   = $instrument_data->resolve_median_insert_size || 0;
        my $upper_bound_on_insert_size= ($sd_above * 5) + $median_insert;
        if($upper_bound_on_insert_size > 0) {
            $self->status_message("Calculated a valid insert size as $upper_bound_on_insert_size.  This will be used when BWA's internal algorithm can't determine an insert size");
        } else {
            $self->status_message("Unable to calculate a valid insert size to run BWA with. Using 600 (hax)");
            $upper_bound_on_insert_size= 600;
        }

        $bwa_sampe_params .= " -a $upper_bound_on_insert_size";
    }

    # store the calculated sampe params
    $self->_bwa_sam_cmd("bwa sampe " . $bwa_sampe_params);
    return $bwa_sampe_params;
}

sub fillmd_for_sam {
    return 0;
}

sub requires_read_group_addition {
    return 0;
}

sub supports_streaming_to_bam {
    return 1;
}

sub multiple_reference_mode {
    my $self = shift;
    return defined $self->reference_build->append_to
        && Genome::Model::Tools::Bwa->supports_multiple_reference($self->aligner_version)
        && $self->instrument_data->is_paired_end;
}

sub accepts_bam_input {
    my $self = shift;
    my $rv = Genome::Model::Tools::Bwa->supports_bam_input($self->aligner_version);
    print STDERR "BWA VERSION " . $self->aligner_version . " accepts bam input: " . $rv . "\n";
    return $rv;
}

sub prepare_reference_sequence_index {
    my $class = shift;
    my $refindex = shift;

    my $staging_dir = $refindex->temp_staging_directory;
    my $staged_fasta_file = sprintf("%s/all_sequences.fa", $staging_dir);

    my $actual_fasta_file = $staged_fasta_file;

    if (-l $staged_fasta_file) {
        $class->status_message(sprintf("Following symlink for fasta file %s", $staged_fasta_file));
        $actual_fasta_file = readlink($staged_fasta_file);
        unless($actual_fasta_file) {
            $class->error_message("Can't read target of symlink $staged_fasta_file");
            return;
        }
    }

    $class->status_message(sprintf("Checking size of fasta file %s", $actual_fasta_file));

    my $fasta_size = -s $actual_fasta_file;
    my $bwa_index_algorithm = ($fasta_size < 11_000_000) ? "is" : "bwtsw";
    my $bwa_version = $refindex->aligner_version;
    my $bwa_path = Genome::Model::Tools::Bwa->path_for_bwa_version($bwa_version);

    $class->status_message(sprintf("Building a BWA index in %s using %s.  The file size is %s; selecting the %s algorithm to build it.", $staging_dir, $staged_fasta_file, $fasta_size, $bwa_index_algorithm));

    # Expected output files from bwa index. Bwa index creates files with the
    # following extensions: amb, ann, bwt, pac, and sa; older versions also
    # create: rbwt, rpac, and rsa.

    # XXX Rethink the separation of functionality between GMT::Bwa and this
    # module. This module should not have flow control based on the value of
    # $bwa_version, but how much functionality should be moved in to GMT::Bwa
    # for other modules to use?
    my @index_extensions = Genome::Model::Tools::Bwa->index_extensions($bwa_version);
    my @output_files = map {sprintf("%s.%s", $staged_fasta_file, $_)} @index_extensions;

    my $bwa_cmd = sprintf('%s index -a %s %s', $bwa_path, $bwa_index_algorithm, $staged_fasta_file);
    my $rv = Genome::Sys->shellcmd(
        cmd          => $bwa_cmd,
        input_files  => [$staged_fasta_file],
        output_files => [@output_files]
    );

    unless($rv) {
        $class->error_message('bwa indexing failed');
        return;
    }

    return 1;
}

1;

