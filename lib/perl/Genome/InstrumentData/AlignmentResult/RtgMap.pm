
package Genome::InstrumentData::AlignmentResult::RtgMap;

use strict;
use warnings;
use File::Basename;
use File::Path;
use Genome;

class Genome::InstrumentData::AlignmentResult::RtgMap{
    is => 'Genome::InstrumentData::AlignmentResult',
    
    has_constant => [
        aligner_name => { value => 'rtg map', is_param=>1 },
    ],
    has => [
        _file_input_option =>   { default_value => 'fastq', is_optional => 1, is_transient => 1 },
    ]
};

sub required_arch_os { 'x86_64' }

sub required_rusage { 
    #"-R 'select[model!=Opteron250 && type==LINUX64 && tmp>90000 && mem>50000] span[hosts=1] rusage[tmp=90000, mem=50000]' -M 50000000 -n 8 -m hmp -q hmp"; #HMP queue params
    "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>90000 && mem>48000] span[hosts=1] rusage[tmp=90000, mem=48000]' -M 48000000 -n 8 -q $ENV{GENOME_LSF_QUEUE_ALIGNMENT_DEFAULT}";
}

sub _decomposed_aligner_params {
    my $self = shift;

    #   -U produce unmapped sam
    #   -Z do not zip sam
    my $aligner_params = ($self->aligner_params || '') . " -U -Z --read-names"; #append core & space

    my $cpu_count = $self->_available_cpu_count;
    $aligner_params .= " -T $cpu_count";

    if (!$self->instrument_data->is_paired_end || $self->force_fragment) {
        $self->debug_message("Running fragment mode.  Params were $aligner_params.  Removing an -E if it exists since that's not applicable to fragment runs.");
        $aligner_params =~ s/-E\s?.*?\s+//g;
    }
    
    return ('rtg_aligner_params' => $aligner_params);
}

sub _run_aligner {
    my $self = shift;
    my @input_pathnames = @_;
    $ENV{'RTG_MEM'} = ($ENV{'TEST_MODE'} ? '1G' : '45G');
    $self->debug_message("RTG Memory request is $ENV{RTG_MEM}");

    # get refseq info
    my $reference_build = $self->reference_build;
    
    my $reference_sdf_path = $self->get_reference_sequence_index->full_consensus_path('sdf');

    unless (-e $reference_sdf_path) {
        $self->error_message("sdf path not found in " . $reference_build->data_directory);
        die $self->error_message;
    }
    
    my $scratch_directory = $self->temp_scratch_directory;
    my $staging_directory = $self->temp_staging_directory;
 
    #   To run RTG, have to first convert ref and inputs to sdf, with 'rtg format', for which you 
    #   have to designate a destination directory

    # disconnect db before long-running action 
    if (Genome::DataSource::GMSchema->has_default_handle) {
        $self->debug_message("Disconnecting GMSchema default handle.");
        Genome::DataSource::GMSchema->disconnect_default_dbh();
    }
    #STEP 1 - convert input to sdf
    my $prechunk_input_sdf = File::Temp::tempnam($scratch_directory, "input-XXX") . ".sdf"; #destination of converted input
    my $rtg_fmt = Genome::Model::Tools::Rtg->path_for_rtg_format($self->aligner_version);
    my $cmd;

    if (@input_pathnames == 1) {
        $self->debug_message("_run_aligner called in single-ended mode.");
        $cmd = sprintf('%s --format=%s --quality-format=sanger -o %s %s',
                        $rtg_fmt,
                        $self->_file_input_option,
                        $prechunk_input_sdf,
                        $input_pathnames[0]);
    } elsif (@input_pathnames == 2) {
        $self->debug_message("_run_aligner called in paired-end mode.");
        $cmd = sprintf('%s --format=%s --quality-format=sanger -o %s -l %s -r %s', #specify paired ends as "l" and "r"
                        $rtg_fmt,
                        $self->_file_input_option,
                        $prechunk_input_sdf,
                        $input_pathnames[0],  
                        $input_pathnames[1]);
    } else {
        $self->error_message("_run_aligner called with " . scalar @input_pathnames . " files.  It should only get 1 or 2!");
        die $self->error_message;
    }
    
    Genome::Sys->shellcmd(
            cmd                 => $cmd, 
            input_files         => \@input_pathnames,
            output_directories  => [$prechunk_input_sdf],
            skip_if_output_is_present => 0,
    );

    #check sdf output was created
    my @idx_files = glob("$prechunk_input_sdf/*");
    if (!@idx_files > 0) {
        $self->error_message(sprintf("rtg formatting of [%s] failed  with %s", join " ", @input_pathnames, $cmd ));
        die $self->error_message;
    }

    my $chunk_path = $scratch_directory . "/chunks";

    $self->debug_message("Chunking....");
    my $chunk_size = ($ENV{'TEST_MODE'} ? 200 : 1000000);
    my $chunk_cmd = sprintf("%s -n %s -o %s %s", Genome::Model::Tools::Rtg->path_for_rtg_sdfsplit($self->aligner_version), $chunk_size, $chunk_path, $prechunk_input_sdf);
    Genome::Sys->shellcmd(
            cmd                 => $chunk_cmd, 
            output_directories  => [$chunk_path],
            skip_if_output_is_present => 0,
    );

    my @chunks;
     # chunk paths all have numeric names 
    for my $chunk (grep {basename($_) =~ m/^\d+$/} glob($chunk_path . "/*")) {
        $self->debug_message("Adding chunk for analysis ... $chunk");
        push @chunks, $chunk;
    }

    $self->debug_message("Removing original unchunked input sdf");
    rmtree($prechunk_input_sdf);

    for my $input_sdf (@chunks) {
        my $output_dir = File::Temp::tempnam($scratch_directory, "output-XXX") . ".sdf";  
        my @output_files; 
        if (@input_pathnames == 1) {
            @output_files = ("$output_dir/alignments.sam", "$output_dir/unmapped.sam");
        } elsif (@input_pathnames == 2) {
            @output_files = ("$output_dir/unmated.sam", "$output_dir/mated.sam", "$output_dir/unmapped.sam");
        } 
        
        #STEP 2 - run rtg map aligner  
        my %aligner_params = $self->_decomposed_aligner_params;
        my $rtg_map = Genome::Model::Tools::Rtg->path_for_rtg_map($self->aligner_version);
        my $rtg_aligner_params = (defined $aligner_params{'rtg_aligner_params'} ? $aligner_params{'rtg_aligner_params'} : "");
        $cmd = sprintf('%s -t %s -i %s -o %s %s ', 
                        $rtg_map,
                        $reference_sdf_path,
                        $input_sdf,
                        $output_dir,
                        $rtg_aligner_params);
        
        #rtg grabs all the CPU it can, so create restriction
            
        Genome::Sys->shellcmd(
            cmd          => $cmd,
            input_files  => [ $input_sdf ],
            output_files => [@output_files, "$output_dir/done"],
            skip_if_output_is_present => 0,
        );

        #STEP 3.1 - Collate, Format and Append sam files 
        my $sam_file = $self->temp_scratch_directory . "/all_sequences.sam";
        my $sam_file_fh = IO::File->new(">> $sam_file");

        foreach my $file_to_append (@output_files)
        {
            my $file_to_append_fh = IO::File->new( $file_to_append);
            while (<$file_to_append_fh>)
            {
                unless ($_ =~ /^@/) # eat headers
                {
                    #fix CIGAR string - RTG puts in "=", replace with "M" - assumes 6th tabbed column    
                    my @line = split ("\t", $_);
                    $line[5]=~s/=/M/g;

                    $sam_file_fh->print (join("\t", @line));
                }
            }
            $file_to_append_fh->close;
        }
        $sam_file_fh->close;

        # confirm that at the end we have a nonzero sam file, this is what'll get turned into a bam and copied out.
        unless (-s $sam_file) {
            die "The sam output file $sam_file is zero length; something went wrong.";
        }

        #STEP 4 - append log files 
        my $log_input = "$output_dir/map.log";
        my $log_output = $self->temp_staging_directory . "/rtg_map.log";
        $cmd = sprintf('cat %s >> %s', $log_input, $log_output);        

        Genome::Sys->shellcmd(
            cmd          => $cmd,
            input_files  => [ $log_input ],
            output_files => [ $log_output ],
            skip_if_output_is_present => 0,
        );

        # cleanup for the next pass, if there is one
        rmtree($output_dir);
        rmtree($input_sdf);

    }

    rmtree($chunk_path);

    return 1;
}

sub aligner_params_for_sam_header {
    my $self = shift;
    my $cmd = Genome::Model::Tools::Rtg->path_for_rtg_map($self->aligner_version);
    my %params = $self->_decomposed_aligner_params;
    my $aln_params = $params{rtg_aligner_params};
    
    return "$cmd $aln_params"; 
}

sub fillmd_for_sam {
    return 1;
}

sub _check_read_count {
    return 1;
}

sub prepare_reference_sequence_index {
    my $class = shift;

    my $refindex = shift;

    my $staging_dir = $refindex->temp_staging_directory;
    my $staged_fasta_file = sprintf("%s/all_sequences.fa", $staging_dir);

    my $actual_fasta_file = $staged_fasta_file;    

    $class->debug_message("Making an RTG index out of $staged_fasta_file");

    my $rtg_path = Genome::Model::Tools::Rtg->path_for_rtg_format($refindex->aligner_version);

    my $cmd = sprintf("%s -o %s/all_sequences.sdf %s", $rtg_path, $staging_dir, $staged_fasta_file);

    my $rv = Genome::Sys->shellcmd(
        cmd=>$cmd,
        output_files => ["$staging_dir/all_sequences.sdf"],
    );

    unless ($rv) {
        $class->error_message("rtg indexing failed");
        return;
    }

    return 1;

}

sub _compute_alignment_metrics {
    return 1;
}

1;

