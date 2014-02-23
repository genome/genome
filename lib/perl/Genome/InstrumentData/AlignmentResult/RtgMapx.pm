package Genome::InstrumentData::AlignmentResult::RtgMapx;

use strict;
use warnings;
use File::Basename;
use File::Path;
use File::Copy;
use Genome;

class Genome::InstrumentData::AlignmentResult::RtgMapx{
    is => 'Genome::InstrumentData::AlignmentResult',
    
    has_constant => [
        aligner_name => { value => 'rtg mapx', is_param=>1 },
    ],
    has => [
        _max_read_id_seen => { default_value => 0, is_optional => 1, is_transient => 1 },
        _file_input_option =>   { default_value => 'fastq', is_optional => 1, is_transient => 1 },
    ]
};

sub required_arch_os { 'x86_64' }

sub required_rusage { 
    "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>90000 && mem>25000] span[hosts=1] rusage[tmp=90000, mem=25000]' -M 25000000 -n 8 -m hmp -q hmp";
}

sub _decomposed_aligner_params {
    my $self = shift;

    # -U = report unmapped reads
    # --read-names print real read names

    my $aligner_params = ($self->aligner_params || '') . " -U --read-names -Z "; 

    my $cpu_count = $self->_available_cpu_count;
    $aligner_params .= " -T $cpu_count";
    
    return ('rtg_aligner_params' => $aligner_params);
}

sub _run_aligner {
    my $self = shift;
    my @input_pathnames = @_;
    
    $ENV{'RTG_MEM'} = ($ENV{'TEST_MODE'} ? '1G' : '23G');
    $self->debug_message("RTG Memory limit is $ENV{RTG_MEM}");

    if (@input_pathnames == 1) {
        $self->debug_message("_run_aligner called in single-ended mode.");
    } elsif (@input_pathnames == 2) {
        $self->debug_message("_run_aligner called in paired-end mode.  We don't actually do paired alignment with MapX though; running two passes.");
    } else {
        $self->error_message("_run_aligner called with " . scalar @input_pathnames . " files.  It should only get 1 or 2!");
        die $self->error_message;
    }

    # get refseq info
    my $reference_build = $self->reference_build;
    
    my $reference_sdf_path = $self->get_reference_sequence_index->full_consensus_path('sdf'); 
    
    unless (-e $reference_sdf_path) {
        $self->error_message("sdf path not found in " . $reference_build->data_directory);
        die $self->error_message;
    }

    my $sam_file = $self->temp_scratch_directory . "/all_sequences.sam";
    my $sam_file_fh = IO::File->new(">>" . $sam_file );
    my $unaligned_file = $self->temp_scratch_directory . "/unaligned.txt";
    my $unaligned_file_fh = IO::File->new(">>" . $unaligned_file); 
        
    my $scratch_directory = $self->temp_scratch_directory;
    my $staging_directory = $self->temp_staging_directory;

    my @chunks;

    foreach my $i(0...$#input_pathnames)
    {

        my $input_pathname = $input_pathnames[$i];

        my $chunk_path = $scratch_directory . "/chunks/chunk-from-" . $i;

        #   To run RTG, have to first convert ref and inputs to sdf, with 'rtg format', 
        #   for which you have to designate a destination directory

        # disconnect db before long-running action 
        if (Genome::DataSource::GMSchema->has_default_handle) {
            $self->debug_message("Disconnecting GMSchema default handle.");
            Genome::DataSource::GMSchema->disconnect_default_dbh();
        }

        #STEP 1 - convert input to sdf
        my $input_sdf = File::Temp::tempnam($scratch_directory, "input-XXX") . ".sdf"; #destination of converted input
        my $rtg_fmt = Genome::Model::Tools::Rtg->path_for_rtg_format($self->aligner_version);
        my $cmd;

        $cmd = sprintf('%s --format=%s --quality-format sanger -o %s %s',
                $rtg_fmt,
                $self->_file_input_option,
                $input_sdf,
                $input_pathname);  

        Genome::Sys->shellcmd(
                cmd                 => $cmd, 
                input_files         => [$input_pathname],
                output_directories  => [$input_sdf],
                skip_if_output_is_present => 0,
                );

        #check sdf output was created
        my @idx_files = glob("$input_sdf/*");
        if (!@idx_files > 0) {
            die("rtg formatting of [$input_pathname] failed  with $cmd");
        }

        # chunk the SDF

        my $chunk_size = ($ENV{'TEST_MODE'} ? 200 : 1000000);;

        $self->debug_message("Chunking....");
        my $chunk_cmd = sprintf("%s -n %s -o %s %s", Genome::Model::Tools::Rtg->path_for_rtg_sdfsplit($self->aligner_version), $chunk_size, $chunk_path, $input_sdf);
        Genome::Sys->shellcmd(
                cmd                 => $chunk_cmd, 
                output_directories  => [$chunk_path],
                skip_if_output_is_present => 0,
        );

        unless (-e $chunk_path . "/done") {
            $self->error_message("Chunk failed! cmd was $chunk_cmd");
            die $self->error_message;
        }
       
        my $log_input = $chunk_path . "/sdfsplit.log"; 
        my $log_output = $self->temp_staging_directory . "sdfsplit.log";
        $cmd = sprintf('cat %s >> %s', $log_input, $log_output);   

        Genome::Sys->shellcmd(
                cmd          => $cmd,
                input_files  => [ $log_input ],
                output_files => [ $log_output ],
                skip_if_output_is_present => 0
        );
    
        # chunk paths all have numeric names 
        for my $chunk (grep {basename($_) =~ m/^\d+$/} glob($chunk_path . "/*")) {
            $self->debug_message("Adding chunk for analysis ... $chunk");
            push @chunks, $chunk;
        }

        $self->debug_message("Removing original unchunked input sdf");
        rmtree($input_sdf);
        
    }

    for my $input_sdf (@chunks) {
        my $output_dir = File::Temp::tempnam($scratch_directory, "output-XXX") . ".sdf";  
        my %output_files = (aligned_file =>"$output_dir/alignments.txt", unaligned_file => "$output_dir/unmapped.txt"); 

        #STEP 2 - run rtg mapx aligner  
        my %aligner_params = $self->_decomposed_aligner_params;
        my $rtg_mapx = Genome::Model::Tools::Rtg->path_for_rtg_mapx($self->aligner_version);
        my $rtg_aligner_params = (defined $aligner_params{'rtg_aligner_params'} ? $aligner_params{'rtg_aligner_params'} : "");
        my $cmd = sprintf('%s -t %s -i %s -o %s %s', 
                $rtg_mapx,
                $reference_sdf_path,
                $input_sdf,
                $output_dir,
                $rtg_aligner_params);

        Genome::Sys->shellcmd(
                cmd          => $cmd,
                input_files  => [ $reference_sdf_path, $input_sdf ],
                output_files => [values (%output_files), "$output_dir/done"],
                skip_if_output_is_present => 0,
                );

        # Copy log files 
        my $log_input = "$output_dir/mapx.log";
        my $log_output = $self->temp_staging_directory . "/rtg_mapx.log";
        $cmd = sprintf('cat %s >> %s', $log_input, $log_output);   

        Genome::Sys->shellcmd(
                cmd          => $cmd,
                input_files  => [ $log_input ],
                output_files => [ $log_output ],
                skip_if_output_is_present => 0
                );

        for (values %output_files) {
            $self->debug_message("Copying $_ into staging...");

            my $output_file = $self->temp_staging_directory . "/" . basename($_);

            Genome::Sys->shellcmd(
                cmd=>sprintf('cat %s >> %s', $_, $output_file),
                input_files => [$_],
                output_files => [$output_file],
                skip_if_output_is_present => 0
            );
        }
    
        $self->debug_message("removing old output to save disk");
        rmtree($output_dir);

    } 
    $self->debug_message("removing chunks to save disk");
    rmtree($self->temp_scratch_directory . "/chunks");

    return 1;
}

#sub input_chunk_size {
#    return 3_000_000;
#}

sub _compute_alignment_metrics 
{
    return 1;
}

sub prepare_scratch_sam_file {
    return 1;
}

sub create_BAM_in_staging_directory {
    return 1;
}

sub postprocess_bam_file {
    return 1;
}

sub _prepare_reference_sequences {
    my $self = shift;
    my $reference_build = $self->reference_build;

    my $ref_basename = File::Basename::fileparse($reference_build->full_consensus_path('fa'));
    my $reference_fasta_path = sprintf("%s/%s", $reference_build->data_directory, $ref_basename);

    unless(-e $reference_fasta_path) {
        $self->error_message("Alignment reference path $reference_fasta_path does not exist");
        die $self->error_message;
    }

    return 1;
}

sub prepare_reference_sequence_index {
    my $class = shift;
    
    $DB::single = 1;

    my $refindex = shift;

    my $staging_dir = $refindex->temp_staging_directory;
    my $staged_fasta_file = sprintf("%s/all_sequences.fa", $staging_dir);

    my $actual_fasta_file = $staged_fasta_file;    

    $class->debug_message("Making an RTG index out of $staged_fasta_file");

    my $rtg_path = Genome::Model::Tools::Rtg->path_for_rtg_format($refindex->aligner_version);

    my $cmd = sprintf("%s --protein -o %s/all_sequences.sdf %s", $rtg_path, $staging_dir, $staged_fasta_file);

    my $rv = Genome::Sys->shellcmd(
        cmd=>$cmd
    );

    unless ($rv) {
        $class->error_message("rtg indexing failed");
        return;
    }

    return 1;
}
