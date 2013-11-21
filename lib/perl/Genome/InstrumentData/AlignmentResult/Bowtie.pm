package Genome::InstrumentData::AlignmentResult::Bowtie;

use strict;
use warnings;
use File::Basename;

use Genome;

class Genome::InstrumentData::AlignmentResult::Bowtie {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => { value => 'Bowtie', is_param=>1 },
    ],
    has_transient_optional => [
         _bwa_sam_cmd => { is=>'Text' }
    ]
};

sub required_arch_os { 'x86_64' }

# fill me in here with what compute resources you need.
sub required_rusage { 
    my $queue = $ENV{GENOME_LSF_QUEUE_ALIGNMENT_DEFAULT};
    $queue = 'alignment-pd' if (Genome::Config->should_use_alignment_pd);
    "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>90000 && mem>10000] span[hosts=1] rusage[tmp=90000, mem=10000]' -M 10000000 -n 4 -m $queue -q $queue";
}


sub _run_aligner {
    my $self = shift;
    my @input_pathnames = @_;

    my $tmp_dir = $self->temp_scratch_directory;

    # get refseq info
    my $reference_build = $self->reference_build;
    # This will search for something along the lines of all_sequences.bowtie. 
    # all_sequences.bowtie should be a symlink to all_sequences.fa
    my $reference_bowtie_index_path = $reference_build->full_consensus_path('bowtie');

    unless (-s $reference_bowtie_index_path) {
        $self->error_message("Bowtie index file not found or is empty at $reference_bowtie_index_path or " . $reference_build->data_directory);
        die;
    }

    my $aligner_params = $self->aligner_params || '';

    my $path_to_bowtie = Genome::Model::Tools::Bowtie->path_for_bowtie_version($self->aligner_version);

    my $output_file = $self->temp_scratch_directory . "/all_sequences.sam";
    my $log_file = $self->temp_staging_directory . "/aligner.log";

    my $temp_aligned_sequences_file = $self->temp_scratch_directory . "/temp_aligned_sequences.sam";
    my $temp_unaligned_fq_file = $self->temp_scratch_directory . "/temp_unaligned_sequences.fq";
	my $temp_unaligned_sam_file = $self->temp_scratch_directory . "/temp_unaligned_sequences.sam";

    if ( @input_pathnames == 1 ) {
        my $cmdline = "$path_to_bowtie $aligner_params --sam-nohead --un $temp_unaligned_fq_file $reference_bowtie_index_path $input_pathnames[0] --sam $temp_aligned_sequences_file >> $log_file && cat $temp_aligned_sequences_file >> $output_file";

        Genome::Sys->shellcmd(
            cmd                         => $cmdline,
            input_files                 => [$reference_bowtie_index_path, $input_pathnames[0]],
            output_files                => [$output_file],
            skip_if_output_is_present   => 0,
        );
    }
    elsif ( @input_pathnames == 2 ) {
        # pass insert size into bowtie
        my $insert_sd   = $self->instrument_data->resolve_sd_insert_size;
        my $insert_size = $self->instrument_data->resolve_median_insert_size;

        if ($insert_size && $insert_sd) {
            $aligner_params .= ' --minins '. ($insert_size - $insert_sd) .' --maxins '. ($insert_size + $insert_sd);
        } 
	
        my $cmdline = "$path_to_bowtie $aligner_params --sam-nohead --un $temp_unaligned_fq_file $reference_bowtie_index_path -1 $input_pathnames[0] -2 $input_pathnames[1] --sam $temp_aligned_sequences_file >>$log_file && cat $temp_aligned_sequences_file >> $output_file";
        # $temp_unaligned_sequences_file still uses original file format.
        
        Genome::Sys->shellcmd(
            cmd                         => $cmdline,
            input_files                 => [$reference_bowtie_index_path, $input_pathnames[0], $input_pathnames[1]],
            output_files                => [$output_file],
            skip_if_output_is_present   => 0,
        );
		
		# Convert unaligned reads from Fastq to SAM, place in $temp_unaligned_sam_file:
		# Note: with 2 input fq files, bowtie outputs 2 also. In this case to .../temp_unaligned_sequences_1.fq and .../temp_unaligned_sequences_2.fq
		# We'll combine these two files so we can convert them to SAM
		print "Done with primary alignment. Appending unaligned reads. \n";
		
		my $temp_fq_1 = $temp_unaligned_fq_file;
		$temp_fq_1 =~ s/temp_unaligned_sequences/temp_unaligned_sequences_1/g;
		my $temp_fq_2 = $temp_unaligned_fq_file;
		$temp_fq_2 =~ s/temp_unaligned_sequences/temp_unaligned_sequences_2/g;

		my $cat_one_to_unaligned_cmd = "cat $temp_fq_1 >> $temp_unaligned_fq_file";
		Genome::Sys->shellcmd(
			cmd                         => $cat_one_to_unaligned_cmd,
			input_files                 => [$temp_fq_1],
			output_files                => [$temp_unaligned_fq_file],
			skip_if_output_is_present   => 0,
		);

		my $cat_two_to_unaligned_cmd = "cat $temp_fq_2 >> $temp_unaligned_fq_file";
		Genome::Sys->shellcmd(
			cmd                         => $cat_two_to_unaligned_cmd,
			input_files                 => [$temp_fq_2],
			output_files                => [$temp_unaligned_fq_file],
			skip_if_output_is_present   => 0,
		);
    }
    else {
        $self->error_message("Input pathnames shouldn't have more than 2...: " . Data::Dumper::Dumper(\@input_pathnames) );
        die $self->error_message;
    }
    
    Genome::Model::Tools::Sam::FastqToSam->execute(
			fastq_file => $temp_unaligned_fq_file,
			sam_file   => $temp_unaligned_sam_file,
	);

    my $cat_unaligned_to_output_cmd = "cat $temp_unaligned_sam_file >> $output_file";
    Genome::Sys->shellcmd(
            cmd                        => $cat_unaligned_to_output_cmd,
            input_files                => [$temp_unaligned_sam_file],
            output_files               => [$output_file],
            skip_if_output_is_present  => 0,
    );   

    return 1;
}

sub aligner_params_for_sam_header {
    my $self = shift;

    return "bowtie " . $self->aligner_params;
    # for bwa this looks like "bwa aln -t4; bwa samse 12345'
}

sub fillmd_for_sam { return 0; }

sub _check_read_count { return 1; }

sub prepare_reference_sequence_index {
    my $class = shift;
    my $refindex = shift;

    my $staging_dir = $refindex->temp_staging_directory;
    my $staged_fasta_file = sprintf("%s/all_sequences.fa", $staging_dir);

    $class->status_message("Doing bowtie indexing.");

    my $bowtie_file_stem = "";
    my $bowtie_extension = "";
    if($refindex->aligner_version =~ /^2/){
      $bowtie_extension = "bt2";
      $bowtie_file_stem = File::Spec->catfile($staging_dir, 'all_sequences.fa');
    }else{
      $bowtie_extension = "ebwt";
      $bowtie_file_stem = File::Spec->catfile($staging_dir, 'all_sequences.bowtie');
    }

    my $bowtie_path = Genome::Model::Tools::Bowtie->path_for_bowtie_version($refindex->aligner_version, "build");

    my $bowtie_cmd = sprintf('%s %s %s', $bowtie_path, $staged_fasta_file, $bowtie_file_stem);
    my $rv = Genome::Sys->shellcmd(
        cmd => $bowtie_cmd,
        input_files => [$staged_fasta_file],
        output_files => ["$bowtie_file_stem.1.$bowtie_extension","$bowtie_file_stem.2.$bowtie_extension","$bowtie_file_stem.3.$bowtie_extension","$bowtie_file_stem.4.$bowtie_extension","$bowtie_file_stem.rev.1.$bowtie_extension","$bowtie_file_stem.rev.2.$bowtie_extension"],    #hardcoding expected names
    );

    unless($rv) {
        $class->error_message('bowtie-build failed.');
        return;
    }

    return 1;
}

