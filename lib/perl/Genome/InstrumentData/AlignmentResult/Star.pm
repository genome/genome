package Genome::InstrumentData::AlignmentResult::Star;

use strict;
use warnings;

use Genome;

my $CPUS = 8;

class Genome::InstrumentData::AlignmentResult::Star {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => {
            value    => 'Star',
            is_param => 1
        },
    ],
};

sub required_arch_os { 'x86_64' }


sub required_rusage {
    my $mem_mb = 32000;
    my $mem_kb = $mem_mb*1024;
    
    my $select  = "select[ncpus >= $CPUS && mem >= $mem_mb] span[hosts=1]";
    my $rusage  = "rusage[mem=$mem_mb]";
    my $options = "-M $mem_kb -n $CPUS";
    my $required_usage = "-R '$select $rusage' $options";
    
    return $required_usage;
}


sub _run_aligner {
    my ($self, @input_pathnames) = @_;

    my $aligner_version   = $self->aligner_version;
    my $aligner_params    = $self->aligner_params;
    my $reference_build   = $self->reference_build;
    my $annotation_build  = $self->annotation_build;

    unless ($annotation_build) {
        die $self->error_message("annotation_build is not valid");
    }

    my $scratch_directory = $self->temp_scratch_directory;
    my $staging_directory = $self->temp_staging_directory;

    my %aligner_params = decomposed_aligner_params($aligner_params);
    my $index_params   = $aligner_params{index_params};
    my $align_params   = $aligner_params{align_params};

    if (Genome::DataSource::GMSchema->has_default_handle) {
        $self->debug_message("Disconnecting GMSchema default handle.");
        Genome::DataSource::GMSchema->disconnect_default_dbh();
    }

    #prepare reference sequence genome index file under directory
    my $annotate_index = Genome::Model::Build::ReferenceSequence::AnnotationIndex->get(
        aligner_name     => $self->aligner_name,
        aligner_version  => $aligner_version,
        aligner_params   => $aligner_params,
        reference_build  => $reference_build,
        annotation_build => $annotation_build,
    );

    unless ($annotate_index) {
        die $self->error_message('Failed to get annotation index for star run');
    }

    my $index_dir = $annotate_index->output_dir;
    my $star_path = Genome::Model::Tools::Star->path_for_star_version($self->aligner_version);
    my $input_fastq_files = join ' ', @input_pathnames;
    
    unless ($align_params and $align_params =~ /runThreadN/) {
        $align_params .= '--runThreadN '. $CPUS;
    }

    my $cmd = "$star_path --genomeDir $index_dir --readFilesIn $input_fastq_files --outFileNamePrefix $scratch_directory/ --outSAMunmapped Within $align_params";
    my $out_sam_file = "$scratch_directory/Aligned.out.sam";

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files  => \@input_pathnames,
        output_files => [$out_sam_file],
    );

    ##Below fillmd and addreadgrouptag should be handled by AlignmentResult base class
    my $sam_path = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);
    my $ref_seq  = $reference_build->full_consensus_path('fa');

    ##To make fillmd run faster sorted bam file is needed, running fillmd with -S sam file is VERY VERY SLOW
    my $tmp_bam = "$scratch_directory/Aligned.out.bam";
    $cmd = "$sam_path view -b -S -o $tmp_bam $out_sam_file";
    
    my $rv = Genome::Sys->shellcmd(
        cmd                          => $cmd,
        input_files                  => [$out_sam_file],
        output_files                 => [$tmp_bam],
        skip_if_output_is_present    => 0,
    );
    die $self->error_message("Fail to run: $cmd") unless $rv == 1;

    my $sort_bam_prefix = "$scratch_directory/Aligned.out.sort";
    my $sort_bam = $sort_bam_prefix . '.bam';
    $cmd = "$sam_path sort -m 402653184 $tmp_bam $sort_bam_prefix";

    $rv = Genome::Sys->shellcmd(
        cmd                          => $cmd,
        input_files                  => [$tmp_bam],
        output_files                 => [$sort_bam],
        skip_if_output_is_present    => 0,
    );
    die $self->error_message("Fail to run: $cmd") unless $rv == 1;

    my $fillmd_sam_file = "$scratch_directory/Aligned.out.fillmd.sam";
    $cmd = "$sam_path fillmd $sort_bam $ref_seq 1> $fillmd_sam_file 2>/dev/null";

    $rv  = Genome::Sys->shellcmd(
        cmd                          => $cmd,
        input_files                  => [$sort_bam, $ref_seq],
        output_files                 => [$fillmd_sam_file],
        skip_if_output_is_present    => 0,
    );
    die $self->error_message("Fail to run: $cmd") unless $rv == 1;

    $self->debug_message('FillMD completed');

    my $sam_fh = Genome::Sys->open_file_for_reading($fillmd_sam_file) or die "failed to open sam file: $fillmd_sam_file for reading\n";

    my $add_rg_cmd = Genome::Model::Tools::Sam::AddReadGroupTag->create(
        input_filehandle  => $sam_fh,
        output_filehandle => $self->_sam_output_fh,
        read_group_tag    => $self->read_and_platform_group_tag_id,
        pass_sam_headers  => 0,
    );

    unless($add_rg_cmd->execute){
        $self->error_message("Adding read group to sam file failed!");
        die $self->error_message;
    }

    $self->debug_message('Read group add completed');

    #promote other misc result files - converted sam will be handled downstream
    map{Genome::Sys->move_file("$scratch_directory/$_", "$staging_directory/$_")}qw(SJ.out.tab Log.progress.out Log.out Log.final.out);

    return 1;
}


sub prepare_annotation_index {
    my ($self, $annot_index) = @_;

    my $ref_build   = $annot_index->reference_build;
    my $annot_build = $annot_index->annotation_build;
    my $out_dir     = $annot_index->temp_staging_directory;

    my %aliger_params = __PACKAGE__->decomposed_aligner_params($annot_index->aligner_params);
    my $params = $aliger_params{index_params};

    my $ref_seq_fasta = $ref_build->full_consensus_path('fa');
    my $gtf_file      = $annot_build->annotation_file('gtf', $ref_build->id);

    my $star_path = Genome::Model::Tools::Star->path_for_star_version($annot_index->aligner_version);

    unless ($params and $params =~ /runThreadN/) {
        $params .= '--runThreadN '. $CPUS;
    }

    unless ($params and $params =~ /sjdbOverhang/) {
        $params .= '--sjdbOverhang 100';
    }

    my $cmd = "$star_path --runMode genomeGenerate --genomeDir $out_dir --genomeFastaFiles $ref_seq_fasta --sjdbGTFfile $gtf_file $params";
    
    my @output_files = map{$out_dir.'/'.$_}qw(SA SAindex);
    my $rv = Genome::Sys->shellcmd(
        cmd          => $cmd,
        input_files  => [$ref_seq_fasta, $gtf_file],
        output_files => [@output_files]
    );

    unless ($rv) {
        $self->error_message('star indexing failed');
        return;
    }

    return 1;
}

sub _check_read_count {
    my $self = shift;
    
    my $fq_rd_ct = $self->_fastq_read_count;
    my $sam_path = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);

    my $cmd = "$sam_path view -F 256 -c " . $self->temp_staging_directory . "/all_sequences.bam";
    my $bam_read_count = `$cmd`;
    my $check = "Read count from bam: $bam_read_count and fastq: $fq_rd_ct";

    unless ($fq_rd_ct == $bam_read_count) {
        $self->error_message("$check does not match.");
        return;
    }
    $self->debug_message("$check matches.");
    return 1;
}


sub decomposed_aligner_params {
    my $params = pop;
    $params ||= '::';

    my @params = split /\:/, $params;

    return (
        index_params => $params[0] || '', 
        align_params => $params[1] || '',
    );
}


sub alignment_bam_file_paths {
    my $self = shift;
    return ($self->output_dir . "/all_sequences.bam");
}

sub aligner_params_for_sam_header {
    my $self = shift;
    return "star " . $self->aligner_params;
}

sub aligner_params_required_for_index {
    return 1;
}

sub fillmd_for_sam {
    return 1;
}

sub requires_read_group_addition {
    return 0;
}

sub supports_streaming_to_bam {
    return 1;
}

sub accepts_bam_input {
    return 0;
}


1;
