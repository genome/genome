package Genome::Model::RnaSeq::Command::ErccSpikeInQc;

use strict;
use warnings;

use Genome;
use Cwd;

class Genome::Model::RnaSeq::Command::ErccSpikeInQc {
    is => 'Command::V2',
    has => [
        build => {
            is => 'Genome::Model::Build::RnaSeq',
            doc => 'The RNA-seq build to evaluate ERCC spike-in QC',
            shell_args_position => 1,
        },
        picard_version => {
            is => 'Text',
            doc => 'The version of Picad to use.',
            example_values => ['1.123'],
        },
        picard_max_records_in_ram => {
            is => 'Integer',
            doc => 'When writing SAM files that need to be sorted, this will specify the number of records stored
    in RAM before spilling to disk. Increasing this number reduces the number of file handles
    needed to sort a SAM file, and increases the amount of RAM needed..',
            example_values => ['4000000'],
        },
        picard_maximum_permgen_memory => {
            is => 'Integer',
            doc => 'the maximum memory (Mbytes) to use for the "permanent generation" of the Java heap (e.g., for interned Strings)',
            example_values => ['256'],
        },
        samtools_version => {
            is => 'Text',
            doc => 'The version of Samtools to use.',
            example_values => ['0.1.19'],
        },
        maximum_memory => {
            is => 'Integer',
            doc => 'The amount of RAM in GB to be used by Samtools and Picard.',
            example_values => ['14'],
        },
        bwa_threads => {
            is => 'Integer',
            doc => 'The number of threads to utilize during BWA alignment.',
            example_values => ['4'],
        },
        bwa_version => {
            is => 'Text',
            doc => 'The version of bwa to use.',
            example_values => ['0.6.2'],
        },
        ercc_mix => {
            is => 'Integer',
            doc => 'The expected ERCC spike-in mix.',
            valid_values => [ '1', '2'],
        },
        ercc_control_analysis_file => {
            is => 'Text',
            doc => 'The control analysis file provided by Agilent for the ERCC spike-in.',
            example_values => ['/gscmnt/gc13001/info/model_data/jwalker_scratch/ERCC/ERCC_Controls_Analysis.txt'],
        },
        ercc_fasta_file => {
            is => 'Text',
            doc => 'The FASTA file provided by Agilent for the ERCC spike-in.',
            example_values => ['/gscmnt/gc13001/info/model_data/jwalker_scratch/ERCC/ERCC92.fa'],
        },
        output_directory => {
            is => 'Text',
            doc => 'The output directory to write summary file and PDF plots.',
        },
    ],
    has_optional => [
        _samtools_path => {},
        _bwa_path => {},
    ],
};

sub help_detail {
    return "Compare the read counts of ERCC spike-in mixes with a known concentration.";
}

sub execute {
    my $self = shift;

    my $samtools_path = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);
    $self->_samtools_path($samtools_path);

    my $bwa_path = Genome::Model::Tools::Bwa->path_for_bwa_version($self->bwa_version);
    $self->_bwa_path($bwa_path);

    my $remapped_bam = $self->_remap_unmapped_reads();
    my $idxstats_hash_ref = $self->_generate_idxstats_hash_ref($remapped_bam);
    my $summary_file = $self->_generate_summary_file($idxstats_hash_ref);
    $self->_generate_r_plots($summary_file);

    return 1;
}

sub _remap_unmapped_reads {
    my $self = shift;
    
    my $alignment_result = $self->build->merged_alignment_result;
    my $unmapped_bam = $self->_extract_unmapped_reads_from_alignment_result($alignment_result);
    my $fasta_index = $self->_bwa_index_reference();
    my $sai_1 = $self->_bwa_aln_read_end($fasta_index,$unmapped_bam,'1');
    my $sai_2 = $self->_bwa_aln_read_end($fasta_index,$unmapped_bam,'2');
    my $remapped_bam = $self->_bwa_sampe($fasta_index,$unmapped_bam,$sai_1,$sai_2);

    $self->_index_bam($remapped_bam);
    
    return $remapped_bam;
}
    
sub _extract_unmapped_reads_from_alignment_result {
    my $self = shift;
    my $alignment_result = shift;

    my $queryname_sorted_bam = Genome::Sys->create_temp_file_path();
    my $queryname_sort_cmd = Genome::Model::Tools::Picard::SortSam->execute(
        input_file => $alignment_result->bam_path,
        output_file => $queryname_sorted_bam,
        maximum_memory => $self->maximum_memory,
        maximum_permgen_memory => $self->picard_maximum_permgen_memory,
        max_records_in_ram => $self->picard_max_records_in_ram,
        use_version => $self->picard_version,
        sort_order => 'queryname',
    );
    unless ($queryname_sort_cmd && $queryname_sort_cmd->result() ) {
        die($self->error_message('Failed to execute Picard SortSam!'));
    }
    
    my $fix_mate_bam = Genome::Sys->create_temp_file_path();
    my $fix_mate_cmd = Genome::Model::Tools::Picard::FixMateInformation->execute(
        input_file => $queryname_sorted_bam,
        output_file => $fix_mate_bam,
        max_records_in_ram => $self->picard_max_records_in_ram,
        maximum_memory => $self->maximum_memory,
        maximum_permgen_memory => $self->picard_maximum_permgen_memory,
        use_version => $self->picard_version,
    );
    unless ($fix_mate_cmd && $fix_mate_cmd->result() ) {
        die($self->error_message('Failed to execute Picard FixMateInformation!'));
    }
    
    my $unmapped_bam_basename = Genome::Sys->create_temp_file_path();
    my $unmapped_bam = $unmapped_bam_basename .'.bam';
    my $samtools_cmd = $self->_samtools_path .' view -u -f 12 '. $fix_mate_bam .' | '. $self->_samtools_path .' sort -n -m '. $self->maximum_memory .'000000000 - '. $unmapped_bam_basename;
    Genome::Sys->shellcmd(
        cmd => $samtools_cmd,
        input_files => [$fix_mate_bam],
        output_files => [$unmapped_bam],
    );
    
    return $unmapped_bam;
}

sub _bwa_index_reference {
    my $self = shift;
    
    my $fasta_index = Genome::Sys->create_temp_file_path();
    Genome::Sys->create_symlink($self->ercc_fasta_file,$fasta_index);

    my $bwa_index_cmd = $self->_bwa_path .' index '. $fasta_index;
    Genome::Sys->shellcmd(
        cmd => $bwa_index_cmd,
        input_files => [$fasta_index],
    );
    return $fasta_index;
}

sub _bwa_aln_read_end {
    my $self = shift;
    my $fasta_index = shift;
    my $unmapped_bam_path = shift;
    my $read_end = shift;

    my $sai_path = Genome::Sys->create_temp_file_path();
    my $aln_cmd = $self->_bwa_path .' aln -t '. $self->bwa_threads .' -b -'. $read_end .' -f '. $sai_path .' '. $fasta_index .' '. $unmapped_bam_path;
    Genome::Sys->shellcmd(
        cmd => $aln_cmd,
        input_files => [$fasta_index,$unmapped_bam_path],
        output_files => [$sai_path],
    );
    return $sai_path;
}

sub _bwa_sampe {
    my $self = shift;
    my $fasta_index = shift;
    my $unmapped_bam = shift;
    my $sai_1 = shift;
    my $sai_2 = shift;
    
    my $remapped_bam_basename = Genome::Sys->create_temp_file_path();
    my $remapped_bam_path = $remapped_bam_basename .'.bam';
    
    my $sampe_cmd = $self->_bwa_path .' sampe '. $fasta_index .' '. $sai_1 .' '. $sai_2 .' '. $unmapped_bam .' '. $unmapped_bam .' | '. $self->_samtools_path .' view -u -S - | '. $self->_samtools_path .' sort -m '. $self->maximum_memory .'000000000 - '. $remapped_bam_basename;
    Genome::Sys->shellcmd(
        cmd => $sampe_cmd,
        input_files => [$fasta_index,$sai_1,$sai_2,$unmapped_bam],
        output_files => [$remapped_bam_path],
    );
    
    return $remapped_bam_path;
}

sub _index_bam {
    my $self = shift;
    my $bam_file = shift;
    
    my $samtools_index_cmd = Genome::Model::Tools::Sam::IndexBam->execute(
        bam_file => $bam_file,
        use_version => $self->samtools_version,
    );
    unless ($samtools_index_cmd && $samtools_index_cmd->result) {
        die('Failed to run samtools index for BAM file: '. $bam_file);
    }
    return 1;
}

sub _generate_idxstats_hash_ref {
    my $self = shift;
    my $bam_file = shift;

    my $samtools_idxstats_path = Genome::Sys->create_temp_file_path();
    my $samtools_idxstats_cmd = Genome::Model::Tools::Sam::Idxstats->execute(
        bam_file => $bam_file,
        output_file => $samtools_idxstats_path,
        use_version => $self->samtools_version,
    );
    unless ($samtools_idxstats_cmd && $samtools_idxstats_cmd->result) {
        die('Failed to run samtools idxstats for BAM file: '. $bam_file);
    }
    
    my $idxstats_hash_ref = $samtools_idxstats_cmd->parse_file_into_hashref($samtools_idxstats_path);
    return $idxstats_hash_ref;
}

sub _generate_r_plots {
    my $self = shift;
    my $summary_file = shift;

    my $r_script_path = $self->__meta__->module_path;
    $r_script_path =~ s/\.pm/\.R/;
    my $cwd = getcwd();
    my $cmd = $r_script_path .' --filename '. $summary_file;
    chdir($self->output_directory);
    Genome::Sys->shellcmd(
        cmd => $cmd,
    );
    chdir($cwd);
    return 1;
}

sub _generate_summary_file {
    my $self = shift;
    my $idxstats_hash_ref = shift;
    
    my @output_headers = (
        'Re-sort ID',
        'ERCC ID',
        'subgroup',
        'ERCC Mix',
        'concentration (attomoles/ul)',
        'count',
    );
    
    my $summary_file = $self->output_directory .'/count_summary.tsv';;
    my $summary_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $summary_file,
        headers => \@output_headers,
        separator => "\t",
    );

    my $concentration_key = 'concentration in Mix '. $self->ercc_mix .' (attomoles/ul)';

    my $ercc_control = $self->_load_ercc_control;
    for my $contig (keys %{$idxstats_hash_ref}) {
        my $ercc_data = $ercc_control->{$contig};
        unless ($ercc_data) {
            die($self->error_message('Missing ERCC control analysis data for reference contig : '. $contig));
        }
        my %data = (
            'Re-sort ID' => $ercc_data->{'Re-sort ID'},
            'ERCC ID' => $ercc_data->{'ERCC ID'},
            'subgroup' => $ercc_data->{'subgroup'},
            'ERCC Mix' => $self->ercc_mix,
            'concentration (attomoles/ul)' => $ercc_data->{$concentration_key},
            'count' => $idxstats_hash_ref->{$contig}{map_read_ct},
        );
        $summary_writer->write_one(\%data);
    }
    $summary_writer->output->close;
    
    return $summary_file;
}


sub _load_ercc_control {
    my $self = shift;

    my %ercc_control;
    my $ercc_control_analysis_reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $self->ercc_control_analysis_file,
        separator => "\t",
    );
    while (my $data = $ercc_control_analysis_reader->next) {
        $ercc_control{$data->{'ERCC ID'}} = $data;
    }
    $ercc_control_analysis_reader->input->close;
    
    return \%ercc_control;
}
