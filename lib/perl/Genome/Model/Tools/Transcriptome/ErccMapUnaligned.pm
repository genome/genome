package Genome::Model::Tools::Transcriptome::ErccMapUnaligned;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Transcriptome::ErccMapUnaligned {
    is => 'Command::V2',
    has_input => [
        bam_file => {
            is => 'Text',
            doc => 'An aligned BAM file to a reference genome without the ERCC transcripts spiked in.',
        },
        ercc_fasta_file => {
            is => 'Text',
            doc => 'The FASTA format sequence for all ERCC transcripts.',
            example_values => ['/gscmnt/gc13001/info/model_data/jwalker_scratch/ERCC/ERCC.fa'],
        },
        ercc_spike_in_file => {
            is => 'Text',
            doc => 'The control analysis file provided by Agilent for the ERCC spike-in.',
            example_values => ['/gscmnt/gc13001/info/model_data/jwalker_scratch/ERCC/ERCC_Controls_Analysis.txt'],
        },
        ercc_spike_in_mix => {
            is => 'Integer',
            doc => 'The expected ERCC spike-in mix.',
            valid_values => [ '1', '2'],
        },
    ],
    has_output => [
        pdf_file => {
            is => 'Text',
            doc => 'The output PDF with histograms and linearity plot.',
            default_value => 'output.pdf',
            is_optional => '1',
        }
    ],
};

sub help_detail {
    return <<EOS
Compare the abundance of an ERCC spike-in with a known concentration.
EOS
}

sub execute {
    my $self = shift;
    my @output_headers = (
        'Re-sort ID',
        'ERCC ID',
        'subgroup',
        'ERCC Mix',
        'concentration (attomoles/ul)',
        'label',
        'count',
    );
    
    my $samtools_version = '0.1.19';
    my $samtools_max_mem = 14000000000;
    my $samtools_path = Genome::Model::Tools::Sam->path_for_samtools_version($samtools_version);
    
    # BWA config
    my $bwa_version = '0.6.2';
    my $bwa_path = Genome::Model::Tools::Bwa->path_for_bwa_version($bwa_version);
    my $ercc_r = $self->__meta__->module_path;
    $ercc_r =~ s/\.pm/\.R/;
    
    # Load ERCC Control Info
    my %ercc_control;
    my $ercc_spike_in_reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $self->ercc_spike_in_file,
        separator => "\t",
    );
    while (my $data = $ercc_spike_in_reader->next) {
        $ercc_control{$data->{'ERCC ID'}} = $data;
    }
    $ercc_spike_in_reader->input->close;

    my $unmapped_bam_basename = Genome::Sys->create_temp_file_path();
    my $unmapped_bam_path = $unmapped_bam_basename .'.bam';

    my $samtools_view_cmd = $samtools_path .' view -u -f 12 '. $self->bam_file .' | '. $samtools_path .' sort -n -m '. $samtools_max_mem .' - '. $unmapped_bam_basename;
    Genome::Sys->shellcmd(
        cmd => $samtools_view_cmd,
        input_files => [$self->bam_file],
        output_files => [$unmapped_bam_path],
    );

    # Create BWA Index
    my $fasta_index_path = Genome::Sys->create_temp_file_path();
    Genome::Sys->create_symlink($self->ercc_fasta_file,$fasta_index_path);
    my $bwa_index_cmd = $bwa_path .' index '. $fasta_index_path;
    Genome::Sys->shellcmd(
        cmd => $bwa_index_cmd,
        input_files => [$fasta_index_path],
    );

    my $sai_read1_path = Genome::Sys->create_temp_file_path();
    my $sai_read2_path = Genome::Sys->create_temp_file_path();

    my $aln_read1_cmd = $bwa_path .' aln -b -1 -f '. $sai_read1_path .' '. $fasta_index_path .' '. $unmapped_bam_path;
    Genome::Sys->shellcmd(
        cmd => $aln_read1_cmd,
        input_files => [$fasta_index_path,$unmapped_bam_path],
        output_files => [$sai_read1_path],
    );
    
    my $aln_read2_cmd = $bwa_path .' aln -b -2 -f '. $sai_read2_path .' '. $fasta_index_path .' '. $unmapped_bam_path;
    Genome::Sys->shellcmd(
        cmd => $aln_read2_cmd,
        input_files => [$fasta_index_path,$unmapped_bam_path],
        output_files => [$sai_read2_path],
    );

    # BWA sampe
    my $remapped_bam_basename = Genome::Sys->create_temp_file_path();
    my $remapped_bam_path = $remapped_bam_basename .'.bam';

    my $sampe_cmd = $bwa_path .' sampe '. $fasta_index_path .' '. $sai_read1_path .' '. $sai_read2_path .' '. $unmapped_bam_path .' '. $unmapped_bam_path .' | '. $samtools_path .' view -u -S - | '. $samtools_path .' sort -m '. $samtools_max_mem .' - '. $remapped_bam_basename;
    Genome::Sys->shellcmd(
        cmd => $sampe_cmd,
        input_files => [$fasta_index_path,$sai_read1_path,$sai_read2_path,$unmapped_bam_path],
        output_files => [$remapped_bam_path],
    );

    # Index BAM
    my $samtools_index_cmd = Genome::Model::Tools::Sam::IndexBam->execute(
        bam_file => $remapped_bam_path,
        use_version => $samtools_version,
    );
    unless ($samtools_index_cmd && $samtools_index_cmd->result) {
        die('Failed to run samtools index for BAM file: '. $remapped_bam_path);
    }

    my $samtools_idxstats_path = Genome::Sys->create_temp_file_path();
    my $samtools_idxstats_cmd = Genome::Model::Tools::Sam::Idxstats->execute(
        bam_file => $remapped_bam_path,
        output_file => $samtools_idxstats_path,
        use_version => $samtools_version,
    );
    unless ($samtools_idxstats_cmd && $samtools_idxstats_cmd->result) {
        die('Failed to run samtools idxstats for BAM file: '. $remapped_bam_path);
    }

    # Parse Idxstats
    my $idxstats_hash_ref = $samtools_idxstats_cmd->parse_file_into_hashref($samtools_idxstats_path);

    my $r_input_file = Genome::Sys->create_temp_file_path();
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $r_input_file,
        headers => \@output_headers,
        separator => "\t",
    );

    for my $chr (keys %{$idxstats_hash_ref}) {
        my $ercc_data = $ercc_control{$chr};
        unless ($ercc_data) {
            die('Missing chromosome: '. $chr);
        }
        my $concentration = $ercc_data->{'concentration in Mix 1 (attomoles/ul)'};
        if ($self->ercc_spike_in_mix == 2) {
            $concentration = $ercc_data->{'concentration in Mix 2 (attomoles/ul)'};
        }
        my %data = (
            'Re-sort ID' => $ercc_data->{'Re-sort ID'},
            'ERCC ID' => $ercc_data->{'ERCC ID'},
            'subgroup' => $ercc_data->{'subgroup'},
            'ERCC Mix' => $self->ercc_spike_in_mix,
            'concentration (attomoles/ul)' => $concentration,
            'label' => 'na',
            'count' => $idxstats_hash_ref->{$chr}{map_read_ct},
        );
        $writer->write_one(\%data);
    }
    $writer->output->close;

    my $cmd = $ercc_r .' --filename '. $r_input_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
    );

    
    return 1;
}


