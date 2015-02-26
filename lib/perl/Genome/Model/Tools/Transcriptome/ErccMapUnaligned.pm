package Genome::Model::Tools::Transcriptome::ErccMapUnaligned;

use strict;
use warnings;

use Genome;
use Path::Class;

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
            example_values => ['/gscmnt/gc13001/info/model_data/jwalker_scratch/ERCC/ERCC92.fa'],
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
        samtools_version => {
            is => 'Text',
            doc => 'The version of samtools to run analysis.',
            example_values => ['0.1.19'],
            default_value => '0.1.19',
            is_optional => '1',
        },
        samtools_max_mem => {
            is => 'Text',
            doc => 'The max memory required by samtools.',
            example_values => ['14000000000'],
            default_value => '14000000000',
            is_optional => '1',
        },
        bwa_version => {
            is => 'Text',
            doc => 'The version of bwa to run analysis.',
            example_values => ['0.6.2'],
            default_value => '0.6.2',
            is_optional => '1',
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
    
    my $ercc_bwa_index = $self->create_ercc_bwa_index();
    my $unmapped_bam = $self->generate_unmapped_bam();

    my $sai_read1 = $self->align_unmapped_reads(
        bam => $unmapped_bam,
        index => $ercc_bwa_index,
        read_pair => 1
    );

    my $sai_read2 = $self->align_unmapped_reads(
        bam => $unmapped_bam,
        index => $ercc_bwa_index,
        read_pair => 2
    );

    my $remapped_bam = $self->generate_remapped_bam(
        index => $ercc_bwa_index,
        sai_read1 => $sai_read1,
        sai_read2 => $sai_read2,
        unmapped_bam => $unmapped_bam,
    );

    $self->index_bam($remapped_bam);
    my $idxstats = $self->generate_idxstats($remapped_bam);
    my $tsv = $self->generate_tsvfile($idxstats);

    $self->run_analysis_script($tsv);

    return 1;
}

sub samtools {
    my $self = shift;
    my $version = $self->samtools_version;
    my $samtools_path =
      Genome::Model::Tools::Sam->path_for_samtools_version($version);
    return $samtools_path;
}

sub bwa {
    my $self = shift;
    my $version = $self->bwa_version;
    my $bwa_path =
      Genome::Model::Tools::Bwa->path_for_bwa_version($version);
    return $bwa_path;
}

sub ERCC_analysis_script {
    my $self = shift;
    my $ercc_r = $self->__meta__->module_path;
    $ercc_r =~ s/\.pm/\.R/;
    return $ercc_r;
}

sub run_analysis_script {
    my ($self, $tsv) = @_;
    my $ercc_r = $self->ERCC_analysis_script;
    my $cmd = "$ercc_r --data $tsv";
    Genome::Sys->shellcmd(cmd => $cmd);
}

sub load_ercc_control_info {
    my $self = shift;
    my %ercc_control;
    my $ercc_spike_in_reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $self->ercc_spike_in_file,
        separator => "\t",
    );
    while (my $data = $ercc_spike_in_reader->next) {
        $ercc_control{$data->{'ERCC ID'}} = $data;
    }
    $ercc_spike_in_reader->input->close;

    return %ercc_control
}

sub generate_unmapped_bam {
    my $self = shift;
    my $unmapped_bam_basename = Genome::Sys->create_temp_file_path();
    my $unmapped_bam_path = $unmapped_bam_basename .'.bam';

    my $max_mem = $self->samtools_max_mem;
    my $samtools_view_cmd = join(' ', 
        $self->samtools, 'view', '-u -f 12', $self->bam_file ,
        '|',
        $self->samtools, 'sort', "-n -m $max_mem", '-',
        $unmapped_bam_basename
    );

    Genome::Sys->shellcmd(
        cmd => $samtools_view_cmd,
        input_files => [$self->bam_file],
        output_files => [$unmapped_bam_path],
    );

    return Path::Class::File->new($unmapped_bam_path);
}

sub create_ercc_bwa_index {
    my $self = shift;

    my $fasta_index_path = Genome::Sys->create_temp_file_path();

    Genome::Sys->create_symlink($self->ercc_fasta_file,$fasta_index_path);
    my $cmd = join(' ', $self->bwa, 'index', $fasta_index_path);

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$fasta_index_path],
    );

    return Path::Class::File->new($fasta_index_path);
}

sub align_unmapped_reads {
    my $self = shift;
    my %args = @_;
    my ($bam, $bwa_index, $read_pair) = @args{'bam', 'index', 'read_pair'};

    my $sai_path = Genome::Sys->create_temp_file_path();

    my $cmd = join(' ',
        $self->bwa,
        "aln",
        "-b -${read_pair}",
        "-f ${sai_path}",
        $bwa_index->stringify,
        $bam->stringify
    );

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => ["$bwa_index", "$bam"],
        output_files => [$sai_path],
    );

    return Path::Class::File->new($sai_path);
}

sub generate_remapped_bam {
    my $self = shift;
    my %args = @_;
    my ($index, $sai_read1_path, $sai_read2_path, $unmapped_bam) =
      @args{'index', 'sai_read1', 'sai_read2', 'unmapped_bam'};

    my $remapped_bam_basename = Genome::Sys->create_temp_file_path();
    my $remapped_bam_path = $remapped_bam_basename .'.bam';

    my $max_mem = $self->samtools_max_mem;
    my $cmd = join(' ',
        $self->bwa, 'sampe',
        "$index",
        "$sai_read1_path", "$sai_read2_path",
        "$unmapped_bam", "$unmapped_bam",
        '|',
        $self->samtools, 'view', '-u -S', '-',
        '|',
        $self->samtools, 'sort', "-m $max_mem", '-',
        $remapped_bam_basename
    );

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [
            "$index",
            "$sai_read1_path",
            "$sai_read2_path",
            "$unmapped_bam"
        ],
        output_files => [$remapped_bam_path],
    );

    return Path::Class::File->new($remapped_bam_path);
}

sub index_bam {
    my ($self, $bam) = @_;
    my $samtools_index_cmd = Genome::Model::Tools::Sam::IndexBam->execute(
        bam_file => "$bam",
        use_version => $self->samtools_version,
    );

    unless ($samtools_index_cmd && $samtools_index_cmd->result) {
        die("Failed to run samtools index for BAM file: $bam");
    }
}

sub generate_idxstats {
    my ($self, $bam) = @_;

    my $samtools_idxstats_path = Genome::Sys->create_temp_file_path();
    my $samtools_idxstats_cmd = Genome::Model::Tools::Sam::Idxstats->execute(
        bam_file => "$bam",
        output_file => $samtools_idxstats_path,
        use_version => $self->samtools_version,
    );

    unless ($samtools_idxstats_cmd && $samtools_idxstats_cmd->result) {
        die "Failed to run samtools idxstats for BAM file: $bam";
    }

    my $idxstats_hash_ref = $samtools_idxstats_cmd->parse_file_into_hashref(
        $samtools_idxstats_path
    );

    return $idxstats_hash_ref;
}

sub generate_tsvfile {
    my ($self, $idxstats) = @_;

    my $r_input_file = Genome::Sys->create_temp_file_path();

    my @headers = (
        'Re-sort ID',
        'ERCC ID',
        'subgroup',
        'ERCC Mix',
        'concentration (attomoles/ul)',
        'label',
        'count',
    );

    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $r_input_file,
        headers => \@headers,
        separator => "\t",
    );

    my %ercc_control = $self->load_ercc_control_info;

    for my $chr (keys %{$idxstats}) {
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
            'count' => $idxstats->{$chr}{map_read_ct},
        );
        $writer->write_one(\%data);
    }
    $writer->output->close;

    return Path::Class::File->new($r_input_file);
}

1;

__END__
