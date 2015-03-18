package Genome::Model::Tools::Transcriptome::ErccMapUnaligned;

use strict;
use warnings;

use Genome;
use Path::Class qw();

class Genome::Model::Tools::Transcriptome::ErccMapUnaligned {
    is => 'Command::V2',
    has_input => [
        model => {
            is => 'Text',
            doc => 'A RNASeq model that has ERCC transcripts spiked in.',
            is_optional => '1,',
        },
        bam_file => {
            is => 'FilePath',
            doc => 'A path to a BAM file that has ERCC transcripts spiked in.',
            is_optional => '1,',
        },
        ercc_fasta_file => {
            is => 'FilePath',
            doc => 'The FASTA format sequence for all ERCC transcripts.',
            default_value => '/gscmnt/gc13001/info/model_data/jwalker_scratch/ERCC/ERCC92.fa',
            is_optional => '1',
        },
        ercc_spike_in_file => {
            is => 'FilePath',
            doc => 'The control analysis file provided by Agilent for the ERCC spike-in.',
            default_value => '/gscmnt/gc13001/info/model_data/jwalker_scratch/ERCC/ERCC_Controls_Analysis.txt',
            is_optional => '1',
        },
        samtools_version => {
            is => 'Text',
            doc => 'The version of samtools to run analysis.',
            example_values => ['1.2'],
            default_value => '1.2',
            is_optional => '1',
        },
        samtools_max_mem => {
            is => 'Integer',
            doc => 'The max memory required by samtools (in bytes).',
            example_values => ['14000000000'],
            default_value => '14000000000',
            is_optional => '1',
        },
        bowtie2_version => {
            is => 'Text',
            doc => 'The version of bwa to run analysis.',
            example_values => ['2.1.0'],
            default_value => '2.1.0',
            is_optional => '1',
        },
    ],
    has_output => [
        pdf_file => {
            is => 'FilePath',
            doc => 'The output PDF with histograms and linearity plot.',
            example_values => ['<model-id>.ERCC.QC.pdf', '<id>.bam.ERCC.QC.pdf'],
            is_optional => '1',
        },
        raw_stats_file => {
            is => 'FilePath',
            doc => 'The raw data (in TSV format) used to create the pdf file.',
            example_values => ['<model-id>.ERCC.QC.tsv', '<id>.bam.ERCC.QC.tsv'],
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

    my $input_bam = $self->get_bam();
    $self->setup_outputs();

    my $ercc_bowtie_index = $self->create_ercc_bowtie_index();

    my $remapped_bam = $self->generate_remapped_bam(
        index => $ercc_bowtie_index,
        input_bam => $input_bam,
    );

    $self->index_bam($remapped_bam);
    my $idxstats = $self->generate_idxstats($remapped_bam);
    my $tsv = $self->generate_tsvfile($idxstats);
    $self->save_tsv_stats($tsv);

    $self->run_analysis_script($tsv);

    return 1;
}

sub get_bam_from_model {
    my $self = shift;
    my $model = Genome::Model->get(id => $self->model);

    my $build = $model->last_succeeded_build;
    unless ($build) {
        die "Couldn't find a recent successful build for model: ",
            $self->model, "\n"
    }
    $self->status_message("Using build: " . $build->id);

    my $bam = $build->merged_alignment_result->bam_path
      or die "Didn't find a bam file associated with build ",
             $self->model, "\n";
    $bam = Path::Class::File->new($bam);
    unless (-e $bam) {
        die "Didn't find bam file: '$bam' on file system!\n";
    }
    $self->status_message("Using BAM: $bam");
    $self->bam_file("$bam");

    return $bam;
}

sub setup_outputs {
    my $self = shift;
    my $dir = Path::Class::Dir->new(); # current working directory
    $self->setup_pdf_file($dir);
    $self->setup_raw_stats_file($dir);
}

sub setup_pdf_file {
    my ($self, $dir) = @_;

    my $name;
    if ($self->model) {
        $name = join('.', $self->model, 'ERCC.QC', 'pdf');
    }
    else {
        my $bam = Path::Class::File->new($self->bam_file);
        my $basename = $bam->basename;
        $name = join('.', $basename, 'ERCC.QC', 'pdf');
    }

    my $f = $dir->file($name);
    $self->pdf_file("$f");
    return $f;
}

sub setup_raw_stats_file {
    my ($self, $dir) = @_;

    my $name;
    if ($self->model) {
        $name = join('.', $self->model, 'ERCC.QC', 'tsv');
    }
    else {
        my $bam = Path::Class::File->new($self->bam_file);
        my $basename = $bam->basename;
        $name = join('.', $basename, 'ERCC.QC', 'tsv');
    }

    my $f = $dir->file($name);
    $self->raw_stats_file("$f");
    return $f;
}

sub get_bam {
    my $self = shift;

    unless ($self->bam_file || $self->model) {
        die "Please specify either a model via '--model' or ",
            "bam file via '--bam_file' to proceed!\n";
    }

    if ($self->bam_file && $self->model) {
        die "Specify only ONE '--model' or '--bam-file', NOT both!\n";
    }

    my $bam;
    if ($self->bam_file) {
        $bam = Path::Class::File->new($self->bam_file);
    }
    else {
        $bam = $self->get_bam_from_model();
    }

    unless (-e $bam) {
        die "Couldn't find bam: '$bam' on file system!\n";
    }

    return $bam;
}

sub _bin_dir {
    return Path::Class::Dir->new("/usr/bin");
}

sub samtools {
    my $self = shift;
    my $version = $self->samtools_version;
#    my $samtools_path =
#      Genome::Model::Tools::Sam->path_for_samtools_version($version);

    my $samtools_path = $self->_bin_dir->file('samtools1.2');

    return $samtools_path;
}

sub bowtie2_align {
    my $self = shift;
    my $bowtie2_align = Genome::Model::Tools::Bowtie->path_for_bowtie_version(
        $self->bowtie2_version,
        'align'
    );
    return Path::Class::File->new($bowtie2_align);
}

sub bowtie2_build {
    my $self = shift;
    my $bowtie2_build = Genome::Model::Tools::Bowtie->path_for_bowtie_version(
        $self->bowtie2_version,
        'build'
    );
    return Path::Class::File->new($bowtie2_build);
}

sub ERCC_analysis_script {
    my $self = shift;
    my $ercc_r = $self->__meta__->module_path;
    $ercc_r =~ s/\.pm$/\.R/;
    return Path::Class::File->new($ercc_r);
}

sub run_analysis_script {
    my ($self, $tsv) = @_;
    my $cmd = join(' ',
        $self->ERCC_analysis_script->stringify,
        "--data $tsv",
        '--output', $self->pdf_file
    );
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

sub create_ercc_bowtie_index {
    my $self = shift;

    my $index_dir = Path::Class::Dir->new(
        Genome::Sys->create_temp_directory('bowtie2-index')
    );
    my $fasta = Path::Class::File->new($self->ercc_fasta_file);

    my $cmd = join(' ',
        "cd $index_dir &&",
        $self->bowtie2_build->stringify, "$fasta", 'ERCC',
        '1>/dev/null'
    );

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$fasta->stringify],
    );

    my $search_pattern = $index_dir->file('ERCC.*');
    my @index_files = glob("$search_pattern");
    unless (@index_files == 6) {
        die "[err] Didn't create the proper bowtie2 index file set in ",
          "$index_dir !\n";
    }

    return $index_dir->file('ERCC');
}

sub generate_remapped_bam {
    my $self = shift;
    my %args = @_;
    my ($index, $input_bam) = @args{'index', 'input_bam'};

    my $remapped_bam_basename = Genome::Sys->create_temp_file_path();
    my $remapped_bam_path =
      Path::Class::File->new($remapped_bam_basename . '.bam');


    my $max_mem = $self->samtools_max_mem;

    # get the unaligned records from the input bam
    my $cmd1 = join(' ',
        $self->samtools, 'view', '-@ 4 -h -b -f 12',
        "$input_bam", q{'*'}
    );

    # gather the unaligned reads
    my $cmd2 = join(' ',
        $self->samtools, 'bam2fq',
        '-s /dev/null',
        '-'
    );

    # align the unaligned reads against the ERCC reference with bowtie2
    my $bowtie_stderr = Path::Class::Dir->new(
        Genome::Sys->create_temp_file_path('bowtie2.stderr')
    );
    my $cmd3 = join(' ',
        $self->bowtie2_align->stringify,
        "-x $index",
        '-U -',
        '--very-fast',
        '--threads 4',
        "2>$bowtie_stderr"
    );

    # sort the remapped alignments and generate the remapped bam
    my $cmd4 = join(' ',
        $self->samtools, 'sort',
        '-@ 4',
        "-m $max_mem",
        '-',
        "$remapped_bam_basename"
    );

    my $stdout =
      Path::Class::Dir->new(
          Genome::Sys->create_temp_file_path('bowtie-remapped.stdout')
      );
    my $stderr = Path::Class::Dir->new(
        Genome::Sys->create_temp_file_path('bowtie-remapped.stderr')
    );

    my $stream = Genome::Sys::ShellPipeline->create(
        pipe_commands => [$cmd1, $cmd2, $cmd3, $cmd4],
        redirects => "> $stdout 2> $stderr",
    );

    $stream->execute
      or die "[err] Trouble executing remapped bam stream command!\n";

    unless (-e $remapped_bam_path) {
        die "[err] Couldn't find the remapped bam: $remapped_bam_path \n!";
    }

    return $remapped_bam_path;
}

sub index_bam {
    my ($self, $bam) = @_;

    my $cmd = join(' ', $self->samtools, 'index', "$bam");
    my $index_file = Path::Class::File->new(
        $bam->stringify . '.bai'
    );

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => ["$bam"],
        output_files => ["$index_file"],
        skip_if_output_is_present => 0
    );
}

sub generate_idxstats {
    my ($self, $bam) = @_;

    my $samtools_idxstats_path = Path::Class::File->new(Genome::Sys->create_temp_file_path());

    my $cmd = join(' ',
        $self->samtools, 'idxstats',
        "$bam",
        '>', "$samtools_idxstats_path"
    );

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => ["$bam"],
        output_files => ["$samtools_idxstats_path"],
        skip_if_output_is_present => 0
    );

    my $idxstats_hash_ref = Genome::Model::Tools::Sam::Idxstats->parse_file_into_hashref(
        "$samtools_idxstats_path"
    );

    return $idxstats_hash_ref;
}

sub generate_tsvfile {
    my ($self, $idxstats) = @_;

    my $r_input_file = Genome::Sys->create_temp_file_path();

    my @headers = (
        q{'Re-sort ID'},
        q{'ERCC ID'},
        q{'subgroup'},
        q{'Mix 1 concentration (attomoles/ul)'},
        q{'Mix 2 concentration (attomoles/ul)'},
        q{'label'},
        q{'count'},
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
        my $mix1_concentration = $ercc_data->{'concentration in Mix 1 (attomoles/ul)'};
        my $mix2_concentration = $ercc_data->{'concentration in Mix 2 (attomoles/ul)'};
        my %data = (
            q{'Re-sort ID'} => $ercc_data->{'Re-sort ID'},
            q{'ERCC ID'} => $ercc_data->{'ERCC ID'},
            q{'subgroup'} => $ercc_data->{'subgroup'},
            q{'Mix 1 concentration (attomoles/ul)'} => $mix1_concentration,
            q{'Mix 2 concentration (attomoles/ul)'} => $mix2_concentration,
            q{'label'} => 'na',
            q{'count'} => $idxstats->{$chr}{map_read_ct},
        );
        $writer->write_one(\%data);
    }
    $writer->output->close;

    return Path::Class::File->new($r_input_file);
}

sub save_tsv_stats {
    my ($self, $tsv) = @_;

    my $dst_file = Path::Class::File->new($self->raw_stats_file);
    if (-e $dst_file) {
        $self->status_message("Removing an earlier version of '$dst_file'");
        $dst_file->remove();
    }

    $self->status_message("Saving raw stats to $dst_file");
    Genome::Sys->copy_file("$tsv", "$dst_file");
}

1;

__END__
