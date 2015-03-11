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
            example_values => ['1.2'],
            default_value => '1.2',
            is_optional => '1',
        },
        samtools_max_mem => {
            is => 'Text',
            doc => 'The max memory required by samtools.',
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
    
    my $ercc_bowtie_index = $self->create_ercc_bowtie_index();

    my $remapped_bam = $self->generate_remapped_bam(
        index => $ercc_bowtie_index,
        input_bam => $self->bam_file,
    );

    $self->index_bam($remapped_bam);
    my $idxstats = $self->generate_idxstats($remapped_bam);
    my $tsv = $self->generate_tsvfile($idxstats);
    $self->copy_to( $tsv, Path::Class::Dir->new() );

    $self->run_analysis_script($tsv);

    return 1;
}

sub samtools {
    my $self = shift;
    my $version = $self->samtools_version;
#    my $samtools_path =
#      Genome::Model::Tools::Sam->path_for_samtools_version($version);

    my $samtools_path = '/usr/bin/samtools1.2';

    return $samtools_path;
}

sub bowtie2_align {
    my $self = shift;
    my $version = $self->bowtie2_version;
    my $bowtie2_dir = Path::Class::Dir->new("/usr/bin");
    my $bowtie2_align = $bowtie2_dir->file('bowtie2-align' . $version);
    unless (-e $bowtie2_align) {
        die "[err] Couldn't find '$bowtie2_align' on the file system!\n";
    }
    return $bowtie2_align;
}

sub bowtie2_build {
    my $self = shift;
    my $version = $self->bowtie2_version;
    my $bowtie2_dir = Path::Class::Dir->new("/usr/bin");
    my $bowtie2_build = $bowtie2_dir->file('bowtie2-build' . $version);
    unless (-e $bowtie2_build) {
        die "[err] Couldn't find '$bowtie2_build' on the file system!\n";
    }
    return $bowtie2_build;
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

sub create_ercc_bowtie_index {
    my $self = shift;

    my $index_dir = Path::Class::Dir->new(
        Genome::Sys->create_temp_directory('bowtie2-index')
    );
    my $fasta = Path::Class::File->new($self->ercc_fasta_file);

    my $cmd = join(' ',
        "cd $index_dir &&",
        $self->bowtie2_build->stringify, "$fasta", 'ERCC');

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

    $input_bam = Path::Class::File->new($input_bam);

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

sub copy_to {
    my $self = shift;
    my ($src, $dest) = @_;
    if (UNIVERSAL::isa($dest, Path::Class::File::)) {
        $dest = $dest->stringify;
        die "Can't copy to file $dest: it is a directory" if -d $dest;
    }
    elsif (UNIVERSAL::isa($dest, Path::Class::Dir::)) {
        $dest = $dest->stringify;
        die "Can't copy to directory $dest: it is a file" if -f $dest;
        die "Can't copy to directory $dest: no such directory"
          unless -d $dest;
    }
    elsif (ref $dest) {
        die "Don't know how to copy files to objects of type '"
            . ref($src) . "'";
    }

    File::Copy::cp($src->stringify, "${dest}")
      or die "[err] Trouble copying from $src to $dest : \n\n $! \n\n";

    return $src->new($dest);
}

1;

__END__
