package Genome::Model::Tools::Tophat::AlignReads;

use strict;
use warnings;

use version;

use Genome;
use Genome::Sys;
use File::Path;

class Genome::Model::Tools::Tophat::AlignReads {
    is  => 'Genome::Model::Tools::Tophat',
    has_input => [
        reference_path => {
            doc => 'The path to the bowtie indexed reference',
        },
        read_1_fastq_list => {
            doc => 'A comma separated list of read 1 fastq files',
        },
        read_2_fastq_list => {
            is_optional => 1,
            doc => 'A comma separated list of read 2 fastq files when paired-end.(NOTICE: order must match read_1_fastq)',
        },
        insert_size => {
            is_optional => 1,
            doc => 'This is the expected (mean) inner distance between mate pairs. For, example, for paired end runs with fragments selected at 300bp, where each end is 50bp, you should set -r to be 200. There is no default, and this parameter is required for paired end runs'
        },
        insert_std_dev => {
            is_optional => 1,
            doc => 'The standard deviation for the distribution on inner distances between mate pairs. The default is 20bp.',
            default_value => 20,
        },
        aligner_params => {
            is_optional => 1,
            doc => 'Additional params passed to the Tophat aligner',
        },
        alignment_directory => {
            doc => 'The output directory where Tophat will write both temporary and final output files',
        },
        bowtie_version => {
            is => 'Text',
            doc => "the version of bowtie for tophat to run internally"
        }
    ],
    has_output => [
        aligner_output_file => {
            calculate_from => ['alignment_directory'],
            calculate => q|
                return $self->alignment_directory .'/tophat.aligner_output';
            |,
        },
        sam_file => {
            calculate_from => ['alignment_directory'],
            calculate => q|
                return $self->alignment_directory .'/accepted_hits.sam';
            |,
        },
        bam_file => {
            calculate_from => ['alignment_directory'],
            calculate => q|
                return $self->alignment_directory .'/accepted_hits.bam';
            |,
        },
        coverage_file => {
            calculate_from => ['alignment_directory'],
            calculate => q|
                return $self->alignment_directory .'/coverage.wig';
            |,
        },
        junctions_file => {
            calculate_from => ['alignment_directory'],
            calculate => q|
                return $self->alignment_directory .'/junctions.bed';
            |,
        },
        insertions_file => {
            calculate_from => ['alignment_directory'],
            calculate => q|
                return $self->alignment_directory .'/insertions.bed';
            |,
        },
        deletions_file => {
            calculate_from => ['alignment_directory'],
            calculate => q|
                return $self->alignment_directory .'/deletions.bed';
            |,
        },
        left_kept_reads_info => {
            calculate_from => ['alignment_directory'],
            calculate => q|
                return $self->alignment_directory .'/left_kept_reads.info';
            |,
        },
        right_kept_reads_info => {
            calculate_from => ['alignment_directory'],
            calculate => q|
                return $self->alignment_directory .'/right_kept_reads.info';
            |,
        },
        tmp_tophat_directory => {
            calculate_from => ['alignment_directory'],
            calculate => q|
                return $self->alignment_directory .'/tmp';
            |,
        },
    ],
};

sub help_synopsis {
    return <<EOS
    A Tophat based utility for aligning reads.
EOS
}

sub help_brief {
    return <<EOS
    A Tophat based utility for aligning reads.
EOS
}

sub help_detail {
    return <<EOS
Provides an interface to the Tophat aligner.

Input fastq files should be in SANGER fastq format.  If Solexa(prior GAPipeline v1.3) or Illumina(GAPipeline v1.3 or later)
fastq files are used as input, one of --solexa-quals or --solexa1.3-quals parameters should be defined respective to the input fastq format.

EOS
}

sub create {
    my $class = shift;
    my $self  = $class->SUPER::create(@_);

    unless ($self) {
        return;
    }

    unless ( $self->use_version ) {
        my $msg = 'use_version is a required parameter to ' . $class;
        $self->delete;
        die($msg);
    }
    unless ( $self->tophat_path ) {
        my $msg =
            'No path found for tophat version '
          . $self->use_version
          . ".  Available versions are:\n";
        $msg .= join( "\n", $self->available_tophat_versions );
        $self->delete;
        die($msg);
    }
    if ($self->read_2_fastq_list) {
        unless (defined($self->insert_size)) {
            die('Failed to provide a insert size for paired-end mates');
        }
    }
    return $self;
}


sub execute {
    my $self = shift;

    my $alignment_directory = $self->alignment_directory;
    $self->debug_message("OUTPUT PATH: $alignment_directory\n");


    # we resolve these first, since we might just print the paths we work with then exit
    my $files_to_align = $self->read_1_fastq_list;
    my $is_paired_end;
    my $insert_size;
    my $insert_sd;
    if ($self->read_2_fastq_list) {
        $files_to_align .= ' '. $self->read_2_fastq_list;
        $insert_sd = $self->insert_std_dev;
        $insert_size = $self->insert_size;
        $is_paired_end = 1;
    } else {
        $is_paired_end = 0;
    }
    $self->debug_message("INPUT PATH(S): $files_to_align\n");

    # prepare the refseq
    my $ref_seq_file =  $self->reference_path;
    unless (-e $ref_seq_file) {
        $self->error_message('Failed to find reference path '. $ref_seq_file);
        die($self->error_message);
    }
    #TODO:  A check for bowtie index files could be added to the reference placeholder...
    my $ref_seq_index_file =  $ref_seq_file .'.1.ebwt';
    unless (-e $ref_seq_index_file) {
        $self->error_message("Reference build index path '$ref_seq_index_file' does not exist.");
        die($self->error_message);
    }
    $self->debug_message("REFSEQ PATH: $ref_seq_file\n");

    # these are general params not infered from the above
    my $aligner_output_file = $self->aligner_output_file;
    my $aligner_params = $self->aligner_params;

    # TODO: Does tophat/bowtie need a param for quality conversion???
    #if ($instrument_data->resolve_quality_converter eq 'sol2phred') {
    #    $aligner_params .= ' --solexa1.3-quals';
    #} elsif ($instrument_data->resolve_quality_converter eq 'sol2sanger') {
    #    $aligner_params .= ' --solexa-quals';
    #} else {
    #    $self->error_message('Failed to resolve fastq quality coversion!');
    #    die($self->error_message);
    #}

    # RESOLVE A STRING OF ALIGNMENT PARAMETERS
    if ($is_paired_end && defined($insert_size) && defined($insert_sd)) {
        #The standard deviation is the above insert size deviation from sls, does below insert size stdev matter?
        $aligner_params .= ' --mate-inner-dist '. $insert_size .' --mate-std-dev '. $insert_sd;
    }

    my $cmdline = $self->tophat_path
    . sprintf(' --output-dir %s %s %s %s ',
              $alignment_directory,
              $aligner_params,
              $ref_seq_file,
              $files_to_align) .' > '. $aligner_output_file .' 2>&1';
    my @input_files = map{ split(',', $_) }  split(' ',$files_to_align);

    $self->debug_message("COMMAND: $cmdline\n");
    my $aligner_output_files = [$self->bam_file];
    if (version->parse($self->use_version) < version->parse('1.1.0')) {
        $aligner_output_files = [$self->sam_file];
    }
    my $bowtie_path = Genome::Model::Tools::Bowtie->path_for_bowtie_version($self->bowtie_version);
    unless($bowtie_path){
        die($self->error_message("unable to find a path for the bowtie version specified!"));
    }

    #put the desired version of bowtie in the path for tophat to run internally
    $bowtie_path =~ s/bowtie$//;
    $ENV{PATH} = $bowtie_path . ":" . $ENV{PATH};

    Genome::Sys->shellcmd(
                                          cmd                         => $cmdline,
                                          input_files                 => \@input_files,
                                          output_files                => $aligner_output_files,
                                          allow_zero_size_output_files => 1,
                                          skip_if_output_is_present   => 1,
                                      );
    if (version->parse($self->use_version) < version->parse('1.1.0')) {
        my $sam_to_bam = Genome::Model::Tools::Sam::SamToBam->execute(
            sam_file    => $self->sam_file,
            bam_file    => $self->bam_file,
            ref_list    => $ref_seq_file .'.fai',
            is_sorted   => 0,
            index_bam   => 1,
            fix_mate    => 1,
            keep_sam    => 1,
        );
        unless ($sam_to_bam) {
            $self->error_message('Error converting SAM file: '. $self->sam_file .' to BAM file '. $self->bam_file);
            die($self->error_message);
        }
    } else {
        unless (Genome::Model::Tools::Sam::IndexBam->execute(
                bam_file => $self->bam_file,
            )) {
            die $self->error_message('Failed to index BAM file '. $self->bam_file);
        }
    }

    unless ($self->verify_aligner_successful_completion) {
        $self->error_message('Failed to verify Tophat successful completion!');
        die($self->error_message);
    }
    #Lets leave these around until we know they are completely unnecessary
    #my @intermediate_files = glob($self->alignment_directory .'/*.fq*');
    #for my $intermediate_file (@intermediate_files) {
        #unless (unlink($intermediate_file)) {
        #    $self->error_message('Failed to remove intermediate tophat file '. $intermediate_file .":  $!");
        #    die($self->error_message);
        #}
    #}
    return 1;
}

sub verify_aligner_successful_completion {
    my $self = shift;

    my @output_files = $self->output_files;
    for my $output_file (@output_files) {
        unless (-e $output_file) {
            $self->error_message("Alignment output file '$output_file' not found.");
            return;
        }
    }
    my $aligner_output_fh = Genome::Sys->open_file_for_reading($self->aligner_output_file);
    unless ($aligner_output_fh) {
        $self->error_message('Failed to open aligner output file '. $self->aligner_output_file .":  $!");
        return;
    }
    while(<$aligner_output_fh>) {
        if (m/^Run complete/) {
            $aligner_output_fh->close();
            return 1;
        }
    }
    $aligner_output_fh->close();
    return;
}

sub output_files {
    my $self = shift;
    my @output_files;
    for my $method (qw/aligner_output_file bam_file junctions_file/) {
        push @output_files, $self->$method;
    }
    if (version->parse($self->use_version) < version->parse('1.1.0')) {
        push @output_files, $self->sam_file;
        push @output_files, $self->coverage_file;
    }
    # NOTE: There may be a gap checking the existence of files between v1.1.0 and v1.3.0
    if (version->parse($self->use_version) >= version->parse('1.3.0')) {
        push @output_files, $self->insertions_file;
        push @output_files, $self->deletions_file;
        push @output_files, $self->left_kept_reads_info;
        #push @output_files, $self->right_kept_reads_info; # This won't exist for single-end reads
    }
    return @output_files;
}

1;
