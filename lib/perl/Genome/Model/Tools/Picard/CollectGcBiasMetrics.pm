package Genome::Model::Tools::Picard::CollectGcBiasMetrics;

use strict;
use warnings;

use Genome;
use File::Basename;


class Genome::Model::Tools::Picard::CollectGcBiasMetrics {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file   => {
            is  => 'String',
            doc => 'The BAM file to run on.',
        },
        output_file  => {
            is  => 'String',
            doc => 'The output metrics file',
        },
        refseq_file  => {
            is  => 'String',
            doc => 'The reference sequence file',
        },
        chart_output => {
            is  => 'String',
            doc => 'The PDF file to render the chart to. Default is GC_bias_chart.pdf in output_file dir',
            is_optional => 1,
        },
        summary_output => {
            is  => 'String',
            doc => 'The text file to write summary metrics to. Default is GC_bias_summary.txt',
            is_optional => 1,
        },
        window_size  => {
            is  => 'Integer',
            doc => 'The size of windows on the genome that are used to bin reads. Default value: 100',
            default_value => 100,
            is_optional   => 1,
        },
        min_genome_fraction => {
            is  => 'Number',
            doc => 'For summary metrics, exclude GC windows that include less than this fraction of the genome. Default value: 1.0E-5.',
            default_value => '1.0E-5',
            is_optional   => 1,
        },
        max_records_in_ram => {
            is => 'Integer',
            doc => 'The number of alignment records to store in RAM before spilling to disk.',
            default_value => 500000,
            is_optional => 1,
        },
        sort_bam => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'Sort the BAM file before cleaning or collectin GC metrics.',
        },
        clean_bam => {
            is => 'Text',
            is_optional => 1,
            valid_values => ['trim','remove','none'],
            default_value => 'trim',
            doc => 'Flag to trim overhanging alignments with picard or remove reads that align beyond the end of chromosomes with bio-samtools',
        },
        clean_bam_summary => {
            is => 'Text',
            is_optional => 1,
            doc => 'A file path to store cleaning output.',
        },
    ],
};

sub help_brief {
    'Tool to collect GC bias metrics from a BAM file.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#CollectGcBiasMetrics
EOS
}

sub execute {
    my $self = shift;

    for my $type qw(input refseq) {
        my $property_name = $type.'_file';
        unless ($self->$property_name and -s $self->$property_name) {
            $self->error_message("$property_name is invalid");
            return;
        }
    }
    my $out_dir = dirname $self->output_file;
    my $input_file = $self->input_file;
    my $basename = basename( $input_file, qw/.bam/ );
    my $tmpdir = $self->temp_directory or die 'no temp_directory';
    if ($self->sort_bam) {
        my $sorted_bam_file = $tmpdir . '/' . $basename . '.sorted.bam';
        ( !-e $sorted_bam_file ) or die "already exists: $sorted_bam_file";
        unless (Genome::Model::Tools::Picard::SortSam->execute(
            input_file => $self->input_file,
            output_file => $sorted_bam_file,
            max_records_in_ram => $self->max_records_in_ram,
            maximum_memory => $self->maximum_memory,
            maximum_permgen_memory => $self->maximum_permgen_memory,
            use_version => $self->use_version,
            temp_directory => $tmpdir,
        )) {
            die('Failed to run G:M:T:Picard::SortSam on BAM file: '. $input_file);
        }
        $input_file = $sorted_bam_file;
    }
    $basename = basename( $input_file, qw/.bam/ );
    my $clean_bam_file = $tmpdir . '/' . $basename . '.cleaned.bam';
    ( !-e $clean_bam_file ) or die "already exists: $clean_bam_file";
    if ($self->clean_bam eq 'remove') {
        unless (Genome::Model::Tools::BioSamtools::CleanBam->execute(
            input_bam_file => $self->input_file,
            output_bam_file => $clean_bam_file,
            summary_output_file => $self->clean_bam_summary,
        )) {
            die('Failed to run G:M:T:BioSamtools::CleanBam on BAM file: '. $input_file);
        }
        $input_file = $clean_bam_file;
    } elsif ($self->clean_bam eq 'trim') {
        unless (Genome::Model::Tools::Picard::CleanSam->execute(
            input_file => $self->input_file,
            output_file => $clean_bam_file,
            log_file => $self->clean_bam_summary,
            use_version => $self->use_version,
            maximum_memory => $self->maximum_memory,
            maximum_permgen_memory => $self->maximum_permgen_memory,
            temp_directory => $tmpdir,
        )) {
            die('Failed to run G:M:T:Picard::CleanSam on BAM file: '. $input_file);
        }
        $input_file = $clean_bam_file;
    }

    # TODO: remove sorted if it was made above? (to save disk)

    my $cmd = $self->picard_path .'/CollectGcBiasMetrics.jar net.sf.picard.analysis.CollectGcBiasMetrics';
    $cmd   .= ' OUTPUT='. $self->output_file  .' INPUT='. $input_file .' REFERENCE_SEQUENCE='. $self->refseq_file;

    my $chart = $self->chart_output   || $out_dir . '/GC_bias_chart.pdf';
    my $sum   = $self->summary_output || $out_dir . '/GC_bias_summary.txt';

    $cmd .= ' CHART_OUTPUT='.$chart .' SUMMARY_OUTPUT='.$sum;

    if ($self->window_size) {
        $cmd .= ' WINDOW_SIZE=' . $self->window_size;
    }
    if ($self->min_genome_fraction) {
        $cmd .= ' MINIMUM_GENOME_FRACTION='. $self->min_genome_fraction;
    }
    if ($self->max_records_in_ram) {
        $cmd .= ' MAX_RECORDS_IN_RAM='. $self->max_records_in_ram;
    }
    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => [$input_file],
        #output_files => [$self->output_file, $chart],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
