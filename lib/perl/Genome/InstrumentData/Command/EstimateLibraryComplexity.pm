package Genome::InstrumentData::Command::EstimateLibraryComplexity;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::EstimateLibraryComplexity {
    is => 'Genome::InstrumentData::Command',
    has => [
        output_file => {
            is_input => 1,
            is => 'Text',
            doc => 'The file path to place the instrument_data library complexity report',
        },
        picard_version => {
            is_input => 1,
            is_optional => 1,
            is => 'Text',
            doc => 'The version of Picard to use(greater than v1.21).',
            default_value => '1.21',
        },
    ],
};

sub help_brief {
    return 'A command to convert fastq files to BAM and run the Picard EstimateLibraryComplexity tool';
}

sub help_detail {
    return <<EOS
    This tool converts Illumina Fastq files into unaligned BAM files.
    Then it runs Picard EstimateLibraryComplexity on the resulting temporary BAM file.
    The duplication ratio and additional metrics from Picard are written to the defined output file per library.
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#EstimateLibraryComplexity
EOS
}

sub execute {
    my $self = shift;


    my $instrument_data = $self->instrument_data;
    
    my $quality_format;
    my @fastq_filenames;
    
    my $quality_converter = $instrument_data->resolve_quality_converter;
    
    if($quality_converter eq 'sol2sanger') {
        $quality_format = 'Solexa';
        @fastq_filenames = $instrument_data->dump_solexa_fastq_files;
    } elsif ($quality_converter eq 'sol2phred') {
        $quality_format = 'Illumina';
        @fastq_filenames = $instrument_data->dump_illumina_fastq_files;
    } elsif ($quality_converter eq 'none') {
        $quality_format = 'Sanger';
        @fastq_filenames = $instrument_data->dump_sanger_fastq_files;
    }
    
    #my @fastq_filenames = $instrument_data->fastq_filenames;
    my $tmp_bam = Genome::Sys->base_temp_directory .'/'. $instrument_data->id .'.bam';
    unless (Genome::Model::Tools::Picard::FastqToSam->execute(
        fastq => Genome::Sys->base_temp_directory.'/'. $instrument_data->read1_fastq_name,
        fastq2 => Genome::Sys->base_temp_directory .'/'. $instrument_data->read2_fastq_name ,
        sample_name => $instrument_data->sample_name,
        library_name => $instrument_data->library_name,
        platform => 'illumina',
        quality_format => $quality_format,
        use_version => $self->picard_version,
        output => $tmp_bam,
    )) {
        die('Failed to generate BAM for instrument data id '. $instrument_data->id);
    }
    unless (Genome::Model::Tools::Picard::EstimateLibraryComplexity->execute(
        use_version => $self->picard_version,
        input_file => [$tmp_bam],
        output_file => $self->output_file,
    )) {
        die('Failed to generate library complexity for instrument data id '. $instrument_data->id);
    }
    return 1;
}
