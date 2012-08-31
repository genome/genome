package Genome::Model::Event::Build::RnaSeq::PrepareReads;

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Build::RnaSeq::PrepareReads {
    is => ['Genome::Model::Event'],
    has => [
    ],
};
sub bsub_rusage {
    return "-R 'select[model!=Opteron250 && type==LINUX64 && mem>16000 && tmp>150000] rusage[tmp=150000, mem=16000]' -M 16000000";
}

sub execute {
    my $self = shift;
    unless (-d $self->build_directory) {
        $self->create_directory($self->build_directory);
        $self->status_message("Created build directory: ".$self->build_directory);
    } else {
        $self->status_message("Build directory exists: ".$self->build_directory);
    }
    my $build_fastq_directory = $self->build->accumulated_fastq_directory;
    unless (-d $build_fastq_directory) {
        Genome::Sys->create_directory($build_fastq_directory);
    }
    my $fastq_directory = $self->fastq_directory;
    unless (-d $fastq_directory) {
        Genome::Sys->create_directory($fastq_directory);
    }
    my $instrument_data = Genome::InstrumentData->get($self->instrument_data_id);

    my $quality_format;
    my $dump_fastq_method;

    if(defined $instrument_data->bam_path && -s $instrument_data->bam_path) {
        $quality_format='Standard';
        $dump_fastq_method = 'dump_sanger_fastq_files';
    } else {
        my $quality_converter = $instrument_data->resolve_quality_converter;
        if($quality_converter eq 'sol2phred') {
            $quality_format = 'Illumina';
            $dump_fastq_method = 'dump_illumina_fastq_files';
        } elsif ($quality_converter eq 'sol2sanger') {
            $quality_format = 'Solexa';
            $dump_fastq_method = 'dump_solexa_fastq_files';
        } elsif ($quality_converter eq 'none') {
            $quality_format = 'Standard';
            $dump_fastq_method = 'dump_sanger_fastq_files';
        } else {
            $self->error_message('Unknown or unsupported quality converter ' . $quality_converter . ' encountered.');
            return;
        }
    }

    unless ($self->instrument_data->$dump_fastq_method(directory => $fastq_directory)) {
        $self->error_message('Failed to dump fastq file for '. $self->instrument_data_id .' to directory '. $fastq_directory);
        return;
    }
    my $picard_version = $self->model->picard_version;
    unless ($picard_version) {
        $picard_version = Genome::Model::Tools::Picard->default_picard_version;
        $self->warning_message('Picard version not defined in processing profile.  Using default picard version: '. $picard_version);
    }
    unless (Genome::Model::Tools::Picard::FastqToSam->execute(
        fastq => $self->fastq_directory .'/'. $instrument_data->read1_fastq_name,
        fastq2 => $self->fastq_directory .'/'. $instrument_data->read2_fastq_name,
        output => $self->fastq_directory .'/s_'. $instrument_data->subset_name .'_sequence.bam',
        quality_format => $quality_format,
        sample_name => $instrument_data->sample_name,
        library_name => $instrument_data->library_name,
        log_file => $self->fastq_directory .'/s_'. $instrument_data->subset_name .'_sequence.log',
        platform => 'illumina',
        platform_unit => $instrument_data->flow_cell_id .'.'. $instrument_data->subset_name,
        read_group_name => $instrument_data->id,
        sort_order => 'queryname',
        use_version => $picard_version,
        maximum_memory => 12,
        maximum_permgen_memory => 256,
        max_records_in_ram => 3000000,
    )) {
        die ('Failed to create per lane, unaligned BAM file: '. $self->fastq_directory .'/s_'. $instrument_data->subset_name .'_sequence.bam');
    }
    return 1;
}

sub fastq_directory {
    my $self = shift;
    return $self->build->accumulated_fastq_directory .'/'. $self->instrument_data_id;
}


1;
