package Genome::Model::RnaSeq::DetectFusionsResult::TophatFusionResult;

use strict;
use warnings;

use Genome;
require File::Which;
use Cwd "chdir";

class Genome::Model::RnaSeq::DetectFusionsResult::TophatFusionResult{
    is => "Genome::Model::RnaSeq::DetectFusionsResult",
    has => [
        known_fusions_result => {
            is => "Genome::Model::RnaSeq::DetectFusionsResult::TophatFusionResult::KnownFusionsFile",
            is_input => 1,
            is_optional => 1,
            id_by => "known_fusions_result_id"
        },
        #TODO - not ideal, hardcoded a default ID for now so we don't have to pass this input all the way up to the RnaSeq model,
        #if this gets used, we'll want to revisit it
        known_fusions_result_id => {
            is_transient => 1,
            is => 'Text',
            example_values => [121542754],
            doc => "id of the known fusions software result"
        },
        detector_params => {
            is_param => 1,
            doc => 'params for tophat-fusion and tophat fusion post, separated by a colon (:)'
        }
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    $self->_prepare_output_directory();
    $self->_prepare_staging_directory();
    $self->_bootstrap_directory_structure();

    my $cmd_path = $self->_path_for_command($self->version, "tophat");
    unless ($cmd_path) {
        die $self->error_message("Failed to find a path for tophat for version " . $self->version . "!");
    }

    my($params, $post_params) = split(/:/, $self->detector_params);
    my $output_directory = $self->temp_staging_directory;
    my $n_threads = $self->_available_cpu_count;

    $self->_put_bowtie_version_in_path();

    my $bowtie_index = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_or_create(
        reference_build_id => $self->alignment_result->reference_build->id,
        aligner_name => 'bowtie',
        aligner_version => $self->alignment_result->bowtie_version,
    );

    #tophat fusion needs to be invoked from within its staging dir
    Cwd::chdir($self->temp_staging_directory);

    my $bowtie_index_path = $bowtie_index->output_dir;

    my @params = ($self->alignment_result->params, $self->alignment_result->inputs);
    my %params_for_alignment = map {my $n = $_->name; $n => $_->$n} grep {$_->name !~ /^instrument_data/ } @params;

    if($params !~/--fusion-search/){
        $params .= " --fusion-search ";
    }
    if($params !~/-p[\s|\d]/){
        $params .= " -p $n_threads ";
    }
    $params_for_alignment{aligner_params} = $params;
    #$params_for_alignment{version} = $self->version if (split(".", $params{version}))[0] < '2';

    my $tophat_fusion_result = Genome::InstrumentData::AlignmentResult::Tophat->get_or_create(%params_for_alignment);
    foreach(glob($tophat_fusion_result->output_dir .'/*')){
        Genome::Sys->create_symlink($_, $output_directory . "/" . (split("/",$_))[-1]);
    }

    #my $tophat_fusion_cmd = "$cmd_path --fusion-search -p $n_threads $params -o $output_directory $bowtie_index_path $fastq1 $fastq2";
    #Genome::Sys->shellcmd(
        #cmd => $tophat_fusion_cmd,
        #input_files => [$fastq1, $fastq2, $output_directory, $bowtie_index_path],
    #);

    my $post_cmd_path = $self->_path_for_command($self->version, "tophat-fusion-post");
    Genome::Sys->shellcmd(
        cmd => "$post_cmd_path -p $n_threads $post_params $bowtie_index_path",
    );

    $self->_promote_data();
    $self->_remove_staging_directory();
    $self->_reallocate_disk_allocation();

    return $self;
}

sub _bootstrap_directory_structure {
    my $self = shift;

    my $index_result = Genome::Model::RnaSeq::DetectFusionsResult::TophatFusionResult::Index->get_or_create(
        reference_build => $self->alignment_result->reference_build
    );

    my @files = glob $index_result->output_dir . "/*";
    for(@files){
        my $filename = (split("/", $_))[-1];
        Genome::Sys->create_symlink($_, $self->temp_staging_dir . "/$filename");
    }

    $self->_create_sample_info_file();

    Genome::Sys->create_symlink($self->known_fusions_result->mcl_file, $self->temp_staging_directory . "/mcl");
}


sub _create_sample_info_file {
    my $self = shift;

    my %sample_data = ();
    my @instrument_data = $self->alignment_result->instrument_data;
    for (@instrument_data) {
        if ($sample_data{$_->sample->name}){
            $sample_data{$_->sample->name}{read_count} += $_->read_count;
        }else{
            $sample_data{$_->sample->name} = {
                fragment_length => $_->library->original_fragment_size,
                read_length => $_->read_length,
                num_reads => $_->read_count
            };
        }
    }

    my $sample_fh = Genome::Sys->open_file_for_writing($self->temp_staging_dir . "/sample_data.txt");

    for(sort keys %sample_data){
        my $sample = $sample_data{$_};
        $sample_fh->say(join("\t", $_, $sample->{fragment_length}, $sample->{read_length}, $sample->{read_count}));
    }
    $sample_fh->close();
}

sub _staging_disk_usage {
    return 60 * 1024 * 1024;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    return 'build_merged_alignments/tophat-fusion/' . $self->id;
}

1;
