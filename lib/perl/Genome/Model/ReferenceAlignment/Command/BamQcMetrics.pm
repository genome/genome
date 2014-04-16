package Genome::Model::ReferenceAlignment::Command::BamQcMetrics;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceAlignment::Command::BamQcMetrics {
    is => 'Genome::Command::Base',
    doc => "Generate tab-delimited (*.tsv) files with BamQc metrics summarized by instrument data id.",
    has => [
        build_id => {
            is => 'Number',
            shell_args_position => 1,
        },
        output_directory => {
            is => 'Text',
            is_optional => 1,
            doc => 'The directory to write output files to.'
        }
    ],
};


sub help_detail {
    return "Generate tab-delimited (*.tsv) files with BamQc metrics summarized by instrument data id for the provided build.";
}


sub execute {
    my $self  = shift;
    my $build = Genome::Model::Build->get($self->build_id);

    die $self->error_message('Please provide valid build id') unless $build;
    die $self->error_message('The provided build '.$build->id. ' is not reference alignment build.')
        unless $build->model->type_name eq 'reference alignment';
    my %instrument_data_bamqc_paths;
    for my $instrument_data ($build->instrument_data) {
        my $instrument_data_id = $instrument_data->id;
        my ($alignment_result) = $build->alignment_results_for_instrument_data($instrument_data);
        unless ($alignment_result) {
            $self->warning_message('Missing alignment result for instrument data: '. $instrument_data->id);
            next;
        }
        my ($bamqc_result) = Genome::InstrumentData::AlignmentResult::Merged::BamQc->get(alignment_result_id => $alignment_result->id);
        unless ($bamqc_result) {
            $self->warning_message('Missing BamQc result for instrument data: '. $instrument_data->id);
            next;
        }
        my $bamqc_path = $bamqc_result->output_dir;
        if ($self->debug) {
            $self->status_message('Found BamQc output: '. $bamqc_path);
        }
        $instrument_data_bamqc_paths{$instrument_data->id} = $bamqc_path;
    }
    my @labels;
    my @directories;
    my @pdf_files;
    my @html_links;
    for my $id (sort keys %instrument_data_bamqc_paths) {
        push @labels, $id;
        push @directories, $instrument_data_bamqc_paths{$id};
        push @pdf_files, glob($instrument_data_bamqc_paths{$id} .'/*.pdf');
        push @html_links, glob($instrument_data_bamqc_paths{$id} .'/*.html');
        push @html_links, glob($instrument_data_bamqc_paths{$id} .'/*_fastqc/*.html');
    }
    my $output_basename = $self->build_id;
    if ($self->output_directory) {
        $output_basename = $self->output_directory .'/'. $self->build_id;
    }
    unless (@labels) {
        $self->warning_message('No BamQc output found for instrument data. Unable to run SummarizeAsText on BamQC output for build: '. $self->build_id);
        return 0;
    }
    unless (Genome::Model::Tools::BamQc::SummarizeAsText->execute(
        labels => join(',',@labels),
        directories => join(',',@directories),
        output_basename => $output_basename,
        labels_are_instrument_data_ids => 1,
    )) {
        die $self->error_message('Failed to run SummarizeAsText on BamQc output for build: '. $self->build_id);
    }

    my $merged_pdf = $output_basename .'.pdf';
    my $cmd = 'pdftk '. join(' ',@pdf_files) .' cat output '. $merged_pdf;
    unless (Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => \@pdf_files,
        output_files => [$merged_pdf],
    )) {
        die('Failed to merge PDF files!');
    }
    
    # TODO: Make a summary html file that provides web links to all instrument data SAMStat and FastQc output.
    my $html_summary = $output_basename.'.html';
    
    return 1;
}

1;

