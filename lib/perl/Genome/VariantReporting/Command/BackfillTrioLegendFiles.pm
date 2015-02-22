package Genome::VariantReporting::Command::BackfillTrioLegendFiles;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Command::BackfillTrioLegendFiles {
    is => 'Command::V2',

    has => [
        trio_process => {
            is => 'Genome::VariantReporting::Process::Trio',
            shell_args_position => 1,
            doc => "The trio process to backfill",
        },
    ],
};

sub execute {
    my $self = shift;

    my $p = $self->trio_process;
    my @merged_results = grep {$_->class =~ /erged/} $p->unique_results;

    # first pass gets first-class merged reports (i.e. input results are not merged)
    for my $merged_result (@merged_results) {
        next if ($merged_result->report_results)[0]->class =~ /erged/;
        $self->backfill($merged_result);
    }

    # first pass gets second-class merged reports (i.e. input results are merged)
    for my $merged_result (@merged_results) {
        next unless ($merged_result->report_results)[0]->class =~ /erged/;
        $self->backfill($merged_result);
    }

    1;
}

sub backfill {
    my $self = shift;
    my $merged_result = shift;

    $merged_result->temp_staging_directory($merged_result->output_dir);
    if (-f $merged_result->legend_path) {
        $self->status_message("MergedReport (%s) ALREADY HAD legend file: %s",
            $merged_result->id, $merged_result->legend_path);
    } else {
        my $new_file = $merged_result->merge_legend_files();
        if($new_file) {
            $self->status_message("MergedReport (%s) NOW HAS legend file: %s",
                $merged_result->id, $new_file);
        } else {
            $self->status_message("MergedReport (%s) doesn't create legend files",
                $merged_result->id);
        }
    }
}

1;
