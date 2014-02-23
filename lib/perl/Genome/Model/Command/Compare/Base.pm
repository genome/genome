package Genome::Model::Command::Compare::Base;
use strict;
use warnings;
use Genome;

class Genome::Model::Command::Compare::Base {
    is => 'Command::V2',
    is_abstract => 1,
    has_input => [
        old_build => {
            is => 'Genome::Model::Build',
            shell_args_position => 1,
            doc => 'one build to be used in the comparison, by convention the one with previous/expected results',
        },
        new_build => {
            is => 'Genome::Model::Build',
            shell_args_position => 2,
            doc => 'the other buid ot be used in the comparision, by convention the one with new/experimental results',
        },
        output_dir => {
            is => 'FilesystemDirectory',
            shell_args_position => 3,
            doc => 'the directory into which results should go'
        },
    ],
    doc => 'compare results between two builds',
};

sub help_detail {
    my $class = shift->class;
    my $command_name = $class->command_name;
    return <<EOS
All comparison modules follow the same API:
    
    genome model PIPELINE compare COMPARISON_TYPE OLD_BUILD NEW_BUILD OUTPUT_DIR

Each pipeline may implement a variety of comparators.  Comparators with the same
name should produce comparable output across pipelines.

The names "old" and "new" are arbitrary, but should be used to in this order to
preserve the sanity of the analyst.  

The comparison tool should produce a directory with similar structure to that of the build director.

The results will have a directory structure similar to that of the input builds, 
except that three files will be made for each pair of comparable input files:
    old-NAME    results present only in the old build
    new-NAME    results present only in the new build
    common-NAME results common to both builds

An additional file is created at the end, counting differences: 
    diff-summary.tsv
EOS
}

sub execute {
    my $self = shift;
    my $old_build = $self->old_build;
    my $new_build = $self->new_build;
    my $output_dir = $self->output_dir;
    $DB::single = 1;
    unless (-d $output_dir) {
        Genome::Sys->create_directory($output_dir);
    }

    $self->debug_message('Compare variants from build '. $self->old_build->__display_name__ .' to build '. $self->new_build->__display_name__ .'\!');

    my $old_dir = $old_build->data_directory;
    my $new_dir = $new_build->data_directory;

    my $summary_path = "$output_dir/diff-summary.tsv";

    my @relative_names = map { /($old_dir)\/(.*)$/ && $2 } glob("$old_dir/effects/*.bed");

    if (-e $summary_path) {
        unlink $summary_path or die "Failed to remove file $summary_path: $!";
    }
    
    my $summary_fh = Genome::Sys->open_file_for_writing($summary_path);
    $summary_fh->print(join("\t",qw/FILE COMMON OLD NEW/),"\n");

    my %counts;
    for my $relative_name (@relative_names) {
        my $old_file = $old_dir . '/' . $relative_name;
        my $new_file = $new_dir . '/' . $relative_name;

        my $base_name = File::Basename::basename($relative_name);
        my $dir_name = File::Basename::dirname($relative_name);
        my $output_subdir_name = $output_dir . '/' . $dir_name;
        Genome::Sys->create_directory($output_subdir_name);

        my $old_only_results_file = $output_subdir_name . '/old-' . $base_name;
        my $new_only_results_file = $output_subdir_name . '/new-' . $base_name;
        my $common_results_file = $output_subdir_name . '/common-' . $base_name;
        
        $self->debug_message("$relative_name: ...");
        
        $self->_compare_files(
            $old_file,
            $new_file,
            $old_only_results_file,
            $new_only_results_file,
            $common_results_file,
        );
    
        my ($old_cnt, $new_cnt, $common_cnt) = (
            Genome::Sys->line_count($old_only_results_file),
            Genome::Sys->line_count($new_only_results_file),
            Genome::Sys->line_count($common_results_file)
        );
        
        $counts{$relative_name} = {
            'old' => $old_cnt, 
            'new' => $new_cnt,
            'common' => $common_cnt,
        };
        $summary_fh->print(join("\t",$relative_name,$common_cnt,$old_cnt,$new_cnt),"\n");
    }

    $summary_fh->close;

    Genome::Sys->shellcmd(cmd => "echo; cat $summary_path | tab2col");
    return 1;
}

1;

