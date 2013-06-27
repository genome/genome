package Genome::Model::SomaticVariation::Command::Compare::Variants;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticVariation::Command::Compare::Variants {
    #is => 'Genome::Model::Compare::Base',
    is => 'Command::V2',
    has_input => [
        from_build => {
            is => 'Genome::Model::Build',
            shell_args_position => 1,
        },
        to_build => {
            is => 'Genome::Model::Build',
            shell_args_position => 2,
        },
        output_dir => {
            is => 'FilesystemDirectory',
            shell_args_position => 3,
        },
    ],
};

sub execute {
    my $self = shift;
    my $from_build = $self->from_build;
    my $to_build = $self->to_build;
    my $output_dir = $self->output_dir;

    unless (-d $output_dir) {
        Genome::Sys->create_directory($output_dir);
    }

    $self->status_message('Compare variants from build '. $self->from_build->__display_name__ .' to build '. $self->to_build->__display_name__ .'\!');

    my $from_dir = $from_build->data_directory;
    my $to_dir = $to_build->data_directory;


    my @relative_names = map { /($from_dir)\/(.*)$/ && $2 } glob("$from_dir/effects/*.bed");

    if (-e "$output_dir/summary.tsv") {
        unlink "$output_dir/summary.tsv";
    }
    
    my $log = Genome::Sys->open_file_for_writing("$output_dir/summary.tsv");
    $log->print(join("\t",qw/FILE COMMON OLD NEW/),"\n");

    my %counts;
    for my $relative_name (@relative_names) {
        my $from_file = $from_dir . '/' . $relative_name;
        my $to_file = $to_dir . '/' . $relative_name;

        my $base_name = File::Basename::basename($relative_name);
        my $dir_name = File::Basename::dirname($relative_name);
        my $output_subdir_name = $output_dir . '/' . $dir_name;
        Genome::Sys->create_directory($output_subdir_name);

        my $miss_a_file = $output_subdir_name . '/old-' . $base_name;
        my $miss_b_file = $output_subdir_name . '/new-' . $base_name;
        my $output_file = $output_subdir_name . '/common-' . $base_name;
        
        $self->status_message("$relative_name: ...");
        Genome::Model::Tools::Joinx::Intersect->execute(
            input_file_a => $from_file,
            input_file_b => $to_file,
            miss_a_file => $miss_a_file,
            miss_b_file => $miss_b_file,
            output_file => $output_file,
            use_version => '1.7',
        );
    
        my ($old_cnt, $new_cnt, $common_cnt) = (
            Genome::Sys->line_count($miss_a_file),
            Genome::Sys->line_count($miss_b_file),
            Genome::Sys->line_count($output_file)
        );
        
        $counts{$relative_name} = {
            'old' => $old_cnt, 
            'new' => $new_cnt,
            'common' => $common_cnt,
        };
        $log->print(join("\t",$relative_name,$common_cnt,$old_cnt,$new_cnt),"\n");
    }

    $log->close;

    
    if (0) {
        Genome::Sys->shellcmd(
            cmd => "wc -l $output_dir/common-* | sed s/common-// > $output_dir/common.counts",
            output_files => ["$output_dir/common.counts"],
        );
        Genome::Sys->shellcmd(
            cmd => "wc -l $output_dir/old-* | sed s/old-// > $output_dir/old.counts",
            output_files => ["$output_dir/old.counts"],
        );
        Genome::Sys->shellcmd(
            cmd => "wc -l $output_dir/new-* | sed s/new-// > $output_dir/new.counts",
            output_files => ["$output_dir/new.counts"],
        );
    }

    #Genome::Sys->shellcmd(cmd => qq{join -1 2 -2 2 $output_dir/old.counts $output_dir/new.counts | perl -nae 'chomp; print join("\t",\@F),"\n"' | tee $output_dir/old_new.counts | tab2col});
    Genome::Sys->shellcmd(cmd => "echo; cat $output_dir/summary.tsv | tab2col");
    return 1;
}

1;

