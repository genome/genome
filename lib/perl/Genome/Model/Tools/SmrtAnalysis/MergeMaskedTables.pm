package Genome::Model::Tools::SmrtAnalysis::MergeMaskedTables;

use strict;
use warnings;

use Genome;

use File::Basename;
use File::Path;
use File::Copy;

class Genome::Model::Tools::SmrtAnalysis::MergeMaskedTables {
    is => 'Genome::Model::Tools::SmrtAnalysis::Base',
    has_input => [
        output_masked_fofn => {  is_output => 1,},
        base_output_directory => { },
    ],
    has_optional_input => [
        input_masked_fofns => { },
        remove_originals => {
            default_value => 1,
        },
    ],
};

sub create {
    my $class = shift;
    my %params = @_;
    my $input_masked_fofns = delete($params{input_masked_fofns});
    my $self = $class->SUPER::create(%params);
    unless ($self) { return; }
    $self->input_masked_fofns($input_masked_fofns);
    return $self;
}

sub execute {
    my $self = shift;

    my $base_output_directory = $self->base_output_directory;
    my $output_directory = $base_output_directory .'/post_control_regions';
    unless (-d $output_directory) {
        unless (Genome::Sys->create_directory($output_directory)) {
            die('Failed to create output directory '. $output_directory);
        }
    }

    my @input_fofns = @{$self->input_masked_fofns};
    my $output_fofn_fh = Genome::Sys->open_file_for_writing($self->output_masked_fofn);
    my @directories;
    for my $input_fofn (@input_fofns) {
        my $input_fofn_fh = Genome::Sys->open_file_for_reading($input_fofn);
        while (my $line = $input_fofn_fh->getline) {
            chomp($line);
            if ($line =~ /^\s*$/) { next; }
            my ($basename,$dirname,$suffix) = File::Basename::fileparse($line,qw/\.rgn\.h5/);
            my $new_file = $output_directory .'/'. $basename . $suffix;
            File::Copy::move($line,$new_file) || die ('Failed to move old file '. $line .' to new location '. $new_file);
            print $output_fofn_fh $new_file ."\n";
            push @directories, $dirname;
        }
        $input_fofn_fh->close;
    }
    $output_fofn_fh->close;
    if ($self->remove_originals) {
        for my $dir (@directories) {
            File::Path::rmtree($dir) || die('Failed to remove intermediate directory '. $dir);
        }
        for my $file (@input_fofns) {
            unlink($file) || die('Failed to remove intermediate file '. $file);
        }
    }
    return 1;
}


1;
