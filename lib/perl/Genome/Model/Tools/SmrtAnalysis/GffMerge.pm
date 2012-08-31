package Genome::Model::Tools::SmrtAnalysis::GffMerge;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SmrtAnalysis::GffMerge {
    is => ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => {
        input_files => { },
        output_file => {  is_output => 1, },
    },
    has_optional_input => [
        remove_originals => {
            default_value => 1,
        },
    ],
};

sub create {
    my $class = shift;
    my %params = @_;
    my $input_files = delete($params{input_files});
    my $self = $class->SUPER::create(%params);
    unless ($self) { return; }
    $self->input_files($input_files);
    return $self;
}

sub execute {
    my $self = shift;
    my @input_files;
    if (ref($self->input_files) eq 'ARRAY') {
        @input_files = @{$self->input_files};
    } else {
        @input_files = split(',',$self->input_files);
    }
    my ($tmp_header_fh,$tmp_header_file) = Genome::Sys->create_temp_file;
    my ($tmp_gff_fh,$tmp_gff_file) = Genome::Sys->create_temp_file;
    my $version_printed = 0;
    for my $input_file (@input_files) {
        my $in_fh = Genome::Sys->open_file_for_reading($input_file);
        while (my $line = $in_fh->getline) {
            chomp($line);
            if ($line =~ /^##/) {
                if ($line =~ /gff-version/) {
                    unless ($version_printed) {
                        print $tmp_header_fh $line ."\n";
                        $version_printed = 1;
                    }
                } else {
                    print $tmp_header_fh $line ."\n";
                }
            } else {
                print $tmp_gff_fh $line ."\n";
            }
        }
    }
    $tmp_header_fh->close;
    $tmp_gff_fh->close;
    unless (Genome::Model::Tools::SmrtAnalysis::CatMerge->execute(
        input_files => [$tmp_header_file,$tmp_gff_file],
        output_file => $self->output_file,
    )) {
        die('Failed to merge the header and GFF file!');
    }
    if ($self->remove_originals) {
        for my $file (@input_files) {
            unlink($file) || die('Failed to remove original file '. $file);
        }
    }
    return 1;
}
