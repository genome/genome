package Genome::Model::Tools::RepeatMasker::MergeOutput;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::RepeatMasker::MergeOutput {
    is => ['Command'],
    has => [
        input_files => { },
        output_file => { },
    ],
};

sub execute {
    my $self = shift;
    my $writer = Genome::Sys->open_file_for_writing($self->output_file);
    unless ($writer) {
        $self->error_message('Failed to open output file '. $self->output_file);
            return;
    }
    for my $input_file (@{$self->input_files}) {
        my $reader = Genome::Sys->open_file_for_reading($input_file);
        unless ($reader) {
            $self->error_message('Failed to read input file '. $input_file);
            return;
        }
        while ( my $line = $reader->getline) {
            chomp($line);
            if ($line =~ /no repetitive sequences detected/) {
                print $writer $line ."\n";
            }
            if ($line =~ /\d+/) {
                print $writer $line ."\n"
            }
        }
        $reader->close;
    }
    $writer->close;
    return 1;
}
