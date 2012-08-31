package Genome::Model::Tools::BioSamtools::ListChromosomes;

use strict;
use warnings;

use Genome;

# The only reason to encapsulate this in a tool is to run with perl v5.10.1 or greater when other processes are running perl5.8

class Genome::Model::Tools::BioSamtools::ListChromosomes {
    is => ['Genome::Model::Tools::BioSamtools'],
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'The path to a SAM or BAM format file with a sequence dictionary.',
        },
        output_file => {
            is => 'Text',
            doc => 'An output file to print one reference chromosome per line.',
            is_optional => 1,
            is_output => 1,
        },
    ],
    has_optional => [
        chromosome_array_ref => {
        },
    ],
};

sub execute {
    my $self = shift;
    my $input_file = $self->input_file;
    my ($basename,$dirname,$suffix) = File::Basename::fileparse($input_file,qw/\.sam \.bam/);
    unless (defined($suffix)) {
        die('Only sam or bam files are supported.  Failed to parse input file suffix: '. $input_file);
    }
    if ($suffix eq '.sam') {
        my $tam = Bio::DB::Tam->open($self->input_file);
        my $header = $tam->header_read();
        my $n_targets = $header->n_targets;
        unless ($n_targets) { die('The sequence dictionary does not have SQ lines per target!'); }
        my $name_arrayref = $header->target_name();
        $self->chromosome_array_ref($name_arrayref);
    } elsif ($suffix eq '.bam') {
        my $bam = Bio::DB::Bam->open($self->input_file);
        my $header = $bam->header_read();
        my $n_targets = $header->n_targets;
        unless ($n_targets) { die('The sequence dictionary does not have SQ lines per target!'); }
        my $name_arrayref = $header->target_name();
        $self->chromosome_array_ref($name_arrayref);
    } else {
        # Should never happen
        die('Please implement logic for parsing '. $suffix . ' format input file '. $self->input_file);
    }
    if (defined($self->output_file)) {
        my $fh = Genome::Sys->open_file_for_writing($self->output_file);
        for my $chr (@{$self->chromosome_array_ref}) {
            print $fh $chr ."\n";
        }
        $fh->close;
    }
    return 1;
}
