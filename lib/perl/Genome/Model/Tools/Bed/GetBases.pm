package Genome::Model::Tools::Bed::GetBases;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::GetBases {
    is => ['Command'],
    has_input => [
        bed_file => {
            is => 'Text',
            doc => 'The BED format file of intervals to pull out bases from reference.',
        },
        reference_sequence_build_id => {
            is => 'Integer',
            doc => 'The Genome Model Build ID of the refrence',
            is_optional => 1,
            default_value => 101947881,
        },
        fasta_file => {
            is => 'Text',
            doc => 'The output fasta file with defined intervals bases.',
        },
    ],
};


sub execute {
    my $self = shift;

    my $build = Genome::Model::Build::ImportedReferenceSequence->get($self->reference_sequence_build_id);
    unless ($build) {
        die('Failed to find ImportedReferenceSequence build for id: '. $self->reference_sequence_build_id);
    }
    my $bed_fh = Genome::Sys->open_file_for_reading($self->bed_file);
    unless ($bed_fh) {
        die;
    }
    my $fasta_fh = Genome::Sys->open_file_for_writing($self->fasta_file);
    unless ($fasta_fh) {
        die;
    }
    while (my $line = $bed_fh->getline) {
        chomp($line);
        if ($line =~ /^#/) { next; }
        my ($chr,$start,$stop,$name) = split("\t",$line);
        $start += 1;
        my $seq = $build->sequence($chr,$start,$stop);
        unless ($name) {
            $name = $chr .':'. $start .'-'. $stop;
        }
        print $fasta_fh '>'. $name ."\n";
        print $fasta_fh $seq ."\n";
    }
    $fasta_fh->close;
    $bed_fh->close;
    return 1;
}
