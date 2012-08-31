package Genome::Model::Tools::RepeatMasker::GenerateTable;

use strict;
use warnings;

use Genome;

use FASTAParse;

class Genome::Model::Tools::RepeatMasker::GenerateTable {
    is => 'Genome::Model::Tools::RepeatMasker::TableI',
    has => [
        fasta_file => {
            is => 'Text',
            doc => 'The masked file output from RepeatMasker(.masked) or original fasta of all sequences(.fa or .fasta).  File is used to calculate total base pair, and NOT the masked base pair.',
        },
        output_file => {
            is => 'Text',
            doc => 'The output from repeat masker(.out)',
        },
    ],
};

sub execute {
    my $self = shift;

    my $fasta_reader = IO::File->new($self->fasta_file,'r');
    unless ($fasta_reader) {
        die('Failed to open fasta file '. $self->fasta_file);
    }
    my $total_bp = 0;
    my $total_count = 0;
    eval {
        local $/ = "\n>";
        while (<$fasta_reader>) {
            if ($_) {
                chomp;
                if ($_ =~ /^>/) { $_ =~ s/\>//g }
                my $myFASTA = FASTAParse->new();
                $myFASTA->load_FASTA( fasta => '>' . $_ );
                my $seqlen = length( $myFASTA->sequence() );
                $total_bp += $seqlen;
                $total_count++;
                # TODO: GC content?
                # TODO: masked bases or just get from alignments below
            }
        }
    };
    $fasta_reader->close;
    if ($@) {die ($@); }
    $self->_total_count($total_count);
    $self->_total_bp($total_bp);
    my $parser = Bio::Tools::RepeatMasker->new(-file => $self->output_file);
    unless ($parser) {
        die ('Failed to create RepeatMasker parser for file: '. $self->output_file);
    }
    my %repeats;
    while (my $result = $parser->next_result) {
        my $tag = $result->hprimary_tag;
        my ($family,$class) = split("/",$tag);
        if ($family =~ /RNA/) {
            $family = 'Small RNA';
        }
        # Take either the hit length or the query length, but we are not accounting insertions/deletions/mismatches
        my $length = (($result->end - $result->start) + 1);
        if ($length < 1) {
            die(Data::Dumper::Dumper($result));
        }
        $repeats{$family}{elements}++;
        $repeats{$family}{base_pair} += $length;
        if ($class) {
            $repeats{$family}{$class}{elements}++;
            $repeats{$family}{$class}{base_pair} += $length;
        }
    }
    $self->print_table_from_hash_ref(\%repeats);
    return 1;
}
