package Genome::Model::GenotypeMicroarray::Command::ConvertGoldSnpBedToGeno;

use Genome;
use Genome::Info::IUB;

class Genome::Model::GenotypeMicroarray::Command::ConvertGoldSnpBedToGeno {
    is => 'Command',
    has => [
        gold_snp_bed => {
            is => 'Text',
        },
        output => {
            is => 'Text',
        },
    ],
};

sub execute {
    my $self = shift;

    my $gold_snp_bed = IO::File->new($self->gold_snp_bed);
    my $output = Genome::Sys->open_file_for_writing($self->output);
    while (my $bed_line = $gold_snp_bed->getline) {
        my ($chr, $start, $end, $call, $score, $depth, $a1t, $a2t) = split("\t", $bed_line);
        my ($ref, $iub) = split('/', $call);
        my $iub_string = Genome::Info::IUB::iub_to_string($iub);
        my ($a1, $a2) = split('', $iub_string);
        if ($a1t ne $a2t) {
            if (lc($a1t) eq 'ref' && lc($a1) ne lc($ref)) {
                ($a1, $a2) = ($a2, $a1);
            } elsif (lc($a2t) eq 'ref' && lc($a2) ne lc($ref)) {
                ($a1, $a2) = ($a2, $a1);
            }
            print STDERR join("\t", $chr, $end, $a1, $a1t, $a2, $a2t, $ref, $iub) . "\n";
        }
        $output->print(join("\t", $chr, $end, "$a1$a2"), "\n");
    }
    $gold_snp_bed->close;
    $output->close;
};

1;
