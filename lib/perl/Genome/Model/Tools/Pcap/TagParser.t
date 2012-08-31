#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use IO::String;
use Test::More;

use_ok('Genome::Model::Tools::Pcap::TagParser');

my $io = IO::String->new();
$io->print
(qq/
CT{
Contig24 comment consed 606 621 050427:142133
COMMENT{
good!
C}
}

CT{
Contig24 oligo consed 606 621 050427:142133
M_BB0392D19.29 ccctgagcgagcagga 60 U
L25990P6000A5 L25990P6000D4
}

CT{
Contig29.2 autoFinishExp autofinish 119 119 060831:122829
C
purpose: weak
0 915 0
dyeTerm customPrimer
fix cons errors: 4.69881 original cons errors: 5.64249
original single subclone bases: 886
primer: ggcaaatatggtgcaataaaac temp: 58 id: Trichinella_spiralis_060315.pcap.scaffold29.ace.AE.1.1
expID_and_template: 1 TPAA-ail08c06
}

CT{
Contig55.555 autoFinishExp autofinish 4000 4500 060831:122829
C
purpose: weak
0 915 5
dyeTerm customPrimer
fix cons errors: 4.69881 original cons errors: 5.64249
original single subclone bases: 886
primer: ggcaaatatggtgcaataaaac temp: 58 id: Trichinella_spiralis_060315.pcap.scaffold29.ace.AE.1.1
expID_and_template: 1 TPAA-ail08c06
COMMENT{
great rxn!
C}
}

/);

$io->seek(0, 0);

my $parser = Genome::Model::Tools::Pcap::TagParser->new();
ok($parser, 'create parser');
foreach my $type (qw/ comment oligo autoFinExp autoFinExp /)
{
    $io->getline;
    my $tag = $parser->parse($io);
    ok(defined $tag, "parsed type: $type");# . Dumper($tag));
}

done_testing();
exit;

