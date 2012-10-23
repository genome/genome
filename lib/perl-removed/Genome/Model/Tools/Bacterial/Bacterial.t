
use strict;
use warnings;

use above "Genome";
use Test::More tests => 4;


BEGIN {
    use_ok('Genome::Model::Tools::Bacterial::DumpSequences');
    use_ok('Genome::Model::Tools::Bacterial::ParseAceFiles');
    use_ok('Genome::Model::Tools::Bacterial::TagOverlaps');
    use_ok('Genome::Model::Tools::Bacterial::AceDumpGenes');
}

