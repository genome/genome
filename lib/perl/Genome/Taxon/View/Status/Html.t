#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above "Genome";

use Genome::Test::Factory::Individual;
use Genome::Test::Factory::Taxon;

use Test::More tests => 5;

use_ok('Genome::Taxon::View::Status::Html') or die "test cannot continue...";

my $subject = Genome::Test::Factory::Taxon->setup_object(
    name => 'Null cat',
    species_latin_name => 'Felis nullus',
    domain => 'Eukaryota',
    estimated_genome_size => '2000000000',
);

for(1..5) {
    Genome::Test::Factory::Individual->setup_object(
        taxon_id => $subject->id,
    );
}

isa_ok($subject, 'Genome::Taxon', "created taxon subject") or die "test cannot continue...";

my $view_obj = $subject->create_view(
    xsl_root => Genome->base_dir . '/xsl',
    rest_variable => '/cgi-bin/rest.cgi',
    toolkit => 'html',
    perspective => 'status',
); 
ok($view_obj, "created a view") or die "test cannot continue...";
isa_ok($view_obj, 'Genome::Taxon::View::Status::Html');

my $html = $view_obj->_generate_content();
ok($html, "view returns HTML") or die "test cannot continue...";

1;
