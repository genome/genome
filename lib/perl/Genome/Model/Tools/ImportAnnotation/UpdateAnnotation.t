
use strict;
use warnings;

use Test::More tests => 2;
use above 'Genome';

BEGIN {
    use_ok("Genome::Model::Tools::ImportAnnotation::UpdateAnnotation");
};


my $u = Genome::Model::Tools::ImportAnnotation::UpdateAnnotation->create( version => '45_36g');

my $path = $u->derive_flatfile_location();
#print $path,"\n";
is($path, '/gsc/var/lib/import/entrez/45_36g/Homo_sapiens.agc', 
    'paths match for 45_36g');
