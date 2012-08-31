package Genome::Taxon::Command;

use strict;
use warnings;

use Genome;
      
class Genome::Taxon::Command {
    is => 'Command::Tree',
    doc => 'work with taxons',
};

use Genome::Command::Crud;
Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::Taxon',
    target_name => 'taxon',
    list => { show => 'id,name,species_latin_name,ncbi_taxon_id,locus_tag,domain', },
    update => { only_if_null => 1, },
    delete => { do_not_init => 1, },
);

1;

