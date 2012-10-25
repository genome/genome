package Genome::DruggableGene::Command::Import::Template;

use strict;
use warnings;
use Genome;

class Genome::DruggableGene::Command::Import::Template {
    has => [
        foo => {
            is => 'Text',
            doc => 'Foo property',
        },
        bar => {
            is => 'Number',
            is_optional => 1,
            doc => 'Bar property, optional',
        },
         citation_base_url => {
            default => '<insert site url here>',
        },
        ciation_site_url => {
            default => '<insert site url here>',
        },
        citation_text => {
            default => '<Insert citation text here>',
        },
    ],
    doc => 'Class documentation',
};

sub help_detail {
    return 'Here is some detailed documentation about this command';
}

sub help_synopsis {
    return 'genome druggable-gene import template --foo blah --bar 4';
}

sub do_stuff {
    my $self = shift;
    return 1;
}

sub execute {
    my $self = shift;
    my $foo = $self->foo; # Get the value from foo, which you provided via the command line
    my $bar = $self->bar; # Get the value from bar, which you MAY have provided via the command line
    
    # The stuff you want executed goes in this method (or other methods you define).
    $self->do_stuff();

    # Retrieve a drug name. You can provide a hash of parameters for your query, like I have below. 
    # Internally, this is converted to SQL that hits the postgres database. You can query on any
    # of the properties that are defined in the Genome/DrugNameReport.pm module, so refer to that
    # for some help.
    my @retrieved_drug_names = Genome::DrugNameReport->get(
        name => 'whatever',
        nomenclature => 'Entrez',
        description => 'whatever',
    );
    
    # This creates a new drug name object with the values you specify. This is NOT committed to the
    # database yet, it only lives locally. A UR::Context->commit is required to have this object 
    # converted into a row in the database. This commit is done automatically once your command
    # successfully returns, though you can do it manually if you wish.
    my $new_drug_name = Genome::DrugNameReport->create(
        name => 'whatever',
        nomenclature => 'whatever',
        source_db_name => 'foo',
        source_db_version => 'bar',
    );

    # Once you are done, return 1 for success. 
    return 1;
}

# After your command has executed and returned, your changes will be committed to the database automatically.
# If you want to do this yourself manually, do 'UR::Context->commit;'
1;

