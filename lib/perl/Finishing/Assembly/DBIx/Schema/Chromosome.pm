package Finishing::Assembly::DBIx::Schema::Chromosome;

use strict;
use warnings;

use base 'DBIx::Class';

use Finfo::Logging 'fatal_msg';

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('chromosome');
__PACKAGE__->add_columns
(
    id => 
    {
        data_type => 'int',
        size => 11,
        is_auto_increment => 1,
        is_nullable => 0,
    },
    name => 
    {
        data_type => 'varchar',
        size => 10,
        is_nullable => 0,
    },
    organism_id =>
    {
        data_type => 'int',
        size => 11,
        is_foreign_key => 1,
        is_nullable => 0,
    },
);
__PACKAGE__->set_primary_key('id');
__PACKAGE__->add_unique_constraint([qw/ organism_id name /]);
__PACKAGE__->belongs_to('organism', 'Finishing::Assembly::DBIx::Schema::Organism', 'organism_id');
__PACKAGE__->has_many
(
    'first_scaffolds',
    'Finishing::Assembly::DBIx::Schema::ChromosomeFirstScaffold', 
    'chromosome_id'
);

sub first_scaffold_for_assembly
{
    my ($self, $assembly) = @_;

    Finfo::Validate->validate
    (
        attr => 'assembly to get first scaffold',
        value => $assembly,
        isa => 'object Finishing::Assembly::Assembly Finishing::Assembly::DBIx::Schema::Assembly',
        msg => 'fatal',
    );

    return $self->first_scaffolds->find
    (
        {
            'assembly_id' => $assembly->id,
        },
    );
}

sub set_first_scaffold
{
    my ($self, $scaffold) = @_;

    Finfo::Validate->validate
    (
        attr => 'scaffold to set first on chromosome',
        value => $scaffold,
        isa => 'object Finishing::Assembly::Scaffold Finishing::Assembly::DBIx::Schema::Scaffold',
        msg => 'fatal',
    );

    my $assembly = $scaffold->assembly;
    my $existing_first_scaffold = $self->first_scaffold_for_assembly($assembly);
    
    $self->fatal_msg
    (
        sprintf
        (
            "First scaffold already exists for assembly (%s) on chromosome (%s)",
            $assembly->name,
            $self->name,
        )
    ) if $existing_first_scaffold;
    
    my $fs = $self->create_related
    (
        'first_scaffolds',
        {
            assembly_id => $assembly->id,
            scaffold_id => $scaffold->id,
        },
    );

    $self->fatal_msg
    (
        sprintf
        (
            "Can't create first scaffold for assembly (%s) on chromosome (%s)",
            $assembly->name,
            $self->name,
        )
    ) unless $fs;

    return $fs;
}

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Schema/Chromosome.pm $
#$Id: Chromosome.pm 31302 2007-12-26 21:55:55Z ebelter $
