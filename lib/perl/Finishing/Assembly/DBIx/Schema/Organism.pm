package Finishing::Assembly::DBIx::Schema::Organism;

use strict;
use warnings;

use base 'DBIx::Class';

use Finfo::Logging 'fatal_msg';
use Finfo::Validate;

__PACKAGE__->load_components(qw/ Core /);
__PACKAGE__->table('organism');
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
      size => 50,
      is_nullable => 0,
    },
    common_name => 
    {
      data_type => 'varchar',
      size => 50,
      is_nullable => 0,
    },
);
__PACKAGE__->set_primary_key('id');
__PACKAGE__->has_many('chromosomes', 'Finishing::Assembly::DBIx::Schema::Chromosome', 'organism_id');
__PACKAGE__->has_many('assemblies', 'Finishing::Assembly::DBIx::Schema::Assembly', 'organism_id');

#- CHROMOSOME-#
sub get_chromosome
{
    my ($self, $first, $second) = @_;

    if (!$second or $first eq 'name') {
        my $name = $second ? $second : $first;

        Finfo::Validate->validate(
            attr => 'chromosome name',
            value => $name,
            isa => 'string',
            msg => 'fatal',
        );
        return $self->chromosomes->find({ name => $name });
    }
    elsif ($first eq 'id') {
        return $self->chromosomes->find({ id => $second });
    }
}

#- ASSEMBLY -#
sub get_assembly
{
    my ($self, $name) = @_;

    Finfo::Validate->validate
    (
        attr => 'assembly name',
        value => $name,
        isa => 'string',
        msg => 'fatal',
    );

    return $self->assemblies->find({ name => $name });
}

sub create_assembly
{
    my ($self, $name) = @_;

    my $existing_assembly = $self->get_assembly($name);

    $self->fatal_msg("Assembly ($name) already exists for " . $self->name) if $existing_assembly;
    
    return $self->_create_assembly($name);
}

sub get_or_create_assembly
{
    my ($self, $name) = @_;

    my $existing_assembly = $self->get_assembly($name);

    return $existing_assembly if $existing_assembly;

    return $self->_create_assembly($name);
}

sub _create_assembly
{
    my ($self, $name) = @_;

    my $assembly = $self->create_related
    (
        'assemblies',
        { name => $name },
    );

    $self->fatal_msg("Can't create assembly ($name) for " . $self->name) unless $assembly;
    
    return $assembly;
}

1;

#$HeadURL$
#$Id$
