package Finishing::Assembly::DBIx::Schema::Assembly;

use strict;
use warnings;

use base 'DBIx::Class';

use Finfo::Logging 'fatal_msg';

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('assembly');
__PACKAGE__->add_columns
(
    id => 
    {
      'data_type' => 'int',
      'size' => 11,
      'is_auto_increment' => 1,
      'is_nullable' => 0,
    },
    organism_id => 
    {
      data_type => 'int',
      size => 11,
      is_auto_increment => 0,
      is_foreign_key => 1,
      is_nullable => 0,
    },
    'name' => 
    {
      'data_type' => 'varchar',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'is_nullable' => 0,
      'size' => '50'
    },
    'file_path' => 
    {
      'data_type' => 'varchar',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'is_nullable' => 1,
      'size' => '150'
    },
);
__PACKAGE__->set_primary_key('id');
__PACKAGE__->add_unique_constraint([qw/ organism_id name /]);
__PACKAGE__->belongs_to('organism', 'Finishing::Assembly::DBIx::Schema::Organism', 'organism_id');
__PACKAGE__->has_many('scaffolds', 'Finishing::Assembly::DBIx::Schema::Scaffold', 'assembly_id');
__PACKAGE__->has_many('contigs', 'Finishing::Assembly::DBIx::Schema::Contig', 'assembly_id');
__PACKAGE__->has_many('assembled_reads', 'Finishing::Assembly::DBIx::Schema::AssembledRead', 'assembly_id');
__PACKAGE__->has_many('tags', 'Finishing::Assembly::DBIx::Schema::AssemblyTag', 'assembly_id');
__PACKAGE__->has_many('improvement_correlations', 'Finishing::Assembly::DBIx::Schema::ImprovementCorrelation', 'assembly_id');

sub get_scaffold
{
    my ($self, $scaffold_num) = @_;

    Finfo::Validate->validate
    (
        attr => 'scaffold number',
        value => $scaffold_num,
        isa => 'int',
        msg => 'fatal',
    );
    
    return $self->scaffolds->find({ scaffold_num => $scaffold_num });
}

sub get_contig
{
    my ($self, $contig_name) = @_;

    $self->fatal_msg
    (
        "Unknown contig name format ($contig_name)"
    ) unless $contig_name =~ /^Contig(\d+)\.(\d+)$/;

    my $scaffold_num = $1;
    my $contig_num = $2;

    my $scaffold = $self->get_scaffold($scaffold_num);
    
    return unless $scaffold;

    return $scaffold->get_contig($contig_num);
}

sub get_assembled_read
{
    my ($self, $read_name) = @_;

    Finfo::Validate->validate
    (
        attr => 'read name',
        value => $read_name,
        isa => 'string',
        msg => 'fatal',
    );

    return $self->assembled_reads->find({ name => $read_name });
}

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Schema/Assembly.pm $
#$Id: Assembly.pm 31277 2007-12-24 17:14:58Z ebelter $
