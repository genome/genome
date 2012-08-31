package Finishing::Assembly::DBIx::Schema::Contig;

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components(qw/ Core /);
__PACKAGE__->table('assembly.contig');
__PACKAGE__->add_columns(
    'id' => {
      'data_type' => 'int',
      'is_auto_increment' => 1,
      'is_nullable' => 0,
      'size' => '11'
    },
    chromosome_id =>
    {
      data_type => 'int',
      size => 11,
      is_foreign_key => 1,
      is_nullable => 1,
    },
    ordered_on_chromosome =>
    {
      data_type => 'int',
      size => 1,
      is_nullable => 1,
    },
    'assembly_id' => {
      'data_type' => 'int',
      'is_foreign_key' => 1,
      'is_nullable' => 0,
      'size' => '11'
    },
    'scaffold_id' => {
      'data_type' => 'int',
      'is_foreign_key' => 1,
      'is_nullable' => 1,
      'size' => '11'
    },
    'contig_num' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'name' => 'contig_num',
      'is_nullable' => 0,
      'size' => '11'
    },
    'length' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'length',
      'is_nullable' => 1,
      'size' => '11'
    },
);
__PACKAGE__->set_primary_key('id');
__PACKAGE__->add_unique_constraint([qw/ scaffold_id contig_num /]);
__PACKAGE__->belongs_to('chromosome', 'Finishing::Assembly::DBIx::Schema::Chromosome', 'chromosome_id');
__PACKAGE__->belongs_to('assembly', 'Finishing::Assembly::DBIx::Schema::Assembly', 'assembly_id');
__PACKAGE__->belongs_to('scaffold', 'Finishing::Assembly::DBIx::Schema::Scaffold', 'scaffold_id');
__PACKAGE__->has_one
(
    'consensus', 'Finishing::Assembly::DBIx::Schema::ConsensusSequence', 'contig_id'
);
__PACKAGE__->has_many('right_gaps', 'Finishing::Assembly::DBIx::Schema::Gap', 'left_contig_id');
__PACKAGE__->has_many('left_gaps', 'Finishing::Assembly::DBIx::Schema::Gap', 'right_contig_id');
__PACKAGE__->has_many
(
    'assembled_reads', 'Finishing::Assembly::DBIx::Schema::AssembledRead', 'contig_id'
);
__PACKAGE__->has_many('tags', 'Finishing::Assembly::DBIx::Schema::ConsensusTag', 'contig_id');
__PACKAGE__->has_many
(
    'left_links',
    'Finishing::Assembly::DBIx::Schema::TemplateLink',
    {
        'foreign.left_contig_id' => 'self.id'
    },
);
__PACKAGE__->has_many
(
    'right_links',
    'Finishing::Assembly::DBIx::Schema::TemplateLink',
    {
        'foreign.right_contig_id' => 'self.id'
    },
);
__PACKAGE__->has_many
(
    'linking_reads', 
    'Finishing::Assembly::DBIx::Schema::AssembledRead',
    'contig_id',
    {
        join => 'mate',
        where => {'mate.contig_id' => \"<> me.contig_id"},
    }
);
__PACKAGE__->has_many
(
    'correlation_contigs', 'Finishing::Assembly::DBIx::Schema::CorrelationContig', 'contig_id',
);
__PACKAGE__->many_to_many('correlations', 'correlation_contigs', 'correlation');

sub name 
{
    my $self = shift;

    return sprintf('Contig%d.%d', $self->scaffold->scaffold_num, $self->contig_num);
}

sub scaffold_position 
{
    my $self = shift;
    my $schema = $self->result_source->schema;
    my $scaffold = $self->scaffold;

    my $rs = $schema->resultset('Contig')->search(
        {
            scaffold_id => $scaffold->id,
            contig_num => {'<' => $self->contig_num},
        },
    );
    my $position = $rs->get_column('length')->sum;
 
    return $position ? $position : 0;
}

sub scaffold_links {
    my $self = shift;
    return $self->left_links(
        {
            'right_contig.scaffold_id' => $self->scaffold_id,
        },
        {
            join => 'right_contig',
        },
    ); 
}

sub ultra_links {
    my $self = shift;
    return $self->left_links(
        {
            'right_contig.scaffold_id' => {'<>', $self->scaffold_id},
        },
        {
            join => 'right_contig',
        },
    ); 
}

sub left_contig
{
    my $self = shift;

    my ($gap) = $self->left_gaps;
    return unless $gap;
    
    my $contig = $gap->left_contig;
    return unless $contig;

    return $contig;
}

sub right_contig
{
    my $self = shift;

    my ($gap) = $self->right_gaps;
    return unless $gap;
    
    my $contig = $gap->right_contig;
    return unless $contig;

    return $contig;
}

#- ASSEMBLED READS -#
sub get_assembled_read
{
    my ($self, $name) = @_;

    Finfo::Validate->validate
    (
        attr => 'assembled read name',
        value => $name,
        isa => 'string',
        msg => 'fatal',
    );

    return $self->assembled_reads->find({ name => $name });
}

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Schema/Contig.pm $
#$Id: Contig.pm 31275 2007-12-24 17:11:33Z ebelter $
