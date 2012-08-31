package Finishing::Assembly::DBIx::Schema::Scaffold;

use strict;
use warnings;

use base 'DBIx::Class';

use Finfo::Logging 'fatal_msg';

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('scaffold');
__PACKAGE__->add_columns
(
    id => 
    {
      'data_type' => 'int',
      'is_auto_increment' => 1,
      'is_foreign_key' => 0,
      'name' => 'id',
      'is_nullable' => 0,
      'size' => '11'
    },
    assembly_id => 
    {
        data_type => 'int',
        is_foreign_key => 1,
        is_nullable => 0,
        size => 11,
    },
    scaffold_num => 
    {
        'data_type' => 'int',
        'is_nullable' => 0,
        'size' => '11'
    },
    'length' => 
    {
        data_type => 'int',
        is_nullable => 1,
        size => 11,
    },
    orientation =>
    {
        data_type => 'enum',
        extra => { list => [qw/ + - /] },
        default_value => '+',
        is_nullable => 1,
    },
);
__PACKAGE__->set_primary_key('id');
__PACKAGE__->add_unique_constraint([qw/ assembly_id scaffold_num /]);
__PACKAGE__->belongs_to('assembly', 'Finishing::Assembly::DBIx::Schema::Assembly','assembly_id');
__PACKAGE__->has_many
(
    'contigs', 
    'Finishing::Assembly::DBIx::Schema::Contig', 
    'scaffold_id',
    {
        proxy => [qw/assembled_reads/],
    }
);
__PACKAGE__->has_one
(
    'chromosomes_first',
    'Finishing::Assembly::DBIx::Schema::ChromosomeFirstScaffold', 
    'scaffold_id'
);

sub get_contig
{
    my ($self, $contig_num) = @_;

    Finfo::Validate->validate
    (
        attr => 'contig number',
        value => $contig_num,
        isa => 'int',
        msg => 'fatal',
    );

    return $self->contigs->find({ contig_num => $contig_num });
}

sub ordered_contigs
{
    my  $self = shift;

    my $contig_rs = ( $self->orientation eq '+' )
    ? $self->contigs->search({},{ order_by => 'contig_num asc' })
    : $self->contigs->search({},{ order_by => 'contig_num desc' });
    
    return ( wantarray ) ? $contig_rs->all : $contig_rs;
}

sub first_contig
{
    my $self = shift;

    my $contig_rs = $self->ordered_contigs;

    return $contig_rs->first;
}

sub ultra_links {
    my $self = shift;
    my $schema = $self->result_source->schema;
    my $ultra_link = $schema->resultset('TemplateLink')->search(
        {
            'left_contig.scaffold_id' => $self->id,
            'right_contig.scaffold_id' => {'<>', $self->id},
        },
        {
            join => ['left_contig','right_contig'],
        },
    );
    return $ultra_link;
}

sub scaffold_links {
    my $self = shift;
    my $schema = $self->result_source->schema;
    my $scaffold_link = $schema->resultset('TemplateLink')->search(
        {
            'left_contig.scaffold_id' => $self->id,
            'right_contig.scaffold_id' => $self->id,
        },
        {
            join => ['left_contig','right_contig'],
        },
    );
    return $scaffold_link;
}

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Schema/Scaffold.pm $
#$Id: Scaffold.pm 31276 2007-12-24 17:12:26Z ebelter $
