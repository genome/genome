package Finishing::Assembly::DBIx::Schema::ImprovementCorrelation;

use base 'DBIx::Class';

use strict;
use warnings;

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('improvement_correlation');
__PACKAGE__->add_columns(
    'id' => {
      'data_type' => 'int',
      'is_auto_increment' => 1,
      'is_nullable' => 0,
      'size' => '11'
    },
    'correlation_type' => {
      'data_type' => 'varchar',
      'is_nullable' => 0,
      'size' => '20'
    },
    'assembly_id' => {
      'data_type' => 'int',
      'is_foreign_key' => 1,
      'is_nullable' => 0,
      'size' => '11'
    },
    'score' => {
      'data_type' => 'float',
      'size' => '8',
      'default_value' => '0',
      'is_nullable' => 1
    },
);
__PACKAGE__->set_primary_key('id');
__PACKAGE__->belongs_to('assembly', 'Finishing::Assembly::DBIx::Schema::Assembly', 'assembly_id');
__PACKAGE__->has_many('correlation_contigs', 'Finishing::Assembly::DBIx::Schema::CorrelationContig', 'correlation_id');
__PACKAGE__->many_to_many('contigs', 'correlation_contigs', 'contig');

sub extended_contigs {
    my $self = shift;

    my $ranges = $self->contigs->search(
        {
        },
        {
            select => [
            'scaffold_id',
            {min => 'contig_num'},
            {max => 'contig_num'},
            ],
            as => ['scaffold_id', 'first_contig','last_contig'],
            group_by => ['scaffold_id'],
        },
    );

    my @where = ();
    while (my $scaffold_range = $ranges->next) {
        push @where, 
        {
            scaffold_id => $scaffold_range->get_column('scaffold_id'),
            contig_num => {
                between => [$scaffold_range->get_column('first_contig')-1,$scaffold_range->get_column('last_contig')+1],
            },
        };
    }
    return $self->result_source->schema->resultset('Contig')->search(
        \@where,
    );
}

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Schema/ImprovementCorrelation.pm $
#$Id: ImprovementCorrelation.pm 31273 2007-12-24 17:06:24Z ebelter $
