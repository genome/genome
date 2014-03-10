#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Site::TGI::Synchronize::DiffLimsAndGenome') or die;

class Iterator {
    has => [
        objects => { is_many => 1, },
        position => { is => 'Integer', default_value => -1, },
    ],
};
sub Iterator::_inc_position {
    my $self = shift;
    return $self->position( $self->position + 1 );
};
sub Iterator::next {
    my $self = shift;
    my $position = $self->_inc_position;
    return ($self->objects)[$position];
};

my $entity_name = 'taxon';
my $lims_class = Genome::Site::TGI::Synchronize::Classes::Dictionary->lims_class_for_entity_name($entity_name);
my @lims_objects = ( 
    $lims_class->__define__(id => -11, name => '__TEST_TAXON__'),
    $lims_class->__define__(id => -13, name => '__TEST_TAXON__'),
);
my $genome_class = $lims_class->genome_class_for_comparison;
my @genome_objects = (
    $genome_class->__define__(id => -12, name => '__TEST_TAXON__'),
    $genome_class->__define__(id => -13, name => '__TEST_TAXON__'),
);
no strict;
*{$lims_class  .'::create_iterator'} = sub{ return Iterator->create(objects => \@lims_objects); };
*{$genome_class.'::create_iterator'} = sub{ return Iterator->create(objects => \@genome_objects); };
use strict;

my $differ = Genome::Site::TGI::Synchronize::DiffLimsAndGenome->create;
ok($differ, 'create differ w/o entity name');
ok(!$differ->execute, 'failed to execute w/o entity_name');

$differ->is_executed(undef);
is($differ->entity_name('taxon'), 'taxon', 'set entity_name');
ok($differ->lims_class, 'lims class');

ok($differ->execute, 'execute');
is_deeply([@{$differ->in_lims_not_genome}], [-11], $entity_name."s in lims not genome");
is_deeply([@{$differ->in_genome_not_lims}], [-12], $entity_name."s in genome not lims");

done_testing();
