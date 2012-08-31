#!/usr/bin/env genome-perl

use strict;
use warnings;
use Data::Dumper;
use Finishing::Assembly::DBIx::Schema;

my $schema = Finishing::Assembly::DBIx::Schema->connect(
    'dbi:mysql:cmap_sorting:mysql1', 'sorting_admin', 's_admin', 
);

$schema->txn_do(
    sub{
        my $tag_i = $schema->resultset('ConsensusTag');
        foreach (1494, 1551, 1608) {
            my $tag = $tag_i->find( {id => $_} );
            $tag->set_columns({text => undef} );
            $tag->update;
            print $tag->text;
        }
    }
);


=pod

=head1 NAME
ScriptTemplate - template for new perl script

=head1 SYNOPSIS

=head1 DESCRIPTION 

=cut

#$HeadURL$
#$Id$


