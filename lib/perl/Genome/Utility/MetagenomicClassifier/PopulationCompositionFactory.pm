package Genome::Utility::MetagenomicClassifier::PopulationCompositionFactory;

use strict;
use warnings;

use base 'Class::Singleton';

require Bio::SeqIO;
use Carp 'confess';
use Data::Dumper 'Dumper';
require Genome::Sys;
require Genome::Utility::MetagenomicClassifier::PopulationComposition;

sub get_composition {
    my ($self, %params) = @_;

    $self = __PACKAGE__->instance unless ref $self;

    for my $req (qw/ classifier fasta_file /) {
        _fatal_message("Required parameter ($req) not found.") unless $params{$req};
    }

    Genome::Sys->validate_file_for_reading($params{fasta_file})
        or _fatal_message( sprintf("Can't validate fasta file (%s)", $params{fasta_file}) );

    my $bioseq_io = Bio::SeqIO->new(
        '-file' => '<'.$params{fasta_file},
        '-format' => 'fasta',
    );
    
    my $population_composition = Genome::Utility::MetagenomicClassifier::PopulationComposition->new(
        confidence_threshold => $params{confidence_threshold},
    );

    while ( my $seq = $bioseq_io->next_seq ) {
        $population_composition->add_classification(
            $params{classifier}->classify($seq)
        );
    }

    return $population_composition;
}

#< PRIVATE >#
sub _fatal_message {
    my ($msg) = @_;

    confess __PACKAGE__." ERROR: $msg\n";
}
 
1;

=pod

=head1 Name

ModuleTemplate

=head1 Synopsis

=head1 Usage

=head1 Methods

=head2 

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head1 See Also

=head1 Disclaimer

Copyright (C) 2005 - 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$

