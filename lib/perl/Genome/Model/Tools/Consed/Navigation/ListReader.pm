package Genome::Model::Tools::Consed::Navigation::ListReader;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Consed::Navigation::ListReader {
    is => 'Genome::Utility::IO::Reader',
};

# Ex of a list
#-contig    -consensus/read        -positions    -comment
#Contig19    (consensus)           10135-10166   base quality below threshold
#Contig22    L24337P6011D9.b1      11132         confirmed by: L24337P6006F12.b1

sub next {
    my $self = shift;

    my $line = $self->getline
        or return;
    chomp $line;
    return if $line eq '';

    my @tokens = split(/\s+/, $line, 4);
    $self->error_msg("Error parsing line:\n$line")
        and return unless @tokens and @tokens == 4;

    my %nav;
    @nav{qw/ contig type pos comment /} = @tokens;
    @nav{qw/ start stop/} = split(/\-/, delete $nav{'pos'});
    $nav{stop} = $nav{start} unless $nav{stop};

    if ( $nav{type} eq '(consensus)' ) {
        $nav{type} = 'CONSENSUS';
    }
    else {
        $nav{'read'} = $nav{type};
        $nav{type} = 'READ';
    }

    return \%nav;
}

1;

=pod

=head1 Name

Genome::Model::Tools::Consed::Navigation::ListReader;

=head1 Synopsis

A stream based reader, parses a consed 'list' file, returning navigation hashrefs.

=head1 Usage

 use Genome::Model::Tools::Consed::Navigation::ListReader;
 
 my $reader = Genome::Model::Tools::Consed::Navigation:ListReader->new(
     input => 'list.txt', # file or object that can 'getline' and 'seek'
 )
     or die;
 
 while ( my $nav = $reader->next ) {
     ...
 }

=head1 Methods

=head2 next

 my $nav = $reader->next;

=over

=item I<Synopsis>   Gets the next navigation form the input.

=item I<Params>     none

=item I<Returns>    navigation (hashref)

=back

=head2 all

 my @navs = $reader->all;

=over

=item I<Synopsis>   Gets all the navigations from the input.

=item I<Params>     none

=item I<Returns>    array (hashrefs)

=back

=head1 See Also

I<Genome::Model::Tools::Consed directory>, I<consed>

=head1 Disclaimer

Copyright (C) 2006-8 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author

B<Eddie Belter> <ebelter@genome.wustl.edu>

=cut

#$HeadURL$
#$Id$
