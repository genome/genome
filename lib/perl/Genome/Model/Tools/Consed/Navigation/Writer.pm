package Genome::Model::Tools::Consed::Navigation::Writer;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Consed::Navigation::Writer {
    is => 'Genome::Utility::IO::Writer',
    has => [
    title => {
        is => 'Text',
        doc => 'Navigation title',
    },
    ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_) 
        or return;
    
    unless ( $self->title ) {
        $self->error_message("Need navigation title");
        $self->delete;
        return;
    }

    $self->print('TITLE: '.$self->title."\n\n");

    return $self;
}

sub write_one {
    my ($self, $nav) = @_;

    return $self->print(
        sprintf(
            "BEGIN_REGION\nTYPE: %s\n%sCONTIG: %s\n%sUNPADDED_CONS_POS: %d %d\nCOMMENT: %s\nEND_REGION\n\n", 
            uc($nav->{type}),
            (( $nav->{acefile} ) ? sprintf("ACEFILE: %s\n", $nav->{acefile}) : ''),
            $nav->{contig},
            (( $nav->{type} eq 'READ' ) ? sprintf("READ: %s\n", $nav->{read}) : ''),
            $nav->{start},
            $nav->{stop},
            $nav->{comment},
        )
    );
}

1;

=pod

=head1 Name

Genome::Model::Tools::Consed::Navigation::Writer

=head1 Usage

 use Genome::Model::Tools::Consed::Navigation::Writer;
 
 my $writer = Genome::Model::Tools::Consed::Navigation::Writer->new(
     output => 'repeats.nav', # REQUIRED, file or an object that can 'print'
     title => 'Repeat Found', # REQUIRED, string
 )
     or die;
  
 $writer->write_one($nav)
    or die;

=head1 Methods

=head2 write_one

 $writer->write_one($nav)
    or die;

=over

=item I<Synopsis>   Writes the navigation to the output.

=item I<Params>     navigation (hashref)

=item I<Returns>    the return value of the write_one method (boolean)

=back

=head1 See Also

I<Genome::Model::Tools::Consed directory>, I<consed>

=head1 Disclaimer

Copyright (C) 2006-8 Washington University Genome Sequencing Center

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author

B<Eddie Belter> <ebelter@genome.wustl.edu>

=cut

#$HeadURL$
#$Id$
