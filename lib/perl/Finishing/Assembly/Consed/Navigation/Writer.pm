package Finishing::Assembly::Consed::Navigation::Writer;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finfo::Writer';

my %title :name(title:o)
    :default('Consed Navigation File')
    :clo('title=s')
    :desc('Title for navigation file');

sub START
{
    my $self = shift;

    $self->io->print
    (
        sprintf
        (
            "TITLE: %s\n\n",
            $self->title
        )
    );

    return 1;
}

sub _write_one
{
    my ($self, $nav) = @_;

    #TODO - accept nav as a hasref
    #$nav = Navigation::Consed->new(%$nav) unless ...;

    return unless Finfo::Validate->validate
    (
        attr => 'navigation',
        value => $nav,
        type => 'inherits_from',
        options => [qw/ Finishing::Assembly::Consed::Navigation /],
        err_cb => $self,
    );

    return $self->io->print
    (
        sprintf
        (
            "BEGIN_REGION\nTYPE: %s\n%sCONTIG: %s\n%sUNPADDED_CONS_POS: %d %d\nCOMMENT: %s\nEND_REGION\n\n", 
            $nav->type,
            ( $nav->acefile ) ? sprintf("ACEFILE: %s\n", $nav->acefile) : '',
            $nav->contig_name,
            ( $nav->type eq 'READ' ) ? sprintf("READ: %s\n", $nav->read_name) : '',
            $nav->start,
            $nav->stop,
            $nav->description
        )
    );
}

1;

=pod

=head1 Name

Finishing::Assembly::Consed::Navigation::Writer

=head1 Usage

 use Finishing::Assembly::Consed::Navigation::Writer;
 
 my $wirter = Finishing::Assembly::Consed::Navigation::Writer->new
 (
     io => 'out.nav', # required, file or IO::* object
     title => '', # optional, sting(new method writes the title to the io)
 )
     or die;
  
 $writer->write_one($nav)
    or die;

 or
 
 $writer->write_many(\@navs)
     or die;

=head1 Methods

=head2 write_one

 Writes a nav to the io
 
=head2 write_many

 Writes many navs to the io 
 
=head1 See Also

=over

=item Finishing::Assembly::Consed::Navigation directory

=item Finfo::Writer (parent class)

=item consed

=back
  
=head1 Disclaimer

Copyright (C) 2006-7 Washington University Genome Sequencing Center

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author

B<Eddie Belter> <ebelter@watson.wustl.edu>

=cut

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Consed/Navigation/Writer.pm $
#$Id: Writer.pm 28692 2007-09-26 22:33:30Z ebelter $
