package Finishing::Assembly::Consed::Navigation::ConvertFromList;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

my %reader :name(reader:r)
    :isa('object Finishing::Assembly::Consed::Navigation::ListReader');
my %writer :name(writer:r)
    :isa('object Finishing::Assembly::Consed::Navigation::Writer');

sub execute
{
    my $self = shift;

    my $count = 0;
    while ( my $nav = $self->reader->next )
    {
        $self->writer->write_one($nav)
            or return; # unlink file?
        $count++;
    }

    $self->info_msg("Done. Wrote $count navs");
    
    return 1;
}
1;

=pod

=head1 Name

Finishing::Assembly::Consed::Navigation::ConvertFromList

=head1 Synopsis

Takes navigations from a consed list reader, then writes them as consed navigations

=head1 Usage

 use Finishing::Assembly::Consed::Navigation::ConvertFromList;
 use Finishing::Assembly::Consed::Navigation::ListReader;
 use Finishing::Assembly::Consed::Navigation::Writer;

 my $reader = Finishing::Assembly::Consed::Navigation::ListReader->new
 (
     io => 'list.txt', # file or IO::* object
     return_as_objs => 1,
 )
     or die;
 
 my $writer = Finishing::Assembly::Consed::Navigation::Writer->new
 (
     io => 'list.txt', # file or IO::* object
     title => 'List Converted to Navs', # optional
 )
     or die;
 
 my $converter = Finishing::Assembly::Consed::Navigation::ConvertFromList->new
 (
     reader => $reader,
     writer => $writer,
 )
     or die;
 
 $converter->execute
     or die;

=head1 Methods

=head2 execute

The main method.

=head1 See Also

=over

=item Finishing::Assembly::Consed::Navigation directory

=item consed

=item list2nav script

=back

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Consed/Navigation/ConvertFromList.pm $
#$Id: ConvertFromList.pm 29586 2007-10-29 15:46:09Z ebelter $

