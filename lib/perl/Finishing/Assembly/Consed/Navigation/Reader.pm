package Finishing::Assembly::Consed::Navigation::Reader;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finfo::Reader';

use Finishing::Assembly::Consed::Navigation;

my %title :name(_title:p);

sub _return_class
{
    return 'Finishing::Assembly::Consed::Navigation';
}   

sub title
{
    return shift->_title;
}

sub START
{
    my $self = shift;

    my $title_line = $self->_getline;
    chomp $title_line;
    my ($title) = $title_line =~ /^TITLE:\s*(.+)/;
    
    $self->fatal_msg("No title in navigation io")
        and return unless $title;
    
    $self->_title($title);
    
    return 1;
}

sub _next
{
    my $self = shift;

    my %nav;
    while ( my $line = $self->_getline )
    {
        chomp $line;
        my ($line_header, $data) = $line =~ /^([\w_]+):?\s?(.+)?$/;
        next unless $line_header;
        last if $line_header eq 'END_REGION';
        next if $line_header eq 'BEGIN_REGION';

        #print "$line_header || $data\n";

        if ($line_header eq 'TYPE')
        {
            $nav{type} = $data;
        }
        elsif ($line =~ 'CONTIG')
        {
            $nav{contig_name} = $data;
        }
        elsif ($line =~ 'READ')
        {
            $nav{read_name} = $data;
        }
        elsif ($line =~ 'UNPADDED_CONS_POS')
        {
            my ($start, $stop) = split(/\s+/, $data);
            $nav{start} = $start;
            $nav{stop} = $stop || $start;
        }
        elsif ($line =~ 'COMMENT')
        {
            $nav{description} = $data;
        }
        elsif ($line =~ 'ACEFILE')
        {
            $nav{acefile} = $data;
        }
    }

    return unless %nav;

    return \%nav;
}

1;

=pod

=head1 Name

 Finishing::Assembly::Consed::Navigation::Reader

=head1 Synopsis

Reads a consed navigation file, returns navigation hashrefs or objects

=head1 Usage

 use Finishing::Assembly::Consed::Navigation::Reader;

 my $reader = Finishing::Assembly::Consed::Navigation::Reader->new
 (
     io => 'cafcop.nav', # can be file, or IO::* object
     return_as_objects => 1, # default is 0 returns hashrefs, 1 returns objects
 )
     or die;
 
 while ( my $nav = $reader->next )
 {
     ...
 }

=head1 Methods

=head2 title

 Return the title string (sans 'TITLE:' part)

=head2 next

 Returns the next navigation objects from the IO
 
=head2 all

 Returns all of the navigation objects from the IO
 
=head1 See Also

=over

=item Finishing::Assembly::Consed::Navigation directory

=item Finfo::Reader (parent class)

=item consed

=back

=head1 Disclaimer

Copyright (C) 2006-7 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details.

=head1 Author

B<Eddie Belter> <ebelter@watson.wustl.edu>

=cut

#$Header$
#$Id: Reader.pm 28691 2007-09-26 22:33:02Z ebelter $
