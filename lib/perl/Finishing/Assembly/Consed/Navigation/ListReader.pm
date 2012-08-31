package Finishing::Assembly::Consed::Navigation::ListReader;

use strict;
use warnings;

use base 'Finfo::Reader';

use Finishing::Assembly::Consed::Navigation;

sub _return_class
{
    return 'Finishing::Assembly::Consed::Navigation';
}

sub _next
{
    my $self = shift;
    
    #Contig19    (consensus)           10135-10166   base quality below threshold
    #Contig22    L24337P6011D9.b1      11132         confirmed by: L24337P6006F12.b1
    my $line = $self->_getline;
    return unless $line;
    chomp $line;
    return if $line =~ /^\s*$/;

    my @tokens = split(/\s+/, $line, 4);

    $self->error_msg("Error parsing line:\n$line")
        and return unless @tokens;

    my %nav;
    @nav{qw/ contig_name type pos description /} = @tokens;

    my $pos = delete $nav{pos};
    ($nav{start}, $nav{stop}) = split(/\-/, $pos);
    $nav{stop} = $nav{start} unless defined $nav{stop};

    if ( $nav{type} eq '(consensus)' )
    {
        $nav{type} = 'CONSENSUS';
    }
    else
    {
        $nav{read_name} = $nav{type};
        $nav{type} = 'READ';
    }

    return \%nav;
}

1;

=pod

=head1 Name

Finishing::Assembly::Navigation::ListReader

=head1 Synopsis

Reads a consed 'list' file, returns navigation hashrefs or objects

=head1 Usage

 use Finishing::Assembly::Consed::Navigation::ListReader
 use IO::File;
 
 my $fh = IO::File->new("< $infile")
     or die "$!\n";
 my $reader = Finishing::Assembly::Consed::Navigation:ListReader->new
 (
     io => $io, # file or IO::* object
     return_as_objects => 1,
 )
     or die;
 
 while ( my $nav = $reader->next )
 {
     ...
 }

=head1 Methods

=head2 next

 Returns the next Navigation objects from the io
 
=head2 all

 Returns all of the Navigation objects from the .IO
 
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

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Consed/Navigation/ListReader.pm $
#$Id: ListReader.pm 28690 2007-09-26 22:32:21Z ebelter $
