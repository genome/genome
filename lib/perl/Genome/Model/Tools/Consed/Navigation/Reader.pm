package Genome::Model::Tools::Consed::Navigation::Reader;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Consed::Navigation::Reader {
    is => 'Genome::Utility::IO::Reader',
};

sub title {
    return $_[0]->{_title};
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    my $title_line = $self->getline;
    chomp $title_line;
    $title_line =~ s#^TITLE:\s*##;
    
    $self->error_msg("No title in navigation io")
        and return unless $title_line;
    
    $self->{_title} = $title_line;
    
    return $self;
}

sub next {
    my $self = shift;

    my %nav;
    while ( my $line = $self->getline ) {
        chomp $line;
        next if $line eq '';
        my ($attr, $value) = split(/\:\s*/, $line);
        next unless $attr;
        last if $attr eq 'END_REGION';
        next if $attr eq 'BEGIN_REGION';

        if ( $attr eq 'UNPADDED_CONS_POS' ) {
            my ($start, $stop) = split(/\s+/, $value);
            $nav{start} = $start;
            $nav{stop} = $stop || $start;
        }
        else {
            $nav{ lc($attr) } = $value;
        }
    }

    return ( %nav ) ? \%nav : undef;
}

1;

=pod

=head1 Name

Genome::Model::Tools::Consed::Navigation::Reader

=head1 Synopsis

A stream based reader, parses a consed navigation file, returning navigation hashrefs.

=head1 Usage

 use Genome::Model::Tools::Consed::Navigation::Reader;

 my $reader = Genome::Model::Tools::Consed::Navigation::Reader->new(
     input => 'cafcop.nav', # REQUIRED, file or an object that can 'getline' and 'seek'
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

=head2 title

 my $title = $reader->title;

=over

=item I<Synopsis>   Gets the title of the navigation input, (sans 'TITLE:')

=item I<Params>     none

=item I<Returns>    string 

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
