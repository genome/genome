package Genome::Model::Tools::Consed::ListToNav;

use strict;
use warnings;

use Genome;

use File::Basename 'basename';

class Genome::Model::Tools::Consed::ListToNav {
    is => 'Command',
    has => [
    list => {
        is => 'Text', 
        doc => 'List input (file from command line) for reading',
    },
    nav => {
        is => 'Text',
        doc => 'Navigation output (file from command line) for writing'
    },
    ],
    has_optional => [
    title => {
        is => 'Text',
        doc => 'Title for nav file',
        default => 'Converted from list',
    },
    ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    $self->{_reader} = Genome::Model::Tools::Consed::Navigation::ListReader->create(
        input => $self->list,
    )
        or return;

    $self->{_writer} = Genome::Model::Tools::Consed::Navigation::Writer->create(
        output => $self->nav,
        title => $self->title,
    )
        or return;

    return $self;
}

sub execute {
    my $self = shift;

    my $count = 0;
    while ( my $nav = $self->{_reader}->next ) {
        $self->{_writer}->write_one($nav)
            or return;
        $count++;
    }

    return 1 if $count;

    $self->error_message('Nothing found in list input');
    return;
}

1;

=pod

=head1 Name

Genome::Model::Tools::Consed::ListToNav

=head1 Synopsis

Converts a consed list input to a consed navigation output.

=head1 Usage

 use Genome::Model::Tools::Consed::ListToNavigation;
 
 my $list2nav = Genome::Model::Tools::Consed::ListToNav->create(
     list => 'repeats.list', # file, or object that can 'getline' and 'seek'
     nav => 'repeats.nav', # file, or object that can 'print'
 )
     or die;
 
 $list2nav->execute
     or die;

=head1 Methods

=head2 execute

 $reader->execute
    or die;

=over

=item I<Synopsis>   Converts the list to nav

=item I<Params>     none

=item I<Returns>    boolean

=back

=head1 See Also

I<Genome::Model::Tools::Consed::Navigation directory>, I<Genome::Model::Tools directory>, I<consed>

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$

