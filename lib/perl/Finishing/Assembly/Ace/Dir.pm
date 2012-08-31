package Finishing::Assembly::Ace::Dir;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;
use File::Basename;
use File::Glob;
use File::Spec;
use Tie::File;

my %dir :name(dir:r)
    :isa(dir_rw)
    :access(ro)
    :clo('acedir=s') 
    :desc('Directory of acefiles');

sub acefile_for_ace
{
    my ($self, $ace) = @_;

    return unless Finfo::Validate->validate
    (
        attr => 'ace',
        value => $ace,
        msg => 'fatal',
    );
    
    return File::Spec->catfile($self->dir, $ace);
}
    
sub get_and_validate_acefile_for_ace
{
    my ($self, $ace) = @_;

    my $acefile = $self->acefile_for_ace($ace);

    return unless $acefile;
    
    return unless Finfo::Validate->validate
    (
        attr => 'acefile',
        value => $acefile,
        isa => 'file_rw',
        msg => 'fatal',
    );

    return $acefile;
}

sub next_iteration_acefile_for_ace
{
    my ($self, $current_ace) = @_;

    return unless Finfo::Validate->validate
    (
        attr => 'ace to get next iteration',
        value => $current_ace,
        isa => 'string',
        msg => 'fatal',
    );
    
    my ($iter) = ( $current_ace =~ s/\.(\d+$)// )
    ? $1
    : 0;

    my $acefile = $self->acefile_for_ace($current_ace);

    return unless $acefile;
    
    my $return_acefile;
    do
    {
        $return_acefile = sprintf('%s.%d', $acefile, ++$iter);
    } until ( not -e $return_acefile );
    
    return unless Finfo::Validate->validate
    (
        attr => 'next iteration acefile',
        value => $return_acefile,
        isa => 'file_w',
        msg => 'fatal',
    );

    return $return_acefile;
}
 
sub all_acefiles
{
    my $self = shift;
    
    my $acedir = $self->dir;
    
    my %acefiles;
    foreach my $acefile ( <$acedir/*.ace*> )
    {
        $acefiles{$acefile} = (stat($acefile))[9];
    }

    return sort { $acefiles{$b} <=> $acefiles{$a} } keys %acefiles;
}

sub all_aces
{
    my $self = shift;
    
    return map { basename($_) } $self->all_acefiles;
}

sub number_of_aces
{
    my $self = shift;

    return scalar($self->all_acefiles);
}

BEGIN
{
    *valid_acefiles = \&acefiles;
    *valid_aces = \&aces;
}

sub acefiles
{
    my $self = shift;

    my @acefiles;
    foreach my $af ( $self->all_acefiles )
    {
        tie(my @lines, 'Tie::File', $af)
            or ( $self->warn_msg("Can't tie acefile($af): $!") and next );
        my $line = $lines[0];
        untie @lines;
        next unless $line and $line =~ /^AS\s\d+\s\d+/;
        push @acefiles, $af;
    }
    
    return @acefiles;
}

sub aces
{
    my $self = shift;
    
    return map { basename($_) } $self->acefiles;
}

sub recent_ace
{
    my $self = shift;

    my @aces = $self->aces;

    return shift @aces;
}

sub recent_acefile
{
    my $self = shift;

    my @acefiles = $self->acefiles;

    return $acefiles[0];
}

sub recent_ace_for_pattern
{
    my ($self, $pattern) = @_;

    my ($ace) = grep { m/$pattern/ } $self->aces;
    
    return $ace;
}

sub _unique_ace_hash : PRIVATE
{
    my $self = shift;

    my %aces;
    for my $ace ( $self->aces )
    {
        my ($base_name) = split(/\.ace/, $ace);
        my $age = $self->age_for_ace($ace, 5);
        if ( exists $aces{$base_name} )
        {
            my ($stored_ace, $stored_age) = @{ $aces{$base_name} };
            #print Dumper([$base_name, $aces{$base_name}]);
            next if $stored_age <= $age;
           $aces{$base_name} = [ $ace, $age ];
        }
        else
        {
           $aces{$base_name} = [ $ace, $age ];
        }
    }

    $self->error_msg( sprintf('No aces found in dir (%s)', $self->dir) )
        and return unless %aces;
    
    return \%aces;
}

sub most_recent_acefiles_from_unique_ace_names
{
    my $self = shift;

    my @aces = $self->most_recent_aces_from_unique_ace_names
        or return;

    return map { $self->acefile_for_ace($_) } @aces;
}

sub most_recent_aces_from_unique_ace_names
{
    my $self = shift;

    my $aces = $self->_unique_ace_hash
        or return;
    
    return sort { $a cmp $b } map { $_->[0] } values %$aces;
}

sub unique_ace_base_names
{
    my $self = shift;

    my $aces = $self->_unique_ace_hash
        or return;

    return sort { $a cmp $b }keys %$aces;
}

sub date_for_ace
{
    my ($self, $ace) = @_;

    my $acefile = $self->get_and_validate_acefile_for_ace($ace)
        or return;
    
    return scalar( localtime( (stat($acefile))[8] ) );
}

sub age_for_ace
{
    my ($self, $ace, $places) = @_;

    my $acefile = $self->get_and_validate_acefile_for_ace($ace)
        or return;
    
    $places = 1 unless defined $places;

    my $sprintf = sprintf('%%.%sf', ( defined $places ) ? $places : 1);
    
    return sprintf($sprintf, -M $acefile);
}

sub owner_for_ace 
{
    my ($self, $ace) = @_;
    
    my $acefile = $self->get_and_validate_acefile_for_ace($ace)
        or return;

    return (getpwuid( (stat($acefile))[4] ))[0];
}

1;

=pod

=head1 Name

Finishing::Assembly::Ace::Dir

=head1 Synopsis

Provides methods to get info about acefiles in a directory.

=head1 Usage

 use Finishing::Assembly::Ace::Dir;

 my $acedir = Finishing::Assembly::Ace::Dir->new(dir => $dir)
    or die;

 my $recent_acefile = $acedir->recent_acefile;

 etc...

=head1 Definitions

=over

=item B<dir>            The directory supplied in the constructor

=item B<acefile>        The full path of the acefile

=item B<ace>            The basename of the acefile

=item B<ace base name>  The basename of the ace, substituting everything past '.ace'

=back

=head1 Methods
 
=head2 dir

 my $dir = $acedir->dir;

=over

=item I<Synopsis>   uses File::Spec->catfile to concatenate the directory and the ace

=item I<Params>     none - read only

=item I<Returns>    dir (string)

=back

=head2 acefile_for_ace

 my $acefile = $acedir->acefile_for_ace($ace);
 
=over

=item I<Synopsis>   uses File::Spec->catfile to concatenate the directory and the ace

=item I<Params>     ace (string) to get acefile for

=item I<Returns>    acefile (string) for ace

=back
 
=head2 get_and_validate_acefile_for_ace

 my $acefile = $acedir->get_and_validate_acefile_for_ace($ace);

=over

=item I<Synopsis>   Gets the full acefile for ace, checking that it exists and is readable

=item I<Params>     ace (string) to get acefile for

=item I<Returns>    acefile (string) for ace

=back 

=head2 next_iteration_acefile_for_ace

 my $acefile = $acedir->next_iteration_acefile_for_ace($ace);

=over

=item I<Synopsis>   Gets the acefile in the next iteration for an ace, ensures the that new acefile does not exist 

=item I<Params>     ace (string) to get acefile for

=item I<Returns>    next iteration acefile (string) for ace

=back

=head2 all_acefiles

 my @acefiles = $acedir->all_acefiles;

=over

=item I<Params>     none

=item I<Returns>    acefiles (array) - all files in dir matching *ace*

=back

=head2 all_aces

 my @ace = $acedir->all_aces;

=over

=item I<Params>     none

=item I<Returns>    aces (array) - basenames for all files in dir matching *ace*

=back 

=head2 acefiles, valid_acefiles

 my @acefiles = $acedir->acefiles;
 
=over

=item I<Params>     none

=item I<Returns>    acefiles (array) - files in dir matching *ace*, with the first line of the file starting with 'AS'

=back

=head2 aces, valid_aces

 my @aces = $acedir->aces;

=over

=item I<Params>     none

=item I<Returns>    aces (array) - basenames of acesfiles in dir matching *ace*, with the first line of the file starting with 'AS'

=back

=head2 recent_ace

 my $ace = $acedir->recent_ace;
 
=over

=item I<Params>     none

=item I<Returns>    ace (string) - the basename of the most recent acefile, with the first line of the file starting with 'AS'
 
=back

=head2 recent_acefile

 my $acefile = $acedir->recent_ace;
 
=over

=item I<Params>     none

=item I<Returns>    acefile (string) - the most recent file (full path), matching *ace*, with the first line of the file starting with 'AS'
    
=back

=head2 most_recent_aces_from_unique_ace_names

 my @aces = $acedir->most_recent_aces_from_unique_ace_names
    or die;
 
=over

=item I<Synopsis>   Grabs the 'base name' of each ace by spliting on '.ace', storing which ace is the most recent

=item I<Params>     none

=item I<Returns>    aces (ary or strings)

=back

=head2 unique_ace_base_names

 my @ace_base_names = $acedir->unique_ace_base_names
    or die;
 
=over

=item I<Synopsis>   Grabs the 'base name' of each ace by spliting on '.ace'

=item I<Params>     none

=item I<Returns>    ace base names (ary or strings)

=back

=head2 date_for_ace

 my $date = $acedir->date_for_ace($ace);
 
=over

=item I<Params>     ace (string)

=item I<Returns>    The date stamp (string) of acefile

=back

=head2 age_for_ace

 my $age = $acedir->age_for_ace($ace, $places);
 
=over

=item I<Params>     ace (string), places (int, default is 1) - the number of places in the return age

=item I<Returns>    acefile's age in days (real)
 
=back

=head2 owner_for_ace

 my $owner = $acedir->owner_for_ace($ace);

=over

=item I<Params>     ace (string)

=item I<Returns>    acefile owner's username (string)

=back

=head1 Disclaimer

Copyright (C) 2006-2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> <ebelter@watson.wustl.edu>

=cut

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Ace/Dir.pm $
#$Id: Dir.pm 30961 2007-12-13 20:47:12Z ebelter $
