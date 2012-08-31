package Genome::Model::Tools::Consed::Directory;

use strict;
use warnings;

use Genome;

use Carp 'confess';
require Cwd;
use Data::Dumper 'Dumper';
require File::Basename;
require File::Spec;

class Genome::Model::Tools::Consed::Directory {
    is => 'UR::Object',
    has => [
        directory => {
            type => 'String',
            doc => 'Consed style directory',
        },
    ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;
    
    eval {
        Genome::Sys->validate_existing_directory($self->directory)
    };
    if($@) {
        $self->warning_message('Directory does not exist: "' . $self->directory . '".');
        return;
    }

    return $self;
}

# Dirs
sub directories {
    return (qw/ edit_dir phd_dir chromat_dir /);
}

sub extended_directories {
    return ( directories(), (qw/ fasta abi_dir reports /) );
}

sub edit_dir {
    return $_[0]->concat_directory('edit_dir');
}

sub phd_dir {
    return $_[0]->concat_directory('phd_dir');
}

sub chromat_dir {
    return $_[0]->concat_directory('chromat_dir');
}

sub fasta_dir {
    return $_[0]->concat_directory('fasta');
}

sub abi_dir {
    return $_[0]->concat_directory('abi_dir');
}

sub reports_dir {
    return $_[0]->concat_directory('reports');
}

sub concat_directory {
    my ($self, $sub_dir) = @_;
    
    $self->error_message("Need sub directory to concatenate")
        and return unless $sub_dir;
    
    return File::Spec->catfile($self->directory, $sub_dir);
}

sub create_consed_directory_structure {
    my $self = shift;
    
    for my $sub_dir ( $self->directories ) {
        $self->create_directory($sub_dir);
    }

    return 1;
}

sub create_extended_directory_structure {
    my $self = shift;
    
    for my $sub_dir ( $self->extended_directories ) {
        $self->create_directory($sub_dir);
    }

    return 1;
}
 
sub create_directory {
    my ($self, $sub_dir) = @_;

    my $dir = $self->concat_directory($sub_dir);
    
    Genome::Sys->create_directory($dir)
        or confess "Can't make directory ($dir)\n";
    
    return 1;
}

#- ACE -#
sub acefile_for_ace
{
    my ($self, $ace) = @_;

    $self->error_message("Need ace to get acefile")
        and return;
    
    return File::Spec->catfile($self->edit_dir, $ace);
}
    
sub get_and_validate_acefile_for_ace
{
    my ($self, $ace) = @_;

    my $acefile = $self->acefile_for_ace($ace)
        or return;

    $self->error_message("Acefile ($acefile) does not exist")
        and return unless -e $acefile;
    $self->error_message("Acefile ($acefile) exists but is empty")
        and return unless -s $acefile;

    return $acefile;
}

sub next_iteration_acefile_for_ace
{
    my ($self, $current_ace) = @_;

    if ( not $current_ace ) {
        $self->error_message('No current ace to get next iteration ace file');
        return;
    }
    
    my ($iter) = ( $current_ace =~ s/\.(\d+$)// )
    ? $1
    : 0;

    my $acefile = $self->acefile_for_ace($current_ace);

    return unless $acefile;
    
    my $return_acefile;
    do {
        $return_acefile = sprintf('%s.%d', $acefile, ++$iter);
    } until ( not -e $return_acefile );
    
    return $return_acefile;
}
 
sub all_acefiles {
    my $self = shift;
    
    my %acefiles;
    for my $acefile ( glob( sprintf('%s/*.ace*', $self->edit_dir) ) ) {
        $acefiles{$acefile} = (stat($acefile))[9];
    }

    return sort { $acefiles{$b} <=> $acefiles{$a} } keys %acefiles;
}

sub all_aces {
    my $self = shift;
    
    return map { File::Basename::basename($_) } $self->all_acefiles;
}

sub number_of_aces {
    my $self = shift;

    return scalar( $self->all_acefiles );
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
    for my $af ( $self->all_acefiles ) {
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

sub _unique_ace_hash {
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

    $self->error_message( sprintf('No aces found in dir (%s)', $self->directory) )
        and return unless %aces;
    
    return \%aces;
}

sub most_recent_acefiles_from_unique_ace_names {
    my $self = shift;

    my @aces = $self->most_recent_aces_from_unique_ace_names
        or return;

    return map { $self->acefile_for_ace($_) } @aces;
}

sub most_recent_aces_from_unique_ace_names {
    my $self = shift;

    my $aces = $self->_unique_ace_hash
        or return;
    
    return sort { $a cmp $b } map { $_->[0] } values %$aces;
}

sub unique_ace_base_names {
    my $self = shift;

    my $aces = $self->_unique_ace_hash
        or return;

    return sort { $a cmp $b }keys %$aces;
}

sub date_for_ace {
    my ($self, $ace) = @_;

    my $acefile = $self->get_and_validate_acefile_for_ace($ace)
        or return;
    
    return scalar( localtime( (stat($acefile))[8] ) );
}

sub age_for_ace {
    my ($self, $ace, $places) = @_;

    my $acefile = $self->get_and_validate_acefile_for_ace($ace)
        or return;
    
    $places = 1 unless defined $places;

    my $sprintf = sprintf('%%.%sf', ( defined $places ) ? $places : 1);
    
    return sprintf($sprintf, -M $acefile);
}

sub owner_for_ace {
    my ($self, $ace) = @_;
    
    my $acefile = $self->get_and_validate_acefile_for_ace($ace)
        or return;

    return (getpwuid( (stat($acefile))[4] ))[0];
}

sub touch_singlets_file_for_acefile
{
    my ($self, $acefile) = @_;

    $self->fatal_msg("Need acefile to touch singlets") unless $acefile;
    
    my $singlets_file = $acefile . '.singlets';

    return 1 if -e $singlets_file;

    system("touch $singlets_file");

    return 1 if -e $singlets_file;
    
    $self->info_msg("Failed to create singlets file for acefile ($acefile)");

    return;
}

1;

=pod

=head1 Name

Genome::Model::Tools::Consed::Directory;
 
=head1 Synopsis

=head1 Usage

=head1 Methods

=head2 edit_dir

 Get the edit_dir.  This is where the read prefixes file will be created.

=head2 chromat_dir

 Get the chromat_dir.  This is where the traces will be moved to.
 
=head2 phd_dir

 Get the phd_dir.  This is where the phds will be created.

=head2 read_prefixes_file

 Get the read prefixes file name for consed.

=head1 Definitions

=over

=item B<directoty>      The directory supplied in the constructor

=item B<acefile>        The full path of the acefile

=item B<ace>            The basename of the acefile

=item B<ace base name>  The basename of the ace, substituting everything past '.ace'

=back

=head1 Methods
 
=head2 directory

 my $dir = $acedir->directory;

=over

=item I<Synopsis>   Get the directory

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

#$HeadURL$
#$Id$
