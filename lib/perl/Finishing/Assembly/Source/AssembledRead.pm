package Finishing::Assembly::Source::AssembledRead;

use strict;
use warnings;

use base 'Finfo::Accessor';

use Data::Dumper;

__PACKAGE__->mk_accessors
(qw/
    name base_string qualities chromat_positions comments
    /);

sub phd_version
{
    my $self = shift;

    my $phd_file = $self->{comment}->{phd_file}
        or return;

    return ($phd_file =~ /\.phd\.(\d+)$/)[0];
}

#- COMMENT ATTRS -#
sub phd_file
{
    my ($self, $phd_file) = @_;

    $self->{comments}->{phd_file} = $phd_file if $phd_file;
    
    return $self->{comments}->{phd_file};
}

sub chromat_file
{
    my ($self, $chromat_file) = @_;

    $self->{comments}->{chromat_file} = $chromat_file if $chromat_file;
    
    return $self->{comments}->{chromat_file};
}

sub abi_thumbprint
{
    my ($self, $abi_thumbprint) = @_;

    $self->{comments}->{abi_thumbprint} = $abi_thumbprint if $abi_thumbprint;
    
    return $self->{comments}->{abi_thumbprint};
}

sub phred_version
{
    my ($self, $phred_version) = @_;

    $self->{comments}->{phred_version} = $phred_version if $phred_version;
    
    return $self->{comments}->{phred_version};
}

sub call_method
{
    my ($self, $call_method) = @_;

    $self->{comments}->{call_method} = $call_method if $call_method;
    
    return $self->{comments}->{call_method};
}

sub quality_levels
{
    my ($self, $quality_levels) = @_;

    $self->{comments}->{quality_levels} = $quality_levels if $quality_levels;
    
    return $self->{comments}->{quality_levels};
}

sub time
{
    my ($self, $time) = @_;

    $self->{comments}->{'time'} = $time if $time;
    
    return $self->{comments}->{'time'};
}

sub trace_array_min_index
{
    my ($self, $trace_array_min_index) = @_;

    $self->{comments}->{trace_array_min_index} = $trace_array_min_index if $trace_array_min_index;
    
    return $self->{comments}->{trace_array_min_index};
}

sub trace_array_max_index
{
    my ($self, $trace_array_max_index) = @_;

    $self->{comments}->{trace_array_max_index} = $trace_array_max_index if $trace_array_max_index;
    
    return $self->{comments}->{trace_array_max_index};
}

sub chem
{    
    my ($self, $chem) = @_;

    $self->{comments}->{chem} = $chem if $chem;
    
    return $self->{comments}->{chem};
}

sub dye
{
    my ($self, $dye) = @_;

    $self->{comments}->{dye} = $dye if $dye;
    
    return $self->{comments}->{dye};
}

#- ARYREF ATTRS -#
#- TAGS -#
sub tags
{
    my ($self, $new_tags) = @_;

    $self->{tags} = $new_tags if $new_tags;

    return $self->{tags} || [];
}

sub add_tags
{
    my ($self, $tags) = @_;

    return push @{ $self->tags }, @$tags;
}

sub remove_tags
{
    my $self = shift;

    return $self->tags([]);
}

#- WR -#
sub wr
{
    my ($self, $new_wr) = @_;

    if ( $new_wr )
    {
        $self->{wr} = $new_wr;
        $self->{template} = undef;
        $self->{lib} = undef;
        $self->{primer} = undef;
    }

    return $self->{wr} || [];
}

sub template
{
    # TODO allow setting?
    my $self = shift;

    foreach my $wr ( @{ $self->wr } )
    {
        if ( $wr =~ /template/s ) 
        {
            $wr =~ /name: (\S*)/s;
            return $self->{template} = $1;
        }
    }

    return;
}

sub lib
{
    # TODO allow setting?
    my $self = shift;

    foreach my $wr ( @{ $self->wr } )
    {
        if ( $wr =~ /template/s ) 
        {
            $wr =~ /lib: (\S*)/s;
            return $self->{lib} = $1;
        }
    }

    return;
}

sub primer
{
    # TODO allow setting?
    my $self = shift;

    foreach my $wr ( @{ $self->wr } )
    {
        if ( $wr =~ /primer/s ) 
        {
            $wr =~ /type: (.+)/s;
            $self->{primer} = $1;
            chomp $self->{primer};
            return $self->{primer};
        }
    }

    return;
}

1;

=pod

=head1 Name

Finishing::Assembly::Source::AssembledRead

=head1 Synopsis

=head1 Usage

=head1 Methods

=head1 See Also

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/Finishing/Assembly/Phd/AssembledRead.pm $
#$Id: AssembledRead.pm 32467 2008-02-20 21:55:49Z ebelter $

