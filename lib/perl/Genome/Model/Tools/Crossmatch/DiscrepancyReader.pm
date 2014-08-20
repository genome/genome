package Genome::Model::Tools::Crossmatch::DiscrepancyReader;

use strict;
use warnings;

class Genome::Model::Tools::Crossmatch::DiscrepancyReader {
    is => 'Genome::Utility::IO::Reader',
    has => [
        _type => { is_optional => 1, default_value => 'list', },
    ]
};

sub create {
    my $class = shift;
    
    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $line;
    while ( $line = $self->getline )
    {
        last if $line =~ /\w\d/;
    }

    my $type;
    if ($line =~ /Discrep(ancy)?\s+(\w+)/) {
        $type = $2;
    }

    unless ( $type )
    {
        $self->input->seek(0, 0); # if no type, assume it's a list
    }
    elsif ( $type eq 'histogram' )
    {
        $self->_type('table');
        $self->getline; # this line is a header
    }
    elsif ( $type eq 'list' )
    {
        # do nothing, current line is a header
    }
    else # unsupported type
    {
        die "Unsupported discrepancy type ($type)";
    }

    return $self;
}

# Discrepancy histogram
# Qual algn  cum    rcum    (%)    unalgn X    N  sub del ins  total (%)   cum  rcum (%)
# 15    120    120    120 (100.00)    11  0    0  10   0   0    10 (8.33)   10   10 (8.33)
#
# Discrep list
# DISCREPANCY   D-3   160  T(-1)    240  atccctTtgctcc
# DISCREPANCY   S     335  C(-1)     63  tttgaaCggcact
# DISCREPANCY   I-2    28  TT(0)    425  ggcacgTTgttggc
# DISCREPANCY   D     588  A(-1)    987  aaatggAtataga
# 
# or
# 
# S   1459565  C(15)  206217  tccaggCtatttt
# D   1470261  T(15)  216914  cggcctTcatatt
# I   1474631  A(15)  221283  ttgtcaAtgaaac
sub next {
    my $self = shift;

    my $line = $self->getline;

    return if not defined $line or $line eq '';

    chomp $line;

    my $method = '_parse_' . $self->_type;

    return $self->$method($line);
}

sub _parse_histogram {
    my ($self, $line) = @_;
    
    my @tokens = split(/\s+/, $line);
    my %discrep = ( type => 'histogram' );
    @discrep {qw/
        quality align align_cum align_rcum align_per unalign X
        N sub del ins discrep_total discrep_per discrep_cum discrep_rcum
        /} = @tokens;

    return \%discrep;
}

sub _parse_list
{
    # DISCREPANCY   D-3   160  T(-1)    240  atccctTtgctcc
    # /^DISCREPANCY\s+([DSI])-?\d\s+(\d+)\s+(\w+)\(\-?\d+\)\s+(\d+)\s+(\w+)$/
    #
    # S   1459565  C(15)  206217  tccaggCtatttt
    # D-2   105  C(15)   7443  ctgctgCtgtatg
    my ($self, $line) = @_;

    $line =~ s/DISCREPANCY//;
    $line =~ s/^\s+//;

    my @tokens = split(/\s+/, $line);
    $self->fatal_msg("Invalid discrepancy line:\n$line")
        and return unless @tokens == 5;

    my ($mut, $num) = split /\-/, $tokens[0];
    $num = 1 unless defined $num;

    my ($base, $bases_in_seq) = split(/\(/, $tokens[2]);

    return 
    {
        mutation => $mut,
        number => $num, 
        query_pos => $tokens[1],
        base => $base,
        subject_pos => $tokens[3],
        sequence => $tokens[4]
    };
}

1;

=pod

=head1 Name

 Genome::Model::Tools::Crossmatch::DiscrepancyStringParser

=head1 Description

 Parses a discrepancy string from a crossmatch output.  Mainly a helper object for Alignment::Crossmatch::Reader.

=head1 Usage

 use Genome::Model::Tools::Crossmatch::DiscrepancyStringParser;

 my $reader = Genome::Model::Tools::Crossmatch::DiscrepancyStringParser->create(
    input => $discrep_io, # required
 )
    or die;

 my @discrepancies = $reader->all
    or die;

=head1 Methods

=head2 next

 my $discrep = $reader->next;

=over

=item I<Synopsis>   Parses, creates and returns a discrepancy from the io

=item I<Params>     none

=item I<Returns>    discrepancy (hashref, scalar)

=back

=head2 all

 my $discrep = $reader->all;

=over

=item I<Synopsis>   Parses, create and returns all of the discrepancies from the io

=item I<Params>     none

=item I<Returns>    all discrepancies (hashrefs, array)

=back

=head1 See Also

=over

=item cross_match

=back

=head1 Disclaimer

Copyright (C) 2006-7 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> <ebelter@watson.wustl.edu>

=cut

