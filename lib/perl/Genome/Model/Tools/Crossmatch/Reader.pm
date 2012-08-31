package Genome::Model::Tools::Crossmatch::Reader;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Crossmatch::Reader {
    is => 'Genome::Utility::IO::Reader',
    has => [
        _header => { is_optional => 1, },
    ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $header;
    my $prev_pos = $self->input->tell;
    my $line = $self->getline;
    until ($self->_line_is_an_alignment_line($line)) {
        $prev_pos = $self->input->tell;
        $header .= $line;
        $line = $self->getline;
    }

    $self->input->seek($prev_pos, 0);

    $self->_header($header);

    return $self;
}

sub next {
    my $self = shift;

    my $alignment;
    my $discrep_io = IO::String->new();
    my $num_of_discreps = 0;
    while ( 1 )
    {
        my $prev_pos = $self->input->tell;
        my $line = $self->getline;
        last unless $line;
        $line =~ s/^\s+//;
        last if $line eq '';
        if ( $self->_line_is_an_alignment_line($line) )
        {
            chomp $line;
            if ( $alignment )
            {
                $self->input->seek($prev_pos, 0);
                last;
            }
            $line =~ s/^ALIGNMENT\s*//;
            $alignment = $self->_create_alignment_ref_from_line($line);
            next;
        }
        $discrep_io->print($line);
        $num_of_discreps++;
    }

    return unless $alignment;
    
    if ( $num_of_discreps )
    {
        $discrep_io->seek(0, 0);
        my $discrep_reader = Genome::Model::Tools::Crossmatch::DiscrepancyReader->new(
            input => $discrep_io,
        )
            or return;
        my @discrepancies = $discrep_reader->all
            or return;
        $alignment->{discrepancies} = \@discrepancies;
    }

    #print Dumper($alignment);

    return $alignment;
}

# Here's some alignments.  They may not have the 'ALIGMENT' at the
# beginning, a star(*) at the end or a 'C' in the middle.  Fun, huh?
# ALIGNMENT    89  8.33 0.00 0.00  Contig169        1   120 (49)  C AluSx_SINE/Alu   (171)   131    12 *
# ALIGNMENT  1159  0.08 0.08 0.00  Contig958        1  1193 (0)    Contig11       15  1208 (8)

sub _line_is_an_alignment_line {
    my ($self, $line) = @_;

    $line =~ s/^\s+//;

    return $line =~/^(ALIGNMENT)?\s*\d+\s+\d+\.\d+\s+\d+\.\d+\s+/;
}

sub _create_alignment_ref_from_line {
    my ($self, $line) = @_;

    my @tokens = split(/\s+/, $line);
    pop @tokens if $tokens[-1] eq '*';

    my %alignment =
    (
        sw_score => $tokens[0],
        per_sub => $tokens[1],
        per_del => $tokens[2],
        per_ins => $tokens[3],
        query_name => $tokens[4],
        query_start => $tokens[5],
        query_stop => $tokens[6],
    );
    
    unless ( defined $tokens[7] ) 
    {
        print Dumper($line);
    }
    
    $tokens[7] =~ s/[\(\)]//g;
    $alignment{bases_after}=$tokens[7];

    $alignment{subject_name} = $tokens[-4];
    if ($tokens[8] =~ /^C$/)
    {
        $tokens[-3] =~ s/[\(\)]//g;
        $alignment{bases_before} = $tokens[-3];
        
        $alignment{subject_start} = $tokens[-2];
        $alignment{subject_stop} = $tokens[-1];
    }
    else
    {
        $alignment{subject_start} = $tokens[-3];
        $alignment{subject_stop} = $tokens[-2];

        $tokens[-1] =~ s/[\(\)]//g;
        $alignment{bases_before} = $tokens[-1];
    }

    return \%alignment;
}

1;

=pod

=head1 Name

 Genome::Model::Tools::Crossmatch::Reader;

=head1 Description

 Crossmatch output file reader.

=head1 Usage

 use Genome::Model::Tools::Crossmatch::Reader;

 my $reader = Genome::Model::Tools::Crossmatch::Reader->new( 
    input => 'cm.out', # required - file or IO::* object
 )
    or die;
 
 my @alignments = $reader->all;

=head1 Methods

=head2 next

 my $alignment = $reader->next;

=over

=item I<Synopsis>   Parses, creates and returns a alignment from the io

=item I<Params>     none

=item I<Returns>    alignment (hashref, scalar)

=back

=head2 all

 my $alignment = $reader->all;

=over

=item I<Synopsis>   Parses, create and returns all of the alignments from the io

=item I<Params>     none

=item I<Returns>    all alignment (hashrefs, array)

=back

=head1 See Also

=over

=item Crossmatch dir

=item Genome::Utility::IO::Reader

=back

=head1 Disclaimer

Copyright (C) 2006-11 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> <ebelter@watson.wustl.edu>

=cut

