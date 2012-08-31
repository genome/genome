package Finishing::Assembly::Ace::Output;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

my %objects :name(objects:r) 
    :ds(aryref);
my %format  :name(format:r) 
    :isa([ 'in_list',  __PACKAGE__->valid_formats ])
    :clo('o|output=s')
    :desc('Output type(s): ' . join(', ', __PACKAGE__->valid_formats),);
my %add_name_to_output  :name(add_name_to_output:o) 
    :default(0)
    :clo('inc-name')
    :desc('Include the acefile base name to the output');

sub valid_formats
{
    return 
    (qw/
        fasta qual raw contignames readnames tagstext status bp_start_end bp_df bp_tags
        /);
}

sub execute
{
    my $self = shift;

    my $method = $self->format;

    return $self->_add_name_to_output( $self->$method );

}

sub _add_name_to_output : PRIVATE
{
    my ($self, $output) = @_;

    $self->error_msg("No output")
        and return unless $output;
    
    return $output unless $self->add_name_to_output;
    
    my $name = $self->add_name_to_output;

    if ( grep { $_ eq $self->format } (qw/ fasta qual /) )
    {
        $output =~ s/>/>$name\./g;
    }
    else
    {
        $output = ">$name\n$output";
    }

    return $output; 
}

sub fasta : PRIVATE
{
    my $self = shift;

    return $self->_bioseqs_to_string('Fasta');
}

sub qual : PRIVATE
{
    my $self = shift;

    return $self->_bioseqs_to_string('qual');
}

sub raw : PRIVATE
{
    my $self = shift;

    return $self->_bioseqs_to_string('raw');
}

sub _bioseqs_to_string : PRIVATE
{
    my ($self, $format) = @_;

    my $io = IO::String->new();

    my $writer = Bio::SeqIO->new('-fh' => $io, '-format' => $format);

    foreach my $bioseq (@{ $self->objects })
    {
        $writer->write_seq($bioseq);
    }

    $io->seek(0, 0);

    return join('', $io->getlines);
}

sub status : PRIVATE
{
    my $self = shift;

    my $io = IO::String->new();

    my $head = sprintf
    (
        "%-20s %-20s\n%-35s\n",
        "Contig",
        "Unpadded Size",
        "------------------------------------------"
    );
    $io->print($head);

    my $total;
    foreach my $ctg (@{ $self->objects })
    { 
        $total += $ctg->{unpad_end};

        $io->print
        (
            sprintf
            (
                "%-20s %-20d\n", $ctg->{name},  $ctg->{unpad_end}
            )
        );
    }

    my $tail = sprintf
    (
        "%-30s\n%-20s %-20d\n\n",
        "-------------------------------------",
        "Total",
        $total
    );
    $io->print($tail);
    
    $io->seek(0, 0);

    return join("", $io->getlines);
}

sub contignames : PRIVATE
{
    my $self = shift;

    return join("\n", @{ $self->objects });
}

sub readnames : PRIVATE
{
    my $self = shift;

    return join("\n", map { $_->name } @{ $self->objects });
}

sub tagstext : PRIVATE
{
    my $self = shift;

    my $text = "*^*^*^*^*^**^*^*^*^*^*^*^*^*^*^*^\n";
    
    foreach my $tag (@{ $self->objects })
    {
        my $comment = $tag->text || 'none';
        $comment =~ s/comment{\n|c}\n//ig; 

        $text .= sprintf
        (
            "Contig: %s\nType: %s\nStart: %d\nStop: %d\nComment:\n%s\n*^*^*^*^*^**^*^*^*^*^*^*^*^*^*^*^\n",
            $tag->parent,
            $tag->type,
            $tag->unpad_start,
            $tag->unpad_stop,
            $comment,
        );
    }

    return $text;
}

sub bp_tags : PRIVATE
{
    my $self = shift;

    my $bp = 0;
    foreach my $tag ( @{ $self->objects } )
    {
        $bp += $tag->unpad_stop - $tag->unpad_start + 1;
    }

    return $bp;
}

sub bp_df : PRIVATE
{
    my $self = shift;

    return $self->bp_tags;
}

sub bp_start_end : PRIVATE
{
    my $self = shift;

    my $bp = 1;
    foreach my $tag (@{ $self->objects })
    {
        next unless defined $tag->text;
        if ( $tag->text =~ /start\s+finished\s+region/i )
        {
            $bp -= $tag->unpad_start;
        }
        elsif ( $tag->text =~ /end\s+finished\s+region/i )
        {
            $bp += $tag->unpad_stop;
        }
        elsif ( $tag->text =~ /COMMENT\{\nend\s+of\s+/i and $tag->unpad_start > 1 )
        {
            $bp = $tag->unpad_stop;
            last;
        }
    }

    return $bp;
}

1;

=pod

=head1 Name

Finishing::Assembly::Ace::Output

=head1 Synopsis

This module takes an aryref of objects and format.  It will then create the output for the objects based on the format type.

=head1 Usage

I<Get fasta sequence from an acefile>

 my $ace_dir = Finishing::Assembly::Ace::Dir->new(dir => '~/seqmgr/M_BB0392D19/edit_dir')
    or die;
    
 my $acefile = $acedir->recent_acefile;
    or die "\n";
    
 my $ace = Finishing::Assembly::Ace->new(input_file => $acefile) # see usage for more details
    or die;

 my $ace_ext = Finishing::Assembly::Ace::Ext->new(ace => $ace) # see usage for more details
    or die;

 my $bioseqs = $ace_ext->contigs_to_bioseqs # see usage for details
    or die; 

 my $aceout = Finishing::Assembly::Ace::Output->new(objects => $bioseqs, type => 'fasta')
    or die;

 my $fasta_seq = $aceout->execute
    or die;

 my $fasta_file = 'ace.fasta';
 my $fasta_fh = IO::File->new("> $fasta_file");
 $fasta_fh->print($output);
 $fasta_fh->close;

Voila!  You just got fasta sequence from an acefile!

=head1 Methods

=head2 valid_formats

 my @formats = $aceout->valid_formats

=over

=item I<Synopsis>    Returns a list of the valid formats for the output.

=item I<Params>      none

=item I<Returns>     output formats (array of strings)

=back

=head2 execute

=over

 my $output = $ace_output->execute
    or die;

=item I<Synopsis>    Creates the output generated from the objects for the format type.

=item I<Params>      none

=item I<Returns>     output (string)

=back

=head1 Disclaimer

Copyright (C) 2006-7 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> <ebelter@watson.wustl.edu>

=cut

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Ace/Output.pm $
#$Id: Output.pm 30495 2007-11-29 17:52:30Z ebelter $
