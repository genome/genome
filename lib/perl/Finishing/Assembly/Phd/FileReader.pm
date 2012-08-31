package Finishing::Assembly::Phd::FileReader;

use strict;
use warnings;

use base 'Finfo::Singleton';

use Data::Dumper;
use Finishing::Assembly::Factory;
use Finishing::Assembly::Source::Tags;

sub execute
{
    my ($self, $io) = @_;

    my %object_builders = 
    ( 
        BEGIN_COMMENT => '_parse_comment',
        BEGIN_DNA => '_parse_DNA',
        "WR{" => '_parse_WR',
        BEGIN_TAG => '_parse_tag',
    );

    my %phd;
    my $header_line = $io->getline;
    $self->fatal_msg("No phd header line") unless $header_line;
    chomp $header_line;
    my @tokens = split(/ /, $header_line);
    $self->fatal_msg("Invalid phd header: $header_line") unless $tokens[0] eq 'BEGIN_SEQUENCE';
    $phd{name} = $tokens[1];

    while ( my $line = $io->getline ) 
    {
        last if $line =~ /^BEGIN_SEQUENCE\s/; #phd.ball compatible
        chomp $line;
        my @tokens = split(/ /, $line);
        my $tag = shift @tokens;
        if ( $tag and my $method = $object_builders{$tag}) 
        {
            $self->$method($io, \%phd);
        }
    }

    return Finishing::Assembly::Factory->connect('source')->create_assembled_read(%phd);
}

# COMMENT EX:
# CHROMAT_FILE: L24337P6000G2.g1
# ABI_THUMBPRINT: 0
# PHRED_VERSION: 0.040406.a_2
# CALL_METHOD: abi
# QUALITY_LEVELS: 99
# TIME: Tue Nov 6 15:08:10 2007
# TRACE_ARRAY_MIN_INDEX: 0
# TRACE_ARRAY_MAX_INDEX: 8736
# CHEM: term
# DYE: big
sub _parse_comment : PRIVATE
{
    my ($self, $io, $phd) = @_;

    while ( my $line = $io->getline ) 
    {
        chomp $line;
        last if $line =~ /END_COMMENT/;
        next if $line =~ /^\s*$/; 
        my @tokens = split(/:\s/, $line);
        $phd->{comments}->{ lc($tokens[0]) } = $tokens[1];
    }
    
    return 1;
}

sub _parse_DNA : PRIVATE
{
    my ($self, $io, $phd) = @_;

    while ( my $line = $io->getline ) 
    {
        chomp $line;
        last if $line =~ /END_DNA/;
        my @pos = split(/ /,$line);
        $phd->{base_string} .= $pos[0];
        push @{ $phd->{qualities} }, $pos[1];
        push @{ $phd->{chromat_positions} }, $pos[2];
    }

    return 1;
}

sub _parse_WR 
{
    my ($self, $io, $phd) = @_;

    my $wr;
    while (my $line = $io->getline ) 
    {
        last if $line =~ /^}/;
        $wr .= $line;
    }
    
    chomp $wr;
    
    return push @{ $phd->{wr} }, $wr;
}

sub _parse_tag 
{
    my ($self, $io, $phd) = @_;
    
    my %tag_data;
    my $text;
    while ( my $line = $io->getline ) 
    {
        chomp $line;
        last if $line =~ /^END_TAG/;
        if ($line =~ /BEGIN_COMMENT/) 
        {
            while ( my $comment_line = $io->getline ) 
            {
                last if $comment_line =~ /END_COMMENT/;
                $text .= $comment_line;
            }
        }
        my @tokens = split(/\:\s/, $line);
        $tag_data{$tokens[0]} = $tokens[1];
    }
    my @pos = split(/\s/, $tag_data{UNPADDED_READ_POS});
    
    my $tag = Finishing::Assembly::Source::ReadTag->new
    (
        parent => $phd->{name},
        type   => $tag_data{TYPE},
        source => $tag_data{SOURCE},
        date   => $self->_convert_from_phd_tag_date( $tag_data{DATE} ),
        unpad_start  => $pos[0],
        unpad_stop   => $pos[1],
    );

    $tag->text($text) if $text;
    
    return push @{ $phd->{tags} }, $tag;
}

sub _convert_from_phd_tag_date : PRIVATE
{
    my ($self, $date) = @_;

    return ( $date =~ /^[89]/ ? ('19' . $date) : ('20' . $date));
}

1;

=pod

=head1 Name

Finishing::Assembly::Phd::FileReader 

=head1 Synopsis

This package is a phd file reader.  Given an IO object, an assemled read object representing the phd infomation in the IO will be returned.  It is a singleton, so "create" by calling instance on the class, then call the execute method.

=head1 Usage

 use Finishing::Assembly::Phd::FileReader;
 use IO::File;

 my $fh = IO::File->new("< read.phd.1")
    or die "$!\n";
 my $reader = Finishing::Assembly::Phd::FileReader->instance;
 my $phd = $reader->execute($fh);
 $fh->close;

 print $phd->name,"\n";

=head1 Methods

 my $phd = Finishing::Assembly::Phd::FileReader->instance->execute($file_handle);
 
=over

=item I<Synopsis>   Parses the phd info in the file handle

=item I<Params>     file handle (IO::File or related object)

=item I<Returns>    Finishing::Assembly::Source::AssembledRead (object)

=back

=cut

#$HeadURL$
#$Id$
