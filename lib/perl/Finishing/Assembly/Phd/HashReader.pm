package Finishing::Assembly::Phd::HashReader;

use strict;
use warnings;

use base 'Finfo::Singleton';

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
        chomp $line;
        my @tokens = split(/ /, $line);
        my $tag = shift @tokens;
        if ( $tag and my $method = $object_builders{$tag}) 
        {
            $self->$method($io, \%phd);
        }
    }

    return \%phd;
}

sub _parse_comment : PRIVATE
{
    my ($self, $io, $phd) = @_;

    my %comments;
    while ( my $line = $io->getline ) 
    {
        chomp $line;
        last if $line =~ /END_COMMENT/;
        next if $line =~ /^\s*$/; 
        my @tokens = split(/: /,$line);
        $phd->{ lc($tokens[0]) } = $tokens[1];
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
    
    push @{ $phd->{wr} }, $wr;

    if ( $wr =~ /template/s ) 
    {
        $wr =~ /name: (\S*)/s;
        $phd->{template} = $1;
    }
    elsif ( $wr =~ /primer/s ) 
    {
        $wr =~ /type: (\S*)/s;
        $phd->{primer} = $1;
    }

    return 1;
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
        my @element = split /: /, $line;
        $tag_data{$element[0]} = $element[1] if defined $element[0];
    }
    my @pos = split / /, $tag_data{UNPADDED_READ_POS};
    
    my $tag = new Finishing::Assembly::Tag
    {
        type   => $tag_data{TYPE},
        source => $tag_data{SOURCE},
        date   => $tag_data{DATE},
        start  => $pos[0],
        stop   => $pos[1],
        parent => $phd->name,
        scope  => 'Phd',
    };
    $tag->text($text) if $text;
    
    push @{ $phd->{tags} }, $tag;

    return 1;
}

1;

=pod

=head1 NAME

PhdReader - Phd file iterator

=head1 SYNOPSIS

    my $reader = new Finishing::Assembly::Phd::Reader(\*STDIN);
    while (my $obj = $reader->nextObject()) {
    }

=head1 DESCRIPTION

Finishing::Assembly::Phd::Reader iterates over a phd file, returning one element at a time.

=head1 METHODS

=item Finishing::Assembly::Phd::Reader::read 

    $phd = $reader->read(\*STDIN);

=cut

#$HeadURL$
#$Id$
