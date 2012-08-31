package Finishing::Assembly::Source::Contig;

use strict;
use warnings;

use base 'Finfo::Accessor';

use Data::Dumper;
use Finishing::Assembly::Source::Tags;
use Finfo::Iterator;

__PACKAGE__->mk_accessors
(qw/
    name base_string qualities complemented 
    align
    /);

#- READS -#
sub assembled_reads
{
	my ($self, $reads) = @_;
	
	if ( $reads )
	{
		$self->{assembled_reads} = $reads;
	}	

	return Finfo::Iterator->new(objects => $self->{assembled_reads});
}

sub get_assembled_read
{
    my ($self, $name) = @_;
	
	my ($read) = $self->assembled_reads->search({'name' => $name});
	
	$self->fatal_msg("Can't find read ($name)") unless $read;

	return $read;
}

#- BASE SEGMENTS -#

sub base_segments{
    my ($self, $base_segments) = @_;
    $self->{base_segments} = {} unless $self->{base_segments};
    $self->{base_segments} = $base_segments if $base_segments;
    return $self->{base_segments};
}


#- TAGS -#

sub tags{
    my ($self, $tags) = @_;
    $self->{tags} = [] unless $self->{tags};
    $self->{tags} = $tags if $tags;
    return $self->{tags}; 
}

sub create_tag
{
    my $self = shift;

    return Finishing::Assembly::Source::Tag->new
    (
        parent => $self->name,
        @_,
    );
}

sub add_tag
{
    my ($self, $tag) = @_;

    my $new_tag = $self->create_tag
    (
        source => $tag->source,
        start => $tag->start,
        stop => $tag->stop,
        date => $tag->date,
        comment => $tag->comment,
        no_trans => $tag->no_trans,
    );

    push @{ $self->tags }, $new_tag;

    return $new_tag;
}

sub add_tags
{
    my ($self, $tags) = @_;

    foreach my $tag ( @$tags )
    {
        $self->add_tag($tag);
    }
    
    return 1;
}

sub contig_num{
    my ($self, $num) = @_;
    if ($num){
        $self->{contig_num} = $num;
        $self->update_name;
    }else{
        unless ($self->{contig_num}){
            my ($match) =  $self->name =~ /Contig(?:(?:\d+)\.)?(\d+)/;
            return undef;
            ($self->{contig_num}) = $match
        }
        return $self->{contig_num};
    }
}

sub scaffold_num{
    my ($self, $num) = @_;
    if ($num){
        $self->{scaffold_num} = $num;
        $self->update_name;
    }else{
        unless ($self->{scaffold_num}){
            my ($match) = $self->name =~ /Contig(\d+)\.(?:\d+)/;
            return undef unless $match;
            ($self->{scaffold_num}) = $match;
        }
        return $self->{scaffold_num};
    }
}

sub update_name{
    my ($self) = @_;
    return undef unless defined $self->scaffold_num or defined $self->contig_num;
    $self->name("Contig".$self->scaffold_num.".".$self->contig_num);
}


1;

#HeadURL$
#$Id: Sources.pm 31125 2007-12-18 20:38:24Z ebelter $
