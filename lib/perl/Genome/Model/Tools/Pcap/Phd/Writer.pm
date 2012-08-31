# 
# Copyright (C) 2004 Washington University Genome Sequencing Center
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#

# Last modified: <Wed, 2006/11/15 16:26:13 ebelter linus108>

package Genome::Model::Tools::Pcap::Phd::Writer;
our $VERSION = 0.01;

=pod

=head1 NAME

PhdWriter - Phd file writer

=head1 SYNOPSIS

    my $writer = new Genome::Model::Tools::Pcap::Phd::Writer();
    $writer->write(\*STDOUT,$phd);

=head1 DESCRIPTION

Genome::Model::Tools::Pcap::Phd::Writer takes a handle to a phd object and writes it to the given file handle

=head1 METHODS

=cut

use strict;
use warnings;
use Carp;

use Genome::Model::Tools::Pcap::Read;
use Genome::Model::Tools::Pcap::Tag;

my $pkg = 'Genome::Model::Tools::Pcap::Phd::Writer';

=pod

=item new 

    my $writer = new Genome::Model::Tools::Pcap::Phd::Writer;

=cut
sub new {
    croak("$pkg:new:no class given, quitting") if @_ < 1;
    my ($caller, $arg) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = {};
    bless ($self, $class);   

    return $self;
}

=pod

=item Genome::Model::Tools::Pcap::Phd::Reader::write 

    $writer->read(\*STDOUT,$phd);

=cut
sub write {
	my ($self,$OUT,$phd) =@_;
	$self->write_sequence_name($OUT,$phd);
	$self->write_comment($OUT,$phd);
	$self->write_DNA($OUT,$phd);
	$self->write_WRs($OUT,$phd);
	$self->write_tags($OUT,$phd);
	return;
}

sub write_sequence_name
{
    my ($self, $OUT, $phd) = @_;
	print $OUT "BEGIN_SEQUENCE ".$phd->name."\n";

}

sub write_comment {
    my ($self, $OUT, $phd) = @_;
	return unless (defined $phd->comments);
    my $comments = $phd->comments;
	print $OUT "\nBEGIN_COMMENT\n\n";
	while(my ($key, $value) = each(%$comments))
	{       
		print $OUT "$key: ";
		print $OUT "$value" if(defined $value);
		print $OUT "\n";
	}
	print $OUT "\nEND_COMMENT\n";    
}

sub write_DNA {
    my ($self, $OUT, $phd) = @_;
	return unless (defined $phd->unpadded_base_string);
	my $bases = $phd->unpadded_base_string;
	my $quality = $phd->unpadded_base_quality;
	my $positions = $phd->unpadded_chromat_positions;
	print $OUT "\nBEGIN_DNA\n";
	
	for(my $i=0;$i<@{$quality};$i++)
	{
		my $base = substr($bases,$i,1);
		print $OUT "$base $$quality[$i] $$positions[$i]\n";
	}	
	print $OUT "END_DNA\n";
	print $OUT "\nEND_SEQUENCE\n";
}


sub write_WRs {
    my ($self, $OUT, $phd) = @_;
    my @wr;
    if (exists $phd->{wr}) {
        @wr = @{$phd->{wr}};
    }
	else
	{
		return;
	}
    foreach my $wr (@wr)
	{
		print $OUT "\nWR{\n";
		print $OUT "$wr";
		print $OUT "}\n";
	}
}

sub write_tags {
    my ($self, $OUT, $phd) = @_;
    return unless (defined $phd->tags);
	my $tags = $phd->tags;
    my $text = undef;
    
	foreach my $tag (@{$tags})
	{
		print $OUT "\nBEGIN_TAG\n"; 
		print $OUT "TYPE: ".$tag->type."\n";
		print $OUT "SOURCE: ".$tag->source."\n";
		print $OUT "UNPADDED_READ_POS: ".$tag->start." ".$tag->stop."\n";
		print $OUT "DATE: ".$tag->date."\n";
		if($tag->text)
		{
			print $OUT "\nBEGIN_COMMENT\n";
			print $OUT $tag->text;
			print $OUT "\nEND_COMMENT\n";		
		}
		print $OUT "END_TAG\n";	
	}    
}

1;

#$HeadURL$
#$Id$
