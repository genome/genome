package Finishing::Assembly::Phd::Writer;

use strict;
use warnings;
use Carp;

my $pkg = 'Finishing::Assembly::Phd::Writer';

sub new {
    croak("$pkg:new:no class given, quitting") if @_ < 1;
    my ($caller) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = {};
    
    bless ($self, $class);   
    return $self;
}


sub write {
	my ($self,$OUT,$phd) =@_;
    my @methods = qw (
        write_sequence_name  write_comment
        write_DNA   write_WRs   write_tags
    );
    map{$self->$_($OUT, $phd)}@methods;
	return;
}


sub write_sequence_name {
    my ($self, $OUT, $phd) = @_;
	print $OUT "BEGIN_SEQUENCE ".$phd->name."\n";

}


sub write_comment {
    my ($self, $OUT, $phd) = @_;
	return unless (defined $phd->comments);
    my $comments = $phd->comments;
	print $OUT "\nBEGIN_COMMENT\n\n";
    #while(my ($key, $value) = each(%$comments)) {       
    for my $key (sort{$a cmp $b}keys %$comments) {
        my $value = $comments->{$key};
		print $OUT "$key: ";
		print $OUT $value if defined $value;
		print $OUT "\n";
	}
	print $OUT "\nEND_COMMENT\n";    
}


sub write_DNA {
    my ($self, $OUT, $phd) = @_;
	return unless (defined $phd->base_string);
    #my $sequence = $phd->sequence;
	my $bases = $phd->base_string;
	my $quality = $phd->qualities;
	my $positions = $phd->chromat_positions;
	print $OUT "\nBEGIN_DNA\n";
	
	for(my $i=0;$i<@{$quality};$i++) {
		my $base = substr($bases,$i,1);
		print $OUT "$base $$quality[$i] $$positions[$i]\n";
	}	
	print $OUT "END_DNA\n";
	print $OUT "\nEND_SEQUENCE\n";
}


sub write_WRs {
    my ($self, $OUT, $phd) = @_;
    my @wr;
    if ($phd->wr) {
        @wr = @{$phd->wr};
    }
	else {
		return;
	}
    foreach my $wr (@wr) {
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
    
	foreach my $tag (@{$tags}) {
		print $OUT "\nBEGIN_TAG\n"; 
		print $OUT "TYPE: ".$tag->type."\n";
		print $OUT "SOURCE: ".$tag->source."\n";
		print $OUT "UNPADDED_READ_POS: ".$tag->start." ".$tag->stop."\n";
		print $OUT "DATE: ".$tag->date."\n";
		if($tag->text) {
			print $OUT "\nBEGIN_COMMENT\n";
			print $OUT $tag->text;
			print $OUT "\nEND_COMMENT\n";		
		}
		print $OUT "END_TAG\n";	
	}    
}

1;
=pod

=head1 NAME

PhdWriter - Phd file writer

=head1 SYNOPSIS

    my $writer = new Finishing::Phd::Writer();
    $writer->write(\*STDOUT,$phd);

=head1 DESCRIPTION

Finishing::Phd::Writer takes a handle to a phd object and writes it to the given file handle

=head1 METHODS

=cut


=pod

=item new 

    my $writer = new Finishing::Phd::Writer;

=cut

=pod

=item Finishing::Phd::Reader::write 

    $writer->read(\*STDOUT,$phd);

=cut


#$HeadURL$
#$Id$
