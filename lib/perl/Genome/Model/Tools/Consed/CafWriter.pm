package Genome::Model::Tools::Consed::CafWriter;

=pod

=head1 NAME
Genome::Model::Tools::Consed::Caf::Writer - Write a caf file one element at a time

=head1 SYNOPSIS

my $writer = new Genome::Model::Tools::Consed::CafWriter(\*STDOUT);

foreach $contig_hashref (@contigs) {
    $writer->write_object($contig_hashref);
}

=head1 DESCRIPTION

Genome::Model::Tools::Consed::CafWriter takes hashrefs representing elements of an assembly and writes them out
to an ace file.

=cut

use strict;
use warnings;
use Carp;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my $pkg = 'Genome::Model::Tools::Consed::CafWriter';

=pod

=item new 

    my $writer = new Genome::Model::Tools::Consed::CafWriter(\*STDOUT);    

=cut
sub new {
    croak("$pkg:new:no class given, quitting") if @_ < 1;
    my ($caller, $output, %args) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = {};
    bless ($self, $class);
    
	$self->{'output'} = $output;
	$self->{contig_names} = [];

	$self->{align_clip} = delete $args{align_clip}||0;
	$self->{qual_clip} = delete $args{qual_clip}||0;
    return $self;
}

=pod

=item output 

$writer->output(\*FH);

$FH = $writer->output();

=cut
sub output {
    my ($self,$output) = @_;
    if (defined $output) {
        $self->{'output'} = $output;
    }
    return $self->{'output'};
}

=pod

=item write_object 

    %contig = ( type => 'contig',
                name => 'Contig1',
                base_count => length $consensus,
                read_count => $read_count,
                base_seg_count => $base_seg_count,
                u_or_c => 'U',
                consensus => $consensus,
                base_qualities => \@base_qualities,
              );
    $writer->write_object(\%contig);

=cut

sub get_map
{
	my ($self, $contig, $read) = @_;
	#find read start and stop in read units
	#and clip the read start and stop to within the boundary of the contig
	my $contig_start = max($read->start_position,1);
	my $contig_stop = min($read->end_position, $contig->length);
	my $read_start = max(1, $contig_start - $read->start_position+1); 
	my $read_stop = min($read->length, $contig_stop - $read->start_position+1);
	my $read_start_save = $read_start;
	my $read_stop_save = $read_stop;
	
	
	if($self->{qual_clip}==1)
	{
		$read_start = max($read_start, $read->qual_clip_start);
		$read_stop = min($read_stop, $read->qual_clip_end);
	
	}
	if($self->{align_clip} == 1)
	{
		$read_start = max($read_start, $read->align_clip_start);
		$read_stop = min($read_stop, $read->align_clip_end);	
	}
	#find the difference between the read start and read stop
	#and adjust the contig start and stop if necessary
	$contig_start += ($read_start - $read_start_save);
	$contig_stop -= ($read_stop_save - $read_stop);
	
	
	if($read->complemented) #swap parent map
	{
		my $temp = $contig_start;
		$contig_start = $contig_stop;
		$contig_stop = $temp;
		#WTF?
		my $old_read_start = $read_start;
		$read_start = 1+ $read->length - $read_stop;
		
		$read_stop = $read->length - ($old_read_start - 1);
	}
	
	return ($contig_start, $contig_stop, $read_start, $read_stop);

}

sub write_contig
{
	my ($self, $contig, @read_order) = @_;
	my %reads = %{$contig->reads};
	my @reads = values %reads;
	 
	#@reads = sort {$a->position <=> $b->position}@reads; 
	foreach my $read_name (@read_order)
	{
		$self->write_read($reads{$read_name});
	}
	$self->write_contig_header($contig, @read_order);


}
sub write_contig_header {
    my ($self, $contig, @read_order) = @_;
    my $output = $self->output;
	
	
	push @{$self->{contig_names}}, $contig->name;
    $output->print("\n");
    $output->print( "Sequence : ",$contig->name,"\n" );
    $output->print("Is_contig\n");
    $output->print("Padded\n");
	my %reads = %{$contig->reads};
    foreach my $read_name (@read_order) {
		my $read = $reads{$read_name};
        my ($contig_start, $contig_stop, $read_start, $read_stop) = $self->get_map($contig, $read);
        $output->print("Assembled_from ",$read->name," ",$contig_start," ",
            $contig_stop," ",$read_start," ",$read_stop,"\n");
    }
	my $orientation = "U";
	$orientation = "C" if ($contig->complemented);
    $output->print("Orientation ",$orientation,"\n");
    foreach my $bs (@{$contig->base_segments}) {
        if($reads{$bs->{read_name}}->complemented) {
        	$output->print("Golden_path $bs->{read_name} $bs->{end_pos} $bs->{start_pos}\n");
		}
		else {
			$output->print("Golden_path $bs->{read_name} $bs->{start_pos} $bs->{end_pos}\n");
		}
    }
    foreach my $tag (@{$contig->tags}) {
        $tag->scope('');
        $self->write_tag($tag);
    }
    $output->print("\n");
    $self->write_sequence($contig->sequence, $contig->name, 23);
}

sub write_read {
    my ($self, $obj) = @_;
    my $output = $self->output;
    $output->print( "Sequence : ",$obj->name,"\n" );
    $output->print("Is_read\n");
    $output->print($obj->sequence->get_sequence_state eq "padded" ? "Padded\n" : "Unpadded\n");
    $output->print("SCF_File " , $obj->name , "\n");
    $output->print("Template ", $obj->template , "\n")if defined $obj->template;
    $output->print("Dye ", $obj->chemistry , "\n")if defined $obj->chemistry;
    $output->print("Dyetype ", $obj->dyetype , "\n") if defined $obj->dyetype;

    $output->print("Strand ",$obj->strand , "\n")if defined $obj->strand;
    foreach my $tag (@{$obj->tags}) {
        $self->write_tag($tag);
    }
	$obj->align_clip_start (0) if $obj->align_clip_start < 0;
    $obj->align_clip_end (0) if $obj->align_clip_end < 0;
    $obj->qual_clip_start (0) if $obj->qual_clip_start < 0;
    $obj->qual_clip_end ( 0) if $obj->qual_clip_end < 0;
	$output->print("Clipping PQUAL ".$obj->qual_clip_start." ".$obj->qual_clip_end." \"\"\n");
    $output->print("Clipping Phrap ".$obj->align_clip_start." ".$obj->align_clip_end." \"\"\n");
    $output->print("Phd_version " , $obj->phd_version , "\n");
    $output->print("Date " , $obj->date , "\n\n");
    $self->write_sequence($obj->sequence, $obj->name);
    $output->print("\n");
}

sub write_tag {
    my ($self, $tag) = @_;
    my $output = $self->output;
    my $type_code = $self->_get_tag_type_code($tag->type);
    my $start = $tag->start || 0;
    my $stop = $tag->stop || 0;
    my $tag_string = "TYPE: ".$tag->type." SOURCE: ".$tag->source;
    $tag_string .= " SCOPE: ".$tag->scope if $tag->scope;
    $tag_string .= " DATE: ".$tag->date;
    my $tag_text = $tag->text;
    if (defined $tag_text) {
        $tag_text =~ s/\n/ | /g;
        $tag_string .= " BEGIN_COMMENT: $tag_text END_COMMENT";

    }
    if (lc $tag->scope eq 'phd') {
        $tag_string = "BEGIN_CONSED_TAG ".$tag_string." END_CONSED_TAG";
    }
    $output->print("Tag $type_code $start $stop \"$tag_string\"\n");
}

sub _get_tag_type_code {
    my ($self, $type) = @_;
    
    if ($type =~ "[Cc]omment")
    {
        return "COMM";
    }
    elsif ($type eq "oligo")
    {
        return "OLIG";
    }
    elsif ($type eq "Compression")
    {
        return "COMP";
    }
    elsif ($type eq "Polymorphism")
    {
        return "POLYM";
    }
    elsif ($type eq "Edit")
    {
        return "EDIT";
    }
    elsif ($type eq "coordinatorComment")
    {
        return "CoordC";
    }
    elsif ($type eq "coordinatorApproval")
    {
        return "CoordA";
    }
    elsif ($type eq "qualityCoreComment")
    {
        return "QCore";
    }
    elsif ($type eq "stolendata")
    {
        return "STOL";
    }
    elsif ($type eq "stolenconsensus")
    {
        return "CSTOL";
    }
    elsif ($type eq "matchElsewhereLowQual")
    {
        return "MatchLQ";
    }
    elsif ($type eq "matchElsewhereHighQual")
    {
        return "MatchHQ";
    }
    else
    {
        return "COMM";
    }
}



sub write_assembly_tags {
    my ($self, $assembly_tags) = @_;
	
	my $contig_names = $self->{contig_names};
    my $output = $self->output;
    $output->print("Sequence : assembly\n");
    $output->print("Is_assembly\n");
    $output->print("Padded\n");
    foreach my $contig_name (@{$contig_names}) {
        $output->print("Subsequence $contig_name\n");
    }
    foreach my $tag (@{$assembly_tags}) {
        $self->write_tag($tag);
    }
}

sub write_sequence {
    my ($self, $seq, $name, $qual_width) = @_;
    $qual_width = 50 if !defined $qual_width;
    my $output = $self->output;
    $output->print("DNA : $name\n");
	my $bases;
	if($seq->get_sequence_state eq "padded")
	{
    	$bases = $seq->padded_base_string;
		$bases =~ tr/*/-/;
	}
	else
	{
		$bases = $seq->unpadded_base_string;
	}
    $bases = uc $bases;
    my $length = $seq->length;
    my $index = 0;
    while ($index < $length) {
        $output->print(substr($bases, $index, 50),"\n");
        $index += 50;
    }
    $output->print("\n");

    my @quality = ();
	if($seq->get_sequence_state eq "padded")
	{
    	@quality = @{$seq->padded_base_quality()} if defined $seq->padded_base_quality();
	}
	else
	{
		@quality = @{$seq->unpadded_base_quality()} if defined $seq->unpadded_base_quality();
	}
	foreach my $qual (@quality)
	{
		$qual = 0 if($qual eq '*');
	}
    if (@quality > 0) {
        $output->print("BaseQuality : $name\n");
        $index = 0;
        my $stop;
        if (scalar @quality != $length) {
            print "$name: " , scalar @quality , " <> $length\n";
        }
        while ($index < $length) {
            $stop = $index + ($qual_width - 1);
            $stop = $length - 1 if $stop >= $length;
            $output->print("@quality[$index..$stop] \n");
            $index += $qual_width;
        }
        $output->print("\n");
        my @positions = ();
		if($seq->get_sequence_state eq "padded")
		{
        	@positions = @{$seq->padded_chromat_positions()} if defined $seq->padded_chromat_positions;
			@positions = @{$seq->replace_chromat_pads(\@positions)} if scalar @positions;
		}
		else
		{
			@positions = @{$seq->unpadded_chromat_positions()} if defined $seq->unpadded_chromat_positions;
		}		
        if (@positions > 0) {
            $output->print("BasePosition : $name\n");
            $index = 0;
            while ($index < $length) {
                $stop = $index + 49;
                $stop = $length - 1 if $stop >= $length;
                $output->print("@positions[$index..$stop] \n");
                $index += 50;
            }
        }
    }
}
1;

