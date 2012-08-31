package Finishing::Assembly::AssembledRead;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::SequencedItem';

my %phd_v :name(phd_v:p); #phd_version, used by sub phd_version below, maybe move this down?

#- OPERATIONS -#

sub complement : CUMULATIVE
{
    my $self = shift;

    unless ($self->qual_clip_start == -1) {
        my $q_start = $self->qual_clip_start;
        $self->qual_clip_start($self->length - $self->qual_clip_stop + 1);
        $self->qual_clip_stop($self->length - $q_start + 1);
    }

    unless ($self->align_clip_start == -1) {
        my $a_start = $self->align_clip_start;
        $self->align_clip_start($self->length - $self->align_clip_stop + 1);
        $self->align_clip_stop($self->length - $a_start + 1);
    }

    my $ReadTags = $self->tags;
    return 1 unless $ReadTags;
		
	for my $tag (@$ReadTags) {
		my $nReadLength = $self->length;
		$tag->start( $nReadLength - $tag->start + 1); 
		$tag->stop($nReadLength - $tag->stop + 1);
		#now swap start and end positions
		my $temp = $tag->stop;
		$tag->stop($tag->start);
		$tag->start($temp); 	
	}
	$self->tags($ReadTags);

    return 1;
}


#- DERIVED NAME ATTRIBUTES -#
sub chemistry 
{
    my $self = shift;

    my $code = $self->code;

    return '454' if $self->chromat_file =~ /\.sff\:/;

    return ( grep { $self->code eq $_ } (qw/ s q r /) )
    ? 'Dye_primer'
    : 'Dye_terminator';            
}

sub dyetype 
{
    my $self = shift;

    return $self->dye;
}

sub code 
{
    my $self = shift;

    return ($self->name =~ /\.(.)\d*$/)[0];
}

sub strand
{
    my $self = shift;

    my $code = $self->code;

    return ( grep { $code eq $_ } (qw/ r y g/) )
    ? 'Reverse'
    : 'Forward';
}

#- POSITION/ALIGNMENT -#
sub start
{
    my $self = shift;
    
	return $self->position(@_);
}

sub align_start
{
    my $self = shift;
    
	return $self->position(@_);
}

sub stop
{
    my $self = shift;

	return $self->position + $self->length - 1;
}

sub align_stop
{
    stop(@_);
}

sub align_end
{
    stop(@_);
}

#- PHD -#
sub unpadded_chromat_positions
{
    my $self = shift;
    #$self->warn_msg("unpadded_chromat_positions is deprecated, use chromat_positions");
    return $self->chromat_positions(@_);
}

sub template  #TODO do the same for lib and primer?
{
    my $self = shift;
    
    if ( my $method = $self->proxy->get_method('template', @_) )
    {
        my $template = $method->();
        return $template if $template;
    }

    # determine template from read name
    my $template;
    if (($template) = $self->name =~ /^(.*?)\./)
    {
    #my ($template) = $self->name =~ /^(.*?)\./;
	if ($template=~/((_\d+)|(_\w\d*)|(_\w\d+\w\d+))$/)
	{
	    my $extra = $1;
	    $template = ~s/$extra//;
	}
    }
    else
    {
	#454 read
	$template = $self->name;
    }

    return $template;
}

#below here check to make sure wrks with object model
=pod

=head1 ace_read

$read->ace_read($ace_read)           

This is a getter/setter that takes or returns a copy of the read hash in the same format as the read hash that is produced by Finishing::Assembly::Ace::Reader/Writer.  This method helps to provide compatibility with the low-level Ace Reader/Writer.

=cut
sub ace_read
{
    my ($self, $ace_read) = @_;
    if(@_>1)
    {
        $self->sequence->padded_base_string( $ace_read->{sequence});
        $self->name ( $ace_read->{name});
        $self->align_clip_start ($ace_read->{align_clip_start});
        $self->align_clip_end ($ace_read->{align_clip_end});
        $self->qual_clip_start ($ace_read->{qual_clip_start});
        $self->qual_clip_end ( $ace_read->{qual_clip_end});
        $self->info_count ($ace_read->{info_count}); 
        $self->chromat_file ( $ace_read->{description}{CHROMAT_FILE});
        $self->phd_file ($ace_read->{description}{PHD_FILE});
		$self->dye($ace_read->{description}{DYE});
        $self->time ($ace_read->{description}{TIME});
		$self->chemistry ($ace_read->{description}{CHEM}) if defined $ace_read->{description}{CHEM};
    }                                         
    
    return { type => "read",
             name => $self->name,
             padded_base_count => $self->length,
             info_count => $self->info_count,
             tag_count => scalar(@{$self->tags}),
             sequence => $self->sequence->padded_base_string,
             qual_clip_start => $self->qual_clip_start,
             qual_clip_end => $self->qual_clip_end,
             align_clip_start => $self->align_clip_start,
             align_clip_end => $self->align_clip_end,
             description => {       
                CHROMAT_FILE => $self->chromat_file,
                PHD_FILE => $self->phd_file,
                TIME => $self->time,
				DYE => $self->dye,
				CHEM => $self->chemistry }};
}

=pod

=head1 ace_read_position

$read->ace_read_position($ace_read_position)           

This is a getter/setter that takes or returns a copy of the read_position hash in the same format as the read_position hash that is produced by Finishing::Assembly::Ace::Reader/Writer.  This method helps to provide compatibility with the low-level Ace Reader/Writer.

=cut

sub ace_read_position
{
    my ($self, $ace_read_position) = @_;
    if(@_>1)
    {
        $self->position ($ace_read_position->{position});
        $self->complemented ($ace_read_position->{u_or_c} =~ /c/i or 0);
    }
    
    return { type => "read_position",
             read_name => $self->name,
             position => $self->position,
             u_or_c => ( $self->complemented ? "C" : "U" )
           };
}

sub scf_file
{
    chromat_file(@_);
}

sub date {
    my ($self, $date) = @_;
    if (defined $date) {
        $date =~ s/(\w*)  (\d*)/$1 $2/;
        $self->time( $date);
    }
    return $self->time;
}

#this is used as a getter for ace2caf, it will need some work to
#also work as a setter
sub phd_version {
    my ($self, $phd_version) = @_;
    if (defined $phd_version) {
        $self->phd_v( $phd_version );
    }
	elsif(!defined $self->phd_v){
		my $phd_file = $self->phd_file;
    	$phd_file =~ /\.phd\.(.*)/;
    	$self->phd_v( $1 );
	}
    return $self->phd_v;
}

sub get_child_position_from_parent_position
{
	my ($self, $parent_pos) = @_;
	return ( $parent_pos + ( 1 - $self->position ) );    #removed -1 at the end
}

sub get_parent_position_from_child_position
{
	my ($self, $child_pos) = @_;
	return ( ( $self->position - 1 ) + $child_pos );    #removed +1 at the end
}

1;

#$HeadURL$
#$Id$
