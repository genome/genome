package Genome::Model::Tools::Pcap::Read;
our $VERSION = 0.01;

use strict;
use warnings;
use Carp;
use Genome::Model::Tools::Pcap::Transform;
use Genome::Model::Tools::Pcap::Tag;

use Genome::Model::Tools::Pcap::SequenceItem;
use Storable;

use base (qw(Genome::Model::Tools::Pcap::SequenceItem));



=pod

=head1 NAME

Read - Class that manages a read's data.

=head1 DESCRIPTION

Genome::Model::Tools::Pcap::Read

=head1 METHODS

=cut


=pod

=head1 new

$read = Genome::Model::Tools::Pcap::Read->new(children => \%children, tags => \@tags, position => $position, length => $length, ace_read => $ace_read, ace_read_position => $ace_read_position);
    
children - optional, a hash containing the items children, indexed by the child names.

tags - optional, an array of tags belonging to the item.

position - optional, the position of the item in the parent string.

length - optional, the length of the item in the parent item in padded base units.

ace_read - optional, will take a read hash as defined by Genome::Model::Tools::Pcap::Ace::Reader, and populate the fields of the read object.

ace_read_position - optional, will populate the read object with the data contained in an read_position.  This hash is produced by Genome::Model::Tools::Pcap::Ace::Reader.

=cut
sub new 
{
    croak("__PACKAGE__:new:no class given, quitting") if @_ < 1;
    my ($caller, %params) = @_;
#    my ($caller, $contig_hash, $reads, $contig_tags, $base_segments) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    
    my $self = $class->SUPER::new(%params);
        
    my $ace_read;
    my $ace_read_position;
	
	if(exists $params{callbacks})
	{
		$self->callbacks($params{callbacks});
	}
	
    if(exists $params{ace_read})
    {
        $ace_read = $params{ace_read};
    }
    if(exists $params{ace_read_position})
    {
        $ace_read_position = $params{ace_read_position}; 
    }    
    
    
    $self->ace_read ($ace_read) if (defined ($ace_read));
    $self->ace_read_position ($ace_read_position) if (defined ($ace_read_position));
    
    return $self;
}

sub align_clip_start 
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
	
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}

sub align_clip_end 
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
	
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}

sub qual_clip_start 
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
	
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}

sub qual_clip_end 
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
	
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}

sub info_count 
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
	
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}

sub chromat_file 
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
	
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}

sub phd_file 
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
	
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}

sub time
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
	
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}

sub dyetype {
    my ($self, $value) = @_;
    return $self->dye($value);
}

sub dye {
    my ($self, $value) = @_;
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;
    
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);    
}

#this is after the second underscore in the read name
sub primer {
     my ($self, $value) = @_;
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;
    
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);    
}

sub chemistry {
    my ($self, $value) = @_;
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;
    
    if (!defined $value and !$self->already_loaded($name)) {
        my $code = $self->code;

        if ($code eq 's' or $code eq 'q' or $code eq 'r') {
            $self->check_and_load_data($name, 'Dye_primer');
        }
        else {
            $self->check_and_load_data($name, 'Dye_terminator');            
        }
    }
    elsif (defined $value) {
        if (lc $value eq 'unknown') {
            $value = 'Unknown_dye';
        }
        elsif (lc $value eq 'term') {
            $value = 'Dye_terminator';
        }
        $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);  
}

sub strand {
    my ($self, $value) = @_;
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;
    if (!defined $value) {
        if (!$self->already_loaded($value) && !defined $self->check_and_load_data($value)) {
            my $code = $self->code;
            if ($code eq 'r' or $code eq 'y' or $code eq 'g') {
                $self->check_and_load_data($name, 'Reverse');
            }
            else {
                $self->check_and_load_data($name,'Forward');
            }
        }
        return $self->check_and_load_data($name);
    }
    else {
        return $self->check_and_load_data($name, $value);
    }
}

#used by phd and phddb

sub comments 
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
	
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}

sub wr 
{
	my ($self, $value) = @_;
	my ($name) = (caller(0))[3] =~ /.+::(.+)/;
	
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}

#need to edit
sub template {
    my ($self, $template) = @_;
    if (defined $template) {
        $self->{template} = $template;
    }
    elsif (!exists $self->{template} or !defined $self->{template}) { 
        # determine template from read name
        my ($template)=($self->name=~/^(.*?)\./);
        if (defined $template && $template=~/((_\d+)|(_\w\d*)|(_\w\d+\w\d+))$/) {
            my $extra=$1;
            $template=~s/$extra//;
        }
        $self->{template} = $template;
    }
    return $self->{template};
}

sub code {
    my ($self, $code) = @_;
    if (!defined $code) {
        if (!exists $self->{code}) {
            $self->name =~ /\.(.)\d*$/;
            $self->{code} = $1;
        }
        return $self->{code};
    }
    else {
        $self->{code} = $code;
    }
}

=pod

=head1 ace_read
    
$read->ace_read($ace_read)           

This is a getter/setter that takes or returns a copy of the read hash in the same format as the read hash that is produced by Genome::Model::Tools::Pcap::Ace::Reader/Writer.  This method helps to provide compatibility with the low-level Ace Reader/Writer.

=cut
sub ace_read
{
    my ($self, $ace_read) = @_;
    if(@_>1)
    {
        $self->padded_base_string( $ace_read->{sequence});
        $self->name ( $ace_read->{name});
        $self->align_clip_start ($ace_read->{align_clip_start});
        $self->align_clip_end ($ace_read->{align_clip_end});
        $self->qual_clip_start ($ace_read->{qual_clip_start});
        $self->qual_clip_end ( $ace_read->{qual_clip_end});
        $self->info_count ($ace_read->{info_count}); 
        $self->chromat_file ( $ace_read->{description}{CHROMAT_FILE});
        $self->phd_file ($ace_read->{description}{PHD_FILE});
		$self->dye ($ace_read->{description}{DYE});
        $self->time ($ace_read->{description}{TIME});
		$self->chemistry ($ace_read->{description}{CHEM}) if defined $ace_read->{description}{CHEM};
    }                                         
    
    return { type => "read",
             name => $self->name,
             padded_base_count => $self->length,
             info_count => $self->info_count,
             tag_count => scalar(@{$self->tags}),
             sequence => $self->padded_base_string,
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

This is a getter/setter that takes or returns a copy of the read_position hash in the same format as the read_position hash that is produced by Genome::Model::Tools::Pcap::Ace::Reader/Writer.  This method helps to provide compatibility with the low-level Ace Reader/Writer.

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

=pod

=head1 base_count

$read->base_count           
 
This is returns the number of bases in unpadded units.   

=cut

sub base_count
{
    my ($self) = @_;

    return length $self->unpadded_base_string;
}

=pod

complemented
    $read->complemented           

    Getter/Setter for complemented boolean value.
=cut

sub complemented
{
	my ($self, $value) = @_;
	
    my ($name) = (caller(0))[3] =~ /.+::(.+)/;
    if(@_>1)
    {   
       return $self->check_and_load_data($name, $value);
    }
    return $self->check_and_load_data($name);
}

sub scf_file
{
    my ($self, $scf_file) = @_;
    
    $self->chromat_file( $scf_file) if defined $scf_file;
    
    return $self->chromat_file;
}

sub date {
    my ($self, $date) = @_;
    if (defined $date) {
        $date =~ s/(\w*)  (\d*)/$1 $2/;
        $self->time( $date);
    }
    return $self->time;
}

sub align_start
{
	my ($self) = @_;
	return $self->position;

}

sub align_end
{
	my ($self) = @_;
	return $self->position + $self->length - 1;
}
#this is used as a getter for ace2caf, it will need some work to
#also work as a setter
sub phd_version {
    my ($self, $phd_version) = @_;
    if (defined $phd_version) {
        $self->{phd_version} = $phd_version;
    }
	elsif(!defined $self->{phd_version}){
		my $phd_file = $self->phd_file;
    	$phd_file =~ /\.phd\.(.*)/;
    	$self->{phd_version} = $1;
	}
    return $self->{phd_version};
}

=pod

=head1 copy_tag

my $read_tag = $read->copy_tag($read_tag);           

Returns a copy of the read tag.

=cut

=pod

=head1 copy

my $read_copy = $read->copy($read);           

Returns a deep copy of the Read.

=cut

1;
