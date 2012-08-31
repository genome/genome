package Genome::Model::Tools::Pcap::Sources::Ace::Read;
our $VERSION = 0.01;

use strict;
use warnings;
use Carp;

use Storable;
use base (qw(Genome::Model::Tools::Pcap::Sources::Ace::SequenceItem));

=cut
Contig Data Structure
Contig:
    sequence (bases and quality)
    children 
    get_child
    padded_base_count
    base_count
    read_count
    complemented
    tags
=cut

sub new 
{
    croak("__PACKAGE__:new:no class given, quitting") if @_ < 1;
    my ($caller, %params) = @_; 
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = $class->SUPER::new(%params);

    return $self;
}

sub padded_base_string
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);      
    my $fh = $self->{fh};
    my $index = $self->{index};
    my $name = $self->{name};
    $fh->seek($index->{reads}{$name}{read}{offset},0);    
    my $sequence;
    my $line = <$fh>;
    while ($line = <$fh>) 
    {
        if ($line =~ /^\s*$/) 
        {
            last;
        }
        chomp $line;
        $sequence .= $line;
    }
	$object->{just_load} = 1;
    $object->padded_base_string($sequence);
	$object->{just_load} = 0;
    return 1;
}

sub padded_base_quality
{
    return 0;
}

sub unpadded_base_string
{
    return 0;
}

sub unpadded_base_quality
{
    return 0;
}

sub has_alignment
{
    return "padded_base_string";
}


sub get_map {
    my ($self) = @_;
    return { self => [ 'self' ],
			 read_position   => [ 'name', 'complemented', 'position' ],
             read => [ 'name', 'padded_base_count', 'length', 'info_count', 'tags' ],             
             qa => [ 'qual_clip_start', 'qual_clip_end', 'align_clip_start', 'align_clip_end'],
             ds => [ 'chromat_file', 'phd_file', 'time', 'chemistry', 'dye' ],
             padded_base_string   => [ 'padded_base_string', 'unpadded_base_string', 'alignment' ],
			 unpadded_base_string => [ 'unpadded_base_string', 'padded_base_string' ]
           };
}

#methods inherited from Item Source

sub position 
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->{index};
    my $reader = $self->{reader};
	my $name = $self->{name};
    $fh->seek($index->{reads}{$name}{read_position}{offset},0);
    my $obj = $reader->next_object;
	$object->{just_load} = 1;
    $object->position($obj->{position});
    $object->complemented($obj->{u_or_c} =~ /c/i or 0) if !$object->already_loaded("complemented");
	$object->{just_load} = 0;
    return 1;    
}


sub length 
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->{index};
    my $reader = $self->{reader};
    my $name = $self->{name};
    $fh->seek($index->{reads}{$name}{read}{offset},0);
    my $line = <$fh>;
    chomp $line;
    my @tokens = <$fh>;
    my %ret_val = (
        padded_base_count => $tokens[1],
        info_count => $tokens[2],
        tag_count => $tokens[3]);
    $object->length($ret_val{padded_base_count});
    $object->info_count($ret_val{info_count}) if !$object->already_loaded("info_count"); 
    return 1;   
}

sub _build_read_tag {
    my ($self, $obj) = @_;
    my $tag = new Genome::Model::Tools::Pcap::Tag(
        type => $obj->{tag_type},
        date => $obj->{date},
        source => $obj->{program},
        parent => $obj->{read_name},
        scope => 'ACE',
        start => $obj->{start_pos},
        stop => $obj->{end_pos},
    );
    return $tag;
}


sub tags
{
    my ($self,$object) = @_; 
    return 1 unless (@_ > 1); 
    my $input = $self->{fh};
    my $index = $self->{index};
    my $reader = $self->{reader};
    my $name = $self->{name};
    my @tags = @{$index->{reads}{$name}{read}{read_tags}};  
    my @read_tags;
    foreach my $tag_index (@tags)
    {
        $input->seek($tag_index->{offset},0);
        my $read_tag = $self->_build_read_tag($reader->next_object);
        push @read_tags, $read_tag;    
    }
	$object->{just_load} = 1;
    $object->tags( \@read_tags );
	$object->{just_load} = 0;
    return 1;
}

sub copy
{
    my ($self) = @_;
    return dclone ($self);
}

sub copy_tag 
{
    my ($self, $tag) = @_;
    return dclone $tag;
}

sub start_position
{
    my ($self) = @_;
    return $self->position;
    return 1;    
}

sub end_position 
{
    my ($self) = @_;
    return $self->position + $self->length;
    return 1;
}

#methods inherited from SequenceItem Source

sub info_count
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->{index};
    my $reader = $self->{reader};
    $fh->seek($index->{reads}{$self->{name}}{read}{offset},0);
    my $line = <$fh>;
    chomp $line;
    my @tokens = split / /,$line;
    my %ret_val = (
        padded_base_count => $tokens[1],
        info_count => $tokens[2],
        tag_count => $tokens[3]);
	$object->{just_load} = 1;
    $object->length($ret_val{padded_base_count}) if !$object->already_loaded("length"); ;
    $object->info_count($ret_val{info_count});
	$object->{just_load} = 0;
    return 1;
}

sub chromat_file 
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->{index};
    my $reader = $self->{reader};
    $fh->seek($index->{reads}{$self->{name}}{read}{ds}{offset},0);
    my $line = <$fh>;
    chomp $line;
     $line =~ s/ (\w+): /|$1|/g; #delimit the key-value pairs
    my @tokens = split(/\|/, $line); 
    shift @tokens; # drop the DS tag
    my %description = @tokens;
    $object->{just_load} = 1;
    $object->chromat_file($description{CHROMAT_FILE});
    $object->phd_file($description{PHD_FILE}) if !$object->already_loaded("phd_file"); ;
    $object->time($description{TIME}) if !$object->already_loaded("time");
    $object->chemistry($description{CHEM}) if !$object->already_loaded("chemistry");
    $object->dyetype($description{DYE}) if !$object->already_loaded("dye");
	$object->{just_load} = 0;
    return 1;    
}

sub phd_file 
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->{index};
    my $reader = $self->{reader};
    $fh->seek($index->{reads}{$self->{name}}{read}{ds}{offset},0);
    my $line = <$fh>;
    chomp $line;
     $line =~ s/ (\w+): /|$1|/g; #delimit the key-value pairs
    my @tokens = split(/\|/, $line); 
    shift @tokens; # drop the DS tag
    my %description = @tokens;
	
    $object->{just_load} = 1;
    $object->chromat_file($description{CHROMAT_FILE}) if !$object->already_loaded("chromat_file");
    $object->phd_file($description{PHD_FILE});
    $object->time($description{TIME}) if !$object->already_loaded("time");
    $object->chemistry($description{CHEM}) if !$object->already_loaded("chemistry");
    $object->dyetype($description{DYE}) if !$object->already_loaded("dye");
	$object->{just_load} = 0;
    return 1;    
}

sub time
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->{index};
    my $reader = $self->{reader};
    $fh->seek($index->{reads}{$self->{name}}{read}{ds}{offset},0);
    my $line = <$fh>;
    chomp $line;
     $line =~ s/ (\w+): /|$1|/g; #delimit the key-value pairs
    my @tokens = split(/\|/, $line); 
    shift @tokens; # drop the DS tag
    my %description = @tokens;
    
    $object->{just_load} = 1;
	$object->chromat_file($description{CHROMAT_FILE}) if !$object->already_loaded("chromat_file");;
    $object->phd_file($description{PHD_FILE}) if !$object->already_loaded("phd_file");
    $object->time($description{TIME});
    $object->chemistry($description{CHEM}) if !$object->already_loaded("chemistry");
    $object->dyetype($description{DYE}) if !$object->already_loaded("dyetype");
	$object->{just_load} = 0;
    return 1;    
}

# aa01a01.x1 -> forward, thermoseq
# aaa01a01.g3 ->
# M_BB0392D19PCR12c14_g14.b1 -> pcr 
# M_BB0392D19PCR12c14_g14e2.b1 -> pcr 
# M_BB0392D19PCR12c14_14.b4e1 -> pcr 
# /_(\w+)(\d+)$/

sub chemistry
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->{index};
    my $reader = $self->{reader};
    $fh->seek($index->{reads}{$self->{name}}{read}{ds}{offset},0);
    my $line = <$fh>;
    chomp $line;
     $line =~ s/ (\w+): /|$1|/g; #delimit the key-value pairs
    my @tokens = split(/\|/, $line); 
    shift @tokens; # drop the DS tag
    my %description = @tokens;
    
    $object->{just_load} = 1;
	$object->chromat_file($description{CHROMAT_FILE}) if !$object->already_loaded("chromat_file");;
    $object->phd_file($description{PHD_FILE}) if !$object->already_loaded("phd_file"); ;
    $object->time($description{TIME}) if !$object->already_loaded("time");
    $object->chemistry($description{CHEM});
    $object->dyetype($description{DYE}) if !$object->already_loaded("dye");
	$object->{just_load} = 0;
    return 1;    




}

sub dyetype
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->{index};
    my $reader = $self->{reader};
    $fh->seek($index->{reads}{$self->{name}}{read}{ds}{offset},0);
    my $line = <$fh>;
    chomp $line;
     $line =~ s/ (\w+): /|$1|/g; #delimit the key-value pairs
    my @tokens = split(/\|/, $line); 
    shift @tokens; # drop the DS tag
    my %description = @tokens;
    
    $object->{just_load} = 1;
	$object->chromat_file($description{CHROMAT_FILE}) if !$object->already_loaded("chromat_file");;
    $object->phd_file($description{PHD_FILE}) if !$object->already_loaded("phd_file"); ;
    $object->time($description{TIME}) if !$object->already_loaded("time");
    $object->chemistry($description{CHEM}) if !$object->already_loaded("chemistry");
    $object->dyetype($description{DYE});
	$object->{just_load} = 0;
    return 1;    


}
sub complemented
{
    my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->{index};
    my $reader = $self->{reader};
    $fh->seek($index->{reads}{$self->{name}}{read_position}{offset},0);
    my $obj = $reader->next_object;
	
	$object->{just_load} = 1;
    $object->position($obj->{position}) if !$object->already_loaded("position");
    $object->complemented($obj->{complemented});
	$object->{just_load} = 0;
    return 1;
}

sub align_clip_start 
{
	my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->{index};
    my $reader = $self->{reader};
    $fh->seek($index->{reads}{$self->{name}}{read}{qa}{offset},0);
    my $line = <$fh>;
    chomp $line;
    my @tokens = split / /,$line;
	shift @tokens;
	
	$object->{just_load} = 1;
    $object->qual_clip_start($tokens[0]) if !$object->already_loaded("qual_clip_start");
    $object->qual_clip_end($tokens[1]) if !$object->already_loaded("qual_clip_end");
    $object->align_clip_start($tokens[2]);
    $object->align_clip_end($tokens[3]) if !$object->already_loaded("align_clip_end");
	$object->{just_load} = 0;
    return 1;	
}

sub align_clip_end 
{
	my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->{index};
    my $reader = $self->{reader};
    $fh->seek($index->{reads}{$self->{name}}{read}{qa}{offset},0);
    my $line = <$fh>;
    chomp $line;
    my @tokens = split / /,$line;
	shift @tokens;
	
	$object->{just_load} = 1;
    $object->qual_clip_start($tokens[0]) if !$object->already_loaded("qual_clip_start");
    $object->qual_clip_end($tokens[1]) if !$object->already_loaded("qual_clip_end");
    $object->align_clip_start($tokens[2]) if !$object->already_loaded("align_clip_start");
    $object->align_clip_end($tokens[3]);
	$object->{just_load} = 0;
    return 1;		
}

sub qual_clip_start 
{
	my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->{index};
    my $reader = $self->{reader};
    $fh->seek($index->{reads}{$self->{name}}{read}{qa}{offset},0);
    my $line = <$fh>;
    chomp $line;
    my @tokens = split / /,$line;
	shift @tokens;
	
	$object->{just_load} = 1;
    $object->qual_clip_start($tokens[0]);
    $object->qual_clip_end($tokens[1]) if !$object->already_loaded("qual_clip_end");
    $object->align_clip_start($tokens[2]) if !$object->already_loaded("align_clip_start");
    $object->align_clip_end($tokens[3]) if !$object->already_loaded("align_clip_end");
	$object->{just_load} = 0;
    return 1;		
}

sub qual_clip_end 
{
	my ($self,$object) = @_;
    return 1 unless (@_ > 1);
    my $fh = $self->{fh};
    my $index = $self->{index};
    my $reader = $self->{reader};
    $fh->seek($index->{reads}{$self->{name}}{read}{qa}{offset},0);
    my $line = <$fh>;
    chomp $line;
    my @tokens = split / /,$line;
	shift @tokens;
	
	$object->{just_load} = 1;
    $object->qual_clip_start($tokens[0]) if !$object->already_loaded("qual_clip_start");
    $object->qual_clip_end($tokens[1]);
    $object->align_clip_start($tokens[2]) if !$object->already_loaded("align_clip_start");
    $object->align_clip_end($tokens[3]) if !$object->already_loaded("align_clip_end");
	$object->{just_load} = 0;
    return 1;		
}

sub _index
{
	my ($self) = @_;
    return $self->{index};

}

1;
