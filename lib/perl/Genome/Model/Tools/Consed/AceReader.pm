package Genome::Model::Tools::Consed::AceReader;

use strict;
use warnings;

use Data::Dumper;

class Genome::Model::Tools::Consed::AceReader {
    has => [
        file => {
            is => 'Text',
            doc => 'Ace file',
        },
        _fh => { is_optional => 1, },
        _previous_contig => { is_optional => 1, },
    ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $fh = eval{ Genome::Sys->open_file_for_reading($self->file); };
    if ( not $fh ) {
        $self->error_message("Failed to open ace file (".$self->file."): $@");
        return;
    }
    $self->_fh($fh);

    return $self;
}

sub next {
    my $self = shift;

    my $ret_val;
    while ( my $line = $self->_fh->getline ) {
        next if $line =~ /^\s*$/;
        chomp $line;
        my @tokens = split(/[ {]/,$line);
        next if not @tokens;
        my $type = shift @tokens;
        my $method = '_build_'.$type;
        return $self->$method($self->_fh, \@tokens);
    }

    return;
}

sub next_contig {
    my $self = shift;

    my $contig = $self->_previous_contig;
    $self->_previous_contig(undef);
    if ( not $contig ) {
        while ( $contig = $self->next ) {
            last if $contig->{type} eq 'contig';
        }
    }
    return if not $contig;

    $contig->{unpadded_consensus} = $contig->{consensus};
    $contig->{unpadded_consensus} =~ s/\*//g;

    CTG: while ( my $obj = $self->next ) {
        if ( $obj->{type} eq 'contig' ) {
            $self->_previous_contig($obj);
            last CTG;
        }
        elsif ( $obj->{type} eq 'base_segment' ) {
            push @{$contig->{base_segments}}, $obj;
        }
        elsif ( $obj->{type} eq 'read_position' ) {
            $contig->{reads}->{ $obj->{read_name} }->{position} = $obj->{position};
            $contig->{reads}->{ $obj->{read_name} }->{u_or_c} = $obj->{u_or_c};
        }
        elsif ( $obj->{type} eq 'read' ) {
            for my $attr (qw/ position u_or_c /) {
                $obj->{$attr} = $contig->{reads}->{ $obj->{name} }->{$attr};
            }
            $contig->{reads}->{ $obj->{name} } = $obj;
            $contig->{reads}->{ $obj->{name} }->{start} = $contig->{reads}->{ $obj->{name} }->{position};
            $contig->{reads}->{ $obj->{name} }->{stop} = $contig->{reads}->{ $obj->{name} }->{start} + length($obj->{sequence}) - 1;
        }
        elsif ( $obj->{type} eq 'read_tag' ) {
            push @{$contig->{reads}->{ $obj->{read_name} }->{tags}}, $obj;
        }
        elsif ( $obj->{type} eq 'contig_tag' ) {
            push @{$self->{_contig_tags}}, $obj;
        }
        elsif ( $obj->{type} eq 'assembly_tag' ) {
            push @{$self->{_assembly_tags}}, $obj;
        }
        else {
            print Dumper({UNKNOWN=>$obj});
        }
    }

    return $contig;
}

sub assembly_tags {
    my $self = shift;
    my $ctg;
    do { $ctg = $self->next_contig } until not defined $ctg;
    return $self->{_assembly_tags} || [];
}

sub contig_tags {
    my $self = shift;
    my $ctg;
    do { $ctg = $self->next_contig } until not defined $ctg;
    return $self->{_contig_tags} || [];
}

sub _build_AS {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val = (
        type => 'assembly',
        contig_count => $token_ary_ref->[0],
        read_count => $token_ary_ref->[1],
    );
    return \%ret_val;
}

sub _build_CO {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val = (
        type => 'contig',
        name => $token_ary_ref->[0],
        base_count => $token_ary_ref->[1],
        read_count => $token_ary_ref->[2],
        base_seg_count => $token_ary_ref->[3],
        u_or_c => $token_ary_ref->[4],
    );

    my $consensus;

    my $line;
    while ($line = <$IN>) {
        if ($line =~ /^\s*$/) {
            last;
        }
        chomp $line;
        $consensus .= $line;
    }
    $ret_val{'consensus'} = $consensus;
    while ($line = <$IN>) {
        if ($line =~ /^BQ/) {
            last;
        }
    }
    my @bq;
    while ($line = <$IN>) {
        if ($line =~ /^\s*$/) {
            last;
        }
        chomp $line;
        $line =~ s/^ //; # get rid of leading space
        push @bq, split(/ /,$line);
    }
    $ret_val{'base_qualities'} = \@bq;
    return \%ret_val;
}

sub _build_AF {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val = (
        type => 'read_position',
        read_name => $token_ary_ref->[0],
        u_or_c => $token_ary_ref->[1],
        position => $token_ary_ref->[2],
    );
    return \%ret_val;
}

sub _build_BS {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val = (
        type => 'base_segment',
        start_pos => $token_ary_ref->[0],
        end_pos => $token_ary_ref->[1],
        read_name => $token_ary_ref->[2],
    );
    return \%ret_val;
}

sub _build_RD {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val = (
        type => 'read',
        name => $token_ary_ref->[0],
        padded_base_count => $token_ary_ref->[1],
        info_count => $token_ary_ref->[2],
        tag_count => $token_ary_ref->[3],
    );
    my $sequence;
    my $line;
    while ($line = <$IN>) {
        if ($line =~ /^\s*$/) {
            last;
        }
        chomp $line;
        $sequence .= $line;
    }
    #my ($seq, $pads) = $self->un_pad_sequence($sequence);
    $ret_val{'sequence'} = $sequence;
    #$ret_val{'pads'} = $pads;
    while ($line = <$IN>) {
        chomp $line;
        if ($line =~ /^QA/) { my @tokens = split(/ /,$line); $ret_val{'qual_clip_start'} = $tokens[1]; $ret_val{'qual_clip_end'} = $tokens[2];
            $ret_val{'align_clip_start'} = $tokens[3];
            $ret_val{'align_clip_end'} = $tokens[4];
        }
        elsif ($line =~ /^DS/) {
            $line =~ s/ (\w+): /|$1|/g; #delimit the key-value pairs
            my @tokens = split(/\|/, $line); 
            shift @tokens; # drop the DS tag
            my %description = @tokens;
            $ret_val{'description'} = \%description;
            last;
        }
    }

    return \%ret_val;
}

sub _build_WA {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val;
    $ret_val{'type'} = 'assembly_tag';
    my $line = <$IN>;
    chomp $line;
    $line =~ s/^\s*// if $line =~ /\w/;
    my @tmp = split(/ /, $line);
    if ( scalar @tmp == 2 ) {
        @ret_val{'tag_type', 'date'} = @tmp;
    } elsif ( scalar @tmp == 3 ) {
        @ret_val{'tag_type', 'program', 'date'} = @tmp;
    }
    my $data;
    while ($line = <$IN>) {
        chomp $line;
        $line =~ s/^\s*// if $line =~ /\w/;
        if ($line =~ /^}/) {
            last;
        }
        $data .= $line;
    }
    $ret_val{'data'} = $data;

    return \%ret_val;
}

sub _build_CT {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val;
    $ret_val{'type'} = 'contig_tag';
    my $line = <$IN>;
    chomp $line;
	$line =~ s/^\s*// if $line =~ /\w/;     
    @ret_val{'contig_name', 'tag_type', 'program', 'start_pos', 'end_pos', 'date', 'no_trans'} = split(/ /, $line);
    while ($line = <$IN>) {
		$line =~ s/^\s*// if $line =~ /\w/;
        if ($line =~ /^}/) {
            last;
        }
        $ret_val{data} .= $line;
    }
    chomp $ret_val{data} if $ret_val{data};
    return \%ret_val;
}

sub _build_RT {
    my ($self, $IN, $token_ary_ref) = @_;
    my %ret_val;
    $ret_val{'type'} = 'read_tag';
    my $line = <$IN>;
    chomp $line;
	$line =~ s/^\s*// if $line =~ /\w/;;
    @ret_val{'read_name', 'tag_type', 'program', 'start_pos', 'end_pos', 'date'} = split(/ /, $line);
    
    while (my $nextline= <$IN>)
    {
        last if $nextline=~/^\s*}\s*\n?$/;
        $ret_val{data}.=$nextline;
    }
    chomp $ret_val{data} if $ret_val{data};
    return \%ret_val;
}

1;

=pod

=head1 NAME

AceReader - Ace file iterator

=head1 SYNOPSIS

    my $reader = new Genome::Model::Tools::Consed::AceReader(input => \*STDIN);
    while (my $obj = $reader->next_object()) {
        if ($obj->{'type'} eq 'contig') {
            ...
        }
        ...
    }

=head1 DESCRIPTION

Genome::Model::Tools::Consed::AceReader iterates over an ace file, returning one element at a time.

=head1 METHODS

=item new 

    my $reader = new Genome::Model::Tools::Consed::AceReader(\*STDIN);

=item next_object 

    my $obj_hashref = $reader->next_object();

    $obj_hashref->{'type'} eq 'contig'

    next_object returns the next object found in the ace file.  The return value is a
    hashref containing a 'type' key, and various other keys depending on the type:

    type eq 'assembly'
        contig_count
        read_count

    type eq 'contig'
        name
        base_count
        read_count
        base_seg_count
        u_or_c
        consensus
        base_qualities

    type eq 'read_position'
        read_name
        u_or_c
        position

    type eq 'base_segment'
        start_pos
        end_pos
        read_name

    type eq 'read'
        name
        padded_base_count
        info_count
        tag_count
        sequence
        qual_clip_start
        qual_clip_end
        align_clip_start
        align_clip_end
        description        - A hashref containing details about the trace
            CHROMAT_FILE
            PHD_FILE
            TIME

    type eq 'assembly_tag'
        tag_type
        program
        date
        data

    type eq 'contig_tag'
        contig_name
        tag_type
        program
        start_pos
        end_pos
        date
        no_trans
        data

    type eq 'read_tag'
        read_name
        tag_type
        program
        start_pos
        end_pos    
        date
        data
        

=cut

