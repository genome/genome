package Genome::Utility::DiffStream;

#:eclark 11/17/2009 Code review.

# One of many reader/writer classes in Genome/Utility.  Need a consistant naming convention and interface for all of them.

use strict;
use warnings;
use Data::Dumper;

#attributes

sub new{
    my $class = shift;
    my $io = shift;
    die "Need IO::Handle" unless $io->isa('IO::Handle');
    my $self = bless({_io => $io }, $class);
    $self->{next_diff} = $self->make_diff($self->{_io}->getline);
    $self->{current_diff_header} = '';
    return $self;
}

sub next_diff{
    my $self = shift;
    my $diff = $self->{next_diff};
    $self->{next_diff} = $self->make_diff($self->{_io}->getline);
    return unless $diff;
    $self->{current_diff_header} = $diff->{header};
    return $diff;
}

sub make_diff{
    my ($self, $line) = @_;
    return unless $line;
    my %diff;

    my ($subject, $pos, $type, $deletion_length, $indel, $pre_diff_sequence, $post_diff_sequence, $read_name) = split(/\s+/, $line);

    $diff{header} = $subject;
    $diff{delete} = uc $indel if $type =~ /^d(el(etion)?)?$/i;
    $diff{insert} = uc $indel if $type =~ /^i(ns(ertion)?)?$/i;
    ($diff{delete}, $diff{insert}) = split(/\//, $indel) if $type =~ /^s(ub(stitution)?)?$/i;
    
    $diff{delete} ||= '';
    $diff{insert} ||= '';
    

    $diff{position} = $pos;
    $diff{type} = $type;

    $diff{deletion_length} = $deletion_length;
    $diff{read_name} = $read_name;
    
    if ( $diff{delete} ){ # in ApplyDiffToFasta deletes start AFTER index, like inserts, also adjuct right flank;
        $diff{position}--;
    }
    
    if ($pre_diff_sequence or $post_diff_sequence){
        $diff{pre_diff_sequence} = uc $pre_diff_sequence unless $pre_diff_sequence eq '-';
        $diff{pre_diff_sequence} ||= '';
        $diff{post_diff_sequence} = uc $post_diff_sequence unless $post_diff_sequence eq '-';
        $diff{post_diff_sequence} ||= '';
    }

    $self->{current_diff_header} = $diff{header};
    return \%diff;
}

sub next_diff_position{
    my $self = shift;
    return unless $self->{next_diff};
    my $pos = $self->{next_diff}->{position};
    my $pre_diff_seq = $self->{next_diff}->{pre_diff_sequence};
    my $post_diff_seq = $self->{next_diff}->{post_diff_sequence};
    my $pre_diff_seq_length = length $pre_diff_seq if $pre_diff_seq;
    my $post_diff_seq_length = length $post_diff_seq if $post_diff_seq;
    return ($pos, $pre_diff_seq_length, $post_diff_seq_length) if $pos;
    return undef;
}

1;

=pod

=head1 DiffStream
diff file input stream used in genome-model-tools apply-diff-to-fasta

=head2 Synopsis

This streams through a diff file, and parses and returns diff objects used in the Tools command ApplyDiffToFasta 

my $ds = Genome::Utility::DiffStream->new( <file_name> )

=head2 Diff File Format
no header
<fasta_header_identifier> <position> <indel_type> <deletion_length> <indel> [<pre_diff_seq> <post_diff_seq>]

The fasta header identifier should match the first \S+ of chars following the '>' on the fasta header line in the reference sequence

valid indel types (case insensitive)
d[el[etion]]
i[ns[ertion]]
s[ub[stitution]]

S SUB del insertion Ins are all valid for this field

The deletion length is always the length of the sequence to be deleted from the reference in deletions and substitutions.  This field should always be 0 for insertions

The post and pre diff sequence are not necessary for apply-diff-to-fasta to run, but are recommended, as they ensure your diff is being applied in the correct position

=head2 Subs

=head3 next_diff
reads and returns the next diff in the file.  Advances one line in the diff file.

=head3 next_diff_position
returns the position of the next diff, but does not advance in the file.

=cut
