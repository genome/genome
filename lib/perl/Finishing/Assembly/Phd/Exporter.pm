package Finishing::Assembly::Phd::Exporter;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;
use Date::Format;
use IO::File;

my %file :name(file:r) :isa(file_rw);
my %read :name(read:r) :isa(object);

my %fh :name(_fh:p) :isa('object IO::File');

sub START
{
    my $self = shift;

    $self->_fh( IO::File->new('>'. $self->file) );

    return 1;
}

sub export { execute(@_); }
sub execute
{
    my $self = shift;
    
    $self->_fh->print
    (
        sprintf
        (
            "BEGIN_SEQUENCE %s\n\n%s\n%s\n%s%sEND_SEQUENCE\n",
            $self->read->name,
            $self->_get_comment_string,
            $self->_get_dna_string,
            $self->_get_wr_string,
            $self->_get_tag_string,
        )
    );
	
    $self->_fh->close;
    
    return 1;
}

sub _get_comment_string
{
    my $self = shift;

    my @fields = 
    (qw/
        chromat_file
        abi_thumbprint
        trace_processor_version
        base_caller_version
        phred_version
        call_method
        quality_levels
        time
        trace_array_min_index
        trace_array_max_index
        trim
        trace_peak_area_ratio
        chem
        dye
    /);

    my %comments = %{ $self->read->comments };
    my $string = "BEGIN_COMMENT\n\n";

    foreach my $field ( @fields )
    {
        my $value = delete $comments{$field};
        next unless defined $value;

        $string .= sprintf
        (
            "%s: %s\n", 
            uc($field), 
            $value,
        );
    }

    $string .= "\nEND_COMMENT\n";

    return $string;
}

sub _get_dna_string
{
    my $self = shift;

    my $read = $self->read;
    my @bases = split(//, $read->unpadded_base_string);
    my $qualities = $read->qualities;
    my $chromats = $read->unpadded_chromat_positions;

    my $string = "BEGIN_DNA\n";
    my $default_qual = 20; # TODO allow as parameter??
    my $trace_width = 0;
    for(my $i = 0; $i <= $#bases; $i++)
    {
        $string .= sprintf
        (
            "%s %d %d\n",
            $bases[$i], 
            ( defined $qualities->[$i] ) ? $qualities->[$i] : $default_qual,
            ( defined $chromats->[$i] ) ? $chromats->[$i] : ($trace_width += 10)
        );
    }

    return $string . "END_DNA\n\n";
}

sub _get_wr_string 
{
    my $self = shift;

    my $string;
    foreach my $wr ( @{ $self->read->wr } )
	{
        chomp $wr;
        $string .= "\nWR{\n$wr\n}\n";
	}

    return ( $string ) ? "$string\n" : '';
}

sub _get_tag_string  : PRIVATE
{
    my $self = shift;

    my $string;
	foreach my $tag ( @{ $self->read->tags } )
	{
        $string .= sprintf
        (
            "\nBEGIN_TAG\nTYPE: %s\nSOURCE: %s\nUNPADDED_READ_POS: %d %d\nDATE: %s\n%sEND_TAG\n",
            $tag->type,
            $tag->source,
            $tag->unpad_start,
            $tag->unpad_stop,
            $self->_convert_to_phd_tag_date( $tag->date ),
            ( $tag->text )
            ? "BEGIN_COMMENT\n" . $tag->text . "\nEND_COMMENT\n"
            : '',
        );
    }    

    return ( $string ) ? "$string\n" : '';
}

sub _convert_to_phd_tag_date : PRIVATE
{
    my ($self, $date) = @_;

    $date =~ s/^\d\d//g;
    
    return $date;
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
