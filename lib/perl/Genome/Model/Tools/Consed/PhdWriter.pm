package Genome::Model::Tools::Consed::PhdWriter;

use strict;
use warnings;

use Data::Dumper;
use Date::Format;

class Genome::Model::Tools::Consed::PhdWriter {
};

sub write {
    my ($self, $file, $read) = @_;
    
    if ( not $file ) {
        $self->error_message('No file to write read as phd');
        return;
    }

    if ( not $read ) {
        $self->error_message('No read to write as phd');
        return;
    }

    my $fh = eval{ Genome::Sys->open_file_for_writing($file); };
    if ( not $fh ) {
        $self->error_message("Failed to open file ($file): $@");
        return;
    }

    $fh->print(
        sprintf(
            "BEGIN_SEQUENCE %s\n\n%s\n%s\n%sEND_SEQUENCE",
            $read->{name},
            $self->_get_comment_string($read),
            $self->_get_dna_string($read),
            $self->_get_tag_string($read),
        )
    );
    if ( my $wr_string = $self->_get_wr_string($read) ) {
        $fh->print("\n$wr_string");
    }
	
    $fh->close;
    
    return 1;
}

sub _get_comment_string {
    my ($self, $read) = @_;

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

    my %comments = %{ $read->{comments} };
    my $string = "BEGIN_COMMENT\n\n";

    foreach my $field ( @fields ) {
        my $value = delete $comments{$field};
        next unless defined $value;

        $string .= sprintf(
            "%s: %s\n", 
            uc($field), 
            $value,
        );
    }

    $string .= "\nEND_COMMENT\n";

    return $string;
}

sub _get_dna_string {
    my ($self, $read) = @_;

    my @bases = split(//, $read->{base_string});
    my $qualities = $read->{qualities};
    my $chromats = $read->{chromat_positions};

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

    return $string . "END_DNA\n";
}

sub _get_wr_string {
    my ($self, $read) = @_;

    return if not $read->{wr};

    my $string;
    foreach my $wr ( @{ $read->{wr} } ) {
        chomp $wr;
        $string .= "WR{\n$wr\n}\n\n";
	}

    return $string;
}

sub _get_tag_string {
    my ($self, $read) = @_;

    return '' if not $read->{tags};

    my $string;
	foreach my $tag ( @{ $read->{tags} } )
	{
        $string .= sprintf
        (
            "\nBEGIN_TAG\nTYPE: %s\nSOURCE: %s\nUNPADDED_READ_POS: %d %d\nDATE: %s\n%sEND_TAG\n",
            $tag->{tag_type},
            $tag->{program},
            $tag->{start_pos},
            $tag->{stop_pos},
            $tag->{date},
            ( $tag->text )
            ? "BEGIN_COMMENT\n" . $tag->text . "\nEND_COMMENT\n"
            : '',
        );
    }    

    return ( $string ) ? "$string\n" : '';
}

sub _convert_to_phd_tag_date {
    my ($self, $date) = @_;

    $date =~ s/^\d\d//g;
    
    return $date;
}

1;

=pod

=head1 NAME

PhdWriter - Phd file writer

=head1 SYNOPSIS

    my $writer = new Genome::Model::Tools::Consed::Phd::Writer();
    $writer->write(\*STDOUT,$phd);

=head1 DESCRIPTION

Genome::Model::Tools::Consed::Phd::Writer takes a handle to a phd object and writes it to the given file handle

=head1 METHODS

=cut


=pod

=item new 

    my $writer = new Genome::Model::Tools::Consed::Phd::Writer;

=cut

=pod

=item Genome::Model::Tools::Consed::Phd::Reader::write 

    $writer->read(\*STDOUT,$phd);

=cut

