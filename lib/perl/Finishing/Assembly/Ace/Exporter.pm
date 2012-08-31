package Finishing::Assembly::Ace::Exporter;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;
use Finfo::ClassUtils 'class';
use Finfo::Iterator;
use IO::File;
use IO::String;
use Tie::File;

my %file :name(file:r) :isa(file_rw);
my %width :name(width:o) :isa('int non_neg') :default(50);
my %overwrite :name(overwrite:o) :isa('boolean') :default(0);

my %fh :name(_fh:p) :isa('object IO::File');
my %ctg_names :name(_contig_names:p) :ds(hashref) :empty_ok(1) :default({});
my %ctg_ct :name(_contig_count:p) :isa('int') :default(0);
my %read_ct :name(_read_count:p) :isa('int') :default(0);
my %at_buffer :name(_at_buffer:p) :isa('object');
my %ct_buffer :name(_ct_buffer:p) :isa('object');

sub START
{
    my $self = shift;

    $self->_fh( IO::File->new('>'. $self->file) );
    $self->_fh->print("\n\n\n");
    $self->_at_buffer( IO::String->new() );
    $self->_ct_buffer( IO::String->new() );

    return 1;
}

sub close
{
    my $self = shift;

    # contig tags
    $self->_ct_buffer->seek(0, 0);
    $self->_fh->print( join('', $self->_ct_buffer->getlines) );
    $self->_ct_buffer->close;

    # assembly tags
    $self->_at_buffer->seek(0, 0);
    $self->_fh->print( join('', $self->_at_buffer->getlines) );
    $self->_at_buffer->close;

    $self->_fh->close;
    
    # assembly header
    tie(my @array, 'Tie::File', $self->file);
    $array[0] = sprintf('AS %d %d', $self->_contig_count, $self->_read_count);
    untie(@array);

    return 1;
}

sub export_assembly_tags
{
    my ($self, $tags) = @_;

    #print Dumper($tags);
    
    return unless $tags and @$tags;

    foreach my $tag ( @$tags )
    {
        $self->_at_buffer->print
        (
            "\nWA{\n%s %s %s%s%s\n}\n",,
            #"CT{\n%s %s %s %d %d %s%s%s\n}\n\n",
            $tag->type,
            $tag->source,
            $self->_convert_to_ace_tag_date( $tag->date ),
            ( $tag->no_trans ) ? ' ' . $tag->no_trans : '',
            $self->_serialize_tag_text($tag) || '',
        );
    }

    return 1;
}

sub export_contigs
{ 
    my ($self, %p) = @_;

    my $ci = delete $p{ci};
    Finfo::Validate->validate
    (
        attr => 'contig iterator to export',
        value => $ci,
        isa => 'object Finishing::Assembly::Iterator',
        msg => 'fatal',
    );
    
    while ( my $contig = $ci->next )
    {
        $self->export_contig
        (
            contig => $contig,
            auto_rename => $p{auto_rename},
        );
    }

    return 1;
}

sub export_contig
{ 
    my ($self, %p) = @_;

    my $contig = delete $p{contig};
    Finfo::Validate->validate
    (
        attr => 'contig to export',
        value => $contig,
        isa => 'object',# Finishing::Assembly::Contig Finishing::Assembly::Ace::Contig',
        msg => 'fatal',
    );

    $self->_increment_contig_count;

    # determine contig name
    my $name = ( $p{new_name} )
    ? delete $p{new_name}
    : ( $p{auto_rename} ) 
    ? sprintf('Contig%s', $self->_contig_count)
    : $contig->name;
   
    # keep track of contig names
    $self->fatal_msg("Contig with name ($name) already written to acefile") if exists $self->_contig_names->{$name};
    $self->_contig_names->{$name} = 1;

    my $ri = $contig->assembled_reads;
    my $read_count = $ri->count;
    $self->_add_to_read_count($read_count);

    my $base_segments = $contig->base_segments;

    $self->_fh->print
    (
        sprintf
        (
            "CO %s %d %d %d %s\n",
            $name,
            $contig->length,
            $ri->count,
            scalar(@$base_segments),
            ( $contig->complemented ? 'C' : 'U' ),
        )
    );
    $self->_write_base_string($contig);
    $self->_fh->print("\n");
    $self->_write_qualities($contig);
    $self->_write_read_positions($ri);
    $ri->reset;
    $self->_write_base_segments($base_segments);
    $self->_write_reads($ri);
    $self->_add_contig_tags_to_buffer($name, $contig->tags);
    
    return 1;
}

sub _increment_contig_count
{
    my $self = shift;

    return $self->_contig_count( $self->_contig_count + 1 );
}

sub _add_to_read_count
{
    my ($self, $count) = @_;

    return $self->_read_count( $self->_read_count + $count );
}

sub _write_base_string
{
    my ($self, $obj) = @_;
	
    my $base_string = $obj->base_string;
    my $width = $self->width;
	my $length = length($base_string);
    for (my $i = 0; $i < $length; $i += $width ) 
    {
        $self->_fh->print
        (
            substr
            (
                $base_string, 
                $i, 
                ( ( $i + $width - 1 ) < $length )
                ? $width 
                : $length - $i,
            ),
            "\n",
        );
    }

    return $self->_fh->print("\n");
}

sub _write_qualities 
{
	my ($self, $obj) = @_;

	$self->_fh->print("BQ\n");

    my $qualities = $obj->qualities;
	my $width = $self->width;
	my $length = scalar(@$qualities) - 1;
    for (my $i = 0; $i <= $length; $i += $width ) 
    {
        my $stop = ( $i + $width - 1 < $length )
        ? $i + $width - 1
        : $length;

        $self->_fh->print(' ', join(' ', @$qualities[$i..$stop]), "\n");
    }

    return $self->_fh->print("\n");
}

sub _write_read_positions
{
    my ($self, $ri) = @_;

    while ( my $read = $ri->next )
    {
        $self->_fh->print
        (
            sprintf("AF %s %s %d\n", $read->name, ( $read->complemented ? 'C' : 'U' ), $read->position) 
        );
    }

    return 1;
}

sub _write_base_segments
{
    my ($self, $base_segments) = @_;

    foreach my $bs ( @$base_segments )
    {
        $self->_fh->print
        ( 
           sprintf("BS %d %d %s\n", $bs->{start}, $bs->{stop}, $bs->{read_name})
        );
    }
    
    return $self->_fh->print("\n");
}

sub _write_reads
{
    my ($self, $ri) = @_;

    while ( my $read = $ri->next )
    {
        my $tag_count = 0;
        my $tags = $read->tags;
        $tag_count = ( $tags ) ? scalar @$tags : 0;
        $self->_fh->print
        (
            sprintf
            (
                "RD %s %d %d %d\n",
                $read->name, 
                $read->length,
                $read->info_count, 
                $tag_count,
            )
        );
        
        $self->_write_base_string($read);
        
        # QA
        $self->_fh->print
        (
            sprintf
            (
                "QA %d %d %d %d\n",
                $read->qual_clip_start,
                $read->qual_clip_stop, 
                $read->align_clip_start,
                $read->align_clip_stop,
            )
        );

        # DS
        $self->_fh->print("DS");
        foreach my $attr (qw/ chromat_file phd_file time /)
        #foreach my $attr (qw/ chromat_file phd_file chem dye time /)
        {
            my $val = $read->$attr;
            next unless $val;
            $self->_fh->print
            (
                sprintf
                (
                    ' %s: %s',
                    uc $attr,
                    $val,
                )
            );
        }

        $self->_fh->print("\n\n");

        # Tags
        $self->_write_read_tags( $read->tags );
    }

    return 1;
}

sub _add_contig_tags_to_buffer
{
    my ($self, $parent, $tags) = @_;

    #print Dumper($tags);
    
    return unless $tags and @$tags;

    foreach my $tag ( @$tags )
    {
        $self->_ct_buffer->print
        (
            sprintf
            (
                "CT{\n%s %s %s %d %d %s%s%s\n}\n\n",
                $parent,
                $tag->type,
                $tag->source,
                $tag->start,
                $tag->stop,
                $self->_convert_to_ace_tag_date( $tag->date ),
                (( $tag->no_trans ) ? ' NoTrans' : ''),
                ($self->_serialize_tag_text($tag) || ''),
            )
        );
    }

    return 1;
}

sub _convert_to_ace_tag_date
{
    my ($self, $date) = @_;

    $date =~ /(19|20)(\d\d)-(\d\d)-(\d\d) (\d\d):(\d\d):(\d\d)/;
    
    return "$2$3$4:$5$6$7";
}

sub _write_read_tags
{
    my ($self, $tags) = @_;

    return unless $tags and @$tags;
    
    foreach my $tag ( @$tags )
    {
        $self->_fh->print
        (
            sprintf
            (
                "RT{\n%s %s %s %d %d %s%s\n}\n\n",
                $tag->parent,
                $tag->type,
                $tag->source,        
                $tag->start,
                $tag->stop,
                $self->_convert_to_ace_tag_date( $tag->date ),
                $self->_serialize_tag_text($tag) || '',
            )
        );
    }

    return 1;
}
   
sub _serialize_tag_text
{
    my ($self, $tag) = @_;

    my $text;

    if ( $tag->text ) 
    {    
        if ( grep { $tag->type eq $_ } (qw/ oligo /) )
        {
            my $serialize_text_method = sprintf
            (
                '_serialize_%s_text', 
                lc($tag->type) . '_tag',
            );

            $text = $self->$serialize_text_method($tag);
        }
        else
        {
            $text .= $tag->text;
            chomp $text;
        }
    }

    if ( my $comment = $tag->comment )
    {
        chomp $comment;
        $text .= sprintf("\nCOMMENT{\n%s\nC}", $comment);
    }
 
    return ( $text ) ? "\n$text" : undef;
}

sub _serialize_oligo_tag_text
{
    my ($self, $tag) = @_;

    # oligo text:
    #
    # C_AB0496C15.2 tgtctggattcaagtcgtgattc 61 C
    # clone

    my $text = sprintf
    (
        "%s %s %s %s\n%s",
        $tag->oligo_name,
        $tag->oligo_seq,
        $tag->oligo_temp,
        $tag->orientation,
        ( $tag->oligo_templates ) ? join(' ', @{ $tag->oligo_templates }) : '',
    );
    
    return $text;
}

1;

=pod

=head1 Name

Finishing::Assembly::Ace::Exporter

=head1 Synopsis

Writes a contig or contigs from an iterator to a file in the ace format.

=head1 Usage

 use Finishing::Assembly::Ace::Exporter;
 use Finishing::Assembly::Factory;

 my $factory =  Finishing::Assembly::Factory->connect->('cmap_user');
 my $organism = $factory->get_organism("pan_troglodytes");
 my $assembly = $organism->get_assembly('2.1_051011');
 my $xporter = Finishing::Assembly::Ace::Exporter->new(file => 'new.ace');

I<write w/ contig iterator>

 $xporter->export_contig_iterator
 (
    contig_iterator => $assembly->contigs,
    auto_rename => 1, # opt, renames contigs to their number as they are written
 );

I<or one contig at a time>

 my $contigs = $assembly->contigs;
 while ( my $contig = $contigs->next )
 {
    $xporter->export_contig_iterator
    (
        contig_iterator => $assembly->contigs,
        # rename contig, optional, use only one:
        new_name => $contigs_new_name,
        auto_rename => 1, # see above
    );
 }

 $xporter->close; # writes the assembly header, contig tags

=head2 export_contigs

 $xporter->export_contigs
 (
    contigs => $contig_iterator,
    auto_rename => 1, # rename contigs, optional
 );

=over

=item I<Synopsis>   goes thru a contig iterator, writing each contig to the file in 'ace' format

=item I<Params>     contig iterator (Finishing::Assembly::Iterator)

=item I<Returns>    true on success

=back

=head2 export_contig

 Finishing::Assembly::Ace::Exporter->instance->export_contig
 (
    contig => $contig,
    # rename contig, optional, use one only:
    new_name => $name,
    auto_rename => 1,
 );

=over

=item I<Synopsis>   writes a contig to the file in 'ace' format

=item I<Params>     contigs (Finishing::Assembly::Contig)

=item I<Returns>    true on success

=back

=head1 See Also

=over

=item B<Finishing::Assembly::Factory>

=back

=head1 Disclaimer

Copyright (C) 2006-2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
