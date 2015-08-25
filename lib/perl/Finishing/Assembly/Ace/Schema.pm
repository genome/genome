package Finishing::Assembly::Ace::Schema;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;
use Finfo::Iterator;
use Finishing::Assembly::Ace::AssembledRead;
use Finishing::Assembly::Ace::Assembly;
use Finishing::Assembly::Ace::Contig;
use Finishing::Assembly::Ace::Exporter;
use Finishing::Assembly::Ace::Tags;
use Finishing::Assembly::Cache;
use IO::File;

require File::Copy;

my %file :name(file:r) :isa(file_r);
my %a_cache :name(_assembly_cache:p) :isa(object);
my %c_cache :name(_contig_cache:p) :isa(object);
my %r_cache :name(_read_cache:p) :isa(object);

sub START
{
    my $self = shift;

    $self->_assembly_cache( Finishing::Assembly::Cache->new() );
    $self->_contig_cache( Finishing::Assembly::Cache->new() );
    $self->_read_cache( Finishing::Assembly::Cache->new() );

    $self->_build_cache;

    return 1;
}

sub connect
{
    my ($class, $file) = @_;

    return $class->new(file => $file);
}

sub disconnect
{
    my $self = shift;

    return 1;
}

sub _fh : PRIVATE
{
    my $self = shift;
    
    my $fh = IO::File->new('<' . $self->file )
        or $self->fatal_msg( sprintf('Can\'t open (%s): %s', $self->file, $!) );

    return $fh;
}

sub _line_from_fh_at
{
    my ($self, $pos) = @_;

    my $fh = $self->_fh;
    $fh->seek($pos, 0);
    my $line = $fh->getline;
    chomp $line;
    $fh->close;

    return $line;
}

sub _read_from_fh
{
    my ($self, $start, $stop) = @_;

    print Dumper([caller, \@_]) unless $stop;
    
    my $fh = $self->_fh;
    $fh->seek($start, 0);
    my $string;
    $fh->read($string, $stop - $start + 1);
    $fh->close;

    return $string;
}

sub tag_factory
{
    return Finishing::Assembly::Ace::TagFactory->instance;
}

sub _build_cache : PRIVATE
{
	my $self = shift;

    my $fh = $self->_fh;
    my $as_line = $fh->getline;
    Finfo::Validate->validate
    (
        attr => sprintf('AS line for acefile (%s)', $self->file),
        value => $as_line,
        #isa => [ 'regex', '^AS \d+ \d+' ], lift this for Newbler format
        isa => [ 'regex', '^AS\s+\d+\s\d+' ],
        msg => 'fatal',
    );
    
    my $assembly_cache = $self->_assembly_cache;
    $assembly_cache->set('assembly', 'tags', []);

    my $contig_cache = $self->_contig_cache;
    my $read_cache = $self->_read_cache;

    my ($contig_name, $read_name);
    my $position = tell($fh);
    
	while ( my $line = $fh->getline )
	{
        $line =~ s/^\s+//;
		my $first_two = substr($line, 0, 2);        
		my $first_three = substr($line, 0, 3);        
        
		if ( $first_three eq 'CO ' )
		{
            #CO Contig1 796 1 1 U
            #CO name length read_count base_segment_count orientation
            my @tokens = split(/\s+/, $line);
            $contig_name = $tokens[1];
            $contig_cache->set($contig_name, 'name', $contig_name);
            $contig_cache->set($contig_name, 'length', $tokens[2]);
            $contig_cache->set($contig_name, 'complemented', ( $tokens[5] eq 'C' ? 1 : 0));
            $contig_cache->set($contig_name, 'base_string_start', $fh->tell);
            $contig_cache->set($contig_name, 'read_names', {});
            $contig_cache->set($contig_name, 'tags', []);
		}
        elsif ( $first_two eq 'BQ' )
        #elsif ( $first_three =~ /^BQ/ )
        {
            $contig_cache->set($contig_name, 'base_string_stop', $position - 1);
            $contig_cache->set($contig_name, 'qualities_start', $fh->tell);
        }
        elsif ( $first_three eq 'AF ' )
        {
            unless ( $contig_cache->get($contig_name, 'qualities_stop') )
            {
                $contig_cache->set($contig_name, 'qualities_stop', $position - 1);
            }
            chomp $line;
            #AF L25990P6001H1.b1 U 1
            #AF read_name orientation position
            my @tokens = split(/\s+/, $line);
            $read_cache->set($tokens[1], 'contig_name', $contig_name);
            $read_cache->set($tokens[1], 'complemented', ( $tokens[2] eq 'C' ? 1 : 0));
            $read_cache->set($tokens[1], 'position', $tokens[3]);
            my $reads = $contig_cache->get($contig_name, 'read_names');
            $reads->{ $tokens[1] } = 1;
        }
        elsif ( $first_three eq 'BS ' )
        {
            unless ( $contig_cache->get($contig_name, 'bs_start') )
            {
                $contig_cache->set($contig_name, 'bs_start', $position);
            }
        }
        elsif ( $first_three eq 'RD ' )
        {
            unless ( $contig_cache->get($contig_name, 'bs_stop') )
            {
                $contig_cache->set($contig_name, 'bs_stop', $position - 1);
            }

            #RD L25990P6002F11.g1 1276 0 2
            #RD read_name length info_count tag_count
            my @tokens = split(/\s+/, $line);
            $read_name = $tokens[1];
            $read_cache->set($read_name, 'length', $tokens[2]);
            $read_cache->set($read_name, 'info_count', $tokens[3]);
            $read_cache->set($read_name, 'tag_count', $tokens[4]);
            $read_cache->set($read_name, 'base_string_start', $fh->tell);
        }
        elsif ( $first_three eq 'QA ' )
        {
            $read_cache->set($read_name, 'base_string_stop', $position - 1);
            $read_cache->set($read_name, 'qa_start', $position);
        }
        elsif ( $first_three eq 'DS ' )
        {
            $read_cache->set($read_name, 'ds_start', $position);
        }
        elsif ( $first_three eq 'RT{' )
        {
            $read_cache->push($read_name, 'tag_file_positions', $fh->tell);
        }
        elsif ( grep { $first_three eq $_ } ( 'WA ', 'WA{' ) )
        {
            $assembly_cache->push
            (
                'assembly', 
                'tags',
                $self->tag_factory->build_assembly_tag($fh),
            );
        }
        elsif ( $first_three eq 'CT{' )
        {
            my $tag = $self->tag_factory->build_contig_tag($fh);
            $contig_cache->push
            (
                $tag->parent, 
                'tags', 
                $tag,
            );
		}

        $position = tell($fh);
	}

    $fh->close;

    return 1;
}

sub _sort_contigs : PRIVATE
{
    $a =~ /(?:Contig)(\d+)(?:\.(\d+))*/;
    my $a_super = $1;
    my $a_reg = ( defined $2 ) ? $2 : -1;
    
    $b =~ /(?:Contig(\d+))(?:\.(\d+))*/;
    my $b_super = $1;
    my $b_reg = ( defined $2 ) ? $2 : -1;

    return $a_super <=> $b_super || $a_reg <=> $b_reg;
}

#- ASSEMBLY -#
sub get_assembly
{
    my $self = shift;

    return Finishing::Assembly::Ace::Assembly->new
    (
        _tags => sub{ $self->_assembly_tags(@_); },
        _scaffolds => sub{ $self->_scaffolds },
        _contig_count => sub{ $self->_contig_count },
        _contig_names => sub{ $self->_contig_names },
        _contigs => sub{ $self->_contigs },
        _get_contig => sub{ $self->_get_contig(@_) },
        _assembled_read_count => sub{ $self->_assembled_read_count },
        _assembled_reads => sub{ $self->_assembled_reads },
        _get_assembled_read => sub{ $self->_get_assembled_read(@_) },
        _destroy => sub{ $self->_flush_assembly; },
    );
}

sub _assembly_tags
{
    my ($self, $new_tags) = @_;

    return $self->_assembly_cache->set('assembly', 'tags', $new_tags) if $new_tags;

    return $self->_assembly_cache->get('assembly', 'tags');
}

sub _contig_names : PRIVATE
{
    my $self = shift;

    return sort { _sort_contigs() } $self->_contig_cache->ids;
}

sub _contig_count : PRIVATE
{
    my $self = shift;

    return scalar $self->_contig_cache->ids;
}

sub _assembled_read_count : PRIVATE
{
    my $self = shift;

    return scalar $self->_read_cache->ids;
}

sub _scaffolds : PRIVATE
{
    my $self = shift;

    #TODO figure out ace scaffolding
    
    return Finfo::Iterator->new
    (
        ids => [ ],
        cb => sub{ },
    );
}

sub _contigs : PRIVATE
{
    my $self = shift;

    return Finfo::Iterator->new
    (
        ids => [ $self->_contig_names ],
        cb => sub{ $self->_get_contig(@_); },
    );
}

sub _flush_assembly
{
    my $self = shift;

    return 1;
}

#_ CONTIG -#
sub _get_contig : PRIVATE
{
    my ($self, $name) = @_;

    $self->fatal_msg("Can't find contig ($name)") unless $self->_contig_cache->id_exists($name);
    
    my $scope = $self->_contig_cache->get($name,'in_scope');
    $scope = 0 if !defined $scope;
    $scope++;
    $self->_contig_cache->set($name, 'in_scope', $scope);
    
    return Finishing::Assembly::Ace::Contig->new
    (
        _name => sub{ return $self->_contig_name($name, @_); },
        _complemented => sub{ $self->_contig_complemented($name, @_); },
        _base_string => sub{ return $self->_contig_base_string($name, @_); },
        _length => sub{ return $self->_contig_length($name); },
        _qualities => sub{ return $self->_contig_qualities($name, @_); },
        _assembled_reads => sub{ return $self->_contig_assembled_reads($name); },
        _get_assembled_read => sub{ return $self->_get_assembled_read(@_); },
        _read_count => sub{ return $self->_contig_read_count($name) },
        _base_segments => sub{ return $self->_contig_base_segments($name, @_); },
        _tags => sub{ return $self->_contig_tags($name, @_) },
        _destroy => sub{ $self->flush_contig($name); },
    );
}

sub _contig_name : PRIVATE
{
    my ($self, $name, $new_name) = @_;

    return $name unless $new_name;

    my $contig_cache = $self->_contig_cache;
    $self->fatal_msg("Contig name ($name) already exists") if $contig_cache->id_exists($name);

    $contig_cache->move($name, $new_name);

    if ( my $tags = $contig_cache->get($new_name, 'tags') )
    {
        foreach my $tag ( @$tags )
        {
            $tag->parent($new_name);
        }
    }

    return $new_name;
}

sub _contig_length : PRIVATE
{
    my ($self, $name) = @_;

    my $contig_cache = $self->_contig_cache;
    my $length = $contig_cache->get($name, 'length');

    unless ( $length )
    {
        my $base_string = $self->_contig_base_string($name);
        $length = CORE::length($base_string);
        $contig_cache->set($name, 'length', $length);
    }

    return $length;
}

sub _contig_complemented : PRIVATE
{
    my ($self, $name, $comp) = @_;

    if ( defined $comp )
    {
        return $self->_contig_cache->set($name, 'complemented', $comp);
    }

    return $self->_contig_cache->get($name, 'complemented');
}

sub _contig_read_count : PRIVATE
{
    my ($self, $name) = @_;

    return scalar( @{ $self->_contig_cache->get($name, 'read_names') } );
}

sub _contig_base_string : PRIVATE
{
    my ($self, $name, $new_base_string) = @_;

    my $contig_cache = $self->_contig_cache;

    if ( defined $new_base_string )
    {
        $contig_cache->set($name, 'base_string', $new_base_string);
        $contig_cache->set($name, 'length', 0); # set as zero, if requested will calculate
        return $new_base_string;
    }

    my $base_string = $contig_cache->get($name, 'base_string');
    unless ( $base_string )
    {
        $base_string = $self->_read_from_fh
        (
            $contig_cache->get($name, 'base_string_start'),
            $contig_cache->get($name, 'base_string_stop') 
        );
        $base_string =~ s/\s+//g;
        $contig_cache->set($name, 'base_string', $base_string);
    }

    return $base_string;
}

sub _contig_qualities : PRIVATE
{
    my ($self, $name, $new_qualities) = @_;

    my $contig_cache = $self->_contig_cache;

    if ( defined $new_qualities )
    {
        $contig_cache->set($name, 'qualities', $new_qualities);
        return $new_qualities;
    }

    my $qualities = $contig_cache->get($name, 'qualities');
    unless ( $qualities )
    {
        my $qual_string = $self->_read_from_fh
        (
            $contig_cache->get($name, 'qualities_start'),
            $contig_cache->get($name, 'qualities_stop') 
        );
        $qual_string =~ s/^\s+//;
        $qual_string =~ s/\n//g; 
        $qual_string =~ s/\s+$//;
        $qualities = [ split(/\s+/, $qual_string) ];
        $contig_cache->set($name, 'qualities', $qualities);
    }

    return $qualities;
}

sub _contig_assembled_reads : PRIVATE
{
    my ($self, $name, $conds, $params) = @_;

    my @names = keys %{ $self->_contig_cache->get($name, 'read_names') };
    if ( $conds ) # sort by
    {
        # this is a search...only sorting by position
        my $read_cache = $self->_read_cache;
        @names = sort 
        {
            $read_cache->get($a, 'position') <=> $read_cache->get($b, 'position')
        } @names;
    }

    return Finfo::Iterator->new
    (
        ids => \@names,
        cb => sub{ $self->_get_assembled_read(@_); },
    );
}

sub _contig_base_segments : PRIVATE
{
    my ($self, $name, $new_base_segments) = @_;

    my $contig_cache = $self->_contig_cache;

    if ( $new_base_segments )
    {
        return $self->_contig_cache->set($name, 'base_segments', $new_base_segments);
    }

    my $base_segments = $contig_cache->get($name, 'base_segments');
    return $base_segments if $base_segments;

    my @base_segments;
    my $base_segment_string = $self->_read_from_fh
    (
        $contig_cache->get($name, 'bs_start'),
        $contig_cache->get($name, 'bs_stop'),
    );
    foreach my $line ( split(/\n/, $base_segment_string) )
    {
        my %bs;
        @bs{qw/ start stop read_name /} = (split(/ /, $line))[1..3];
        push @base_segments, \%bs;
    }

    return \@base_segments;
}

sub _contig_tags : PRIVATE
{
    my ($self, $contig_name, $tags) = @_;

    if ( $tags )
    {
        return $self->_contig_cache->set($contig_name, 'tags', $tags);
    }

    return $self->_contig_cache->get($contig_name, 'tags');
}

sub flush_contig {
    my ($self, $name) = @_;

    $self->fatal_msg("Can't find contig ($name)") unless $self->_contig_cache->id_exists($name);

    return $self->_flush_contig($name);
}

sub _flush_contig : PRIVATE {
    my ($self, $name) = @_;

    #$self->warn_msg("Flushing contig $name");

    my $scope = $self->_contig_cache->get($name,'in_scope');
    $scope--;
    $self->_contig_cache->set($name, 'in_scope', $scope);
    
    return 1 unless $scope == 0;
    foreach my $attr (qw/ base_string qualities base_segments /) {
        $self->_contig_cache->undef_attr($name, $attr);
    }

    # Flush reads unless they are in scope
    #for my $read_name ( keys %{ $self->_contig_cache->get($name, 'read_names') } ) {
    #    next if $self->_read_cache->get($name, 'in_scope');
    #    $self->_flush_assembled_read($read_name);
    #}
    
    return 1;
}

#- ASSEMBLED READ -#
sub _assembled_reads : PRIVATE
{
    my $self = shift;

    return Finfo::Iterator->new
    (
        ids => [ $self->_read_cache->ids ],
        cb => sub{ $self->_get_assembled_read(@_); },
    );
}

sub _get_assembled_read : PRIVATE
{
    my ($self, $name) = @_;

    $self->fatal_msg("Can't find read ($name)") unless $self->_read_cache->id_exists($name);

    my $scope = $self->_read_cache->get($name, 'in_scope');
    $scope = 0 if !defined $scope;
    $scope++;
    $self->_read_cache->set($name, 'in_scope', $scope);
    
    return Finishing::Assembly::Ace::AssembledRead->new
    (
        _name => sub{ return $self->_assembled_read_name($name, @_); },
        _rename => sub{ return $self->_rename_assembled_read($name, @_); },
        _tags => sub{ $self->_assembled_read_tags($name, @_); },
        _position => sub{ return $self->_assembled_read_position($name, @_); },
        _complemented => sub{ return $self->_assembled_read_complemented($name, @_); },
        _info_count => sub{ return $self->_assembled_read_info_count($name); },
        # assembled read's attributes
        _base_string => sub{ return $self->_assembled_read_base_string($name, @_); },
        _length => sub{ return $self->_assembled_read_length($name); },
        # qa
        _qual_clip_start => sub{ return $self->_assembled_read_qa_attribute($name, 'qual_clip_start', @_); },
        _qual_clip_stop => sub{ return $self->_assembled_read_qa_attribute($name, 'qual_clip_stop', @_); },
        _align_clip_start => sub{ return $self->_assembled_read_qa_attribute($name, 'align_clip_start', @_); },
        _align_clip_stop => sub{ return $self->_assembled_read_qa_attribute($name, 'align_clip_stop', @_); },
        # ds
        _time => sub{ return $self->_assembled_read_ds_attribute($name, 'time', @_); },
        _chromat_file => sub{ return $self->_assembled_read_ds_attribute($name, 'chromat_file', @_); },
        _phd_file => sub{ return $self->_assembled_read_ds_attribute($name, 'phd_file', @_); },
        _chem => sub{ return $self->_assembled_read_ds_attribute($name, 'chem', @_); },
        _dye => sub{ return $self->_assembled_read_ds_attribute($name, 'dye', @_); },
        # other
        _destroy => sub{ 
                    if(!defined $self) 
                    {
                        #$DB::single = 1;
                        return 1;
                        #print "$name not defined\n";                        
                    }
                    else
                    {
                        #print "$name\n";
                    }
                    return $self->_flush_assembled_read($name); },
    );
}

sub _assembled_read_name : PRIVATE
{
    my ($self, $name, $new_name) = @_;

    if ( $new_name )
    {
        $self->fatal_msg("Use method 'rename' to set a read's name");
    }

    return $name;
}

sub _rename_assembled_read
{
    my ($self, $name, %p) = @_;

    my $new_name = delete $p{new_name};
    $self->fatal_msg("Need new name to rename assembled read") unless $new_name;

    $self->_read_cache->move($name, $new_name);
    my $contig_name = $self->_read_cache->get($new_name, 'contig_name');
    my $read_names = $self->_contig_cache->get($contig_name, 'read_names');
    delete $read_names->{$name};
    $read_names->{$new_name} = 1;

    my $new_read = $self->_get_assembled_read($new_name);
    my $tags = $new_read->tags;
    if ( @$tags )
    {
        my @updated_tags;
        foreach my $tag ( @$tags )
        {
            $tag->parent($new_name);
            push @updated_tags, $tag;
        }
        $new_read->tags(\@updated_tags);
    }
    
    if ( $p{update_phd} )
    {
        my (undef, $iteration) = split(/\.phd\./, $new_read->phd_file);
        $new_read->phd_file( sprintf('%s.phd.%d', $new_name, $iteration) )
    }
    
    $new_read->chromat_file($new_name) if $p{update_chromat};

    return $new_read;
}

sub _assembled_read_length : PRIVATE
{
    my ($self, $name) = @_;

    my $read_cache = $self->_read_cache;
    my $length = $read_cache->get($name, 'length');

    unless ( $length )
    {
        my $base_string = $self->_assembled_read_base_string($name);
        $length = CORE::length($base_string);
        $read_cache->set($name, 'length', $length);
    }

    return $length;
}

sub _assembled_read_base_string : PRIVATE
{
    my ($self, $name, $new_base_string) = @_;

    my $read_cache = $self->_read_cache;

    if ( defined $new_base_string )
    {
        $read_cache->set($name, 'base_string', $new_base_string);
        $read_cache->set($name, 'length', 0); # set as zero, if requested will calculate
        return $new_base_string;
    }

    my $base_string = $read_cache->get($name, 'base_string');
    unless ( $base_string )
    {

        $base_string = $self->_read_from_fh
        (
            $read_cache->get($name, 'base_string_start'),
            $read_cache->get($name, 'base_string_stop') 
        );
        $base_string =~ s/\s+//g;
        $read_cache->set($name, 'base_string', $base_string);
    }

    return $base_string;
}

sub _assembled_read_position : PRIVATE
{
    my ($self, $name, $new_position) = @_;

    if ( defined $new_position )
    {
        return $self->_read_cache->set($name, 'position', $new_position);
    }
    
    return $self->_read_cache->get($name, 'position');
}

sub _assembled_read_complemented : PRIVATE
{
    my ($self, $name, $new_complemented) = @_;

    if ( defined $new_complemented )
    {
        return $self->_read_cache->set($name, 'complemented', $new_complemented);
    }

    return $self->_read_cache->get($name, 'complemented');
}

sub _assembled_read_tag_count : PRIVATE
{
    my ($self, $name) = @_;

    return $self->_read_cache->get($name, 'tag_count');
}

sub _assembled_read_info_count : PRIVATE
{
    my ($self, $name) = @_;

    return $self->_read_cache->get($name, 'info_count');
}

sub _assembled_read_qa_attribute : PRIVATE
{
    my ($self, $name, $attr, $new_value) = @_;

    my $read_cache = $self->_read_cache;
    my $qa = $read_cache->get($name, 'qa');
    unless ( $qa )
    {
        my $line = $self->_line_from_fh_at( $read_cache->get($name, 'qa_start') );
        my %qa;
        @qa{qw/ qual_clip_start qual_clip_stop align_clip_start align_clip_stop /} = (split(/ /,$line))[1..4];
        $qa = \%qa;
        $read_cache->set($name, 'qa', $qa);
    }

    if ( defined $new_value )
    {
        $qa->{$attr} = $new_value;
        return $new_value;
    }

    return $qa->{$attr};
}

sub _assembled_read_ds_attribute : PRIVATE
{
    my ($self, $name, $attr, $new_value) = @_;

    my $read_cache = $self->_read_cache;
    my $ds = $read_cache->get($name, 'ds');
    unless ( $ds )
    {
        my $line = $self->_line_from_fh_at( $read_cache->get($name, 'ds_start') );
        $line =~ s/ (\w+): /|$1|/g; #delimit the key-value pairs
        my @tokens = split(/\|/, $line); 
        shift @tokens; # drop the DS tag
        chomp @tokens;
        my %desc = @tokens;
        while ( my ($key, $val) = each %desc )
        {
            $ds->{ lc $key } = $val;
        }
        $read_cache->set($name, 'ds', $ds);
    }

    if ( defined $new_value )
    {
        $ds->{$attr} = $new_value;
        return $new_value;
    }

    return $ds->{$attr};
}

sub _assembled_read_tags : PRIVATE
{
    my ($self, $name, $tags) = @_;

    my $read_cache = $self->_read_cache;

    if ( $tags )
    {
        Finfo::Validate->validate
        (
            attr => "read's ($name) tags",
            value => $tags,
            ds => 'aryref',
            isa => 'object',
            empty_ok => 1,
            msg => 'fatal',
        );
        
        $read_cache->set($name, 'tag_count', scalar @$tags); 
        return $read_cache->set($name, 'tags', $tags); 
    }

    $tags = $read_cache->get($name, 'tags');
    return $tags if $tags;

    my @tags;
    if ( my $tag_file_positions = $read_cache->get($name, 'tag_file_positions') )
    {
        my $fh = $self->_fh;
        foreach my $position ( @$tag_file_positions )
        {
            $fh->seek($position, 0);
            push @tags, $self->tag_factory->build_read_tag($fh);
        }
        $fh->close;
    }

    return \@tags;
}

sub _assembled_read_infos : PRIVATE
{
    my ($self, $name, $infos) = @_;

    my $read_cache = $self->_read_cache;

    if ( $infos )
    {
        Finfo::Validate->validate
        (
            attr => "read's ($name) tags",
            value => $infos,
            ds => 'aryref',
            isa => 'object',
            empty_ok => 1,
            msg => 'fatal',
        );
        
        $read_cache->set($name, 'info_count', scalar @$infos); 
        return $read_cache->set($name, 'infos', $infos); 
    }

    $infos = $read_cache->get($name, 'infos');
    return $infos if $infos;

    my @infos;
    if ( my $info_file_positions = $read_cache->get($name, 'info_file_positions') )
    {
        my $fh = $self->_fh;
        foreach my $position ( @$info_file_positions )
        {
            $fh->seek($position, 0);
            push @infos, $self->tag_factory->build_read_info($fh);
        }
        $fh->close;
    }

    return \@infos;
}

sub flush_assembled_read {
    my ($self, $name) = @_;

    $self->fatal_msg("No read name to flush") unless $name;

    return $self->_flush_assembled_read($name);
}

sub _flush_assembled_read {
    my ($self, $name) = @_;

    #$self->warn_msg("Flushing read ($name)");
    my $scope = $self->_read_cache->get($name,'in_scope');
    $scope--;
    $self->_read_cache->set($name, 'in_scope', $scope);
    return 1 unless ($scope == 0);
    
    # don't flush read if contig still in scope
    my $contig_name = $self->_read_cache->get($name, 'contig_name');
    return 1 if $self->_contig_cache->get($contig_name, 'in_scope');
    #print "flushing assembled read.\n";
    for my $attr (qw/ base_string length base_string position complemented
        qa_attribute ds_attribute tags infos /) {
        
        $self->_read_cache->undef_attr($name, $attr);
    }

    return 1;
}

1;

=pod

=head1 Name

Finishing::Assembly::Ace::Schema

=head1 Synopsis

=head1 Usage

=head1 Methods

=head1 See Also

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$

