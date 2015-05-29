package Finishing::Assembly::Assembly;

use strict;
use warnings;

use base 'Finishing::Assembly::Item';

use Bio::Seq::Quality;

#- SCAFFOLDS -#
sub scaffold_count
{
    return shift->scaffolds->count;
}

#- CONTIGS -#
sub contigs
{
    my ($self, $conds, $params) = @_;

    my $method = $self->proxy->get_method('contigs');

    my $iterator = $method->();

    return $iterator unless $conds or $params;

    return $iterator->search($conds, $params);
}

sub get_bioseqs_from_all_contigs
{
    my $self = shift;

    my $positions = { map({ name => $_ }, $self->get_contig_names) };

    return $self->get_bioseqs_for_positions($positions);
}

sub get_bioseqs_from_tags
{
    my ($self, $tags) = @_;

    Finfo::Validate->validate
    (
        attr => 'tags to get bioseqs',
        value => $tags,
        ds => 'aryref',
        isa => 'object',
        msg => 'fatal',
    );
    
    my $positions = 
    {
        map(
        {
            name => $_->parent,
            start => $_->start,
            stop => $_->stop,
            type => 'padded',
        }, @$tags)
    };

    return $self->get_bioseqs_for_positions($positions);
}

sub get_bioseqs_from_contigs_for_positions
{
    my ($self, $positions) = @_;

    Finfo::Validate->validate
    (
        attr => 'contig positions',
        value => $positions,
        ds => 'hashref',
        msg => 'fatal',
    );

    my @bioseqs;
    foreach my $p ( @$positions )
    {
        my $contig = $self->get_contig( $p->{name} );
        push @bioseqs, $contig->to_bioseq(%$p);
    }

    return \@bioseqs;
}

sub contig_count{
    return shift->contigs->count;
}

#- READS -#
sub reads
{
    my $self = shift;

    return $self->assembled_reads;
}

sub read_count{
    return shift->assembled_reads->count;
}

sub assembled_read_count{
    return shift->assembled_reads->count;
}

#- TAGS -#
sub consensus_tags
{
    my ($self, %p) = @_;

    my $contig_conditions = delete $p{contig_conditions} || [];
    Finfo::Validate->validate
    (
        attr => 'contig search conditions',
        value => $contig_conditions,
        ds => 'aryref',
        empty_ok => 1,
        msg => 'fatal',
    );
    
    my $tag_conditions = delete $p{tag_conditions} || [];
    Finfo::Validate->validate
    (
        attr => 'tag search conditions',
        value => $tag_conditions,
        ds => 'aryref',
        empty_ok => 1,
        msg => 'fatal',
    );

    my $contig_iterator = $self->get_contig_iterator(@$contig_conditions);

    my @tags;
    while ( my $contig = $contig_iterator->next )
    {
        push @tags, $contig->tags(@$tag_conditions);
    }

    return @tags;
}

sub add_consensus_tags
{
    my ($self, $tags) = @_;

    Finfo::Validate->validate
    (
        attr => 'tags to add to assembly',
        value => $tags,
        ds => 'aryref',
        isa => 'object',
    );
    
    my %contigs_tags;
    foreach my $tag ( @$tags )
    {
        my $name = $tag->parent;
        my $contig = $self->get_contig($name);
        $self->fatal_msg("Can't find contig ($name)") unless $contig;
        $contig->add_tags($tag);
    }
    
    return 1;
}

sub remove_all_consensus_tags
{
    my $self = shift;

    my $ci = $self->get_contig_iterator(@_);
    while ( my $contig = $ci->next )
    {
        $contig->remove_all_tags;
    }
    
    return 1;
}

sub remove_duplicate_consensus_tags
{
    my $self = shift;

    my $ci = $self->get_contig_iterator(@_);
    while ( my $contig = $ci->next )
    {
        $contig->remove_duplicate_tags;
    }
    
    return 1;
}

sub remove_tags
{
    my ($self, %p) = @_;

    my $contig_conditions = delete $p{contig_conditions} || [];
    Finfo::Validate->validate
    (
        attr => 'contig search conditions',
        value => $contig_conditions,
        ds => 'aryref',
        empty_ok => 1,
        msg => 'fatal',
    );
    
    my $tag_conditions = delete $p{tag_conditions} || [];
    Finfo::Validate->validate
    (
        attr => 'tag search conditions',
        value => $tag_conditions,
        ds => 'aryref',
        empty_ok => 1,
        msg => 'fatal',
    );

    my $contig_iterator = $self->get_contig_iterator(@$contig_conditions);
    while ( my $contig = $contig_iterator->next )
    {
        $contig->remove_tags(@$tag_conditions);
    }

    return 1;
}

#- FUNCTIONS -#
#- Syncing PHD/READ TIMES -#
sub set_read_times_to_phd_times
{
    # TODO 
    my ($self, $phd_obj) = @_;
    
    Finfo::Validate->validate
    (
        attr => 'phd object',
        value => $phd_obj,
        isa => 'object',
        msg => 'fatal',
    );

    my $ci = $self->get_contig_iterator;
    while ( my $contig = $ci->next )
    {
        my $ri = $contig->get_assembled_read_iterator;
        while ( my $read = $ri->next )
        {
            my $phd = $phd_obj->get_phd($read->name . '.phd.1');
            next unless $phd;

            $read->time( $phd->time ); 
        }
    }

    return 1;
}

#- Merges/Joins -#
sub merge_scaffolds
{
    my ($self, %p) = @_;
    
}

sub join_contigs_in_scaffold
{
    my ($self, $scaffold) = @_;
}

sub join_contigs
{
    my ($self, %p) = @_;

    my $phd = delete $p{phd} #GSC::IO::Assembly::PhdDB->new;
        or $self->fatal_msg("Need phd to join contigs");
    my $file = delete $p{file};
    my $fh;
    if ( $file )
    {
        $fh = IO::File->new(">mergestats$ARGV[2]")
            or $self->fatal;
    }
    my $left_contig = delete $p{left_contig}
        or $self->fatal_msg("Need left contig to join");
    my $right_contig = delete $p{right_contig}
        or $self->fatal_msg("Need left contig to join");
    my $ct = GSC::IO::Assembly::ContigTools->new();

    my %merge_p = 
    ( 
        cutoffs => 
        { 
            hqlength => 1000,
            hq_percent_identity => 90,
        },
    );
    $merge_p{statsfh} = $fh if $fh;
    
    my $new_contig = $ct->merge($left_contig, $right_contig, $phd, %merge_p);

    $fh->print( sprintf("%s %s\n\n\n", $left_contig->name, $right_contig->name) );
    $fh->print("Finished detecting merges\n");
    $fh->close;

    return $new_contig;
}

###################################
#TODO below here, not updated for object model
sub join_scaffolds{

    my ($self, %p) = @_;
    my @missing_params;
    my $left_scaffold = delete $p{left_scaffold} or push @missing_params, 'left_scaffold';
    my $right_scaffold = delete $p{right_scaffold} or push @missing_params, 'right_scaffold';
    my $scaffold_num = delete $p{scaffold_num} or push @missing_params, 'scaffold_num';
    $self->fatal_msg('Missing required params for join_scaffolds:'.join(' ',@missing_params)) if @missing_params;

    my $left_contig = $left_scaffold->get_contig_iterator->last;
    my $right_contig = $right_scaffold->get_contig_iterator->first;
    my $right_neighbor = $left_contig->right_contig;
    my $left_neighbor = $right_contig->left_contig;
    $self->fatal_msg("Contigs joined to each other are not neighbors according to agp information!") unless $right_neighbor->id == $right_contig->id and $left_neighbor->id == $left_contig->id;

    $self->fatal_msg("Scaffolds to be merged have different orientations!") unless $right_scaffold->orientation == $left_scaffold->orientation;


    my $next_number = $left_contig->contig_num + 1;
    $self->schema->txn_do(
        sub{
            $self->schema->txn_do(
                sub{
                    my @contigs = $right_scaffold->contigs;
                    @contigs = sort { $a->contig_num <=> $b->contig_num } @contigs;
                    foreach my $contig ($right_scaffold->contigs){
                        $contig->set_contig_num($next_number);
                        $contig->set_scaffold($left_scaffold);
                        $next_number++;
                    }
                }
            );

            $self->schema->txn_do(
                sub{
                    $self->fatal_msg("Don't want to delete a scaffold with contigs in it!") unless $right_scaffold->contig_count == 0;
                    $right_scaffold->delete;
                }
            );
        }
    );
    $self->info_msg("Scaffold merge complete!");
    return 1;
}
1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Assembly.pm $
#$Id: Assembly.pm 31361 2007-12-28 18:16:32Z adukes $
