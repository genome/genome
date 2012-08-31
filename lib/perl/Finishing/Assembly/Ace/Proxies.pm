####################
#- ASSEMBLY PROXY -#
####################

package Finishing::Assembly::Ace::AssemblyProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::Proxy';

use Storable;

sub methods_for_source_method
{
    return 
    (
        name => undef,
        # contig
        contig_count => undef,
        contigs => 'contig_iterator',
        get_contig => 'contig',
        delete_contig => undef,
        # read
        read_count => undef,
        assembled_reads => 'assembled_read_iterator',
        get_assembled_read => 'assembled_read',
        tags => 'assembly_tag',
    );
}

sub get_longest_contig{
    my $self = shift;
    my ($contig )= sort {$b->length <=> $a->length} $self->source->contigs->all;
    $self->_construct_object('contig', $contig);
}

##########################
#- SEQUENCED ITEM PROXY -#
##########################

package Finishing::Assembly::Ace::SequencedItemProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::Proxy';

use Finishing::Assembly::Sequence;

my %seq :name(_seq:p) :isa('object');

sub _sequence
{
    my $self = shift;

    my $caller = caller(1);

    $self->fatal_msg
    (
        "Private method (_sequence) can not be accessed from an outside caller ($caller)"
    ) unless grep 
    {
        $caller eq "Finishing::Assembly::$_" 
    } (qw/ AssembledRead Contig SequencedItem /);
    
    #return $self->_seq or $self->_seq
    return $self->_seq
    (
        Finishing::Assembly::Sequence->new
        (
            base_string => $self->source->base_string,
            qualities => $self->source->qualities,
        )
    )
}

sub base_string
{
    my ($self, $bases) = @_;

    if ( $bases )
    {
        $self->undef_attribute('_seq');
        return $self->source->base_string($bases);
    }

    return $self->source->base_string;
}

sub qualities
{
    my ($self, $qualities) = @_;

    if ( $qualities )
    {
        $self->undef_attribute('_seq');
        return $self->source->qualities($qualities);
    }

    return $self->source->qualities;
}

##################
#- CONTIG PROXY -#
##################

package Finishing::Assembly::Ace::ContigProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::Ace::SequencedItemProxy';

sub methods_for_source_method
{
    return 
    (
        name => undef,
        complemented => undef,
        assembled_reads => 'assembled_read_iterator',
        get_assembled_read => 'assembled_read',
        read_count => undef,
        base_segments => undef,
        tags => 'consensus_tag',
    );
}

sub contig_num{
    my $self = shift;
    my ($num) = $self->source->name =~ /Contig(?:\d+\.)?(\d+)/;
    return $num;
}

################
#- READ PROXY -#
################

package Finishing::Assembly::Ace::AssembledReadProxy;

use base 'Finishing::Assembly::Ace::SequencedItemProxy';

sub methods_for_source_method
{
    return 
    (
        name => undef,
        position => undef,
        complemented => undef,
        tags => 'read_tag',
        qual_clip_start => undef,
        qual_clip_stop => undef,
        align_clip_start => undef,
        align_clip_stop => undef,
        'time' => undef,
        chromat_file => undef,
        phd_file => undef,
        chem => undef,
        dye => undef,
        info_count => undef,
        'length' => undef,
    );
}

sub qualities
{
    return [];
}

sub direction{
    my $self = shift;
    my ($template, $direction) = $self->source->name =~ /(.+)\.([^0-9]+)\d+$/;

    if ($direction){
        if ($direction eq 'b') {
            $direction = 'f';
        }elsif ($direction eq 'g') {
            $direction = 'r';
        }
    }else{
        $direction = 'u';
    }
    return $direction;
}

#################
#- TAG PROXIES -#
#################

package Finishing::Assembly::Ace::TagProxy;

use base 'Finishing::Assembly::Proxy';

sub methods_for_source_method : CUMULATIVE
{
    return 
    (
        type => undef,
        source => undef,
        date => undef,
        text => undef,
        comment => undef,
    );
}

#-#

package Finishing::Assembly::Ace::AssemblyTagProxy;

use base 'Finishing::Assembly::Ace::TagProxy';

#-#

package Finishing::Assembly::Ace::SequenceTagProxy;

use base 'Finishing::Assembly::Ace::TagProxy';

sub methods_for_source_method : CUMULATIVE
{
    return
    (
        parent => undef,
        start => undef,
        stop => undef,
        unpad_start => undef,
        unpad_stop => undef,
        scope => undef,
    );
}

#-#

package Finishing::Assembly::Ace::ConsensusTagProxy;

use base 'Finishing::Assembly::Ace::SequenceTagProxy';

my %is_oligo :name(is_oligo:p) :isa(boolean);
my %is_auto_finish_exp :name(is_auto_finish_exp:p) :isa(boolean);

sub START{
    my $self = shift;
    $self->is_oligo(1) if $self->source->type eq 'oligo';
    $self->is_auto_finish_exp(1) if $self->source->type eq 'autoFinishExp';
    return 1;
}

sub methods_for_source_method : CUMULATIVE
{
    my $self=shift;
    my %source_methods =  ( no_trans => undef );
    my %oligo_methods = (
        oligo_name  => undef,
        oligo_seq  => undef,
        oligo_temp  => undef,
        oligo_templates  => undef,
        complemented => undef,
        orientation => undef,
    );
    my %auto_exp_methods = (

        orientation => undef,
        num1 => undef,
        num2 => undef,
        num3 => undef,
        chem => undef,
        primer_type => undef,
        purpose  => undef,
        fix_cons_errors => undef,
        original_cons_errors => undef,
        original_single_subclone_bases => undef, 
        primer => undef,
        temp => undef,
        id => undef,
        exp_id_and_template => undef,
        oligo_name => undef,
        oligo_seq => undef,
        oligo_temp => undef,
    );
    return (%source_methods, %oligo_methods) if $self->is_oligo;
    return (%source_methods, %auto_exp_methods) if $self->is_auto_finish_exp;
    return %source_methods;
}

#-#

package Finishing::Assembly::Ace::ReadTagProxy;

use base 'Finishing::Assembly::Ace::SequenceTagProxy';

1;

#$HeadURL$
#$Id$
