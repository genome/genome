package Finishing::Assembly::DBIx::Proxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::Proxy';

use Data::Dumper;
use Finfo::ClassUtils 'class';

sub methods_for_source_method : CUMULATIVE
{
    return ( id => undef );
}

sub update
{
    return shift->source->update;
}

sub update_attribute
{
    my ($self, $attribute, $value) = @_;

    $self->source->$attribute($value);

    return $self->source->update;
}

##UNUSED METHOD, LEFT IN CASE SOMETHING BREAKS 12/31/07 
#sub get_attribute
#{
#    my ($self, $attribute) = @_;

#    return $self->source->$attribute;
#}


sub delete
{
    return shift->source->delete;
}

#- SCHEMA/RESULTSET -#
sub schema
{
    my $self = shift;

    return $self->source->result_source->schema;
}

#- METHODS CALLED FROM CHILD PROXIES -#
sub related : RESTRICTED
{
    my ($self, %p) = @_;

    my $method_params = $p{params};
    my ($method) = sprintf('_%s_related', shift(@$method_params));
    
    return $self->$method($method_params, $p{args});
}

sub _find_related : PRIVATE
{
    my ($self, $method_params, $args) = @_;

    my ($related_method, $object_type) = @$method_params;
    my %params = @$args;
    Finfo::Validate->validate
    (
        attr => "params for get $object_type",
        value => \%params,
        ds => 'hashref',
        msg => 'fatal',
        'caller' => [ caller ],
        #caller_level => 2,
    );

    my $object = $self->source->$related_method->find(\%params);
    $self->fatal_msg
    (
        sprintf
        (
            "Can't find %s with params (%s) for %s",
            $object_type,
            join('', keys %params),
            class( $self->source ),
        ),
        { caller_level => 2 },
    ) unless $object;

    return $self->_construct_object($object_type, $object);
}

sub _search_related : RESTRICTED
{
    my ($self, $method_params, $args) = @_;

    my ($related_method, $object_type, $return_iterator) = @$method_params;
    my %params = @$args;
    if ( %params )
    {
        Finfo::Validate->validate
        (
            attr => "params for search $object_type",
            value => \%params,
            ds => 'hashref',
            msg => 'fatal',
            'caller' => [ caller ],
        );
    }

    my $rs = $self->source->$related_method->search(\%params);
    $self->fatal_msg
    (
        sprintf
        (
            "Can't find %s with search params (%s) for %s",
            $object_type,
            join('', keys %params),
            class( $self->source ),
        ),
    ) unless $rs and $rs->count > 0;

    return ( $return_iterator ) 
    ? $self->_construct_object($object_type . '_iterator', $rs)
    : $self->_construct_objects($object_type, $rs->all);
}

sub _create_related : RESTRICTED
{
    my ($self, $method_params, $args) = @_;

    my ($related_method, $object_type) = @$method_params;
    my %params = @$args;
    Finfo::Validate->validate
    (
        attr => "params for create $object_type",
        value => \%params,
        ds => 'hashref',
        msg => 'fatal',
        'caller' => [ caller ],
    );
    
    my $existing_object = $self->source->$related_method->find(\%params);
    $self->fatal_msg
    (
        sprintf
        (
            "Can't create %s for %s with params (%s), it already exists",
            $object_type,
            class( $self->source ),
            join('', keys %params),
        ),
    ) if $existing_object;

    my $object = $self->source->create_related($related_method, \%params);
    $self->fatal_msg
    (
        sprintf
        (
            "Error creating %s for %s with params (%s)",
            $object_type,
            class( $self->source ),
            join('', keys %params),
        ),
    ) unless $object;

    return $self->_construct_object($object_type, $object);
}

sub _add_to_related : RESTRICTED
{
    my ($self, $method_params, $args) = @_;

    my ($relationship, $object_type, @required_params) = @$method_params;
    my %params = @$args;
    Finfo::Validate->validate
    (
        attr => "params for create $object_type",
        value => \%params,
        ds => 'hashref',
        msg => 'fatal',
        'caller' => [ caller ],
    );
    
    my $add_to_method = "add_to_$relationship";
    
    #my @missing_params;  #TODO necessary?
    #foreach $req (@required_params){
    #    push @missing_params, $_ unless grep { $_ eq $req } keys %$params;
    #}
    #$self->fatal_msg("Missing the following required params for $add_to_method: ".join(" ", @missing_params)) if @missing_params;
    

    my $object = $self->source->$add_to_method(\%params);
    $self->fatal_msg
    (
        sprintf
        (
            "Error adding %s to %s with params (%s)",
            $object_type,
            class( $self->source ),
            join('', keys %params),
        ),
    ) unless $object;

    return $self->_construct_object($object_type, $object);
}

sub _find_or_create_related : RESTRICTED
{
    my ($self, $method_params, $args) = @_;

    my ($related_method, $object_type) = @$method_params;
    my %params = @$args;
    Finfo::Validate->validate
    (
        attr => "params for get or create $object_type",
        value => \%params,
        ds => 'hashref',
        msg => 'fatal',
        'caller' => [ caller ],
    );

    my $object = $self->source->$related_method->find(\%params);
    
    return $self->_construct_object($object_type, $object) if $object;

    $object = $self->source->create_related($related_method, \%params);
    $self->fatal_msg
    (
        sprintf
        (
            "Error creating %s for %s with params (%s)",
            $object_type,
            class( $self->source ),
            join('', keys %params),
        ),
    ) unless $object;

    return $self->_construct_object($object_type, $object);
}

#- NAMES -#
sub _names_related : RESTRICTED
{
    my ($self, $method_params) = @_;

    my ($related_method) = $method_params->[0];

    my @names;
    my $rs = $self->source->$related_method;
    while ( my $object = $rs->next )
    {
        push @names, $object->name;
    }
    
    return @names;
}

#- COUNT -#
sub _count_related : RESTRICTED
{
    my ($self, $method_params) = @_;
   
    my ($related_method) = $method_params->[0];

    return $self->source->$related_method->count;
}

###########################################################################################

package Finishing::Assembly::DBIx::OrganismProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::DBIx::Proxy'; 

sub methods_for_source_method : CUMULATIVE
{
    return  
    (
        name => undef,
        # chr
        get_chromosome => 'chromosome',
        chromosomes => 'chromosome_iterator',
        # assem
        get_assembly => 'assembly',
        assemblies => 'assembly_iterator',
    );
}

sub methods : CUMULATIVE
{
    return (qw/ related /);
}

sub methods_for_related : RESTRICTED
{
    return 
    (
        # assem
        create_assembly => [qw/ create assemblies assembly /],
        get_or_create_assembly => [qw/ find_or_create assemblies assembly /],
        get_assembly_names => [qw/ names assemblies /],
        # chr
        create_chromosome => [qw/ create chromosomes chromosome /],
        get_or_create_chromosome => [qw/ find_or_create chromosomes chromosome /],
    );
}

###########################################################################################

package Finishing::Assembly::DBIx::AssemblyProxy;

use strict;
use warnings;
no warnings 'reserved';

use base (qw/Finishing::Assembly::DBIx::Proxy/); 

sub methods_for_source_method : CUMULATIVE
{
    return 
    (
        name => undef,
        file_path => undef,
        organism => 'organism',
        # scaff
        scaffolds => 'scaffold_iterator',
        get_scaffold => 'scaffold',
        # contig
        contigs => 'contig_iterator',
        get_contig => 'contig',
        # reads
        assembled_reads => 'assembled_read_iterator',
        get_assembled_read => 'assembled_read',
        # imp cor
        improvement_correlations => 'improvement_correlation_iterator',
        tags => 'assembly_tag',
    );
}

sub methods : CUMULATIVE
{
    return (qw/ related /);
}

sub methods_for_related : RESTRICTED
{
    return 
    (
        # scaff
        create_scaffold => [qw/ create scaffolds scaffold /],
        get_or_create_scaffold => [qw/ find_or_create scaffolds scaffold /],
        #add_to
        add_to_tags => [qw/ add_to tags assembly_tag /],
               
    );
}

sub add_tags
{
    my ($self, @tags) = @_;
    $self->add_tag($_) foreach @tags;
}

sub add_tag
{
    my ($self, $tag) = @_;
    return $self->add_to_tags( #TODO id needed?
        
            tag_type        => $tag->type,
            source         => $tag->source,
            creation_date   => $tag->date, #TODO date_conversion
            data            => $tag->data,
        
    );
}


#############################################################################################

package Finishing::Assembly::DBIx::SequencedItemProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::DBIx::Proxy'; 

#- BASES/QUALs -#
sub base_string
{
    my ($self, $base_string) = @_;

    my $seq_method = $self->_sequence_method;
    my $seq = $self->source->$seq_method;
    return unless $seq;

    if ( $base_string )
    {
        $seq->bases($base_string);
        $seq->update;
    }

    return $seq->bases;
}

sub qualities
{
    my ($self, $new_quals) = @_;

    my $seq_method = $self->_sequence_method;
    my $seq = $self->source->$seq_method;
    return  unless $seq;

    if ( $new_quals )
    {
        $seq->qualities( join(' ', @$new_quals) );
        $seq->update;
        return $new_quals;
    }

    # FIXME temp 'til pads are removed from db
    #return [ split(/\s+/, $seq->qualities) ];
    
    my $qual_string =  $seq->qualities;
    $qual_string =~ s/ \*//g;
    return [ split(/\s+/, $qual_string) ];
}

sub create_sequence{
    my ($self, %p) = @_;
    my $seq_method = $self->_sequence_method;
    my $seq = $self->source->create_related($seq_method, \%p);
    $self->fatal_msg("Couldn't create $seq_method with params %p") unless $seq;
    return $seq;
}

sub length
{
    my $self = shift;

    return $self->source->length;
}

###########################################################################################

package Finishing::Assembly::DBIx::AssembledReadProxy;

use strict;
use warnings;
no warnings 'reserved';

use base (qw/Finishing::Assembly::DBIx::SequencedItemProxy/); 

sub methods_for_source_method : CUMULATIVE 
{
    return 
    (
        name => undef,
        template_name => undef,
        direction => undef,
        length => undef,
        start_position => undef,
        stop_position => undef,
        complemented => undef, 
        # org
        organism_id => undef,
        get_organism => 'organism',
        # assem
        assembly_id => undef,
        assembly => 'assembly',
        # scaff
        scaffold_id => undef,
        scaffold => 'scaffold',
        # ctg
        contig_id => undef,
        contig => 'contig',
        # tags
        tags => 'read_tag',
    );
}

sub methods : CUMULATIVE
{
    return (qw/ related /);
}

sub methods_for_related : RESTRICTED
{
    return 
    (
        add_to_tags => [qw/ add_to tags read_tag /],
    )
}

sub _sequence_method
{
    return 'sequence';
}

#- ALIASES -#
sub orientation
{
    my $self = shift;

    return ( $self->source->complemented ) ? 'C' : 'U';
}

sub position
{
    my $self = shift;

    return $self->source->start_position(@_);
}

sub add_tag{
    my ($self, $tag) = @_;
    return $self->add_to_tags( #TODO id needed?
        
            tag_type        => $tag->type,
            source         => $tag->source,
            start_position  => $tag->start,
            stop_position   => $tag->stop,
            creation_date   => $tag->date, #TODO date_conversion
            text            => $tag->text,
            comment         => $tag->comment,
        
    );
}


###########################################################################################

package Finishing::Assembly::DBIx::ContigProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::DBIx::SequencedItemProxy'; 

use Finishing::Assembly::Sequence;

my %_sequence :name(_sequence:p) :isa('object Finishing::Assembly::Sequence');

sub methods_for_source_method : CUMULATIVE
{
    return 
    (
        name => undef,
        contig_num => undef,
        length => undef,
        # chr
        chromosome_id => undef, 
        chromosome => 'chromosome',
        ordered_on_chromosome=> undef, 
        # assem
        assembly_id => undef,
        assembly => 'assembly',
        # scaff
        scaffold_id => undef,
        scaffold => 'scaffold',
        # ctgs
        left_contig => 'contig',
        right_contig => 'contig',
        # tags
        tags => 'consensus_tag',
        # projects
        project => 'project',
        #reads
        assembled_reads => 'assembled_read_iterator',
        get_assembled_read => 'assembled_read',
        #improv correl
        improvement_correlations => 'improvement_correlation_iterator',
        #gaps
        left_gaps => 'gap_iterator',
        right_gaps => 'gap_iterator',
    );
}    


sub methods : CUMULATIVE
{
    return (qw/ related /);
}

sub methods_for_related : RESTRICTED 
{
    return 
    (
        #tags
        add_to_tags => [qw/ add_to tags consensus_tag /],
        add_to_assembled_reads => [qw/ add_to assembled_reads assembled_read/]
    );
}

sub _sequence_method
{
    return 'consensus';
}

#- COMPLEMENTED -#
sub complemented
{
    return 0;
}

#- NAMING -#
sub scaffold_num
{
    my $self = shift;

    return $self->source->scaffold->scaffold_num;
}

#- SCAFFOLD -#
sub set_scaffold{
    my ($self, $scaffold) = @_;
    $self->fatal_msg("Need a scaffold object!") unless ref $scaffold =~ /Finishing::Assembly::Scaffold/;
    my $id = $scaffold->id;
    $self->fatal_msg("Need a scaffold object with an id!") unless $id;
    $self->_set_scaffold_id($id);
}

sub _set_scaffold_id :RESTRICTED
{
    my ($self, $id) = @_;

    $self->update_attribute('scaffold_id', $id);
}

#- CHROMOSOME -#
sub set_chromosome 
{
    my ($self, $chr, $ordered) = @_;

    $self->update_attribute('chromosome_id', $chr->id) if $chr;
    $self->update_attribute('ordered_on_chromosome', $ordered) if defined $ordered;

    return $self->_construct_object('chromosome', $self->source->chromosome);
}

#- GAPS/NEIGHBOR CTGS -#
sub create_left_gap
{
    my ($self, %params) = @_;

    %params = $self->dummy_gap_params if $params{dummy};

    my $gap = $self->source->create_related('left_gaps', \%params);

    return $self->_construct_object('gap', $gap);
}

sub create_right_gap
{
    my ($self, %params) = @_;

    %params = $self->dummy_gap_params if $params{dummy};
    
    my $gap = $self->source->create_related('right_gaps', \%params);

    return $self->_construct_object('gap', $gap);
}

sub dummy_gap_params
{
    return (
        length => 0,
        type => 'dummy',
        linkage => 0, 
    );
}

sub set_right_contig
{
    my $self = shift;

    return $self->_set_neighbor_contig('right', @_);
}

sub set_left_contig
{
    my $self = shift;

    return $self->_set_neighbor_contig('left', @_);
}

sub _set_neighbor_contig : PRIVATE
{
    my ($self, $direction, $contig) = @_;

    my $gaps = $direction . "_gaps";
    my $contig_method = $direction . "_contig";

    my $gapped = 0;
    foreach my $gap ($self->source->$gaps->all)
    {
        #dealing with the dbix objects here so we can set directly
        $gap->$contig_method($contig->source);
        $gap->update;
        $gapped++;
    }
    $self->error_msg("No gaps to set  $direction neighbor contig on! create $direction gap first!") unless $gapped;

    return $self->_construct_object('contig', $self->source->$contig_method);
}

#- READS -#
sub add_assembled_read
{
    my ($self, $read) = @_;
    my $source = $self->source->add_to_assembled_reads(
        {
            name            => $read->name,
            #contig_id       => $self->id, #TODO, is this needed
            template_name   => $read->template,
            direction       => $read->direction, 
            complemented    => $read->complemented,
            length          => $read->length,
            start_position  => $read->position, 
            stop_position   => $read->position + $read->length - 1,
            assembly_id     => $self->source->assembly_id,
        }
    );
    return $self->_construct_object('assembled_read', $source);
}

#- TAGS -#
sub add_tag
{
    my ($self, $tag) = @_;
    
    my $date    = sprintf("TO_DATE('%s', 'YYYY-MM-DD HH24:MI:SS')", $tag->date);
    my $notrans = $tag->no_trans ? 1 : 0;

    my $arg = {
        tag_type       => $tag->type,
        no_trans       => $notrans,
        source         => $tag->source,
        start_position => $tag->start,
        stop_position  => $tag->stop,
        creation_date  => \$date, #TODO date conversion
    };
    
    $arg->{text}        = $tag->text     if $tag->text;
    $arg->{tag_comment} = $tag->comment  if $tag->comment;

    my $source = $self->source->add_to_tags($arg);
    return $self->_construct_object('consensus_tag', $source);
}

###########################################################################################

package Finishing::Assembly::DBIx::ChromosomeProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::DBIx::Proxy'; 

sub methods_for_source_method : CUMULATIVE
{
    return
    (
        name => undef,
        first_scaffold_for_assembly => 'chromosome_first_scaffold',
        set_first_scaffold => 'chromosome_first_scaffold',
        organism_id => undef,
        organism => 'organism',
        first_scaffolds => 'chromosome_first_scaffold_iterator',
    );
}

sub methods : CUMULATIVE
{
    return (qw/ related /);
}

sub methods_for_related : RESTRICTED
{
    return 
    (
    );
}

###########################################################################################

package Finishing::Assembly::DBIx::ChromosomeFirstScaffoldProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::DBIx::Proxy'; 

sub methods_for_source_method : CUMULATIVE
{
    return  
    (
        name => undef,
        chromosome => 'chromosome',
        assembly => 'assembly',
        scaffold => 'scaffold',
    );
}

###########################################################################################

package Finishing::Assembly::DBIx::ScaffoldProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::DBIx::Proxy'; 

sub methods_for_source_method : CUMULATIVE
{
    return 
    (
        scaffold_num => undef,
        orientation => undef,
        'length' => undef,
        assembly => 'assembly',
        chromosome => 'chromosome',
        get_contig => 'contig',
        first_contig => 'contig',
        contigs => 'contig_iterator',
        ordered_contigs => 'contig_iterator',
    );
}

sub methods : CUMULATIVE
{
    return (qw/ related /);
}

sub methods_for_related : RESTRICTED
{
    return 
    (
        create_contig =>  [qw/ create contigs contig /],
        get_or_create_contig =>  [qw/ find_or_create contigs contig /],
    );
}

sub add_contig
{
    my ($self, $contig, $contig_num) = @_;
    $contig_num = $contig->contig_num unless defined $contig_num;
    my $source = $self->source->add_to_contigs(
        {
            #scaffold_id => $self->source->id,
            assembly_id => $self->source->assembly->id,
            contig_num => $contig_num,
            length => $contig->length,
        }
    );
    return $self->_construct_object('contig', $source);
}

sub dump{
    my ($self,$object) = @_;
    my $assembly = $self->assembly;
    my $string = "Assembly:".$assembly->name."\n";
    $string .= "Scaffold:".$self->name."orientation".$self->orientation."\n";
    my $ci = $self->source->contigs;
    while (my $contig = $ci->next){
        $string .="Contig".$contig->contig_num.".".$contig->scaffold->scaffold_num.":id-".$contig->id.", length:".$contig->length.", read count:".$contig->assembled_reads->count."\n";
    }
    return $string;
}

###########################################################################################

package Finishing::Assembly::DBIx::ImprovementCorrelationProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::DBIx::Proxy'; 

sub methods_for_source_method : CUMULATIVE
{
    return  
    (
        score => undef,
        # ctgs
        contigs => 'contig_iterator',
    );
}

sub methods : CUMULATIVE
{
    return (qw/ related /);
}

sub methods_for_related : RESTRICTED
{
    return 
    (
    );
}

###########################################################################################

package Finishing::Assembly::DBIx::ProjectProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::DBIx::Proxy'; 

sub methods_for_source_method : CUMULATIVE
{
    return
    (
        name => undef,
        owner => undef,
        priority => undef,
        directory => undef,
        project_type => undef,
        organism => 'organism',
        contigs => 'contig_iterator',
        created => undef,
        modified => undef,
    );
}

sub methods : CUMULATIVE
{
    return (qw/ related /);
}

sub methods_for_related : RESTRICTED
{
    return 
    (
        contig_count => [qw/ count contigs /],
    );
}

#sub name
#{
#    my $self = shift;

#    return $self->source->id;
#}

sub organism_name
{
    my $self = shift;

    return $self->source->organism->name;
}

sub organism_nickname
{
    my $self = shift;

    return $self->source->organism->nickname;
}

sub dump{
    my ($self,$object) = @_;
    my $assembly = $self->assembly;
    my $string = "Assembly:".$assembly->name."\n";
    my $contigs = $self->contig_iterator;
    while (my $contig = $contigs->next){
        $string .=$contig->name.":id-".$contig->id.",length:".$contig->length."read count:".$contig->read_count."\n";
    }
    return $string;
}

###########################################################################################

package Finishing::Assembly::DBIx::GapProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::DBIx::Proxy'; 

sub methods_for_source_method : CUMULATIVE
{
    return  
    (
        right_contig_id => undef,
        right_contig => 'contig',
        left_contig_id => undef,
        left_contig => 'contig',
        length => undef,
        type => undef,
        linkage => undef,
    );
}

sub set_left_contig
{
    my ($self, $new_contig) = @_;

    $self->fatal_msg("need new contig!") unless $new_contig;

    $self->update_attribute('left_contig_id', $new_contig->id) and return 1;
}

sub set_right_contig
{
    my ($self, $new_contig) = @_;

    $self->fatal_msg("need new contig!") unless $new_contig;

    $self->update_attribute('right_contig_id', $new_contig->id) and return 1;
}

###########################################################################################

package Finishing::Assembly::DBIx::TagProxy;

use base 'Finishing::Assembly::DBIx::Proxy';

sub methods_for_source_method : CUMULATIVE
{
    return 
    (
        source => undef,
        text => undef,
        comment => undef, 
        creation_date => undef,
    );
}


sub date
{
    my $self = shift;
    return $self->source->creation_date;
}


###########################################################################################

package Finishing::Assembly::DBIx::SequenceTagProxy;

use base 'Finishing::Assembly::DBIx::TagProxy';

sub methods_for_source_method : CUMULATIVE
{
    return
    (
        tag_type => undef,
        start_position => undef,
        stop_position => undef,
    );
}

sub start
{
    my $self = shift;

    return $self->source->start_position(@_);
}

sub stop
{
    my $self = shift;

    return $self->source->stop_position(@_);
}

sub type
{
    my $self = shift;

    return $self->source->tag_type(@_);
}

###########################################################################################

package Finishing::Assembly::DBIx::ConsensusTagProxy;

use base 'Finishing::Assembly::DBIx::SequenceTagProxy';

my %is_oligo :name(is_oligo:p) :isa(boolean);
my %is_auto_finish_exp :name(is_auto_finish_exp:p) :isa(boolean);

sub START{
    my $self = shift;
    $self->is_oligo(1) if $self->source->tag_type eq 'oligo';
    $self->is_auto_finish_exp(1) if $self->source->tag_type eq 'autoFinishExp';
    return 1;
}

sub methods_for_source_method : CUMULATIVE
{
    my $self = shift;
    my %source_methods = (
        contig_id => undef,
        contig => 'contig',
        no_trans => undef,
    );
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

sub parent{
    return shift->contig->name;
}

###########################################################################################

package Finishing::Assembly::DBIx::ReadTagProxy;

use base 'Finishing::Assembly::DBIx::SequenceTagProxy';

sub methods_for_source_method : CUMULATIVE
{
    return 
    (
        read_id => undef,
        assembled_read => 'assembled_read',
    );
}

sub parent{
    return shift->assembled_read->name;
}

###########################################################################################
package Finishing::Assembly::DBIx::AssemblyTagProxy;

use base 'Finishing::Assembly::DBIx::TagProxy';

sub methods_for_source_method : CUMULATIVE
{
    return 
    (
        assembly_id => undef, 
        assembly => 'assembly',
    );
}

sub parent{
    return shift->assembly->name;
}

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/DBIxProxies.pm $
#$Id: DBIxProxies.pm 31442 2008-01-03 23:47:59Z adukes $
