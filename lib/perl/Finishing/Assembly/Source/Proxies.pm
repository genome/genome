package Finishing::Assembly::Source::ContigProxy;

use base 'Finishing::Assembly::Proxy';

sub methods_for_source_method
{
    return 
    (
        name => undef,
        assembled_reads => 'assembled_read_iterator',
        base_string => undef,
        qualities => undef,
        tags => 'consensus_tag',
		base_segments => undef,
		get_assembled_read => 'assembled_read',
        align => undef,
        complemented => undef,
        create_tag => 'consensus_tag',
        add_tags => 'consensus_tag',
        contig_num => undef,
        scaffold_num => undef,
    );
}

##############################################################

package Finishing::Assembly::Source::AssembledReadProxy;

use base 'Finishing::Assembly::Proxy';

sub methods_for_source_method
{
	return 
    (
        name => undef,
        base_string => undef,
        qualities => undef,
        chromat_positions => undef,
        position => undef,
        complemented => undef,
        tags => 'read_tag',
        remove_tags => 'read_tag',
        add_tags => 'read_tag',
        comments => undef,
        'time' => undef,
        chromat_file => undef,
        phd_file=> undef,
        chem => undef,
        dye => undef,
        'length' => undef,
        wr => undef,
        template => undef,
        primer => undef,
        lib => undef,
    );
}

##############################################################

package Finishing::Assembly::Source::ProjectProxy;

use base 'Finishing::Assembly::Proxy';

sub methods_for_source_method : CUMULATIVE
{
    return 
    (
        name => undef,
        directory => undef,
    );
}

##############################################################

#################
#- TAG PROXIES -#
#################

package Finishing::Assembly::Source::TagProxy;

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

package Finishing::Assembly::Source::AssemblyTagProxy;

use base 'Finishing::Assembly::Source::TagProxy';

#-#

package Finishing::Assembly::Source::SequenceTagProxy;

use base 'Finishing::Assembly::Source::TagProxy';

sub methods_for_source_method : CUMULATIVE
{
    return
    (
        parent => undef,
        scope => undef,
        start => undef,
        stop => undef,
        unpad_start => undef,
        unpad_stop => undef,
    );
}

#-#

package Finishing::Assembly::Source::ConsensusTagProxy;

use base 'Finishing::Assembly::Source::SequenceTagProxy';

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

#vvv AUTOEXP and OLIGO methods go here vvv
#-#

package Finishing::Assembly::Source::ReadTagProxy;

use base 'Finishing::Assembly::Source::SequenceTagProxy';

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Proxies.pm $
#$Id: Proxies.pm 31129 2007-12-18 20:41:16Z ebelter $
