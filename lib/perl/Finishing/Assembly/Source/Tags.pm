package Finishing::Assembly::Source::Tag;
{

    use strict;
    use warnings;
    
    use base 'Finfo::Accessor';

    __PACKAGE__->mk_accessors(qw/ type source date text comment scope /);
}

##################################################################

package Finishing::Assembly::Source::AssemblyTag;
{

    use strict;
    use warnings;
    
    use base 'Finishing::Assembly::Source::Tag';
}

##################################################################

package Finishing::Assembly::Source::SequenceTag;
{

    use strict;
    use warnings;
    
    use base 'Finishing::Assembly::Source::Tag';

    __PACKAGE__->mk_accessors(qw/ parent start stop unpad_start unpad_stop /);
}

##################################################################

package Finishing::Assembly::Source::ConsensusTag;
{
    use base 'Finishing::Assembly::Source::SequenceTag';
    
    __PACKAGE__->mk_accessors(qw/ no_trans /);
}

##################################################################

package Finishing::Assembly::Source::AutoFinishExpTag;
{
    use strict;
    use warnings;
              
    use base 'Finishing::Assembly::Source::ConsensusTag';

    __PACKAGE__->mk_accessors
    (qw/ 
        orientation num1 num2 num3 chem primer_type purpose 
        fix_cons_errors original_cons_errors original_single_subclone_bases 
        primer temp id exp_id_and_template
        oligo_name oligo_seq oligo_temp
        /);
}

##################################################################

package Finishing::Assembly::Source::OligoTag;
{
    use strict;
    use warnings;

    use base 'Finishing::Assembly::Source::ConsensusTag';

    __PACKAGE__->mk_accessors(qw/ oligo_name oligo_seq oligo_temp oligo_templates complemented /);

    sub orientation{
        my $self = shift;
        $self->fatal_msg("orientation is not a setter, use complemented()") if @_;
        return $self->complemented? 'C' : 'U';
    }
}

##################################################################

package Finishing::Assembly::Source::ReadTag;
{
    use strict;
    use warnings;

    use base 'Finishing::Assembly::Source::SequenceTag';
}

1;

=pod
=cut

#$HeadURL$
#$Id$
