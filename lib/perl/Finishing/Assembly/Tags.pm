package Finishing::Assembly::Tag;
{
    use strict;
    use warnings;
    no warnings 'reserved';

    use base 'Finishing::Assembly::Item';
}

######################################################################

package Finishing::Assembly::SequenceTag;
{
    use strict;
    use warnings;
    no warnings 'reserved';

    use base 'Finishing::Assembly::Tag';

    my %unpad_start :name(unpad_start:o);
    my %unpad_stop :name(unpad_stop:o);

    sub length
    {
        my $self = shift;

        return $self->stop - $self->start + 1;
    }
}

######################################################################

package Finishing::Assembly::ReadTag;
{
    use base 'Finishing::Assembly::SequenceTag';

    sub read_name
    {
        assembled_read_name();
    }

    sub assembled_read_name
    {
        return shift->parent;
    }

    1;
}

######################################################################

package Finishing::Assembly::ConsensusTag;

use base 'Finishing::Assembly::SequenceTag';

sub contig_name
{
    return shift->parent;
}

######################################################################
#TODO removed these consensus tag types so this information can be handled at proxy level,
#may need to move some of these methods down to proxy
######################################################################
=cut
package Finishing::Assembly::OligoTag;
{
    use strict;
    use warnings;
    no warnings 'reserved';

    use base 'Finishing::Assembly::ConsensusTag';

    sub oligo_num
    {
        my $self = shift;

        my ($num) = $self->oligo_name =~ /\.(\d+)$/;

        return $num;
    }

    1;
}

######################################################################

package Finishing::Assembly::AutoFinishExpTag;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::ConsensusTag';

sub ace
{
    my $self = shift;

    return unless defined $self->id;

    my ($ace) = $self->id =~ /^(.+)\.\d+$/;

    return $ace;
}

sub exp_ids_and_templates
{
    my $self = shift;

    return $self->_exp_ids_and_templates if $self->_exp_ids_and_templates;

    my %ids_and_templates;
    foreach my $rxn ( split /, /, $self->exp_id_and_template )
    {
        my ($id, $template)  = split / /, $rxn;
        $ids_and_templates{$id} = $template;
    }

    return $self->_exp_ids_and_templates(\%ids_and_templates);
}

sub exp_ids
{
    my $self = shift;

    my $ids_and_templates = $self->exp_ids_and_templates;

    return unless defined $ids_and_templates and %$ids_and_templates;
    
    return keys %$ids_and_templates;
}

sub template_for_id
{
    my ($self, $id) = @_;

    die "Need exp id to get template\n" unless defined $id;
    
    my $ids_and_temps = $self->exp_ids_and_templates;

    return unless defined $ids_and_temps and %$ids_and_temps;
    
    return $ids_and_temps->{$id};
}
=cut
######################################################################

package Finishing::Assembly::AssemblyTag;

use base 'Finishing::Assembly::Tag';

######################################################################

package Finishing::Assembly::TagInfo;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finfo::Singleton';

require Date::Format;

sub timestamp
{
    # 050427:142133
    return Date::Format::time2str('%y%m%d:%H%M%S', time);
}

sub tag_sources
{
    return sort # list not complete, no consed constraint on source
    (qw/
        consed consed_duplicated_from_consed
        /);
}

sub tag_types
{
    return sort # should be complete
    (qw/
        Annotation G_dropout QA SingleSubclone Transposon autoFinishExp becomeConsensus
        changedGenotype chimera cloneEnd cloningVector comment compression consedFixedGoldenPath
        contigEndPair contigName coordinatorApproval coordinatorComment coordinatorComment dataNeeded
        doFinish doNotDoPCR doNotFinish ecoliContamination edit editable heterozygoteAC heterozygoteAG
        heterozygoteAT heterozygoteCG heterozygoteCT heterozygoteGT heterozygoteIndel homozygoteAA
        homozygoteCC homozygoteGG homozygoteIndel homozygoteTT ignoreMatches ignoreMismatches indelSite
        markedHighQuality markedLowQuality matchElsewhereHighQual matchElsewhereLowQual oligo
        oligo3PrimeEnd oligoPrimerSequence phred30 phred30conUnedited polyPhredRank1 polyPhredRank2
        polyPhredRank3 polyPhredRank4 polyPhredRank5 polyPhredRank6 polymorphism polymorphismConfirmed
        polymorphismDeletion polymorphismInsertion polymorphismSubstitution prefinish qualityCoreComment
        repeat sequencingVector significantDiscrepancy startNumberingConsensus stolenconsensus
        stolendata structuralDeletion tagsOverlap vector annotation-AmbiguousBase
        annotation-DigestComments annotation-EndOfClone annotation-GSSmRNAOnly
        annotation-NonRepetitiveButUnresolvedRegion annotation-OtherComments
        annotation-PCROnly annotation-PCROnlyGenomic annotation-RepresubmitInfo
        annotation-SingleSubclone annotation-Shatter annotation-TransposonExcisedFromFinishedRegion
        annotation-TransposonInFinishedRegion annotation-TransposonInVector
        annotation-UnresolvedDuplication annotation-UnresolvedInverted
        annotation-UnresolvedMono annotation-UnresolvedSsr annotation-UnresolvedTandem
        /);
}

sub is_tag_source_valid
{
    my $self = shift;

    return $self->_is_valid('source', shift)
}
 
sub is_tag_type_valid
{
    my $self = shift;

    return $self->_is_valid('type', shift)
}
    
sub _is_valid : PRIVATE
{
    my ($self, $attr, $value) = @_;

    return unless Finfo::Validate->validate
    (
        attr => 'valid attribute type',
        value => $attr,
        type => 'in_list',
        options => [qw/ type source /],
        err_cb => $self,
    );
 
    return unless Finfo::Validate->validate
    (
        attr => 'tag type',
        value => $value,
        type => 'defined',
        err_cb => $self,
    );
 
    my $valid_attr_method = sprintf('tag_%ss', $attr);

    $self->info_msg($value);
    
    return grep { $value eq $_ } $self->$valid_attr_method;
}

1;

#$HeadURL$
#$Id$
