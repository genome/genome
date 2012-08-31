package Finishing::Assembly::Ace::TagFactory;
{
    use strict;
    use warnings;

    use base 'Finfo::Singleton';

    use Data::Dumper;

    sub build_assembly_tag 
    {
        my ($self, $io) = @_;

        my $line = $io->getline;
        chomp $line;
        $line =~ s/^\s*// if $line =~ /\w/;;

        my %tag;
        @tag{qw/ type source date /} = split(/ /, $line);
        $tag{date} = $self->_convert_from_ace_tag_date($tag{date});
        $self->_parse_tag_text_and_comment($io, \%tag);

        return Finishing::Assembly::Ace::AssemblyTag->new(%tag);
    }

    sub build_contig_tag
    {
        my ($self, $io) = @_;

        my $line = $io->getline;
        chomp $line;
        $line =~ s/^\s*// if $line =~ /\w/;     

        my %tag;
        @tag{qw/ parent type source start stop date no_trans /} = split(/ /, $line);
        $tag{date} = $self->_convert_from_ace_tag_date($tag{date});
        $self->_parse_tag_text_and_comment($io, \%tag);

        my $class = 'Finishing::Assembly::Ace::';
        my %tag_types_and_class_additions = 
        (
            oligo => 'OligoTag',
            autoFinishExp => 'AutoFinishExpTag',
        );
        if ( my $addition = $tag_types_and_class_additions{ $tag{type} } )
        {
            $class .= $addition;
        }
        else
        {
            $class .= 'ConsensusTag';
        }

        return $class->new(%tag);
    }
    
    sub build_read_tag 
    {
        my ($self, $io) = @_;

        my $line = $io->getline;
        chomp $line;
        $line =~ s/^\s*// if $line =~ /\w/;;

        my %tag;
        @tag{qw/ parent type source start stop date /} = split(/ /, $line);
        $tag{date} = $self->_convert_from_ace_tag_date($tag{date});
        $self->_parse_tag_text_and_comment($io, \%tag);

        return Finishing::Assembly::Ace::ReadTag->new(%tag);
    }

    sub _convert_from_ace_tag_date
    {
        my ($self, $date) = @_;

        $date =~ /(\d\d)(\d\d)(\d\d):(\d\d)(\d\d)(\d\d)/;
        $date = ($1 > 80 ? "19$1" : "20$1") . "-$2-$3 $4:$5:$6";

        return $date;
    }

    sub _parse_tag_text_and_comment
    {
        my ($self, $io, $tag) = @_;

        my $attr = 'text';
        while ( my $line = $io->getline ) 
        {
            last if $line =~ /^\s*}/;
            next if $line =~ /^C}/;

            $line =~ s/^\s*// if $line =~ /\w/;

            if ( $line =~ /COMMENT{/ )
            {
                $attr = 'comment';
                next;
            }

            $tag->{$attr} .= $line;
        }

        chomp $tag->{text} if $tag->{text};
        chomp $tag->{comment} if $tag->{comment};

        if ( $tag->{type} eq 'oligo' )
        {
            my $text_parse_method = '_parse_oligo_tag_text';
            $self->$text_parse_method($tag);
        }
        elsif ( $tag->{type} eq 'autoFinishExp' )
        {
            my $text_parse_method = '_parse_autoFinishExp_tag_text';
            $self->$text_parse_method($tag);
        }

        return 1;
    }

    sub _parse_autoFinishExp_tag_text : PRIVATE
    {
        my ($self, $tag) = @_;

        #   CT{
        #   Contig29.2 autoFinishExp autofinish 119 119 060831:122829
        # 0 C
        # 1 purpose: weak
        # 2 0 915 0
        # 3 dyeTerm customPrimer
        # 4 fix cons errors: 4.69881 original cons errors: 5.64249
        # 5 original single subclone bases: 886
        # 6 primer: ggcaaatatggtgcaataaaac temp: 58 id: Trichinella_spiralis_060315.pcap.scaffold29.ace.AE.1.1
        # 7 expID_and_template: 1 TPAA-ail08c06
        #   }

        #my (@lines) = split /\n/, $attrs->{text};

        my $patterns = 
        {
            orientation => '(\w)',
            purpose => 'purpose: (.+)',
            fix_cons_errors => 'fix cons errors: (\d+\.\d+)',
            original_cons_errors => 'original cons errors: (\d+\.\d+)',
            original_single_subclone_bases => 'original single subclone bases: (\d+)',
            oligo_seq => 'primer: (\w+)',
            oligo_temp => 'temp: (\d+)',
            oligo_name => 'id: (.+)',
            exp_id_and_template => 'expID_and_template: (.+)',
        };

        my @lines = split /\n/, $tag->{text};

        if ( $lines[2] =~ /(\d+)\s+(\d+)\s+(\d+)/ )
        {
            $tag->{num1} = "$1";
            $tag->{num2} = "$2";
            $tag->{num3} = "$3";
        }

        if ( $lines[3] =~ /\s?(\w+)\s+(\w+)\n/ )
        {
            $tag->{chem} = "$1";
            $tag->{primer_type} = "$2";
        }

        foreach my $key ( keys %$patterns )
        {
            my $pattern = $patterns->{$key};

            $tag->{$key} = "$1" if $tag->{text} =~ /$pattern/;
        }

        return 1;
    }

    sub _parse_oligo_tag_text : PRIVATE
    {
        my ($self, $tag) = @_;

        # CT{
        # Contig24 oligo consed 606 621 050427:142133
        # M_BB0392D19.29 ccctgagcgagcagga 60 U
        # L25990P6000A5 L25990P6000D4
        # }

        my (@lines) = split(/\n/, $tag->{text});

        my ($name, $seq, $temp, $ori) = split /\s+/, $lines[0];

        $tag->{oligo_name} = $name;
        $tag->{oligo_seq} = $seq;
        $tag->{oligo_temp} = $temp;
        $tag->{complemented} = ( $ori eq 'C' ) ? 1 : 0;

        $tag->{oligo_templates} = [ split(/\s+/, $lines[1]) ] if $lines[1] and $lines[1] !~ /^\s*$/;

        return 1;
    }

    sub build_tag_from_reader_hashref
    {
        my ($self, $tag) = @_;

        my %types_and_classes =
        (
            assembly_tag => 'Finishing::Assembly::Ace::AssemblyTag',
            contig_tag => 'Finishing::Assembly::Ace::ConsensusTag',
            read_tag => 'Finishing::Assembly::Ace::ReadTag',
        );

        my %tag_types_and_class_additions = 
        (
            oligo => 'OligoTag',
            autoFinishExp => 'AutoFinishExpTag',
        );

        my $class = $types_and_classes{ delete $tag->{object_type} };
        if ( my $addition = $tag_types_and_class_additions{ $tag->{type} } )
        {
            $class =~ s/::[^:]+$/::$addition/;
        }

        return $class->new(%$tag);
    }

    1;
}

##################################################################

package Finishing::Assembly::Ace::Tag;
{
    use strict;
    use warnings;

    use base 'Finfo::Accessor';

    __PACKAGE__->mk_accessors(qw/ type source date text comment /);

    sub new
    {
        my ($class, %p) = @_;

        return bless \%p, $class;
    }

    1;
}

##################################################################

package Finishing::Assembly::Ace::AssemblyTag;
{
    use strict;
    use warnings;
    
    use base 'Finishing::Assembly::Ace::Tag';
    
    1;
}

##################################################################

package Finishing::Assembly::Ace::SequenceTag;
{
    use strict;
    use warnings;
    
    use base 'Finishing::Assembly::Ace::Tag';
    
    __PACKAGE__->mk_accessors(qw/ parent start stop /);
    #unpad_start unpad_stop 

    sub scope
    {
        return 'ace';
    }
    
    1;
}

##################################################################

package Finishing::Assembly::Ace::ConsensusTag;
{
    use base 'Finishing::Assembly::Ace::SequenceTag';

    __PACKAGE__->mk_accessors(qw/ no_trans /);

    1;
}

##################################################################

package Finishing::Assembly::Ace::AutoFinishExpTag;
{
    use strict;
    use warnings;
              
    use base 'Finishing::Assembly::Ace::ConsensusTag';

    __PACKAGE__->mk_accessors
    (qw/ 
        orientation num1 num2 num3 chem primer_type purpose fix_cons_errors
        original_cons_errors original_single_subclone_bases primer temp id
        exp_id_and_template oligo_name oligo_seq oligo_temp 
        /);

    1;
}

##################################################################

package Finishing::Assembly::Ace::OligoTag;
{
    use strict;
    use warnings;

    use base 'Finishing::Assembly::Ace::ConsensusTag';

    __PACKAGE__->mk_accessors(qw/ oligo_name oligo_seq oligo_temp oligo_templates complemented /);

    sub orientation{
        my $self = shift;
        $self->fatal_msg("orientation is not a setter, use complemented()") if @_;
        return $self->complemented? 'C' : 'U';
    }

    1;
}

##################################################################

package Finishing::Assembly::Ace::ReadTag;
{

    use base 'Finishing::Assembly::Ace::SequenceTag';

    1;
}

=pod

=cut

1;

#$HeadURL$
#$Id$
