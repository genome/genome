package Genome::Model::Tools::Pcap::TagParser;

use strict;
use warnings;

use Genome::Model::Tools::Pcap::Tag;
use Genome::Model::Tools::Pcap::Tag::AutoFinishExp;
use Genome::Model::Tools::Pcap::Tag::Oligo;

our $VERSION = 0.01;

sub new
{
    return bless {}, shift;
}

sub parse
{
    my ($self, $io) = @_;
    
    my %tag_atts; # tag attributes

    $io->getline;
    my $line = $io->getline;
    
    chomp $line;
    
	$line =~ s/^\s*// if $line =~ /\w/;     

    @tag_atts{'parent', 'type', 'source', 'start', 'stop', 'date', 'no_trans'} = split ' ', $line;

    my $comment_flag = 0;
    while (my $line = $io->getline) 
    {
        last if $line =~ /^\s*}/;

		$line =~ s/^\s*// if $line =~ /\w/;

        $tag_atts{text} .= $line;
        
        if ( $line =~ /COMMENT{/ )
        {
            $comment_flag = 1;
            next;
        }
        elsif ( $line =~ /C\}/ )
        {
            $comment_flag = 0;
        }

        $tag_atts{comment} .= $line if $comment_flag;
    }

    $tag_atts{scope} = 'ACE',
    
    return $self->_create_tag(\%tag_atts);
}

sub _create_tag
{
    my ($self, $atts) = @_;
    
    my $class = 'Genome::Model::Tools::Pcap::Tag';

    if ( $atts->{type} eq 'oligo' )
    {
        $class .= '::Oligo';
        $self->_parse_oligo_text($atts);
    }
    elsif ($atts->{type} eq 'autoFinishExp')
    {
        $class .= '::AutoFinishExp';
        $self->_parse_autoFinishExp_text($atts);
    }

    return $class->new(%$atts);
}

sub _parse_autoFinishExp_text
{
    my ($self, $atts) = @_;
    
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

    #my (@lines) = split /\n/, $atts->{text};

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

    my @lines = split /\n/, $atts->{text};
    
    if ( $lines[2] =~ /(\d+)\s+(\d+)\s+(\d+)/ )
    {
        $atts->{num1} = "$1";
        $atts->{num2} = "$2";
        $atts->{num3} = "$3";
    }

    if ( $lines[3] =~ /\s?(\w+)\s+(\w+)\n/ )
    {
        $atts->{chem} = "$1";
        $atts->{primer_type} = "$2";
    }

    foreach my $key ( keys %$patterns )
    {
        my $pattern = $patterns->{$key};
        
        $atts->{$key} = $1 if $atts->{text} =~ /$pattern/;
    }

    return 1;
}

sub _parse_oligo_text
{
    my ($self, $atts) = @_;
    
    # CT{
    # Contig24 oligo consed 606 621 050427:142133
    # M_BB0392D19.29 ccctgagcgagcagga 60 U
    # L25990P6000A5 L25990P6000D4
    # }

    my (@lines) = split "\n", $atts->{text};

    my ($name, $seq, $temp, $u_or_c) = split /\s+/, $lines[0];

    $atts->{oligo_name} = $name;
    $atts->{oligo_seq} = $seq;
    $atts->{oligo_temp} = $temp;
    $atts->{oligo_u_or_c} = $u_or_c;

    if ($name =~ /\.(\d+)$/)
    {
        $atts->{oligo_num} = "$1";
    }
    
    $atts->{oligo_templates} = [ split /\s+/, $lines[1] ] if defined $lines[1] 
        and $lines[1] !~ /^\s*\n?$/;
    
    return 1;
}

sub serialize
{
    my ($self, $tag) = @_;

    die "not ready";
    
    return;
}

1;

=pod

=head1 Name

 Genome::Model::Tools::Pcap::TagParser

  > Parses tags from an acefice, or formats a tag for an acefile (not ready)
  
=head1 Synopsis

 my $tp = Genome::Model::Tools::Pcap::TagParser->new();

 $tp->parse($io);

=head1 Methods

=head2 parse

 Takes an io stream (inherits from IO::Handle) and parses the tag info.  The IO 
 stream must be where a tag begins.  Returns a Genome::Model::Tools::Pcap::Tag.

   # CT{
   # Contig24 oligo consed 606 621 050427:142133
   # M_BB0392D19.29 ccctgagcgagcagga 60 U
   # L25990P6000A5 L25990P6000D4
   # }

=head2 serialize

 NOT READY!
 formats a tag to be written to an acefile

=head1 Disclaimer

 Copyright (C) 2006 Washington University Genome Sequencing Center

 This module is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 License for more details.

=head1 Author(s)

 Eddie Belter <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
