package Genome::Model::Tools::Pcap::Tag::Oligo;

use strict;
use warnings;

use base qw(Genome::Model::Tools::Pcap::Tag);

Genome::Model::Tools::Pcap::Tag::Oligo->mk_accessors
(qw/
    oligo_name
    oligo_num
    oligo_seq
    oligo_temp
    oligo_templates
    /);

our $VERSION = 0.01;

1;

=pod

=head1 Name

 Genome::Model::Tools::Pcap::Tag::Oligo

  > Inherits from Genome::Model::Tools::Pcap::Tag;
 
=head2 Contig Tag Format

 CT{
 Contig24 oligo consed 606 621 050427:142133
 M_BB0392D19.29 ccctgagcgagcagga 60 U
 L25990P6000A5 L25990P6000D4
 }

=head1 Methods

 Getters Only!! Add to TagParse.pm and Ace.pm to make these Setters to
  produce the correct output in an ace file

=head2 oligo_name

=head2 oligo_num

=head2 oligo_seq

=head2 oligo_temp

=head2 oligo_templates

=head1 Disclaimer

 Copyright (C) 2006 Washington University Genome Sequencing Center

 This module is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 License for more details.

=head1 Author(s)

 Lynn Carmicheal <lcarmich@watson.wustl.edu>
 Jon Schindler <jschindl@watson.wustl.edu>
 Eddie Belter <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
