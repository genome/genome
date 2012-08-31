package Genome::Utility::MetagenomicClassifier::Rdp::Reader;

use strict;
use warnings;

use Genome;

use Bio::Taxon;
use Data::Dumper 'Dumper';

class Genome::Utility::MetagenomicClassifier::Rdp::Reader {
    is => 'Genome::Utility::IO::Reader',
};

sub next {
    my $self = shift;

    my $line = $self->getline
        or return;

    chomp $line;
    my @tokens = split(/[:;]\s*/, $line);

    my %classification;
    $classification{name} = shift @tokens;
    $classification{classifier} = 'rdp';

    my $complemented = shift @tokens;
    $classification{complemented} = ( $complemented eq '-' ) ? 1 : 0;

    my $root = shift @tokens;
    unless ( $root eq 'Root' ) {
        $self->error_message("Mal formed line for rdp output: $line");
        die;
    }

    #< TAXON >#
    # current rdp files do not output the kingdom, since it is essentially the same as domain
    my $taxon = Bio::Taxon->new(
        '-id' => $root,
        '-rank' => 'root',
    );
    $taxon->add_tag_value('confidence', shift(@tokens));
    $classification{taxon} = $taxon; # this is the root taxox
    my @ranks = grep { $_ ne 'kingdom' } Genome::Utility::MetagenomicClassifier->taxonomic_ranks;
    my $prev_taxon = $classification{taxon}; # set as prev taxon to add descendants
    while ( @tokens ) {
        my $id = shift @tokens;
        my $confidence = shift @tokens;
        $confidence = $confidence / 100 if $confidence =~ s/\%//; # convert confidence % to decimal
        my $taxon = Bio::Taxon->new( # now taxon 
            '-id' => $id,
            '-rank' => shift(@ranks),
        );
        $taxon->add_tag_value('confidence', $confidence);
        $prev_taxon->add_Descendent($taxon);
        $prev_taxon = $taxon;
    }

    return Genome::Utility::MetagenomicClassifier::SequenceClassification->new(%classification);
}

1;

=pod

=head1 Tests

=head1 Disclaimer

 Copyright (C) 2006 Washington University Genome Sequencing Center

 This script is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 License for more details.

=head1 Author(s)

 Eddie Belter <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$

