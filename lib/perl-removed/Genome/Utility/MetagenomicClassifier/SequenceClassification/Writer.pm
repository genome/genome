package Genome::Utility::MetagenomicClassifier::SequenceClassification::Writer;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Utility::MetagenomicClassifier::SequenceClassification::Writer {
    is => 'Genome::Utility::IO::Writer',
    has => [
        format => {
            is => 'Text',
            is_optional => 1,
            valid_values => [qw/ hmp_fix_ranks hmp_all_ranks/],
            default_value => 'hmp_fix_ranks',
            doc => <<DOC,
The format of the output. Default is hmp_fix_ranks.
  hmp_fix_ranks => name;complemented('-' or ' ');taxon:confidence;[taxon:confidence;]
    prints only root, domain, phylum, class, order, family, genus from classification
  hmp_all_ranks => name;complemented('-' or ' ');taxon:confidence;[taxon:confidence;]
    prints ALL taxa in classification
DOC
        },
    ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    # Format
    my $valid_formats = $self->__meta__->property_meta_for_name('format')->valid_values;
    unless ( grep { $self->format eq $_ } @$valid_formats ) {
        $self->error_message('Invalid format ('.$self->format.'). Valid formats are: '.join(', ', @$valid_formats));
        $self->delete;
        return;
    }

    return $self;
}

sub write_one {
    my ($self, $classification) = @_;

    my $format = '_string_for_'.$self->format.'_format';
    my $string = $self->$format($classification) # any error should be reported in format method
        or return;
    $self->output->print($string);
    
    return 1;
}

sub _string_for_hmp_fix_ranks_format {
    my ($self, $classification) = @_;

    my $string = join(
        ';',
        $classification->get_name,
        ($classification->is_complemented ? '-' : ' '),
        '',
    );
    my $prev_taxon_id;
    my $prev_conf;
    for my $rank (qw/ root domain phylum class order family genus /) {
        my ($taxon_id, $conf) = $classification->get_taxon_name_and_confidence_for_rank($rank);
        unless ( $taxon_id ) {
            $taxon_id = $prev_taxon_id.'_no_'.$rank;
            $conf = $prev_conf || '';
        }
        else {
            $prev_taxon_id = $taxon_id;
            $prev_conf = $conf;
        }
        $string .= $taxon_id.':'.$conf.';';
    }
    $string .= "\n";

    return $string;
}

sub _string_for_hmp_all_ranks_format {
    my ($self, $classification) = @_;

    my $string = join(
        ';',
        $classification->get_name,
        ($classification->is_complemented ? '-' : ' '),
        '',
    );
    for my $taxon ( $classification->get_taxa ) {
        if ( $taxon->id ) {
            $string .= $taxon->id.':'.($taxon->get_tag_values('confidence'))[0].';';
        }
        else {
            $string .= ';';
        }
        
    }
    $string .= "\n";

    return $string;
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
