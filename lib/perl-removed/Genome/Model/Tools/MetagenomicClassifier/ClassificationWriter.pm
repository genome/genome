package Genome::Model::Tools::MetagenomicClassifier::ClassificationWriter;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::MetagenomicClassifier::ClassificationWriter {
    has => [
        file => {
            is => 'Text',
            doc => 'File to write classifications.',
        },
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
        _io => {},
    ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $valid_formats = $self->__meta__->property_meta_for_name('format')->valid_values;
    unless ( grep { $self->format eq $_ } @$valid_formats ) {
        $self->error_message('Invalid format ('.$self->format.'). Valid formats are: '.join(', ', @$valid_formats));
        return;
    }

    if ( not $self->file ) {
        $self->error_message('No file given to classification writer');
        return;
    }

    my $fh = eval { Genome::Sys->open_file_for_writing($self->file) };
    if ( not $fh ) {
        $self->error_message('Failed to open file ('.$self->file.') for writing: '.$@);
        return;
    }
    $fh->autoflush(1);
    $self->_io($fh);

    return $self;
}

sub write {
    my ($self, $classification) = @_;

    my $format = '_string_for_'.$self->format.'_format';
    my $string = $self->$format($classification);
    return if not $string;

    return $self->_io->print($string);
}

sub _string_for_hmp_fix_ranks_format {
    my ($self, $classification) = @_;

    my $string = join(
        ';',
        $classification->{id},
        ($classification->{complemented} ? '-' : ' '),
        '',
    );
    my $prev_taxon_id;
    my $prev_conf = '';
    for my $rank (qw/ root domain phylum class order family genus /) {
        my ($taxon) = grep { $_->{rank} eq $rank } @{$classification->{taxa}};
        my $taxon_id = $classification->{$rank}->{id};
        my $conf = $classification->{$rank}->{confidence};
        unless ( $taxon_id ) {
            $taxon_id = $prev_taxon_id.'_no_'.$rank;
            $conf = $prev_conf;
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
        $classification->{id},
        ($classification->{complemented} ? '-' : ' '),
        '',
    );

    for my $rank (qw/ root domain phylum class order family genus /) {
        if ( $classification->{$rank}->{id} ) {
            $string .= $classification->{$rank}->{id}.':'.$classification->{$rank}->{confidence}.';';
        }
        else {
            $string .= ';';
        }
        
    }

    if ( $classification->{species}->{id} ) { # old behavior
        $string .= $classification->{species}->{id}.':'.$classification->{species}->{confidence}.';';
    }

    $string .= "\n";

    return $string;
}

1;

