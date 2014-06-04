package Genome::InstrumentData::Imported::Reimport;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Imported::Reimport{ 
};

sub attribute_label_for_reimported_from { 'reimported_from' }

sub attribute_labels_to_ignore_when_reimporting {
    (qw/ bam_path genotype_file genotype_file_name ignored import_date import_format user_name /);
}

1;

