package Genome::Notable;
use strict;
use warnings;
use Genome;

class Genome::Notable {
    is_abstract => 1,
    has => [
        notes => {
            is => 'Genome::MiscNote',
            is_many => 1,
            reverse_as => 'subject',
            specify_by => 'header_text',
            where => ['-order_by' => 'entry_date'],
        },
    ],
};

1;

