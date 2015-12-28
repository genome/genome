package Genome::Role::Notable;
use strict;
use warnings;
use Genome;
use UR::Role;

role Genome::Role::Notable {
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

