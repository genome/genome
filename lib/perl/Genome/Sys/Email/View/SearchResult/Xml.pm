package Genome::Sys::Email::View::SearchResult::Xml;

use strict;
use warnings;

use Genome;

class Genome::Sys::Email::View::SearchResult::Xml {
    is => 'Genome::View::SearchResult::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            default => [
                'id',
                'subject',
                'blurb',
                'body',
                'list_name',
                'month',
                'message_id',
                'mail_server_path',
                'mail_list_path',
            ]
        }
    ]
};


1;
