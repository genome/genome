package Genome::Sys::User::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::Sys::User::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'user'
        },
        display_type => {
            is  => 'Text',
            default => 'User',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_sys_user_16',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub { return join ('?id=', '/view/genome/sys/user/status.html',$_[0]->email()); },
        },
        display_label1 => {
            is  => 'Text',
            default => 'send mail',
        },
        display_url1 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => sub { return 'mailto:' . $_[0]->email(); },
        },
        display_label2 => {
            is  => 'Text',
            default => 'wiki',
        },
        display_url2 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => sub { return $ENV{GENOME_SYS_SERVICES_WIKI_URL} . 'User:' . $_[0]->username(); },
        },
        display_label3 => {
            is  => 'Text',
        },
        display_url3 => {
            is  => 'Text',
        },
        display_title => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub { return shift->name },
        },
        default_aspects => {
            is => 'ARRAY',
            default => [
                {
                    name => 'name',
                    position => 'content',
                },
                {
                    name => 'email',
                    position => 'content',
                },
            ],
        }
    ]
};

1;
