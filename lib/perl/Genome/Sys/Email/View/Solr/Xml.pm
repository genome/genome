package Genome::Sys::Email::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::Sys::Email::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'mail'
        },
        display_type => {
            is => 'Text',
            default => 'Email',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_sys_email_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate =>  q{
                return $subject->message_url();
            }
        },
        display_label1 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => q{
                return $subject->list_name() . ' archive';
            }
        },
        display_url1 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => q{
                return $subject->list_archive_url();
            }
        },
        display_label2 => {
            is  => 'Text',
            default => 'All lists',
        },
        display_url2 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => q{
                return $subject->mail_list_path();
            }
        },
        display_label3 => {
            is  => 'Text',
            default => '',
        },
        display_url3 => {
            is  => 'Text',
            default => '',
        },
        default_aspects => {
            is => 'ARRAY',
            default => [
                {
                    name => 'subject',
                    position => 'title',
                },
                {
                    name => 'body',
                    position => 'content',
                },
                {
                    name => 'subject',
                    position => 'display_title',
                }
            ],
        }
    ]
};


sub _reconstitute_from_doc {
    my $class = shift;
    my $solr_doc = shift;
    
    unless($solr_doc->isa('WebService::Solr::Document')) {
        $class->error_message('content_doc must be a WebService::Solr::Document');
        return;
    }
    
    my $subject_class_name = $solr_doc->value_for('class');
    my $subject_id = $solr_doc->value_for('object_id');
    unless($subject_class_name eq 'Genome::Sys::Email') {
        $class->error_message('content_doc for this view must point to a Genome::Sys::Email');
        return;
    }
    
    my $mail = $subject_class_name->get($subject_id);
    unless($mail) {
        $class->error_message('Could not get Genome::Sys::Email object.');
    }
    
    unless($mail->is_initialized) {
        $mail->initialize($solr_doc);
    }
    
    return $class->SUPER::_reconstitute_from_doc($solr_doc, @_);
}

1;
