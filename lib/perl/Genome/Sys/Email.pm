package Genome::Sys::Email;

use strict;
use warnings;

use Genome;
use Email::Simple;

use LWP::UserAgent;

class Genome::Sys::Email {
    is => ['UR::Value','Genome::Searchable'],
    doc => 'Represents a mailing list message',
    has => [
        list_name => {
            is => 'Text',
            calculate_from => ['id'],
            calculate => q{ ($self->_parse_id($id))[0]; },
            doc => 'The name of the mailing list to which the message was posted',
        },
        month  => {
            is => 'Text',
            calculate_from => ['id'],
            calculate => q{ ($self->_parse_id($id))[1]; },
            doc => 'The month of the archives of the mailing list in which the message can be found',
        },
        message_id => {
            is => 'Text',
            calculate_from => ['id'],
            calculate => q{ ($self->_parse_id($id))[2]; },
            doc => 'The id of the message in the mailing list system',
        },    
        body => {
            is => 'Text',
            doc => 'The body of the message',
            is_transient => 1,
        },
        subject => {
            is => 'Text',
            doc => 'The subject of the message',
            is_transient => 1,
        },
        mail_server_path => {
            is => 'Text',
            calculate => q{ 'http://gscsmtp.wustl.edu/pipermail' },
        },
        mail_list_path => {
            is => 'Text',
            calculate => q{ 'http://gscsmtp.wustl.edu/cgi-bin/mailman/listinfo' },
        },
    ],
};

__PACKAGE__->add_observer(
    aspect   => 'load',
    callback => sub {
        my ($self) = @_;

        $self->initialize();
    }
);


sub get {
    my $class = shift;
    my @params = @_;
    
    my $id; #want to verify the ID matches our expected format
    if(scalar(@params) eq 1) {
        $id = $params[0];
    } else {
       my %params = @params;
       if(exists $params{id}) {
           $id = $params{id};
       }
    }
    
    unless($class->_parse_id($id)) {
        return;
    }
    
    return $class->SUPER::get(@_);
}


sub initialize {
    my $self = shift;
    my $source = shift;
    
    unless($source) {
        return $self->_initialize_from_server;
    }
    
    if (ref $source eq 'Email::Simple') {
        my $source_id = $source->header('X-Genome-Search-ID');
        unless($source_id eq $self->id) {
            Carp::confess('Source object identifier ' . $source_id . ' does not appear to match this identifier (' . $self->id . ').');
        }

        my $body = $source->body();
        $body =~ s/[\x01-\x08\x0B-\x0C\x0E-\x1F\x7F-\x84\x86-\x9F]//go ;
        $self->body($body);

        my $subject = $source->header('Subject');
        $subject =~ s/[\x01-\x08\x0B-\x0C\x0E-\x1F\x7F-\x84\x86-\x9F]//go ;
        $self->subject($subject);
        
    } elsif (ref $source eq 'WebService::Solr::Document') {
        my $source_id = $source->value_for('object_id');
        unless($source_id eq $self->id) {
            Carp::confess('Source object identifier ' . $source_id . ' does not appear to match this identifier (' . $self->id . ').');
        }
       
        my $body = $source->value_for('content');
        $body =~ s/[\x01-\x08\x0B-\x0C\x0E-\x1F\x7F-\x84\x86-\x9F]//go ;
        $self->body($body);

        my $subject = $source->value_for('title'); 
        $subject =~ s/[\x01-\x08\x0B-\x0C\x0E-\x1F\x7F-\x84\x86-\x9F]//go ;
        $self->subject($subject);

    } else {
        Carp::confess('Invalid source object ' . $source . ' passed to initialize().');
    }
    
    return 1;
}

sub _initialize_from_server {
    my $self = shift;
    
    my $ua = LWP::UserAgent->new();
    
    #Just scrape http://gscsmtp.wustl.edu/pipermail/[list]/[year-month]/[message_id].html to extract subject and body 
    
    my $response = $ua->get($self->message_url());
    
    unless ($response->is_success) {
        Carp::confess('Failed to load message from the mail server.');
    }
    
    my $content = $response->content;
    
    my ($subject) = $content =~ m!<TITLE>(.*)</TITLE>!msi;
    my ($body) = $content =~ m/<!--beginarticle-->(.*)<!--endarticle-->/msi;
    
    chomp $subject;
    chomp $body;
    
    $subject =~ s!\n!!g; #Remove newlines
    $body =~ s!<.*?>!!g; #Remove HTML tags
    $body =~ s!^\n!!gms; #Remove initial newlines
    $body =~ s!\n$!!gms; #Remove final newlines
    
    $self->subject($subject);
    $self->body($body);
    
    return 1;
}

sub _parse_id {
    my $class = shift;
    my $id = shift;
    
    my ($list_name, $year_month, $message_id) = split(/\//,$id);
    
    unless($list_name and $year_month and $message_id) {
        return;
    }
    
    return wantarray?
        ($list_name, $year_month, $message_id) :
        $message_id;
}

sub blurb {
    my $self = shift;
    my ($query) = @_; #optionally try to highlight a word in the text
    
    my $text = $self->body;

    $text =~ s/\s{2}/ /g;
    $text =~ s/-------- Original Message --------(.|\n)*//g;
    $text =~ s/\n/ /g;

    my $summarystart = 0;

    if($query) {
        #find summary region around query in result
        my $querypos = index($query, $text);

        if($querypos -75 > $summarystart) {
          $summarystart = $querypos - 75;
        }


    }

    my $summary = substr($text,$summarystart,150);

    if (length($text) > length($summary)) {
        $summary .= ' ...';
    }

    if ($summarystart > 0) {
        $summary = '... ' . $summary;
    }

    return $summary;
}

sub message_url {
    my ($self) = @_;

    my $url = join('/', $self->list_archive_url, $self->message_id . '.html');
    return $url; 
}


sub list_archive_url {
    my ($self) = @_;

    my $url = join('/', $self->mail_server_path, $self->list_name, $self->month);
    return $url; 
}




1;

