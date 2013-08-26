package Genome::Site::TGI::PSEErrorLogEntry;
use strict;
use warnings;
use Genome;
use Workflow;

class Genome::Site::TGI::PSEErrorLogEntry {
    table_name => 'pse_error_log_entry',
    schema_name => 'public',
    data_source => 'Genome::DataSource::OldPostgres',
    id_generator => '-uuid',
    id_by => [
        id => { is => 'Text' },
    ],
    has => [
        entry_date         => { is => 'TIMESTAMP' },
        auto_truncate_message => { is => 'Boolean', default => '1', is_transient => 1},
    ],
    has_optional => [
        message     => { is => 'Text', len => 1000, default => '' },
        file        => { is => 'Text' },
        subroutine  => { is => 'Text' },

        pse => { is => 'GSC::PSE', id_by => 'pse_id' },
        pse_id => { is => 'NUMBER', implied_by => 'pse'},
        status  => {is => 'Text' },
        username        => { is => 'Text' },
        sudo_username   => { is => 'Text' },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless $self;

    $self->entry_date(Workflow::Time->now) unless $self->entry_date;
    $self->username(Genome::Sys->username) unless $self->username;
    $self->sudo_username(Genome::Sys->sudo_username) unless $self->sudo_username;

    $self->_auto_truncate_message if $self->auto_truncate_message;

    return $self;
}

sub _auto_truncate_message {
    my $self = shift;
    my $message_max_length = $self->class->__meta__->property('message')->data_length;

    if ($message_max_length && length($self->message) > $message_max_length) {
        $self->message(substr($self->message, 0, $message_max_length));
    }
    return 1;
}

1;
