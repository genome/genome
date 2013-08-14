package Genome::MiscNote;

use strict;
use warnings;

use Genome;
class Genome::MiscNote {
    type_name => 'misc note',
    table_name => 'subject.misc_note',
    id_generator => '-uuid',
    id_by => [
        id => {
            is => 'Text',
            len => 32,
        },
    ],
    has => [
        subject_class_name => { is => 'Text' },
        subject_id         => { is => 'Text' },
        header_text        => { is => 'Text' },
        subject            => { is => 'UR::Object', id_class_by => 'subject_class_name', id_by => 'subject_id' },
        editor_id          => { is => 'Text' },
        entry_date         => { is => 'DateTime' },
        entry_date_sort    => { 
            is => 'Number',
            calculate_from => ['entry_date'],
            calculate => q{ (my $sort_date = $entry_date) =~ s/[-: ]//g; return $sort_date; } # for easier XSLT sorting
        },
        auto_truncate_body_text => { is => 'Boolean', default => '1', is_transient => 1},
    ],
    has_optional => [
        body_text          => { is => 'VARCHAR2', len => 4000, default => '' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless $self;

    unless ($self->entry_date) {
        $self->entry_date($UR::Context::current->now);
    }

    unless ($self->editor_id) {
        $self->editor_id(Genome::Sys->username);
    }

    my $body_text = $self->body_text;
    my $sudo_username = Genome::Sys->sudo_username;
    if ($sudo_username) {
        $self->body_text($sudo_username . ' is running as ' . $self->editor_id . '. ' . $body_text);
    }

    $self->_auto_truncate_body_text if $self->auto_truncate_body_text;

    return $self;
}

sub _auto_truncate_body_text {
    my $self = shift;

    my $body_text_max_length = $self->class->__meta__->property('body_text')->data_length;

    my $body_text = $self->body_text;
    if ($body_text_max_length and length($body_text) > $body_text_max_length) {
        $body_text = substr($body_text, 0, $body_text_max_length);
        $self->body_text($body_text);
    }

    return 1;
}

1;
