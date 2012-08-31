package Genome::Utility::SeqcleanReport::Writer;

#REVIEW fdu 11/20/2009
#This can be changed to use G::U::IO::Writer as base class
#and place it under new namespace G::U::IO::SeqcleanReport::Writer

use strict;
use warnings;

use Genome;

my @header_fields = qw(accession pc_undetermined start end length trash_code comments);

class Genome::Utility::SeqcleanReport::Writer {
    is => 'UR::Object',
    has => [
            file => {
                     is => 'string',
                 },
            header => {
                       is => 'Boolean',
                       default_value => '0',
                   }
            ],
    has_optional => [
                     _file_handle => {
                                      is => 'IO::File',
                                  },
                 ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    my $fh = IO::File->new($self->file,'w');
    unless ($fh) {
        $self->error_message('Failed to create output file handle');
        return;
    }
    $self->_file_handle($fh);

    if ($self->header) {
        print $self->_file_handle, join("\t",@header_fields);
    }
    return $self;
}

sub write_record {
    my $self = shift;
    my $record = shift;
    for (my $i = 0; $i < scalar(@header_fields); $i++) {
        my $header_field = $header_fields[$i];
        unless (defined $$record{$header_field}) {
            die "Value not found for field '$header_field' in record";
        }
        $self->_file_handle->print($$record{$header_field});
        if ($i < (scalar(@header_fields) - 1)) {
            $self->_file_handle->print("\t");
        }
    }
    $self->_file_handle->print("\n");
    return 1;
}

sub close {
    my $self = shift;
    $self->_file_handle->close;
}



1;
