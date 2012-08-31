# FIXME ebelter
# Needed?? nothig in Genome uses this.
# If needed: convert to G:U:IO:Writer
#
package Genome::Utility::PSL::Writer;

use strict;
use warnings;

use Genome;

my @header_fields = qw(matches misMatches repMatches nCount qNumInsert qBaseInsert tNumInsert tBaseInsert strand qName qSize qStart qEnd tName tSize tStart tEnd blockCount blockSizes qStarts tStarts qSeq tSeq);

class Genome::Utility::PSL::Writer {
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
