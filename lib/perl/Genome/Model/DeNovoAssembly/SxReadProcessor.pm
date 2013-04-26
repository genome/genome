package Genome::Model::DeNovoAssembly::SxReadProcessor;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
use Parse::RecDescent;

class Genome::Model::DeNovoAssembly::SxReadProcessor { 
    has => [
        read_processor => {
            is => 'Text',
            doc => '',
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            doc => 'Instrument data to be processed with SX.',
        },
    ],
    has_transient_optional => [
        _default_read_processing => { is => 'HASH', },
        _read_processings => { is => 'ARRAY', },
    ],
};

our $parser; # singleton
sub parser {
    my $self = shift;

    return $parser if $parser;

    my $grammar = q{

        read_processor: specification(s)
        { $item[1]; }

        specification: default "(" processor ", coverage " coverage ")"
        { $return = { condition => $item{default}, processor => $item{processor}, coverage => $item{coverage}, } }
        | default "(" processor ")"
        { $return = { condition => $item{default}, processor => $item{processor}, } }
        | condition(s) "(" processor ", coverage " coverage ")"
        { $return = { condition => $item{'condition(s)'}, processor => $item{processor}, coverage => $item{coverage}, } }
        | condition(s) "(" processor ")"
        { $return = { condition => $item{'condition(s)'}, processor => $item{processor}, } }

        default: /DEFAULT/
        { $item[1]; }

        condition: /([\w\d\.\=\<\>\*]+)/
        { $item[1]; }

        processor: /([\w\d\-\.\/\_\:\= ]+)/
        { $item[1]; }

        coverage: /(\d+)X/
        { $item[1] =~ s/X//; $item[1]; }
    };

    $parser = Parse::RecDescent->new($grammar);
    if ( not $parser ) {
        $self->error_message("Failed to create parser from sx read processor grammar!");
        return;
    }

    return $parser
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;

    my $parser = $self->parser;
    my $read_processor = $self->read_processor;
    my $parsed_read_processings;
    if ( $read_processor =~ /^DEF/ ) {
        $parsed_read_processings = $parser->read_processor($read_processor);
        if ( not $parsed_read_processings ) {
            return UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ read_processor /],
                desc => 'Failed to parse read processor specification!',
            );
        }
        $parsed_read_processings = [ $parsed_read_processings ] if not ref $parsed_read_processings eq 'ARRAY';
    }
    else {
        $read_processor =~ s/^gmt sx //;
        $parsed_read_processings = [{
                condition => 'DEFAULT',
                processor => $read_processor,
            },];
    }
    #print Dumper($parsed_read_processors);

    my $default_processing_count = grep { $_->{condition} eq 'DEFAULT' } @$parsed_read_processings;
    if ( $default_processing_count != 1 ) {
        return UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ read_processor /],
            desc => ( $default_processing_count ? 'Multiple' : 'No' ).' DEFAULT read processor(s) in specification!',
        );
    }
    $self->_default_read_processing(shift @$parsed_read_processings); # always the first
    if ( not $self->_default_read_processing->{condition} eq 'DEFAULT' ) {
        return UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ read_processor /],
            desc => 'Could not find DEFAULT read processor from specification!',
        );
    }
    $self->_read_processings($parsed_read_processings);

    return;
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my @errors = $self->__errors__;
    if ( @errors ) {
        for my $error ( @errors ) {
            $self->error_message( $error->__display_name__ );
        }
        return;
    }

    return $self;
}

sub determine_processing_for_instrument_data {
    my ($self, $instrument_data) = @_;

    my $processing;

    return $processing;
}

1;

