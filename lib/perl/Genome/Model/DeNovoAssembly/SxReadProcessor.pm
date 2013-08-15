package Genome::Model::DeNovoAssembly::SxReadProcessor;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
use Parse::RecDescent;

class Genome::Model::DeNovoAssembly::SxReadProcessor { 
    has_optional => [
        processor => {
            is => 'Text',
            doc => 'An SX command or instrument data conditions and matching SX commnds.',
        },
    ],
    has_transient_optional => [
        default_processing => { is => 'HASH', },
        additional_processings => { is => 'ARRAY', },
        _instrument_data => { is => 'ARRAY', },
        _sx_result_params => { is => 'ARRAY', },
        _merged_sx_result_params => { is => 'ARRAY', },
    ],
};

our $parser; # singleton
sub parser {
    my $self = shift;

    return $parser if $parser;

    my $grammar = q{

        evaluate_processor: specification(s)
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

        processor: /([\w\d\-\.\/\_\:\=| ]+)/
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

    # Parse processor
    my $parser = $self->parser;
    my $processor = $self->processor;
    my $parsed_processings;
    if ( not $processor ) {
        $parsed_processings = [{
                condition => 'DEFAULT',
                processor => '',
            },];
    }
    elsif ( $processor =~ /^DEF/ ) {
        $parsed_processings = $parser->evaluate_processor($processor);
        if ( not $parsed_processings ) {
            return UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ processor /],
                desc => 'Failed to parse read processor specification! '.$processor,
            );
        }
        $parsed_processings = [ $parsed_processings ] if not ref $parsed_processings eq 'ARRAY';
    }
    else {
        $processor =~ s/^gmt sx //;
        $parsed_processings = [{
                condition => 'DEFAULT',
                processor => $processor,
            },];
    }

    # Check default processor, set processors
    my $default_processing_count = grep { $_->{condition} eq 'DEFAULT' } @$parsed_processings;
    if ( $default_processing_count != 1 ) {
        return UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ processor /],
            desc => ( $default_processing_count ? 'Multiple' : 'No' ).' DEFAULT read processor(s) in specification! '.$processor,
        );
    }
    $self->default_processing(shift @$parsed_processings); # always the first
    if ( not $self->default_processing->{condition} eq 'DEFAULT' ) {
        return UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ processor /],
            desc => 'Could not find DEFAULT read processor from specification! '.$processor,
        );
    }
    $self->additional_processings($parsed_processings);

    return if not $processor; # do not validate SX command

    # Validate SX commands
    for my $processing ( $self->default_processing, @{$self->additional_processings} ) {
        $processing->{processor}= $self->default_processing->{processor} if $processing->{processor} eq 'DEFAULT';
        my @processor_parts = split(/\s+\|\s+/, $processing->{processor});
        if ( not @processor_parts ) {
            return UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ processor /],
                desc => "Could not find read processors in string: ".$processing->{processor},
            );
        }

        for my $processor_part ( @processor_parts ) {
            my $processor_is_ok = eval{ Genome::Model::Tools::Sx::Validate->validate_command('gmt sx '.$processor_part); };
            if ( not $processor_is_ok ) {
                return UR::Object::Tag->create(
                    type => 'invalid',
                    properties => [qw/ processor /],
                    desc => "Cannot validate processor part ($processor_part)! See above error(s)",
                );
            }
        }
    }

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

sub determine_processing {
    my ($self, @instrument_data) = @_;

    Carp::confess('No instrument data given to determine_processing!') if not @instrument_data;

    # Determine processing for instrument data
    my %processings_and_instrument_data;
    for my $instrument_data ( @instrument_data ) {
        my $processing = $self->determine_processing_for_instrument_data($instrument_data);
        if ( not $processings_and_instrument_data{ $processing->{condition} } ) {
            $processings_and_instrument_data{ $processing->{condition} } = $processing;
        }
        push @{$processings_and_instrument_data{ $processing->{condition} }->{instrument_data}}, $instrument_data;
    }

    # Determine SX results and merged SX results
    my (%sx_result_params, @merged_sx_result_params);
    for my $processing ( values %processings_and_instrument_data ) {
        for my $instrument_data ( @{$processing->{instrument_data}} ) {
            $sx_result_params{ $instrument_data->id } =  {
                instrument_data_id => $instrument_data->id,
                read_processor => $processing->{processor},
                output_file_count => ( $instrument_data->is_paired_end ? 2 : 1 ),
                output_file_type => 'sanger',
                test_name => ($ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef),
            };
        }

        if ( not defined $processing->{coverage} ) {
            next;
        }

        my %output_file_counts = map { my $count = $_->is_paired_end ? 2 : 1; $count => 1; } @{$processing->{instrument_data}};
        my @output_file_counts = keys %output_file_counts;
        if ( @output_file_counts > 2 ) {
            $self->error_message("Instrument data for merged result must be all paired or all unpaired!");
            return;
        }

        push @merged_sx_result_params, {
            instrument_data_id => [ map { $_->id } @{$processing->{instrument_data}} ],
            read_processor => $processing->{processor},
            coverage => $processing->{coverage},
            output_file_count => $output_file_counts[0],
            output_file_type => 'sanger',
            test_name => ($ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef),
        };
    }

    $self->_instrument_data(\@instrument_data);
    $self->_sx_result_params(\%sx_result_params);
    $self->_merged_sx_result_params(\@merged_sx_result_params);

    return 1;
}

sub sx_result_params_for_instrument_data {
    my ($self, $instrument_data) = @_;

    Carp::confess('No instrument data given to sx_result_params_for_instrument_data!') if not $instrument_data;
    Carp::confess('Need to call "determine_processing" with all instrument data before before calling sx_result_params_for_instrument_data.') if not $self->_instrument_data;

    my $sx_result_params = $self->_sx_result_params;
    if ( not $sx_result_params ) {
        $self->error_message('Need to set sx params with all instrument data by using set_sx_result_params_for_instrument_data" before getting sx result params for instrument data.');
        return;
    }

    return $sx_result_params->{ $instrument_data->id };
}

sub merged_sx_result_params_for_instrument_data {
    my ($self, $instrument_data) = @_;

    Carp::confess('No instrument data given to merged_sx_result_params_for_instrument_data!') if not $instrument_data;
    Carp::confess('Need to call "determine_processing" with all instrument data before before calling merged_sx_result_params_for_instrument_data.') if not $self->_instrument_data;

    for my $merged_sx_result_params ( @{$self->_merged_sx_result_params} ) {
        next if not grep { $instrument_data->id eq $_ } @{$merged_sx_result_params->{instrument_data_id}};
        return $merged_sx_result_params;
    }

    return;
}

sub final_sx_result_params {
    my $self = shift;

    Carp::confess('Need to call "determine_processing" with all instrument data before before calling sx_result_params_for_instrument_data.') if not $self->_instrument_data;

    my %final_sx_result_params;
    for my $instrument_data ( @{$self->_instrument_data} ) {
        my $sx_result_params = $self->merged_sx_result_params_for_instrument_data($instrument_data);
        $sx_result_params = $self->sx_result_params_for_instrument_data($instrument_data) if not $sx_result_params;
        my $id = ( ref $sx_result_params->{instrument_data_id} ) 
        ? join(' ', sort { $a cmp $b } @{$sx_result_params->{instrument_data_id}})
        : $sx_result_params->{instrument_data_id};
        $final_sx_result_params{$id} = $sx_result_params;
    }

    return map { $final_sx_result_params{$_} } sort keys %final_sx_result_params;
}

sub determine_processing_for_instrument_data {
    my ($self, $instrument_data) = @_;

    Carp::confess('No instrument data given to determine_processing_for_instrument_data!') if not $instrument_data;

    my @matched_processings;
    if ( @{$self->additional_processings} ) {
        for my $processing ( @{$self->additional_processings} ) {
            push @matched_processings, $processing if $self->_does_instrument_data_match_condition($instrument_data, @{$processing->{condition}});
        }
    }

    # Bad - found more than one
    if ( @matched_processings > 1 ) {
        $self->error_message(
            'Found multiple processings for instrument data! '.$instrument_data->id."\n".Dumper(\@matched_processings)
        );
        return;
    }

    # Use the one we found, or if none, the default
    my %processing = ( @matched_processings == 1  ? %{$matched_processings[0]} : %{$self->default_processing} );
    return \%processing;
}

sub _does_instrument_data_match_condition {
    my ($self, $instrument_data, @conditions) = @_;

    Carp::confess('No instrument data given to _does_instrument_data_match_condition!') if not $instrument_data;
    Carp::confess('No condition given to _does_instrument_data_match_condition!') if not @conditions;

    my $processing;
    for my $word ( @conditions ) {
        for my $property (qw/ original_est_fragment_size read_length /) {
            if ( $word eq $property ) {
                $word = $instrument_data->$word;
            }
        }
    }

    my $full_condition = join(" ", @conditions);
    return eval $full_condition;
}

1;

