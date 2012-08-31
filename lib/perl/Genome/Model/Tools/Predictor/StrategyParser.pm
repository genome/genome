package Genome::Model::Tools::Predictor::StrategyParser;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Predictor::StrategyParser {
    has => [
        strategy => {
            is => 'Text',
            doc => 'Strategy to be parsed',
        },
    ],
    has_transient => [
        _parsed_result => {
            doc => 'Data structure resulting from parsing of strategy',
        },
    ],
};

sub predictor_names {
    my $self = shift;
    my %results = %{$self->get_result};
    my @names = sort keys %results;
    return @names;
}

sub predictor_classes {
    my $self = shift;
    my %results = %{$self->get_result};
    my @classes;
    for my $predictor (sort keys %results) {
        push @classes, $results{$predictor}{class};
    }
    return @classes;
}

sub get_result {
    my $self = shift;
    return $self->parse;
}

sub parse {
    my $self = shift;
    return $self->_parsed_result if $self->_parsed_result;

    my @phrases = split(/\s+union\s+/, $self->strategy);
    my %predictor_data;
    for my $phrase (@phrases) {
        # Each phrase should look something like: predictor_name [params] [version]
        # Params and version are optional, but the [ ] are not.
        my ($name, $params, $version) = $phrase =~ /^(.+)\s+\[(.*)\]\s+\[(.*)\]$/; 
        $params ||= '';
        $version ||= '';

        my $class = join('::', $self->base_predictor_class, Genome::Utility::Text::string_to_camel_case($name));
        $predictor_data{$name}{class} = $class;
        $predictor_data{$name}{parameters} = $params;
        $predictor_data{$name}{version} = $version;
    }

    $self->_parsed_result(\%predictor_data);
    return $self->_parsed_result;
}

sub base_predictor_class {
    return 'Genome::Model::Tools::Predictor';
}

1;

