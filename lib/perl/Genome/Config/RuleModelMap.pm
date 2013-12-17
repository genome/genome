package Genome::Config::RuleModelMap;

use strict;
use warnings;

use Genome;

class Genome::Config::RuleModelMap {
    is => 'UR::Object',
    is_transactional => 0,
    has => [
        rules => {
            is => 'Genome::Config::Rule',
            is_many => 1,
        },
        models => {
            is => 'HASH',
        },
        config => {
            is => 'Genome::Config::Profile::Item'
        }
    ]
};

sub match {
    my $self = shift;
    my $instrument_data = shift;

    die('You must supply an instrument data instance to match against!')
        unless $instrument_data;

    for my $rule ($self->rules) {
        my $result = $self->_evaluate_rule($rule, $instrument_data);
        return unless $result;
    }

    $self->config->concretize();
    return 1;
}

sub _evaluate_rule {
    my $self = shift;
    my $rule = shift;
    my $instrument_data = shift;

    my $actual_value = _evaluate_method_chain($instrument_data, @{$rule->method_chain});
    return $actual_value eq $rule->expected_value;
}

sub _evaluate_method_chain {
    return unless defined($_[0]);

    if (@_ >= 2) {
        my $obj = shift;
        my $meth = shift;
        return _evaluate_method_chain($obj->$meth, @_);
    } else {
       return shift;
    }
}

1;
