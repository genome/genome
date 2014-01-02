package Genome::Config::Profile;

use strict;
use warnings;

use Genome;
use UR::Util;

use feature qw(switch);

class Genome::Config::Profile {
    is => 'UR::Object',
    is_transactional => 0,
    has => [
        config_rule_maps => {
            is_many => 1,
            is => 'Genome::Config::RuleModelMap'
        },
    ]
};


sub create_from_analysis_project {
    my $class = shift;
    my $analysis_project = shift;

    my @config_rule_maps = map {
        Genome::Config::Translator->get_rule_model_map_from_config($_)
    } $analysis_project->config_items;

    return $class->create(
        config_rule_maps => \@config_rule_maps
    );
}

sub get_config {
    my $self = shift;
    my $inst_data = shift;

    my @model_hashes = map { $_->models }
        grep { $_->match($inst_data) }
        $self->config_rule_maps;

    return $self->_merge_model_hashes(@model_hashes)
}

sub _merge_model_hashes {
    my $self = shift;
    my @source_hashes = @_;

    my $destination_hash = {};

    for my $source (@source_hashes) {
        for my $key (keys %$source) {
            $destination_hash->{$key} ||= [];
            given(ref $source->{$key}) {
                when('ARRAY') { push @{$destination_hash->{$key}}, @{$source->{$key}}; }
                when('HASH')  { push @{$destination_hash->{$key}}, $source->{$key}; }
                default       { die $self->error_message('unexpected value for key %s', $key); }
            }
        }
    }

    return UR::Util::deep_copy($destination_hash);
}

1;
