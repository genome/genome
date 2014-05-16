package Genome::Annotation::Plan;

use strict;
use warnings;
use Genome;
use YAML;
use Params::Validate qw(validate_pos);
use Set::Scalar;
use Data::Compare qw(Compare);
use JSON;

my $_JSON_CODEC = new JSON->allow_nonref;

class Genome::Annotation::Plan {
    is => 'Genome::Annotation::Plan::Base',
    has => [
        expert_plans => {
            is => 'Genome::Annotation::Plan::ExpertPlan',
            is_many => 1,
        },
        reporter_plans => {
            is => 'Genome::Annotation::Plan::ReporterPlan',
            is_many => 1,
        },
    ],
};

sub get_plan {
    my ($self, $category, $name) = validate_pos(@_, 1, 1, 1);

    my $accessor = sprintf("%s_plans", $category);
    my @category_plans = $self->$accessor;
    for my $plan (@category_plans) {
        if ($plan->name eq $name) {
            return $plan;
        }
    }
    die "Couldn't find plan of category ($category) and name ($name)";
}

sub object {
    # return the master workflow generator
    return;
}

sub children {
    my $self = shift;
    return ('experts' => [$self->expert_plans], 'reporters' => [$self->reporter_plans]);
}

sub write_to_file {
    my $self = shift;
    my $filename = shift;

    YAML::DumpFile($filename, $self->as_hashref);
}

sub create_from_file {
    my $class = shift;
    my $file = shift;

    my ($hashref, undef, undef) = YAML::LoadFile($file);

    my $self = $class->create_from_hashref($hashref);

    my $understood_hashref = $self->as_hashref;
    unless(Compare($understood_hashref, $hashref)) {
        die $self->error_message("Problems encountered when loading the YAML. File contains invalid information for the plan. Parsed file contents: (%s) Understood file contents: (%s)",
            Data::Dumper::Dumper($hashref), Data::Dumper::Dumper($understood_hashref) );
    }

    return $self;
}

sub as_hashref {
    my $self = shift;

    my $hashref = $self->SUPER::as_hashref;
    my $result = $hashref->{root};
    delete $result->{params};

    return $result;
}

sub create_from_hashref {
    my $class = shift;
    my $hashref = shift;

    my $self = $class->SUPER::create(name => 'root',
        params => {});
    my @expert_plans;
    for my $expert_name (keys %{$hashref->{experts}}) {
        push @expert_plans, Genome::Annotation::Plan::ExpertPlan->create(
            name => $expert_name,
            params => $hashref->{experts}->{$expert_name},
        );
    }
    $self->expert_plans(\@expert_plans);

    my @reporter_plans;
    for my $reporter_name (keys %{$hashref->{reporters}}) {
        push @reporter_plans, Genome::Annotation::Plan::ReporterPlan->create_from_hashref(
            $reporter_name, $hashref->{reporters}->{$reporter_name},
        );
    }
    $self->reporter_plans(\@reporter_plans);

    return $self;
}

sub as_json {
    my $self = shift;

    return $_JSON_CODEC->canonical->encode($self->as_hashref);
}

sub create_from_json {
    my $class = shift;
    my $json = shift;

    my $hashref = $_JSON_CODEC->decode($json);
    return $class->create_from_hashref($hashref);
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__;

    my $have = Set::Scalar->new(map {$_->name} $self->expert_plans);
    my $total_needed = Set::Scalar->new();
    for my $reporter_plan ($self->reporter_plans) {
        my $needed = Set::Scalar->new();
        for my $plan ($reporter_plan->interpreter_plans, $reporter_plan->filter_plans) {
            # This may die because the module may not exist
            my $object;
            eval {
                $object = $plan->object;
            };

            if ($@) {
                push @errors, UR::Object::Tag->create(
                    type => 'error',
                    properties => [],
                    desc => $@,
                );
            } else {
                $needed->insert($object->requires_experts);
            }

        }
        $total_needed += $needed;

        if (my $still_needed = $needed - $have) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => [$still_needed->members],
                desc => sprintf("Experts required by reporter (%s) but not provided: (%s)",
                    $reporter_plan->name, join(",", $still_needed->members)),
            );
        }
    }

    if (my $not_needed = $have - $total_needed) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => [$not_needed->members],
            desc => sprintf("Experts provided by plan but not required by any reporters: (%s)",
                join(",", $not_needed->members)),
        );
    }

    return @errors;
}

1;
