package Genome::VariantReporting::Framework::Plan::MasterPlan;

use strict;
use warnings;
use Genome;
use YAML;
use Params::Validate qw(validate_pos);
use Set::Scalar;
use Data::Compare qw(Compare);
use JSON;

my $_JSON_CODEC = new JSON->allow_nonref;

class Genome::VariantReporting::Framework::Plan::MasterPlan {
    is => 'Genome::VariantReporting::Framework::Plan::Base',
    has => [
        expert_plans => {
            is => 'Genome::VariantReporting::Framework::Plan::ExpertPlan',
            is_many => 1,
        },
        report_plans => {
            is => 'Genome::VariantReporting::Framework::Plan::ReportPlan',
            is_many => 1,
        },
        needs_translation => {
            is => 'Boolean',
        },
    ],
};

sub category {
    return 'master';
}

sub validate_translation_provider {
    my $self = shift;
    my $provider = shift;
    my @errors = $self->__translation_errors__($provider->translations);
    if (@errors) {
        $self->print_errors(@errors);
        die "Plan is incompatible with translations file.";
    }
    return;
}

sub __translation_errors__ {
    my $self = shift;
    my $provider = shift;

    my @errors;
    for my $plan ($self->expert_plans, $self->report_plans) {
        push @errors, $plan->__translation_errors__($provider);
    }
    return @errors;
}

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
    return ('experts' => [$self->expert_plans], 'reports' => [$self->report_plans]);
}

sub write_to_file {
    my $self = shift;
    my $filename = shift;

    YAML::DumpFile($filename, $self->as_hashref);
}

sub create_from_file {
    my $class = shift;
    my $file = shift;

    Genome::Sys->validate_file_for_reading($file);
    my ($hashref, undef, undef) = YAML::LoadFile($file);
    $hashref->{needs_translation} = 1 unless exists $hashref->{needs_translation};

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
    $hashref->{root}->{needs_translation} = $self->needs_translation;
    my $result = $hashref->{root};
    delete $result->{params};

    return $result;
}

sub create_from_hashref {
    my $class = shift;
    my $hashref = shift;

    # TODO make a copy of the hashref so we don't change the original via perl autovivify when we access filters (or anything else that is optional)
    # we will need to specifically initialize these optional things to empty refs
    my %params_for_create = (
        name => 'root',
        params => {},
    );
    if (exists $hashref->{needs_translation}) {
        $params_for_create{needs_translation} = $hashref->{needs_translation};
    }
    my $self = $class->SUPER::create(%params_for_create);
    my @expert_plans;
    while (my ($expert_name, $adaptor_params) = each %{$hashref->{experts}}) {
        push @expert_plans, Genome::VariantReporting::Framework::Plan::ExpertPlan->create(
            name => $expert_name,
            adaptor_params => $adaptor_params,
        );
    }
    $self->expert_plans(\@expert_plans);

    my @report_plans;
    while (my ($report_name, $report_params) = each %{$hashref->{reports}}) {
        push @report_plans, Genome::VariantReporting::Framework::Plan::ReportPlan->create_from_hashref(
            $report_name, $report_params,
        );
    }
    $self->report_plans(\@report_plans);

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

sub __plan_errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__plan_errors__;
    my $have = Set::Scalar->new(map {$_->name} $self->expert_plans);
    for my $report_plan ($self->report_plans) {
        my $needed = Set::Scalar->new($report_plan->requires_annotations);

        if (my $still_needed = $needed - $have) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => [$still_needed->members],
                desc => sprintf("Annotations required by report (%s) but not provided by any experts: (%s)",
                    $report_plan->name, join(",", $still_needed->members)),
            );
        }
    }

    return @errors;
}

sub __class_errors__ {
    return;
}

sub __object_errors__ {
    return;
}

sub translate {
    my $self = shift;
    return unless ($self->needs_translation);

    for my $plan ($self->expert_plans, $self->report_plans) {
        $plan->translate(@_);
    }
    $self->needs_translation(0);
    return;
}

1;
