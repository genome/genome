package Genome::Model::Command::Define::SomaticValidation;

use strict;
use warnings;

use Genome;

use Cwd;
use File::Basename;

class Genome::Model::Command::Define::SomaticValidation {
    is => 'Command::V2',
    has_optional_input => [
        name => {
            is => 'Text',
            doc => 'A name for the model',
        },
        variants => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            doc => 'One or more BED files (or database ids) of the variants to validate',
            is_many => 1,
            shell_args_position => 3,
        },
        design => {
            is => 'Genome::FeatureList',
            doc => 'BED file (or database id) of the designs for the probes',
            shell_args_position => 1,
        },
        target => {
            is => 'Genome::FeatureList',
            doc => 'BED file (or database id) of the target region set',
            shell_args_position => 2,
        },
        tumor_sample => {
            is => 'Genome::Sample',
            doc => 'If there are no variants, specify the "tumor" sample directly',
        },
        normal_sample => {
            is => 'Genome::Sample',
            doc => 'If there are no variants, specify the "normal" sample directly',
        },
        sample_list_file => {
            is => 'Text',
            doc => 'A file of samples for which to make models.  Each line should have the samples that should be paired together with the control to the left.  If more than two are on each line, each sample on the line will be paired with every sample to its right in turn. (Provide this in lieu of variants or the --tumor-sample and --normal-sample.)'
        },
        region_of_interest_set => {
            is => 'Genome::FeatureList',
            doc => 'Specify this if reference coverage should be run on a different set than the target',
        },
        processing_profile => {
            is => 'Genome::ProcessingProfile::SomaticValidation',
            doc => 'Processing profile for the model',
        },
        groups => {
            is => 'Genome::ModelGroup',
            is_many => 1,
            doc => 'Group(s) to which to add the newly created model(s)',
        },
        auto_assign_inst_data => {
            is => 'Boolean',
            default => 1,
            doc => 'Automatically assign instrument data using the cron',
        },
        force => {
            is => 'Boolean',
            default => 0,
            doc => 'force model creation to occur, even if samples don\'t match',
        },
    ],
    has_transient_optional_input => [
        reference_sequence_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'Reference used for the variant calls',
        },
    ],
    has_transient_optional_output => [
        result_model_ids => {
            is => 'Number',
            via => 'result_models',
            to => 'id',
            is_many => 1,
            doc => 'ID of the model created by this command',
        },
        result_models => {
            is => 'Genome::Model::SomaticValidation',
            is_many => 1,
            doc => 'model created by this command',
        },
    ],
    doc => 'define a new somatic validation model',
};


sub help_detail {
    return <<'EOHELP'
To set up the model to run the validation process, three pieces of information are needed: the design (as sent to the vendor), the target set (as received from the vendor), and the variants to be validated. Each of these constituent parts are tracked individually by the analysis system, and this model takes the individual pieces and links them together.

First, the individual pieces need to be added to the system. For the designs we send to the vendor and targets we get back from the vendor, the files are stored as feature lists. For the lists of variants, we track them as detect variants results, either directly from the Somatic Variation pipeline or from manual curation. Then the parts are assembled with this command. The two main commands to add the individual pieces are:

`genome feature-list create` to create the feature lists, once for the design, and once for the target set.

`genome model somatic-validation manual-result` to record the manually curated results, if necessary. (One per file of variants.)
EOHELP
;
}

sub execute {
    my $self = shift;

    my @variants = $self->variants;
    if(scalar(grep { $_ } ($self->tumor_sample, scalar(@variants), $self->sample_list_file)) > 1) {
        die $self->error_message('Please supply only one of the following: variants, a tumor/normal sample pair, a sample list file');
    }

    if($self->target) {
        my $t = $self->target;
        if(!$t->content_type) {
            $t->content_type('validation');
        } elsif($t->content_type ne 'validation') {
            die $self->error_message('Specified target set has content-type ' . $t->content_type . ' when validation expected.');
        }
    }

    $self->resolve_reference_sequence_build;

    my @m;
    if($self->tumor_sample) {
        push @m, $self->_define_model(
            tumor_sample => $self->tumor_sample,
            normal_sample => $self->normal_sample,
        );
    } elsif(@variants) {
        my %variants_for_samples;
        for my $v (@variants) {
            my ($tumor_sample, $normal_sample) = $self->resolve_samples_for_variant_list($v);
            my $normal_sample_id = $normal_sample? $normal_sample->id : ''; #normal sample will be undef for tumor-only case
            $variants_for_samples{$tumor_sample->id}{$normal_sample_id} ||= [];
            push @{$variants_for_samples{$tumor_sample->id}{$normal_sample_id}}, $v;
        }

        for my $tumor_sample_id (keys %variants_for_samples) {
            for my $normal_sample_id (keys %{$variants_for_samples{$tumor_sample_id}}) {
                my $v = $variants_for_samples{$tumor_sample_id}{$normal_sample_id};
                push @m, $self->_define_model(
                    tumor_sample => Genome::Sample->get($tumor_sample_id),
                    ($normal_sample_id? (normal_sample => Genome::Sample->get($normal_sample_id)) : ()),
                    variants => $v,
                );
            }
        }
    } elsif($self->sample_list_file) {
        my @lines = Genome::Sys->read_file($self->sample_list_file);
        chomp @lines;

        #first check that all our samples really exist
        my %all_names;
        for my $line (@lines) {
            my @names = split(/\t/, $line);
            for my $name (@names) {
                $all_names{$name} = undef;
            }
        }

        my @samples = Genome::Sample->get(name => [keys %all_names]);
        for my $s (@samples) {
            $all_names{$s->name} = $s;
        }

        local %Command::V2::ALTERNATE_FROM_CLASS = (
            'Genome::Sample' => {
                'Genome::Library' => ['sample'],
                'Genome::InstrumentData' => ['sample'],
            },
        );

        my @bad_names;
        for my $name (keys %all_names) {
            unless(defined $all_names{$name}) {
                my $sample = eval { return $self->resolve_param_value_from_cmdline_text({name => 'sample', class => 'Genome::Sample', value => [$name]}); };
                if($sample) {
                    $all_names{$name} = $sample;
                } else {
                    push @bad_names, $name;
                }
            }
        }

        if(@bad_names) {
            die $self->error_message('Could not find samples for the following names: ' . join(", ", @bad_names));
        }

        #all samples found--now make the pairs
        for my $line (@lines) {
            my @names = split(/\t/, $line);

            if(scalar(@names) == 1) {
                #crete a single-bam model and move on
                push @m, $self->_define_model( tumor_sample => $all_names{$names[0]} );
            } else {
                for my $i (0..$#names-1) {
                    for my $j ($i+1..$#names) {
                        push @m, $self->_define_model(
                            normal_sample => $all_names{$names[$i]},
                            tumor_sample => $all_names{$names[$j]},
                        );
                    }
                }
            }
        }
    } else {
        die $self->error_message('Please supply one of the following: variants, a tumor/normal sample pair, a sample list file');
    }


    $self->result_models(\@m);
    $self->status_message('New models: ' . join(',', map($_->id, @m)));
    if($self->groups) {
        for my $group ($self->groups) {
            $group->assign_models(@m);
        }
    }
    return scalar @m;
}

sub _define_model {
    my $self = shift;
    my %params = @_;

    my $tumor_sample = $params{tumor_sample};

    if($tumor_sample and $tumor_sample->common_name and $tumor_sample->common_name eq 'normal') {
        $self->warning_message('Sample specified as tumor ' . $tumor_sample->__display_name__ . ' indicates it is a normal!');
    }

    my $normal_sample = $params{normal_sample};
    my $variants = $params{variants};
    my $processing_profile = $self->resolve_processing_profile($tumor_sample, $normal_sample);
    my $subject = $self->_resolve_subject_from_samples($tumor_sample, $normal_sample);

    #for now, use defaults and don't allow customization at definition time
    my $annotation_build = Genome::Model::ImportedAnnotation->annotation_build_for_reference($self->reference_sequence_build);
    my $previously_discovered_build = Genome::Model::ImportedVariationList->dbsnp_build_for_reference($self->reference_sequence_build);

    my $variant_results_by_type = $self->resolve_variant_list_types(@$variants);

    my @params;
    push @params, processing_profile => $processing_profile;
    push @params, subject => $subject;
    push @params, name => $self->name
        if defined $self->name;
    push @params, reference_sequence_build => $self->reference_sequence_build
        if defined $self->reference_sequence_build;
    push @params, design_set => $self->design
        if defined $self->design;
    push @params, target_region_set => $self->target
        if defined $self->target;
    push @params, snv_variant_list => $variant_results_by_type->{snv}
        if defined $variant_results_by_type->{snv};
    push @params, indel_variant_list => $variant_results_by_type->{indel}
        if defined $variant_results_by_type->{indel};
    push @params, sv_variant_list => $variant_results_by_type->{sv}
        if defined $variant_results_by_type->{sv};
    push @params, tumor_sample => $tumor_sample
        if defined $tumor_sample;
    push @params, normal_sample => $normal_sample
        if defined $normal_sample;
    push @params, annotation_build => $annotation_build
        if defined $annotation_build;
    push @params, previously_discovered_variations_build => $previously_discovered_build
        if defined $previously_discovered_build;
    push @params, auto_assign_inst_data => $self->auto_assign_inst_data 
        if defined $self->auto_assign_inst_data;

    if($self->region_of_interest_set) {
        push @params, region_of_interest_set => $self->region_of_interest_set;
    } elsif($self->target) {
        push @params, region_of_interest_set => $self->target;
    }

    my $m = Genome::Model->create(@params);
    return unless $m;

    $self->status_message('Successfully defined model: '.$m->__display_name__);
    return $m;
}

sub resolve_samples_for_variant_list {
    my $self = shift;
    my $potential_source = shift;

    my ($tumor_sample, $control_sample);

    if($potential_source->can('sample') and $potential_source->can('control_sample')) {
        #this is the expected case for manual results
        $tumor_sample = $potential_source->sample;
        $control_sample = $potential_source->control_sample;
    } elsif($potential_source->can('users')) {
        my @users = $potential_source->users;
        my @user_objects = map($_->user, @users);
        my @candidate_objects = grep($_->isa('Genome::Model::Build::SomaticVariation'), @user_objects);
        if(@candidate_objects) {
            for my $c (@candidate_objects) {
                if($c->can('normal_model') and $c->can('tumor_model')) {
                    #this is the expected case for results directly from somatic variation models
                    $tumor_sample = $c->tumor_model->subject;
                    $control_sample = $c->normal_model->subject;
                }
            }
        }
    }


    unless($tumor_sample or $control_sample) {
        $self->error_message('At least one sample is required to define a model. None found for ' . $potential_source->__display_name__);
        return;
    }

    return ($tumor_sample, $control_sample) if wantarray;
    die('Two items to return in a scalar context') if defined wantarray;
    return;
}

sub _resolve_subject_from_samples {
    my $self = shift;
    my $tumor_sample = shift;
    my $control_sample = shift;

    my $subject;
    if($tumor_sample) {
        if($control_sample and $tumor_sample->source ne $control_sample->source) {
            unless($self->force){
                my $problem = 'Tumor (' . $tumor_sample->name . ') and control (' . $control_sample->name . ') samples do not appear to have come from the same individual.';
                my $answer = $self->_ask_user_question(
                    $problem . ' Continue anyway?',
                    300,
                    "y.*|n.*",
                    "no",
                    "[y]es/[n]o",
                    );
                unless($answer and $answer =~ /^y/) {
                    $self->error_message($problem);
                    return;
                }
            }
        }

        $subject = $tumor_sample->source;
    } elsif($control_sample) {
        $subject = $control_sample->source;
    }

    return $subject;
}

sub resolve_processing_profile {
    my $self = shift;

    my $tumor_sample = shift;
    my $normal_sample = shift;

    return $self->processing_profile if $self->processing_profile;

    my $pp;
    if($tumor_sample and not $normal_sample) {
        $pp = Genome::Model::SomaticValidation->default_single_bam_profile();
    } else {
        #Nov 2011 default Somatic Validation
        $pp = Genome::Model::SomaticValidation->default_profile();
    }

    return $pp;
}


sub resolve_reference_sequence_build {
    my $self = shift;

    return $self->reference_sequence_build if $self->reference_sequence_build;

    my $rsb;
    for my $potential_indicator ($self->region_of_interest_set, $self->design, $self->target, $self->variants) {
        next unless $potential_indicator;
        if($potential_indicator->can('reference')) {
            $rsb = $potential_indicator->reference;
            last;
        } elsif($potential_indicator->can('reference_build')) {
            $rsb = $potential_indicator->reference_build;
        }
    }

    $self->reference_sequence_build($rsb);
    return $rsb;
}

sub resolve_variant_list_types {
    my $self = shift;
    my @variant_results = @_;

    my $results_by_type = {};
    for my $variant_result (@variant_results) {
        my $variant_type;
        if($variant_result->can('variant_type')) {
            $variant_type = $variant_result->variant_type;
        } elsif ($variant_result->can('_variant_type')) {
            $variant_type = $variant_result->_variant_type;
        } else {
            my @files = glob($variant_result->output_dir . '/*.hq');
            if(scalar @files > 1) {
                $self->error_message('Multiple .hq files found in result ' . $variant_result->id . ' at ' . $variant_result->output_dir);
                return;
            }
            unless(scalar @files) {
                $self->error_message('No .hq file found in result ' . $variant_result->id . ' at ' . $variant_result->output_dir);
                return;
            }

            ($variant_type) = $files[0] =~ /(\w+)\.hq/; #lame solution
        }

        $variant_type =~ s/s$//;
        if(exists $results_by_type->{$variant_type}) {
            $self->error_message('Multiple variant results have same type!');
            return;
        }
        $results_by_type->{$variant_type} = $variant_result;
    }

    return $results_by_type;
}

1;

