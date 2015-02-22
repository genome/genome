package Genome::Model::Comparison;
use strict;
use warnings;
use Genome;

class Genome::Model::Comparison {
    is => 'Genome::Model',
    has_input => [
        from_models => {
            is => 'Genome::Model',
            is_many => 1,
            doc => 'the models built in a prior way for which a new model will be made and tested (i.e. "groups.name=cancer-test1")',
        },
        changes => {
            is => 'Text', #'UR::BoolExpr' with extensions
            is_many => 1,
            doc => 'changes to the "from" models which are being tested (param=value or input=value)',
        },
    ],
    has_param => [
        aspects => {
            is => 'Text', # used to find a module to perform the comparison
            is_many => 1,
            is_optional => 1,
            doc => 'the aspects of each build pair to compare (i.e. "variants", "alignments")',
        },
    ],
    doc => "pipeline to compare models across processing changes (IN DEVELOPMENT)",
};

sub define_by { return 'Genome::Model::Command::Define::BaseMinimal'; }

sub _help_synopsis {
    return <<EOS;
 genome model define comparison \
    --name test-bwasw1-on-somatic \
    --from test-cancer-aml31 \
    --changes "tumor_model.read_aligner_name=bwa-sw
                and normal_model.read_aligner_name=bwa-sw
                and tumor_model.read_aligner_version=0.6.1
                and normal_model.read_aligner_version=0.6.1 " \
    --aspects alignments,variants,metrics

 genome model define comparison \
    --name test-clinseq-noexome-dgidb \
    --from id:2890260793/2890224790 \
    --changes "exome_model=''" \
    #--aspects drug-gene-interactions,metrics,run-time
    --processing-profile "compare clinseq default"
EOS
}

sub _help_detail {
    return <<EOS
Compare a set of models (specifiable by expression) to a newlly built set with a specific set of changes to test.

The changes are in the form key=value, where a key is a processing-profile parameter on the "from" model.

The aspects to be compared are determined by the processing profile used by *this* model.

WARNING: This is in development and currently does no work.
EOS
}

sub create {
    my $class = shift;
    my $subject = Genome::Taxon->get(name => "unknown");
    my $bx = $class->define_boolexpr(@_)->add_filter(subject => $subject);
    return $class->SUPER::create($bx);
}

sub _execute_build {
    my $self = shift;
    my $build = shift;

    #
    # put the "from" models into a group owned by this model
    #

    my @from_builds = sort $build->from_builds;

    my @from_models = sort {$a->id cmp $b->id } map { $_->model } @from_builds;

    # Get From group
    my $from_group_name = $self->name . '.from';
    my $from_group = Genome::ModelGroup->get(name => $from_group_name);

    # Get To group
    my $to_group_name = $self->name . '.to';
    my $to_group = Genome::ModelGroup->get(name => $to_group_name);

    # TODO: Rather that relying on the above names there should be a linkage or sorts to
    # show the connections between these from/to models and groups

    if ($from_group) {
        my @members = sort $from_group->models;
        unless ("@members" eq "@from_models") {
            $self->_rename_model_group($from_group);
            $from_group = undef;
            $self->_rename_model_group($to_group);
            $to_group = undef;
        }
    }
    unless ($from_group) {
        $from_group = Genome::ModelGroup->create(name => $from_group_name);
        for my $from_model (sort { $a cmp $b } @from_models) {
            $from_group->assign_models($from_model);
        }
    }
    $self->status_message('Found "from" model group : '. $from_group->__display_name__);

    #
    # make a set of "to" models with the changes and build them
    #

    # TODO: # pull the logic from this copy command and call it on each $from_model
    # then make an entity representing a pair (UR::Value::Pair?)
    # and set these on the build in some way.
    # Right now we rely on a call to ->members to always sort the same.

    # TODO: when <changes> is ambiguous (in-clause instead of single value),
    # break each down into N model group copies.  If multiple are ambiguous
    # (M) there will be a matrix of model groups created of M dimensions
    # for all combinations.

    # TODO: when changes has a field name which is indirect through an input model,
    # make a copy of the input model and use it as an input so the change is "true"
    # for the new model

    unless ($to_group) {
        my @changes = $build->changes;
        Genome::ModelGroup::Command::Copy->execute(
            from => $from_group,
            to => $to_group_name,
            changes => \@changes,
        );
        $to_group = Genome::ModelGroup->get(name => $to_group_name);
        unless ($to_group) {
            die $self->error_message("Failed to create model group $to_group_name!");
        }
    }
    $self->status_message('Found "to" model group : '. $to_group->__display_name__);
    my @to_models = $to_group->models;

    #
    # go through each of the build pairs and compare aspects
    #

    my @aspects = $build->processing_profile->aspects;
    for my $aspect (@aspects) {
        my $aspect_dir = $build->data_directory .'/'. $aspect;
        Genome::Sys->create_directory($aspect_dir);
    }

    for (my $n = 0; $n < scalar(@to_models); $n++) {
        my $to_model = $to_models[$n];
        my $to_build = $to_model->last_complete_build;

        $DB::single=1;
        unless ($to_build) {
            my $existing_build = $to_model->current_build;
            if ($existing_build && $existing_build->status eq 'Running') {
                $self->status_message('Please wait for existing current build '. $existing_build->__display_name .' to complete running.');
            } elsif ( $existing_build && $existing_build->status eq 'Unstartable') {
                $self->status_message('Existing current build '. $existing_build->__display_name__ .' is '. $existing_build->status );
            } else {
                if ($existing_build) {
                    $self->status_message('Existing current build '. $existing_build->__display_name__ .' is '. $existing_build->status );
                }
                $self->status_message('Starting a build for model '. $to_model->__display_name__);
                my $build = $to_model->add_build;
                unless ($build->start) {
                    die('Failed to start build: '. $build->id);
                }
            }
            next;
        }
        my $from_build = $from_builds[$n];
        my $from_model = $from_build->model;

        $self->status_message("Compare build " . $from_build->__display_name__ . " to " . $to_build->__display_name__);
        for my $aspect (@aspects) {
            my $aspect_camel_case = join('', map { ucfirst(lc($_)) } split('-',$aspect)); # TODO: switch to UR::Util method
            my $aspect_dir = $build->data_directory .'/'. $aspect;
            my $compare_aspect_dir = $aspect_dir .'/'. $n;

            Genome::Sys->create_directory($compare_aspect_dir);

            $self->status_message('Comparing '. $n .' : '. $from_build->model->name .' to '. $to_build->model->name);
            Genome::Sys->shellcmd( cmd => 'touch '.$compare_aspect_dir .'/from_'. $from_build->id );
            Genome::Sys->shellcmd( cmd => 'touch '. $compare_aspect_dir .'/to_'. $to_build->id );

            my $from_model = $from_build->model;
            my $from_model_class = $from_model->class;
            my @classes_to_check = ($from_model_class,$from_model_class->inheritance);
            my @comparison_classes;
            for my $model_class (@classes_to_check) {
                my $comp_class = $model_class . '::Command::Compare::' .$aspect_camel_case;
                my $comp_class_meta = UR::Object::Type->get($comp_class);
                if ($comp_class_meta) {
                    push @comparison_classes, $comp_class;
                    $build->status_message("doing comparison for $comp_class on pair $n");
                    Genome::Sys->shellcmd(cmd => 'touch '. $compare_aspect_dir .'/'. $comp_class);
                    last;
                }
            }

            unless (@comparison_classes) {
                die "No comparisons defined for $aspect!";
            }

            for my $comparison_class (@comparison_classes) {
                unless ($comparison_class->execute(
                    from_build => $from_build,
                    to_build => $to_build,
                    output_dir => $compare_aspect_dir,
                )) {
                    die('Failed to execute '. $comparison_class);
                }
            }

            # TODO: figure out how to pair @from_sr and @to_sr
            #my %comparisons;
            #my @from_sr = $from_build->software_results;
            #my @to_sr = $to_build->software_results;
            #for my $from_sr (@from_sr) {
            #    my $sr_class = $from_sr->class;
            #    my $compare_class = $sr_class . "::Command::Compare::" . ucfirst(lc($aspect));
            #    if (UR::Object::Type->get($sr_class)) {
            #        my $a = $comparisons{$compare_class} ||= [];
            #        push @$a, $from_sr;
            #    }
            #}

            # each ::Compare::X module should produce a software result for the build pair
            # each comparison should have an output directory linked under the build directory
        }
    }
    $DB::single=1;
    return 1;
}

sub _rename_model_group {
    my ($self, $group) = @_;
    Carp::confess("TODO: implement me to get stale subordinate models out-of-the-way!");
}

sub _doc_examples {
    return <<EOS
 More examples.  Edit these.

 # this takes two clinseq models and performs a null comparison
 # that profile is mostly for testing, or to make it easy to build things
 # and decide later what to compare when you copy the model
 # (this runs now)
 genome model define comparison \
    --from id:2890260793/2890224790 \
    --changes "exome_model=''" \
    --processing-profile name="compare nothing" \
    --name my-comparison1

 # see how drug-gene-interactions change when exome data is not used
 # (presumes Genome::Model::ClinSeq::Comparison::DrugGeneInteractions exists and is put into a processing profile)
 genome model define comparison \
    --from id:2890260793/2890224790 \
    --changes "exome_model=''" \
    --processing-profile name="compare dgidb results" \
    --name my-comparison2

 # see how switching breakdancer version and subsampling to 1/2 depth affects metrics
 # (presumes Genome::Model::Comparison::Metrics exists)
 # (presumes we hae a subsample_reads parameter which can take N% or Nx)
 genome model define comparison \
    --from groups.name=testdata-aml31-somatic-variation \
    --changes "sv_detection_strategy=~s/breakdancer 1.3/breakdancer 1.3.7/g" \
    --changes "subsample_reads=50%" \
    --processing-profile name="compare metrics" \
    --name my-comparison3

 # see how varying bwa -q from nothing to 0, 2, 5, and 10 compares
 # because the changes has an in-clause, we will do multiple output model groups
 # because aspects are listed explicitly it will dynamically will get/create a processing profile to compare those things
 # presumes Genome::Model::SomaticVariation::Comparison::Variants exists
 # presumes Genome::Model::SomaticVariation::Comparison::Alignments (probably delegates to refalign comparison)
 # presumes Genome::Model::Comparison::Metrics (on the base class: usable across model types)
 genome model define comparison \
    --from groups.name=testdata-cancer-aml31 \
    --changes "aligner_params in ['', '-q 0', '-q 2', '-q 5', '-q 10']"
    --aspects metrics,alignments,variants
    --name my-comparison4

EOS
}

1;

