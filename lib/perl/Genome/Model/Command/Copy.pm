package Genome::Model::Command::Copy;

use strict;
use warnings;

use Genome;

use Regexp::Common;

class Genome::Model::Command::Copy {
    class_name => __PACKAGE__,
    is => 'Command::V2',
    has_input => [
        model => {
            is => 'Genome::Model',
            shell_args_position => 1,
            doc => 'The source model to copy.'
        },
        overrides => {
            is_many => 1,
            is_optional => 1,
            shell_args_position => 3,
            doc => 'Properties to override in the new model in the form key_name=value. Properties include name, subject, processing_profile, instrument_data, auto_build_alignments, auto_assign_inst_data and all model inputs.'
        },
    ],
    has_param => [
        recurse => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'also copy all inputs recursively which are models',
        },
        start => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'start the new model (with recurse down-stream models will be set to auto build)',
        },
    ],
    has_output => [
        _new_model => { is_optional => 1, },
        copied_model_pairs => {
            is => 'Genome::Model::Pair',
            is_transient => 1,
            is_input => 1,
            is_output => 1,
            is_optional => 1,
            is_many => 1,
        },
    ],
    doc => 'create a new genome model based on an existing one'
};

sub sub_command_sort_position { 2 }

sub help_synopsis {
    return <<"EOS"

Copy model 123456789 to new model using default name and overriding the subject

 genome model copy 123456789 subject=name=H_SAMPLE

Copy model 123456789 to new model with name "Copy of my Awesome Model", overriding processing profile and auto build alignments:

 genome model copy 123456789 "name=Copy of my Awesome Model" processing_profile="use this processing profile instead" auto_build_alignments=0

Copying with the "recurse" option will copy input models, and their input models, etc.:

 genome model copy name=apipe-test-clinseq-v4 name=apipe-test-clinseq-v5 --recurse

Recursion happens automatically when changes are made to underlying models, but only where necessary.
This copy notes that the PP must change for the refalign models under the wgs_model, and as such requires 3 new models to achieve the change:

 > genome model copy 2888708572 name=my-all1-bwamem wgs_model.tumor_model.processing_profile=2828673 wgs_model.normal_model.processing_profile=2828673
 NEW MODEL: my-all1-bwamem (2892350615)
 NEW MODEL: my-all1-bwamem.wgs_model (2892350616)
 NEW MODEL: my-all1-bwamem.wgs_model.normal_model (2892350617)
 NEW MODEL: my-all1-bwamem.wgs_model.tumor_model (2892350618)

When a model is shared as an input in several places during recursion, it is only copied once:

 > genome model copy name=tst1-clinseq name=tst1-clinseq-copy1
 NEW MODEL: tst1-clinseq-copy2 (403f62e0936c4f7ca97bb94e5409b7b0)
 NEW MODEL: tst1-clinseq-copy2.tumor_rnaseq_model (6f31ee70ebcd40868ef828ae7ffd76d1)
 NEW MODEL: tst1-clinseq-copy2.exome_model (fc7c6c2a4cfb471ab07b9b597d650fd6)
 NEW MODEL: tst1-clinseq-copy2.exome_model.tumor_model (2f3566a600dd48a6b43a1892ad33944f)
 NEW MODEL: tst1-clinseq-copy2.exome_model.tumor_model.genotype_microarray (2b93995ae05b4bc5a97580be4931ce2f)
 NEW MODEL: tst1-clinseq-copy2.exome_model.normal_model (a024d98c577b4c5fbc8feed296dcc9cb)
 NEW MODEL: tst1-clinseq-copy2.exome_model.normal_model.genotype_microarray (00db05ec23384eb48d13adf0ecddfb73)
 NEW MODEL: tst1-clinseq-copy2.wgs_model (a5ddcdf51a184d4d90bc395a0a471f4d)
 NEW MODEL: tst1-clinseq-copy2.wgs_model.tumor_model (ec13a9043b68434aabce500934e8f1a1)
 RE-USING tst1-clinseq-copy2.exome_model.tumor_model.genotype_microarray as tst1-clinseq-copy2.wgs_model.tumor_model.genotype_microarray
 NEW MODEL: tst1-clinseq-copy2.wgs_model.normal_model (375bc615cfe24cad9dce14a3a3f85e65)
 RE-USING tst1-clinseq-copy2.exome_model.normal_model.genotype_microarray as tst1-clinseq-copy2.wgs_model.normal_model.genotype_microarray

EOS
}

sub help_detail {
    return <<"EOS"
Copy a model to a new model overriding some of the original model's properties. This will copy the model's definition only, and not any builds.

Override properties are actual model properties. If you want to override reference_sequence_build, you use that, and not reference_sequence_build_id. Give override as space separated bare arguments. Use the format: property=value. The value can be an id or filter string.  If the property can have many values, use a filter string, or multiple key=value pairs. To set something to undefined, use 'property='.

These properties can be overridden on all models:

 name
 processing_profile
 subject

All model input specific to that type can also be overriden.
See model input names and values with:

 genome model input show MODELID

Examples:

 name=\$NEW_NAME
 processing_profile=\$ID
 processing_profile=name=\$NAME
 instrument_data=library_id=\$LIB_ID,is_paired_end=1
 db_dnp_build=id=\$DB_SNP_ID
 db_snp_build=model_id=\$MODEL_ID,status=Succeeded

If no name is given in the list of overrides, a new model name will be auto-generated.

Commas in the filter string will be interpreted as key-value seperators.

To set any property to undefined have nothing after the equal sign.  To set no instrument_data, for instance:
 'instrument_data='

For properties which have multiple values, they can be specified multiple times:

 instrument_data=\$ID1 instrument_data=\$ID2

To specify a new subject, you can use the ID, or any expression to match the subject, such as name or comon_name:

 subject=\$ID
 subject=name=\$NAME
 subject=common_name=\$NAME

When the "recurse" option is set, inputs which are also models will also be copied.
The new models will be named after the primary model with the suffix "." followed by the name of the input.

 genome model copy name=apipe-test-clinseq-v4 name=apipe-test-clinseq-v5 --recurse

 'model' may require verification...
 Resolving parameter 'model' from command argument 'name=apipe-test-clinseq-v4'... found 1
 NEW MODEL: apipe-test-clinseq-v5 (2890809239)
 NEW MODEL: apipe-test-clinseq-v5.tumor_rnaseq_model (2890809240)
 NEW MODEL: apipe-test-clinseq-v5.exome_model (2890809241)
 NEW MODEL: apipe-test-clinseq-v5.exome_model.normal_model (2890809242)
 NEW MODEL: apipe-test-clinseq-v5.exome_model.normal_model.genotype_microarray (2890809243)
 NEW MODEL: apipe-test-clinseq-v5.exome_model.tumor_model (2890809244)
 NEW MODEL: apipe-test-clinseq-v5.exome_model.tumor_model.genotype_microarray (2890809245)
 NEW MODEL: apipe-test-clinseq-v5.wgs_model (2890809246)
 NEW MODEL: apipe-test-clinseq-v5.wgs_model.normal_model (2890809247)
 NEW MODEL: apipe-test-clinseq-v5.wgs_model.normal_model.genotype_microarray (2890809248)
 NEW MODEL: apipe-test-clinseq-v5.wgs_model.tumor_model (2890809249)
 NEW MODEL: apipe-test-clinseq-v5.wgs_model.tumor_model.genotype_microarray (2890809250)

EOS
}

sub copied_model_pairs_as_hash {
    my $self = shift;
    my %hash;
    for my $pair ($self->copied_model_pairs) {
        my $list = $hash{$pair->first_id} ||= [];
        push @$list, $pair->second_id;
    }
    return \%hash;
}

sub execute {
    my $self = shift;

    my $model = $self->model;
    if ( not $model ) {
        $self->error_message('No model to copy');
        return;
    }
    #$self->status_message('Copy model: '.$model->__display_name__);

    my %overrides = $self->params_from_param_strings($model->class, $self->overrides);
    if (not %overrides) {
        die "eror parsing overrides!";
    }

    my %indirect_overrides;
    my @errors;
    for my $name (keys %overrides) {
        if ($name =~ /^(.*?)\.(.*)$/) {
            my $first = $1;
            my $rest = $2;
            unless ($model->can($first)) {
                push @errors, "No property $first on model " . $model->__display_name__ . ".  Cannot set indirect property $name!";
                next;
            }
            $indirect_overrides{$first}{$rest} = delete $overrides{$name};
        }
    }
    if (@errors) {
        die $self->error_message(join("\n",@errors));
    }

    my $new_model = $model->copy(%overrides);
    return if not $new_model;
    $self->_new_model($new_model);

    $self->status_message("NEW MODEL: ".$new_model->__display_name__);

    $self->add_copied_model_pair(
        Genome::Model::Pair->get(first => $model, second => $new_model)
    );
    my @copied_model_pairs = $self->copied_model_pairs;
    my $copied_model_pairs = $self->copied_model_pairs_as_hash();

    if ($self->recurse or %indirect_overrides) {
        my @input_assoc = $new_model->input_associations();
        my $n = 0;
        for my $assoc (@input_assoc) {
            my $input_name = $assoc->name;
            my $input_value = $assoc->value;
            unless ($input_value->isa("Genome::Model")) {
                next;
            }

            my $new_input_model_name = $new_model->name . "." . $assoc->name;

            #TODO make this an object which parses from this text blob instead of a text blob
            my @overrides;
            my $more_overrides = delete $indirect_overrides{$input_name};
            if ($more_overrides) {
                for my $key (keys %$more_overrides) {
                    my $value = $more_overrides->{$key};
                    push @overrides, "$key=$value";
                }
            }

            if ($self->recurse or @overrides) {
                my $new_input_model;
                # if we already copied this model during recursive calls to this command,
                # the pairing of the original and new model will be available
                my $previously_copied_model_ids = $copied_model_pairs->{$input_value->id};
                if ($previously_copied_model_ids) {
                    for my $previously_copied_model_id (@$previously_copied_model_ids) {
                        my $previously_copied_model = Genome::Model->get($previously_copied_model_id);
                        die "no model found: $previously_copied_model_id!" unless $previously_copied_model;
                        my $bx = $self->_overrides_as_boolexpr($previously_copied_model->class,@overrides);
                        if ($bx->evaluate($previously_copied_model)) {
                            $self->status_message("RE-USING " . $previously_copied_model->name . " as " . $new_input_model_name);
                            $new_input_model = $previously_copied_model;
                            last;
                        }
                    }
                    unless ($new_input_model) {
                        $self->status_message("Making an ADDITIONAL copy of " . $input_value->name . " since previous copies do not match (@$previously_copied_model_ids)");
                    }
                }

                # if we did not already copy the model, or we need this copy to be different
                # we make a new copy
                unless ($new_input_model) {
                    my @new_copied_model_pairs = eval {
                        my $cmd = __PACKAGE__->create(
                            model => $input_value,
                            recurse => $self->recurse,
                            start => $self->start,
                            overrides => [ "name=$new_input_model_name", @overrides ],
                            copied_model_pairs => [@copied_model_pairs],
                        );
                        $cmd->execute();
                        return $cmd->copied_model_pairs;
                    };

                    my $err = $@ || '';

                    $new_input_model = Genome::Model->get(name => $new_input_model_name);
                    if ($err or not $new_input_model) {
                        $new_model->delete;
                        die "failed to generate model $new_input_model_name for $input_name replacing $input_value: $err!";
                    }

                    my ($both,$old,$new) = UR::Util::intersect_lists([@copied_model_pairs],[@new_copied_model_pairs]);
                    if (@$old or not @$new) {
                        die "New list of copied model pairs should be a superset of the old list got (old,new,both)!: " . Data::Dumper::Dumper($old,$new,$both);
                    }

                    unless (@overrides) {
                        # If there were overrides, this new model cannot safely be used
                        # in place of other copies of the original model.
                        # An improvement on this logic would see if the overrides happened to match.  In lieu of that this is conservative.
                        $self->copied_model_pairs(\@new_copied_model_pairs);
                        @copied_model_pairs = $self->copied_model_pairs;
                        $copied_model_pairs = $self->copied_model_pairs_as_hash();
                    }
                }

                $new_model->$input_name($new_input_model);
                $n++;
            }
        }
        if (%indirect_overrides) {
            die "indirect overrides for values that are not set on the original model are not supported yet, sadly: " . Data::Dumper::Dumper(\%indirect_overrides);
        }

        if ($self->start) {
            if ($n) {
                # this model will build when input models complete
                $new_model->auto_build(1);
            }
            else {
                # no input models: build now
                my $build = $new_model->add_build();
                $self->status_message("...starting...");
                $build->start();
            }
        }
    }

    return 1;
}

sub _overrides_as_boolexpr {
    my $self = shift;
    my $class = shift;
    my $bx;
    if (@_) {
        my @params = $self->params_from_param_strings($class,@_);
        $bx = $class->define_boolexpr(@params);
    }
    else {
        $bx = $class->define_boolexpr();
    }
    return $bx;
}

sub params_from_param_strings {
    # preprocess command line inputs for Model->copy function.
    my ($self, $class, @param_strings) = @_;

    Carp::confess('No param strings to convert to params') if not @param_strings;

    my %params;
    my $meta = $class->__meta__;
    for my $param_string ( @param_strings ) {
        my ($key, $value) = split('=', $param_string, 2);
        if ($key =~ /\./) {
            # defer processing
            # this will be handled in recursive calls to copy()
            $params{$key} = $value;
            next;
        }

        my $property = $meta->property_meta_for_name($key);
        if ( not $property ) {
            $class->error_message("Failed to find model property: $key");
            return;
        }

        if ( my ($unallowed) = grep { $property->$_ } (qw/ is_calculated is_constant is_transient /) ){
            $class->error_message("Property ($key) cannot be given on the command line because it is '$unallowed'");
            return;
        }

        if ( not defined $value or $value eq '' ) {
            $params{$key} = undef;
            next;
        }

        my @values = $value;
        my $data_type = $property->data_type;
        if (defined($data_type) and not grep { $data_type =~ /^$_$/i } (qw/ boolean integer number string text ur::value /) ) { # hacky...if u kno a better way...
            my $filter = ( $value =~ /^$RE{num}{int}$/ || $value =~ /^[[:xdigit:]]{32}$/) ? 'id='.$value : $value;
            my $bx = eval { UR::BoolExpr->resolve_for_string($data_type, $filter); };
            if ( not $bx ) {
                $class->error_message("Failed to create expression for $key ($data_type) from '$value'");
                return;
            }
            @values = $data_type->get($bx);
            if ( not @values ) {
                $DB::single = 1;
                die $class->error_message("Failed to get $key ($data_type) for $value");
            }
        }

        if ( $property->is_many ) {
            push @{$params{$key}}, @values;
        }
        elsif ( @values > 1 or exists $params{$key} ) {
            $class->error_message(
                "Singular property ($key) cannot have more than one value (".join(', ', grep { defined } (@values, $params{$key})).')'
            );
            return;
        }
        else {
            $params{$key} = $values[0];
        }
    }

    return %params;
}


1;

