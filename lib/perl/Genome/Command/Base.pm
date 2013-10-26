package Genome::Command::Base;

use strict;
use warnings;

use Genome;

use Data::Dumper;
use IO::File;
require File::Basename;

class Genome::Command::Base {
    is => ['Genome::Command::ColorMixin', 'Command'],
    is_abstract => 1,
    has => [
        show_display_command_summary_report => {
            is => 'Boolean',
            default => 1,
            doc => 'Enable/disable summary report.',
            is_transient => 1, # hide this from docs
        },
    ],
    has_optional => [
        _total_command_count => {
            is => 'Integer',
            default => 0,
        },
        _command_errors => {
            is => 'HASH',
            doc => 'values can be an array ref if multiple errors for on command',
            default => {},
        },
    ],
    attributes_have => [
        require_user_verify => {
            is => 'Boolean',
            is_optional => 1,
            # 0 = prevent verify, 1 = force verify, undef = allow auto verify
        },
    ],
};

our %ALTERNATE_FROM_CLASS = (
    # find_class => via_class => via_class_methods
    # first method is the default method
    # the default method is used automatically if not the paramater
    # data type so it should be the most verbose option
    'Genome::Sample' => {
        'Genome::PopulationGroup' => ['samples'],
    },
    'Genome::InstrumentData' => {
        'Genome::Model' => ['instrument_data'],
        'Genome::Model::Build' => ['instrument_data'],
        'Genome::ModelGroup' => ['instrument_data'],
    },
    'Genome::Model' => {
        'Genome::Model::Build' => ['model'],
        'Genome::ModelGroup' => ['models'],
    },
    'Genome::Model::Build' => {
        'Genome::Model' => ['builds'],
    },
);
# This will prevent infinite loops during recursion.
our %SEEN_FROM_CLASS = ();
our @error_tags;
our $MESSAGE;

sub resolve_param_value_from_cmdline_text {
    my ($self, $param_info) = @_;
    my $param_name  = $param_info->{name};
    my $param_class = $param_info->{class};
    my @param_args  = @{$param_info->{value}};
    my $param_str   = join(',', @param_args);

    if (ref($param_class) eq 'ARRAY') {
        my @param_class = @$param_class;
        if (@param_class > 1) {
            die 'Multiple data types on command arguments are not supported.';
        } else {
            $param_class = $param_class[0];
        }
    }

    my $param_resolve_message = "Resolving parameter '$param_name' from command argument '$param_str'...";
    my $pmeta = $self->__meta__->property($param_name);
    my $require_user_verify = $pmeta->{'require_user_verify'};

    my @results;
    my $bx = eval { UR::BoolExpr->resolve_for_string($param_class, $param_str) };
    my $bx_error = $@;
    if ($bx) {
        @results = $param_class->get($bx);
        if (@results > 1 && !defined($require_user_verify)) {
            $require_user_verify = 1;
        }
    } else {
        for my $arg (@param_args) {
            %SEEN_FROM_CLASS = ();

            # call resolve_param_value_from_text without a via_method to "bootstrap" recursion
            my @arg_results = $self->resolve_param_value_from_text($arg, $param_class);

            if (@arg_results != 1 && !defined($require_user_verify)) {
                $require_user_verify = 1;
            }

            push @results, @arg_results;
        }
    }
    if (@results) {
        # the ALTERNATE_FROM_CLASS stuff leads to non $param_class objects in results
        @results = $self->_unique_elements(@results);
        @results = grep { $_->isa($param_class) } @results;

        $self->status_message($param_resolve_message . " found " . @results);
    }
    else {
        if ($bx_error) {
            $self->status_message($bx_error);
        }
        $self->status_message($param_resolve_message . " none found.");
    }

    return unless (@results);

    my $limit_results_method = "_limit_results_for_$param_name";
    if ( $self->can($limit_results_method) ) {
        @results = $self->$limit_results_method(@results);
        return unless (@results);
    }
    @results = $self->_unique_elements(@results);
    if ($require_user_verify) {
        if (!$pmeta->{'is_many'} && @results > 1) {
            $MESSAGE .= "\n" if ($MESSAGE);
            $MESSAGE .= "'$param_name' expects only one result.";
        }
        @results = $self->_get_user_verification_for_param_value($param_name, @results);
    }
    while (!$pmeta->{'is_many'} && @results > 1) {
        $MESSAGE .= "\n" if ($MESSAGE);
        $MESSAGE .= "'$param_name' expects only one result, not many!";
        @results = $self->_get_user_verification_for_param_value($param_name, @results);
    }

    if (wantarray) {
        return @results;
    }
    elsif (not defined wantarray) {
        return;
    }
    elsif (@results > 1) {
        Carp::confess("Multiple matches found!");
    }
    else {
        return $results[0];
    }
}

sub resolve_param_value_from_text {
    my ($self, $param_arg, $param_class, $via_method) = @_;

    unless ($param_class) {
        $param_class = $self->class;
    }

    $SEEN_FROM_CLASS{$param_class} = 1;

    my @results_by_string;
    if ($param_class->can('_resolve_param_value_from_text_by_name_or_id')) {
        @results_by_string = $param_class->_resolve_param_value_from_text_by_name_or_id($param_arg);
    }
    else {
        @results_by_string = $self->_resolve_param_value_from_text_by_name_or_id($param_class, $param_arg);
    }

    my @results = @results_by_string;
    unless (@results) {
        @results = $self->_resolve_param_value_via_related_class_method($param_class, $param_arg, $via_method);
    }

    if ($via_method) {
        @results = map { $_->$via_method } @results;
    }

    if (wantarray) {
        return @results;
    }
    elsif (not defined wantarray) {
        return;
    }
    elsif (@results > 1) {
        Carp::confess("Multiple matches found!");
    }
    else {
        return $results[0];
    }
}

sub _resolve_param_value_via_related_class_method {
    my ($self, $param_class, $param_arg, $via_method) = @_;
    my @results;
    my $via_class;
    if (exists($ALTERNATE_FROM_CLASS{$param_class})) {
        $via_class = $param_class;
    }
    else {
        for my $class (keys %ALTERNATE_FROM_CLASS) {
            if ($param_class->isa($class)) {
                if ($via_class) {
                    $self->error_message("Found additional via_class $class but already found $via_class!");
                }
                $via_class = $class;
            }
        }
    }
    if ($via_class) {
        my @from_classes = sort keys %{$ALTERNATE_FROM_CLASS{$via_class}};
        while (@from_classes && !@results) {
            my $from_class  = shift @from_classes;
            my @methods = @{$ALTERNATE_FROM_CLASS{$via_class}{$from_class}};
            my $method;
            if (@methods > 1 && !$via_method && !$ENV{GENOME_NO_REQUIRE_USER_VERIFY}) {
                $self->status_message("Trying to find $via_class via $from_class...\n");
                my $method_choices;
                for (my $i = 0; $i < @methods; $i++) {
                    $method_choices .= ($i + 1) . ": " . $methods[$i];
                    $method_choices .= " [default]" if ($i == 0);
                    $method_choices .= "\n";
                }
                $method_choices .= (scalar(@methods) + 1) . ": none\n";
                $method_choices .= "Which method would you like to use?";
                my $response = $self->_ask_user_question($method_choices, 0, '\d+', 1, '#');
                if ($response =~ /^\d+$/) {
                    $response--;
                    if ($response == @methods) {
                        $method = undef;
                    }
                    elsif ($response >= 0 && $response <= $#methods) {
                        $method = $methods[$response];
                    }
                    else {
                        $self->error_message("Response was out of bounds, exiting...");
                        exit;
                    }
                    $ALTERNATE_FROM_CLASS{$via_class}{$from_class} = [$method];
                }
                elsif (!$response) {
                    $self->status_message("Exiting...");
                }
            }
            else {
                $method = $methods[0];
            }
            unless($SEEN_FROM_CLASS{$from_class}) {
                #$self->debug_message("Trying to find $via_class via $from_class->$method...");
                @results = eval {$self->resolve_param_value_from_text($param_arg, $from_class, $method)};
            }
        } # END for my $from_class (@from_classes)
    } # END if ($via_class)
    return @results;
}

sub _resolve_param_value_from_text_by_bool_expr {
    my ($self, $param_class, $arg) = @_;

    my @results;
    my $bx = eval {
        UR::BoolExpr->resolve_for_string($param_class, $arg);
    };
    if ($bx) {
        @results = $param_class->get($bx);
    }
    else {
        die "Not a valid BoolExpr";
    }
    #$self->debug_message("B: $param_class '$arg' " . scalar(@results));

    return @results;
}

sub _resolve_param_value_from_text_by_name_or_id {
    my ($self, $param_class, $str) = @_;
    my (@results);

    my $class_meta = $param_class->__meta__;
    my @id_property_names = $class_meta->id_property_names;
    if (@id_property_names == 0) {
        die "Failed to determine id property names for class $param_class.";
    }

    my $first_type = $class_meta->property_meta_for_name($id_property_names[0])->data_type || '';
    if (@id_property_names > 1 or $first_type eq 'Text' or $str =~ /^-?\d+$/) { # try to get by ID
        @results = $param_class->get($str);
    }
    if (!@results && $param_class->can('name')) {
        @results = $param_class->get(name => $str);
        unless (@results) {
            @results = $param_class->get("name like" => "$str");
        }
    }
    #$self->debug_message("S: $param_class '$str' " . scalar(@results));

    return @results;
}

sub _get_user_verification_for_param_value {
    my ($self, $param_name, @list) = @_;

    my $n_list = scalar(@list);
    if ($n_list > 200 && !$ENV{GENOME_NO_REQUIRE_USER_VERIFY}) {
        my $response = $self->_ask_user_question("Would you [v]iew all $n_list item(s) for '$param_name', (p)roceed, or e(x)it?", 0, '[v]|p|x', 'v');
        if(!$response || $response eq 'x') {
            $self->status_message("Exiting...");
            exit;
        }
        return @list if($response eq 'p');
    }

    my @new_list;
    while (!@new_list) {
        @new_list = $self->_get_user_verification_for_param_value_drilldown($param_name, @list);
    }

    my @ids = map { $_->id } @new_list;
    $self->status_message("The IDs for your selection are:\n" . join(',', @ids) . "\n\n");
    return @new_list;
}
sub _get_user_verification_for_param_value_drilldown {
    my ($self, $param_name, @results) = @_;
    my $n_results = scalar(@results);
    my $pad = length($n_results);

    # Allow an environment variable to be set to disable the require_user_verify attribute
    return @results if ($ENV{GENOME_NO_REQUIRE_USER_VERIFY});
    return if (@results == 0);

    my @dnames = map {$_->__display_name__} grep { $_->can('__display_name__') } @results;
    my $max_dname_length = @dnames ? length((sort { length($b) <=> length($a) } @dnames)[0]) : 0;
    my @statuses = map {$_->status || 'missing_status'} grep { $_->can('status') } @results;
    my $max_status_length = @statuses ? length((sort { length($b) <=> length($a) } @statuses)[0]) : 0;
    @results = sort {$a->__display_name__ cmp $b->__display_name__} @results;
    @results = sort {$a->class cmp $b->class} @results;
    my @classes = $self->_unique_elements(map {$_->class} @results);

    my $response;
    my @caller = caller(1);
    while (!$response) {
        $self->status_message("\n");
        # TODO: Replace this with lister?
        for (my $i = 1; $i <= $n_results; $i++) {
            my $param = $results[$i - 1];
            my $num = $self->_pad_string($i, $pad);
            my $msg = "$num:";
            $msg .= ' ' . $self->_pad_string($param->__display_name__, $max_dname_length, 'suffix');
            my $status = ' ';
            if ($param->can('status')) {
                $status = $param->status || 'missing_status';
            }
            $msg .= "\t" . $self->_pad_string($status, $max_status_length, 'suffix');
            $msg .= "\t" . $param->class if (@classes > 1);
            $self->status_message($msg);
        }
        if ($MESSAGE) {
            $MESSAGE = "\n" . '*'x80 . "\n" . $MESSAGE . "\n" . '*'x80 . "\n";
            $self->status_message($MESSAGE);
            $MESSAGE = '';
        }
        my $pretty_values = '(c)ontinue, (h)elp, e(x)it';
        my $valid_values = '\*|c|h|x|[-+]?[\d\-\., ]+';
        if ($caller[3] =~ /_trim_list_from_response/) {
            $pretty_values .= ', (b)ack';
            $valid_values .= '|b';
        }
        $response = $self->_ask_user_question("Please confirm the above items for '$param_name' or modify your selection.", 0, $valid_values, 'h', $pretty_values.', or specify item numbers to use');
        if (lc($response) eq 'h' || !$self->_validate_user_response_for_param_value_verification($response)) {
            $MESSAGE .= "\n" if ($MESSAGE);
            $MESSAGE .=
            "Help:\n".
            "* Specify which elements to keep by listing them, e.g. '1,3,12' would keep\n".
            "  items 1, 3, and 12.\n".
            "* Begin list with a minus to remove elements, e.g. '-1,3,9' would remove\n".
            "  items 1, 3, and 9.\n".
            "* Ranges can be used, e.g. '-11-17, 5' would remove items 11 through 17 and\n".
            "  remove item 5.";
            $response = '';
        }
    }
    if (lc($response) eq 'x') {
        $self->status_message("Exiting...");
        exit;
    }
    elsif (lc($response) eq 'b') {
        return;
    }
    elsif (lc($response) eq 'c' | $response eq '*') {
        return @results;
    }
    elsif ($response =~ /^[-+]?[\d\-\., ]+$/) {
        @results = $self->_trim_list_from_response($response, $param_name, @results);
        return @results;
    }
    else {
        die $self->error_message("Conditional exception, should not have been reached!");
    }
}

sub _validate_user_response_for_param_value_verification {
    my ($self, $response_text) = @_;
    $response_text = substr($response_text, 1) if ($response_text =~ /^[+-]/);
    my @response = split(/[\s\,]/, $response_text);
    for my $response (@response) {
        if ($response =~ /^[xbc*]$/) {
            return 1;
        }
        if ($response !~ /^(\d+)([-\.]+(\d+))?$/) {
            $MESSAGE .= "\n" if ($MESSAGE);
            $MESSAGE .= "ERROR: Invalid list provided ($response)";
            return 0;
        }
        if ($3 && $1 && $3 < $1) {
            $MESSAGE .= "\n" if ($MESSAGE);
            $MESSAGE .= "ERROR: Inverted range provided ($1-$3)";
            return 0;
        }
    }
    return 1;
}

sub _trim_list_from_response {
    my ($self, $response_text, $param_name, @list) = @_;

    my $method;
    if ($response_text =~ /^[+-]/) {
        $method = substr($response_text, 0, 1);
        $response_text = substr($response_text, 1);
    }
    else {
        $method = '+';
    }

    my @response = split(/[\s\,]/, $response_text);
    my %indices;
    @indices{0..$#list} = 0..$#list if ($method eq '-');

    for my $response (@response) {
        $response =~ /^(\d+)([-\.]+(\d+))?$/;
        my $low = $1; $low--;
        my $high = $3 || $1; $high--;
        die if ($high < $low);
        if ($method eq '+') {
            @indices{$low..$high} = $low..$high;
        }
        else {
            delete @indices{$low..$high};
        }
    }
    #$self->debug_message("Indices: " . join(',', sort(keys %indices)));
    my @new_list = $self->_get_user_verification_for_param_value_drilldown($param_name, @list[sort keys %indices]);
    unless (@new_list) {
        @new_list = $self->_get_user_verification_for_param_value_drilldown($param_name, @list);
    }
    return @new_list;
}

sub _pad_string {
    my ($self, $str, $width, $pos) = @_;
    my $padding = $width - length($str);
    $padding = 0 if ($padding < 0);
    if ($pos && $pos eq 'suffix') {
        return $str . ' 'x$padding;
    }
    else {
        return ' 'x$padding . $str;
    }
}

sub _can_interact_with_user {
    my $self = shift;
    if ( -t STDERR ) {
        return 1;
    }
    else {
        return 0;
    }
}

sub _shell_args_property_meta
{
    my $self = shift;
    my $class_meta = $self->__meta__;

    # Find which property metas match the rules.  We have to do it this way
    # because just calling 'get_all_property_metas()' will product multiple matches
    # if a property is overridden in a child class
    my $rule = UR::Object::Property->define_boolexpr(@_);
    my %seen;
    my (@positional,@required,@optional);
    foreach my $property_meta ( $class_meta->get_all_property_metas() ) {
        my $property_name = $property_meta->property_name;

        next if $seen{$property_name}++;
        next unless $rule->evaluate($property_meta);

        next if $property_name eq 'id';
        next if $property_name eq 'result';
        next if $property_name eq 'is_executed';
        next if $property_name eq 'original_command_line';
        next if $property_name =~ /^_/;

        next if $property_meta->implied_by;
        next if $property_meta->is_calculated;
        # Kept commented out from UR's Command.pm, I believe is_output is a workflow property
        # and not something we need to exclude (counter to the old comment below).
        #next if $property_meta->{is_output}; # TODO: This was breaking the G::M::T::Annotate::TranscriptVariants annotator. This should probably still be here but temporarily roll back
        next if $property_meta->is_transient;
        next if $property_meta->is_constant;
        if (($property_meta->is_delegated) || (defined($property_meta->data_type) and $property_meta->data_type =~ /::/)) {
            next unless($self->can('resolve_param_value_from_cmdline_text'));
        }
        else {
            next unless($property_meta->is_mutable);
        }
        if ($property_meta->{shell_args_position}) {
            push @positional, $property_meta;
        }
        elsif ($property_meta->is_optional) {
            push @optional, $property_meta;
        }
        else {
            push @required, $property_meta;
        }
    }

    my @result;
    @result = (
        (sort { $a->property_name cmp $b->property_name } @required),
        (sort { $a->property_name cmp $b->property_name } @optional),
        (sort { $a->{shell_args_position} <=> $b->{shell_args_position} } @positional),
    );

    return @result;
}

sub _missing_parameters {
    my ($self, $params) = @_;

    my $class_meta = $self->__meta__;

    my @property_names;
    if (my $has = $class_meta->{has}) {
        @property_names = $self->_unique_elements(keys %$has);
    }
    my @property_metas = map { $class_meta->property_meta_for_name($_); } @property_names;

    my @missing_property_values;
    for my $property_meta (@property_metas) {
        my $pn = $property_meta->property_name;

        next if $property_meta->is_optional;
        next if $property_meta->implied_by;
        next if defined $property_meta->default_value;
        next if defined $params->{$pn};

        push @missing_property_values, $pn;
    }

    @missing_property_values = map { $_ =~ s/_/-/g; "--$_" } @missing_property_values;
    if (@missing_property_values) {
        $self->status_message('');
        $self->error_message("Missing required parameter(s): " . join(', ', @missing_property_values) . ".");
    }
    return @missing_property_values;
}

sub __errors__ {
    my ($self,@property_names) = @_;

    return (@error_tags, $self->SUPER::__errors__);
}

sub _can_resolve_type {
    my ($self, $type) = @_;

    return 0 unless($type);

    my $non_classes = 0;
    if (ref($type) ne 'ARRAY') {
        $non_classes = $type !~ m/::/;
    } else {
        $non_classes = scalar grep { ! m/::/ } @$type;
    }
    return $non_classes == 0;
}

sub _params_to_resolve {
    my ($self, $params) = @_;
    my @params_to_resolve;
    if ($params) {
        my $cmeta = $self->__meta__;
        my @params_will_require_verification;
        my @params_may_require_verification;

        for my $param_name (keys %$params) {
            next if ($param_name eq 'help');
            my $pmeta = $cmeta->property($param_name);
            unless ($pmeta) {
                # This message was a die after a next, so I guess it isn't supposed to be fatal?
                $self->warning_message("No metadata for property '$param_name'");
                next;
            }

            my $param_type = $pmeta->data_type;
            next unless($self->_can_resolve_type($param_type));

            my $param_arg = $params->{$param_name};
            if (my $arg_type = ref($param_arg)) {
                next if $arg_type eq $param_type; # param is already the right type
                if ($arg_type ne 'ARRAY') {
                    $self->error_message("no handler for property '$param_name' with argument type " . ref($param_arg));
                    next;
                }
            } else {
                $param_arg = [$param_arg];
            }
            next unless (@$param_arg);

            my $resolve_info = {
                name => $param_name,
                class => $param_type,
                value => $param_arg,
            };
            push(@params_to_resolve, $resolve_info);

            my $require_user_verify = $pmeta->{'require_user_verify'};
            if ( defined($require_user_verify) ) {
                push @params_will_require_verification, "'$param_name'" if ($require_user_verify);
            } else {
                push @params_may_require_verification, "'$param_name'";
            }
        }

        my @adverbs = ('will', 'may');
        my @params_adverb_require_verification = (
            \@params_will_require_verification,
            \@params_may_require_verification,
        );
        for (my $i = 0; $i < @adverbs; $i++) {
            my $adverb = $adverbs[$i];
            my @param_adverb_require_verification = @{$params_adverb_require_verification[$i]};
            next unless (@param_adverb_require_verification);

            if (@param_adverb_require_verification > 1) {
                $param_adverb_require_verification[-1] = 'and ' . $param_adverb_require_verification[-1];
            }
            my $param_str = join(', ', @param_adverb_require_verification);
            $self->status_message($param_str . " $adverb require verification...");
        }
    }
    return @params_to_resolve;
}

sub resolve_class_and_params_for_argv {
    my $self = shift;
    my ($class, $params, $error_tag_list) = $self->SUPER::resolve_class_and_params_for_argv(@_);

    if ($self ne $class or ( $error_tag_list and @$error_tag_list)) {
        return ($class, $params, $error_tag_list);
    }
    unless (@_ && scalar($self->_missing_parameters($params)) == 0) {
        return ($class, $params, $error_tag_list);
    }

    local $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;

    my @params_to_resolve = $self->_params_to_resolve($params);
    for my $p (@params_to_resolve) {
        my $param_arg_str = join(',', @{$p->{value}});
        my $pmeta = $self->__meta__->property($p->{name});

        my @params;
        eval {
            @params = $self->resolve_param_value_from_cmdline_text($p);
        };

        if ($@) {
            push @error_tags, UR::Object::Tag->create(
                type => 'invalid',
                properties => [$p->{name}],
                desc => "Errors while resolving from $param_arg_str: $@",
            );
        }
        if (@params and $params[0]) {
            if ($pmeta->{'is_many'}) {
                $params->{$p->{name}} = \@params;
            }
            else {
                $params->{$p->{name}} = $params[0];
            }
        }
        else {
            push @error_tags, UR::Object::Tag->create(
                type => 'invalid',
                properties => [$p->{name}],
                desc => "Problem resolving from $param_arg_str.",
            );
            $self->error_message();
        }
    }

    if (@error_tags) {
        return ($class, undef);
    }
    else {
        return ($class, $params);
    }
}

sub terminal_input_filehandle {
    my $self = shift;

    my $fh = IO::File->new('/dev/tty', 'r');
    unless ($fh) {
        Carp::carp("Couldn't open /dev/tty for terminal input: $!\n    Using STDIN...");
        $fh = *STDIN;
    }
    return $fh;
}


sub _ask_user_question {
    my $self = shift;
    my $question = shift;
    my $timeout = shift;
    my $valid_values = shift || "yes|no";
    my $default_value = shift || undef;
    my $pretty_valid_values = shift || $valid_values;
    $valid_values = lc($valid_values);
    my $input;
    $timeout = 60 unless(defined($timeout));

    local $SIG{ALRM} = sub { print STDERR "Exiting, failed to reply to question '$question' within '$timeout' seconds.\n"; exit; };
    print STDERR "\n$question\n";
    print STDERR "Reply with $pretty_valid_values: ";

    unless ($self->_can_interact_with_user) {
        print STDERR "\n";
        die $self->error_message("Attempting to ask user question but cannot interact with user!");
    }

    my $terminal = $self->terminal_input_filehandle();

    alarm($timeout) if ($timeout);
    chomp($input = $terminal->getline());
    alarm(0) if ($timeout);

    print STDERR "\n";

    if(lc($input) =~ /^$valid_values$/) {
        return lc($input);
    }
    elsif ($default_value) {
        return $default_value;
    }
    else {
        $self->error_message("'$input' is an invalid answer to question '$question'\n\n");
        return;
    }
}

sub _unique_elements {
    my ($self, @list) = @_;
    my %seen = ();
    my @unique = grep { ! $seen{$_} ++ } @list;
    return @unique;
}

sub display_command_summary_report {
    my $self = shift;

    return unless $self->show_display_command_summary_report;

    my $total_count = $self->_total_command_count;
    my %command_errors = %{$self->_command_errors};

    if (keys %command_errors) {
        $self->status_message("\n\nErrors Summary:");
        for my $key (keys %command_errors) {
            my $errors = $command_errors{$key};
            $errors = [$errors] unless (ref($errors) and ref($errors) eq 'ARRAY');
            my @errors = @{$errors};
            printf "%s: \n" , $self->_color($key, 'red');
            for my $error (@errors) {
                $error = $self->truncate_error_message($error);
                print "\t- $error\n";
            }
        }
    }

    if ($total_count > 1) {
        my $error_count = scalar(keys %command_errors);

        my $successful_count = ($total_count - $error_count);
        my $success_msg = " Successful: $successful_count";
        my $error_msg   = "     Errors: $error_count";
        my $total_msg   = "      Total: $total_count";

        $success_msg = $self->_color($success_msg, 'green') if $successful_count;
        $error_msg   = $self->_color($error_msg, 'red') if $error_count;

        $self->status_message("\n\nCommand Summary:");
        $self->status_message($success_msg);
        $self->status_message($error_msg);
        $self->status_message($total_msg);
    }
}

sub append_error {
    my $self = shift;
    my $key = shift || die;
    my $error = shift || die;

    my $command_errors = $self->_command_errors;
    push @{$command_errors->{$key}}, $error;
    $self->_command_errors($command_errors);

    return 1;
}

sub truncate_error_message {
    my $self = shift;
    my $error = shift || die;

    # truncate errors so they are actually a summary
    ($error) = split("\n", $error);

    # meant to truncate a callstack as this is meant for user/high-level
    $error =~ s/\ at\ \/.*//;

    return $error;
}

1;
