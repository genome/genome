package Finfo::OldValidate;

use strict;
use warnings;
no warnings 'reserved';

use base qw(Class::Accessor);

use Carp;
use Data::Dumper;
use File::Basename;
use Regexp::Common;
use Scalar::Util 'blessed';

__PACKAGE__->mk_ro_accessors
(qw/
    err_cb type attr value ref_type options custom_validation caller 
    /);

sub new { validate(@_) }

sub validate
{
    my ($class, %p) = @_;

    $p{'caller'} = [ caller ] unless $p{'caller'};

    my $self = bless \%p, $class;

    $self->_communicate_error("No attr given to validate")
        and return unless defined $self->attr;
 
    $self->_communicate_error( sprintf("No type given to validate attribute (%s)", $self->attr) )
        and return unless defined $self->type;
 
    unless ( UNIVERSAL::can($self, $self->type) )
    {
        $self->_communicate_error(sprintf('Invalid validation type (%s)', $self->type) );
        return;
    }
    
    if ( $self->{obj} )
    {
        warn "obj param is deprecated, please set obj to err_cb\n" unless $class::_OBJ_WARNED;
        $class::_Obj_WARNED = 1;
        $self->{err_cb} = delete $self->{obj};
    }

    if ( $self->{valid_options} )
    {
        warn "valid_options param is deprecated, please use options\n" unless $class::_VALID_OPTS_WARNED;
        $class::_VALID_OPTS_WARNED = 1;
        $self->{options} = delete $self->{valid_options};
    }

    if ( $self->options or $self->type_requires_options( $self->type ) )
    {
        return unless $self->_validate_options;
    }
    
    my $method = $self->_test_type( $self->type );

    return unless $method;

    return $self->$method;
}

sub _test_type
{
    my ($self, $type) = @_;

    return unless $type;

    return if $type =~ /^_/;
    
    return UNIVERSAL::can($self, $type);
}

sub type_requires_options
{
    my ($self, $type) = @_;
    
    return grep { $type eq $_ } 
    (qw/
        int_gt int_gte int_lt int_lte int_between
        integer_gt integer_gte integer_lt integer_lte integer_between
        real_gt real_gte real_lt real_lte real_between
        float_gt float_gte float_lt float_lte float_between
        inherits_from inherits_from_aryref
        in_list in_list_aryref
        /);
}

sub _communicate_error
{
    my ($self, $error) = @_;

    confess "No error defined to communicate\n" unless defined $error;
    
    if ( not defined $self->err_cb )
    {
        my $caller = $self->caller;
        confess sprintf
        (
            "\n[ FATAL ] %s\n\tin package %s at line", 
            $error, 
            $caller->[0],
            $caller->[2],
        );
    }
    elsif ( ref $self->err_cb eq 'CODE' )
    {
        $self->err_cb->($error);
    }
    elsif ( blessed($self->err_cb) and $self->err_cb->can('fatal_msg'))
    {
        $self->err_cb->fatal_msg
        (
            $error,
            { 'caller' => $self->caller },
        );
    }
    else
    {
        confess sprintf
        (
            "\n[ FATAL ] Unsupported err_cb (%s)\n[ FATAL] %s\n", 
            $self->err_cb, 
            $error
        );
    }

    return 1;
}

sub _validate_options
{
    my $self = shift;

    $self->_communicate_error
    (
        sprintf
        (
            "Options are required for attribute (%s) type (%s)", 
            $self->attr,
            $self->type,
        )
    )
        and return unless $self->options;

    $self->_communicate_error( sprintf("Options for %s is not an aryref", $self->attr) )
        and return unless ref($self->options) eq 'ARRAY';

    $self->_communicate_error( sprintf("Options for %s is an empty aryref", $self->attr) )
        and return unless @{ $self->options };

    return 1;
}

sub custom
{
    my $self = shift;

    return unless $self->code("Custom validation code", $self->custom_validation);

    return $self->custom_validation(@_);
}

sub defined
{
    my $self = shift;

    $self->_communicate_error($self->attr . ' is not defined')
        and return unless defined $self->value;
    
    return 1;
}

sub string
{
    my $self = shift;

    return unless $self->defined;

    $self->_communicate_error($self->attr . ' is not a string')
        and return if ref($self->value); 

    if ( $self->options )
    {
        my ($type, @opts) = @{ $self->options };
        __PACKAGE__->validate
        (
            attr => "length of " . $self->attr,
            value => length($self->value),
            type => $type,
            options => ( @opts ) ? \@opts : undef,
            err_cb => $self->err_cb,
        );
    }
    
    return 1;
}

# Numbers
# integer
*integer = sub{ shift->int };
# positive, non negative, negative
*pos_int = sub{ shift->_pos('int') };
*pos_integer = \&pos_int;
*positive_integer = \&pos_int;
*non_neg_int = sub{ shift->_non_neg('int') };
*non_neg_integer = \&non_neg_int;
*non_negative_integer = \&non_neg_int;
*neg_int = sub{ shift->_neg('int') };
*neg_integer = \&neg_int;
*negitive_integer = \&neg_int;

# real
*float = sub{ shift->real };
# positive, non negative, negative
*pos_real = sub{ shift->_pos('real') };
*pos_float = \&pos_real;
*positive_real = \&pos_real;
*positive_float = \&pos_real;
*non_neg_real = sub{ shift->_non_neg('real') };
*non_neg_float = \&non_neg_real;
*non_negative_real = \&non_neg_real;
*non_negative_float = \&non_neg_real;
*neg_real = sub{ shift->_neg('real') };
*neg_float = \&neg_real;
*negitive_real = \&neg_real;
*negitive_float = \&neg_real;

# greater than
*int_gt = sub{ shift->_gt('int') };
*integer_gt = \&int_gt;
*real_gt = sub{ shift->_gt('real') };
*float_gt = \&real_gt;
# greater than or equal to
*int_gte = sub{ shift->_gte('int') };
*integer_gte = \&int_gte;
*real_gte = sub{ shift->_gte('real') };
*float_gte = \&real_gte;

# less than
*int_lt = sub{ shift->_lt('int') };
*integer_lt = \&int_lt;
*real_lt = sub{ shift->_lt('real') };
*float_lt = \&real_lt;
# less than or equal to
*int_lte = sub{ shift->_lte('int') };
*integer_lte = \&int_lte;
*real_lte = sub{ shift->_lte('real') };
*float_lte = \&real_lte;

# between
*int_between = sub{ shift->_between('int') };
*integer_between = \&int_between;
*real_between = sub{ shift->_between('real') };
*float_between = \&real_between;

sub int
{
    my $self = shift;

    return unless $self->defined;

    $self->_communicate_error
    (
        sprintf('%s (%s) needs to be an integer', $self->attr, $self->value) 
    )
        and return unless $self->value =~ /^$RE{num}{int}$/;

    return 1;
}

sub real
{
    my $self = shift;

    return unless $self->defined;

    $self->_communicate_error
    (
        sprintf("%s (%s) must be a real number", $self->attr, $self->value)
    ) 
        and return unless $self->value =~ /^$RE{num}{real}$/;

    return 1;
}

sub _pos
{
    my ($self, $type) = @_;

    return unless $self->$type;

    $self->_communicate_error
    (
        sprintf('%s (%s) needs to be a positive %s', $self->attr, $self->value, $type) 
    )
        and return unless $self->value > 0;
    
    return 1;
}

sub _non_neg
{
    my ($self, $type) = @_;

    return unless $self->$type;

    $self->_communicate_error
    (
        sprintf('%s (%s) needs to be a non negative %s', $self->attr, $self->value, $type) 
    )
        and return unless $self->value >= 0;
    
    return 1;
}

sub _neg
{
    my ($self, $type) = @_;

    return unless $self->$type;

    $self->_communicate_error
    (
        sprintf('%s (%s) needs to be a negative %s', $self->attr, $self->value, $type) 
    )
        and return unless $self->value < 0;
    
    return 1;
}

sub _gt
{
    my ($self, $type) = @_;

    return unless $self->$type;

    my ($min) = $self->options->[0];
    $self->_communicate_error($self->attr . " needs a min in options to check $type gt")
        and return unless defined $min and $min =~ /$RE{num}{$type}/;

    $self->_communicate_error
    (
        sprintf('%s (%s) needs to be greater than %s', $self->attr, $self->value, $min) 
    )
        and return unless $self->value > $min;

    return 1;
}

sub _gte
{
    my ($self, $type) = @_;

    return unless $self->$type;

    my ($min) = $self->options->[0];
    $self->_communicate_error($self->attr . " needs a min in options to check $type gte")
        and return unless defined $min and $min =~ /$RE{num}{$type}/;

    $self->_communicate_error
    (
        sprintf('%s (%s) needs to be greater than or equal to %s', $self->attr, $self->value, $min) 
    )
        and return unless $self->value >= $min;

    return 1;
}

sub _lt
{
    my ($self, $type) = @_;

    return unless $self->$type;

    my ($max) = $self->options->[0];
    $self->_communicate_error($self->attr . " needs a max in options to check $type between")
        and return unless defined $max and $max =~ /$RE{num}{$type}/;

    $self->_communicate_error
    (
        sprintf('%s (%s) needs to be less than %s', $self->attr, $self->value, $max) 
    )
        and return unless $self->value < $max;

    return 1;
}

sub _lte
{
    my ($self, $type) = @_;

    return unless $self->$type;

    my ($max) = $self->options->[0];
    $self->_communicate_error($self->attr . " needs a max in options to check $type between")
        and return unless defined $max and $max =~ /$RE{num}{$type}/;
    
    $self->_communicate_error
    (
        sprintf('%s (%s) needs to be less than or equal to %s', $self->attr, $self->value, $max) 
    )
        and return unless $self->value <= $max;

    return 1;
}

sub _between
{
    my ($self, $type) = @_;

    return unless $self->$type;

    my ($min, $max) = @{ $self->options };

    $self->_communicate_error
    (
        $self->attr . " needs a min in options to check $type between"
    )
        and return unless defined $min and $min =~ /$RE{num}{$type}/;

    $self->_communicate_error
    (
        $self->attr . " needs a max in options to check $type between"
    )
        and return unless defined $max and $max =~ /$RE{num}{$type}/;
    

    $self->_communicate_error
    (
        $self->attr . " max is less than min for $type between"
    )
        and return unless $max > $min;
    
    $self->_communicate_error
    (
        sprintf('%s (%s) needs to be between %s and %s', $self->attr, $self->value, $min, $max) 
    )
        and return unless $self->value >= $min and $self->value <= $max;

    return 1;
}

# ARYREF
sub aryref
{  
    my $self = shift;

    return unless $self->defined;

    $self->_communicate_error($self->attr . " is not an array reference")
        and return unless ref($self->value) eq 'ARRAY';
    
    return $self->_check_aryref_against_ref_type;
}

sub non_empty_aryref
{
    my $self = shift;

    return unless $self->aryref;

    $self->_communicate_error($self->attr . " is an empty array reference")
        and return unless @{ $self->value };
    
    return $self->_check_aryref_against_ref_type;
}

# HASHREF
sub hashref
{
    my $self = shift;

    return unless $self->defined;

    $self->_communicate_error($self->attr . " is not an hash reference")
        and return unless ref($self->value) eq 'HASH';

    return 1;
}

sub non_empty_hashref
{
    my $self = shift;

    return unless $self->hashref;

    $self->_communicate_error($self->attr . " is an empty hash reference")
        and return unless %{ $self->value };

    return $self->_check_hashref_against_ref_type;
}

# REF_TYPE
sub _test_ref_type
{
    my $self = shift;

    return 1 unless $self->ref_type;
    
    return $self->test_type( $self->ref_type );
}

sub _check_aryref_against_ref_type
{
    my $self = shift;

    return 1 unless $self->ref_type;
    
    unless ( $self->_test_type($self->ref_type) )
    {
        $self->_communicate_error("Invalid ref type for " . $self->attr);
        return;
    }

    my $aryref = $self->value;
    foreach my $val ( @$aryref )
    {
        return unless __PACKAGE__->validate
        (
            attr => sprintf('value in aryref (%s)', $self->attr),
            value => $val,
            type => $self->ref_type,
            options => $self->options,
            err_cb => $self->err_cb,
        );
    }

    return 1;
}

sub _check_hashref_against_ref_type
{
    my $self = shift;

    return 1 unless $self->ref_type;
    
    unless ( $self->_test_type($self->ref_type) )
    {
        $self->_communicate_error("Invalid ref type for " . $self->attr);
        return;
    }

    my $hashref = $self->value;
    foreach my $key ( keys %$hashref )
    {
        return unless __PACKAGE__->validate
        (
            attr => sprintf('%s with key %s', $self->attr, $key),
            value => $hashref->{$key},
            type => $self->ref_type,
            options => $self->options,
            ref_type => $self->ref_type,
            err_cb => $self->err_cb,
        );
    }

    return 1;
}

# FILE
sub file
{
    my $self = shift;

    return unless $self->defined;

    $self->_communicate_error
    (
        sprintf("%s (%s) is a directory (should be a file)", $self->attr, $self->value)
    )
        and return if -d $self->value;

    $self->_communicate_error
    (
        sprintf("%s (%s) does not exist", $self->attr, $self->value)
    )
        and return unless -e $self->value;
 
    return 1;
}

sub input_file
{
    my $self = shift;

    return unless $self->file;

    $self->_communicate_error
    (
        sprintf("%s (%s) is empty or does not exist", $self->attr, $self->value)
    )
        and return unless -s $self->value;
    
    $self->_communicate_error
    (
        sprintf("%s (%s) is not readable", $self->attr, $self->value)
    )
        and return unless -r $self->value;
    
    return 1;
}

sub output_file
{
    my $self = shift;
    
    return unless $self->defined;

    $self->_communicate_error( sprintf("%s (%s) exists", $self->attr, $self->value) )
        and return if -e $self->value;
    
    my ($file, $dir) = fileparse($self->value);
    
    $self->_communicate_error("Directory ($dir) for file ($file) is not writable")
        and return unless -w $dir;
 
    return 1;
}

sub rw_file
{
    my $self = shift;

    return unless $self->defined;

    $self->_communicate_error
    (
        sprintf("%s (%s) is a directory (should be a file)", $self->attr, $self->value)
    )
        and return if -d $self->value;

    my ($file, $dir) = fileparse($self->value);
    
    $self->_communicate_error("Directory ($dir) for file ($file) is not readable")
        and return unless -r $dir;

    $self->_communicate_error("Directory ($dir) for file ($file) is not writable")
        and return unless -w $dir;
 
    return 1;

    
    return 1;
}

*_executable = \&exe;
sub exe
{
    my ($self, $attr, $exe) = @_;

    return unless $self->input_file;
    
    $self->communicate_error( sprintf("File (%s) is not executable", $self->attr) )
        and return unless -x $self->value;
    
    return 1;
}

# PATH/DIR
*dir = \&path;
*input_dir = \&input_path;
*output_dir = \&output_path;
*io_dir = \&io_path;

sub path
{
    my $self = shift;

    $self->_communicate_error
    (
        sprintf('%s (%s) does not exist', $self->attr, $self->value)
    )
        and return unless -e $self->value;
    
    $self->_communicate_error
    (
        sprintf('%s (%s) is not a directory', $self->attr, $self->value)
    )
        and return unless -d $self->value;
    
    return 1;
}

sub input_path
{
    my $self = shift;

    return unless $self->path;

    $self->_communicate_error
    (
        sprintf('Can\'t read from %s (%s)', $self->attr, $self->value)
    )
        and return unless -r $self->value;

    return 1;
}

sub output_path
{
    my $self = shift;

    return unless $self->path;

    $self->_communicate_error
    (
        sprintf('Can\'t write to %s (%s)', $self->attr, $self->value)
    )
        and return unless -r $self->value;

    return 1;
}

sub io_path
{
    my $self = shift;
    
    return unless $self->input_path;
    
    return unless $self->output_path;

    return 1;
}

# LIST
sub in_list
{
    my $self = shift;

    return unless $self->defined;

    $self->_communicate_error
    (
        sprintf
        (
            "%s (%s) is not in the list of options (%s)", 
            $self->attr,
            $self->value,
            join(', ', @{ $self->options })
        )
    )
        and return unless grep { $self->value eq $_ } @{ $self->options };

    return 1;
}

sub in_list_aryref
{
    my $self = shift;

    return unless $self->non_empty_aryref;

    foreach my $value ( @{ $self->value } )
    {
        $self->_communicate_error
        (
            sprintf
            (
                "%s (%s) is not in the list of options (%s)", 
                $self->attr,
                $value,
                join(', ', @{ $self->options })
            )
        )
            and return unless grep { $value eq $_ } @{ $self->options };
    }

    return 1;
}

# OBJECTS/INHERITANCE
sub object
{
    my $self = shift;

    return unless $self->defined;

    my $class = blessed($self->value);

    $self->_communicate_error($self->attr . " is not a blessed reference")
        and return unless defined $class;

    return $class;
}

sub inherits_from
{
    my $self = shift;

    my $class = $self->object;

    return unless defined $class;

    $self->_communicate_error
    (
        "$class does not inherit from " . join(', ', @{ $self->options } )
    )
        and return unless grep { UNIVERSAL::isa($self->value, $_) } @{ $self->options };

    return 1;
}

sub inherits_from_aryref
{
    my $self = shift;

    return unless $self->non_empty_aryref;

    foreach my $obj ( @{ $self->value } )
    {
        my $class = blessed($obj);

        $self->_communicate_error($self->attr . " is not a blessed reference")
            and return unless defined $class;

        $self->_communicate_error
        (
            "$class does not inherit from " . join(', ', @{ $self->options } )
        )
            and return unless grep { UNIVERSAL::isa($obj, $_) } @{ $self->options };
    }

    return 1;
}

# REGEX
sub regex
{
    my $self = shift;

    return unless $self->defined;

    my $val = $self->value;

    $self->_communicate_error
    (
        sprintf
        (
            "%s (%s) does not match the given patterns (%s)",
            $self->attr,
            $val,
            join(', ', @{ $self->options }),
        )
    )
        and return unless grep { $val =~ /$_/i } @{ $self->options };

    return 1;
}

sub in_list_regex
{
    warn "in_list_regex validation type is deprecated, please use regex\n";
    return shift->regex(@_);
}

sub regex_aryref
{
    my $self = shift;

    return unless $self->non_empty_aryref;

    foreach my $value ( @{ $self->value } )
    {
        $self->_commucate_error($self->attr . " has an undefined value in the aryref")
            and return unless defined $value;

        $self->_communicate_error
        (
            sprintf
            (
                "%s (%s) does not match the given patterns (%s)",
                $self->attr,
                $value,
                join(', ', @{ $self->options }),
            )
        )
            and return unless grep { $value =~ /$_/i } @{ $self->options };
    }

    return 1;
}

sub in_list_regex_aryref
{
    warn "in_list_regex_aryref validation type is deprecated, please use regex_aryref\n";
    return shift->regex_aryref(@_);
}

sub code
{
    my $self = shift;

    return unless $self->defined;

    $self->_communicate_error($self->attr . " is not a code ref")
        and return unless ref($self->value) eq 'CODE';
    
    return 1;
}

sub is_validation_type
{
    my $self = shift;

    return 1 if $self->_test_type( $self->value );

    $self->_communicate_error( sprintf('Invalid %s (%s)', $self->attr, $self->value) );
    
    return;
}

sub y_or_n
{
    my $self = shift;

    return unless $self->defined;

    $self->_communicate_error($self->attr . " needs to be [y]es or [n]o")
        and return unless $self->value =~ /^[ny]/i;
    
    return 1;
}

sub not_blank
{
    my $self = shift;

    return unless $self->defined;

    $self->_communicate_error($self->attr . " is blank")
        and return if $self->value eq '';
    
    return 1;
}

sub project_name
{
    my $self = shift;

    return unless $self->defined;

    use DBI;
    my $dbh = DBI->connect
    (
        'dbi:Oracle:gscprod',
        'gscuser',
        'g_user',
        { AutoCommit => 0, RaiseError => 1 },
    );
    
    $self->_communicate_error("Can't connect: $DBI::errstr")
        and return unless $dbh;

    my $name = $self->value;

    my $sth = $dbh->prepare("select * from projects where name = '$name'");
    $sth->execute;
    my $aryref = $sth->fetchall_arrayref;
    $dbh->disconnect;
    
    $self->_communicate_error( sprintf("Invalid %s (%s)", $self->attr, $name) )
        and return unless @$aryref;
    
    my %data;
    @data{ @{ $sth->{NAME_lc} } } = @{ $aryref->[0] };

    return \%data;
}
    
1;

=pod

=head1 Name

 Finfo::Validate

=head1 Synopsis

 Package with generic validation methods.  Evaluates one variable at a time, 
  executing a callback on failure.  First tries to execute the error callback (err_cb),
  then if an object is given, an attempt will be made to set the obj's 
  error_msg.  If these are not successful, the error will be confess'd.
 
=head1 Usage

 Required
  attr              The name of the value, used in error communication
  type              Type of validation to use
  
 Optional
  value             Value to be checked ( may be undef to check if defined )
  object            The object to set the error_msg on failure
  err_cb            Executed with error msg as the first param on failure
  ref_type          string, checks the values of an ary of hash ref against this type
  options           Aryref of valid values or other options, if applicable
  custom_validation Code ref that checks the incoming value for the attribute
  
 Examples 
 
 True:
 
 Finfo::Validate->validtate
 (
    obj => $self,
    attr => 'answer',
    value => 'yep',
    type => 'in_list',  # should use type y_or_n for this
    options => [qw/ yes yeah yo yep /], 
 )
    or $self->gui_popup( $self->error_msg); # would not happen

 False, sets error_msg in $self:

 Finfo::Validate->validtate
 (
    obj => $self,
    attr => 'answer',
    value => 'cool',
    type => 'in_list',
    options => [qw/ yes yeah yo yep /], 
  ) 
    or $self->gui_popup( $self->error_msg);

=head1 Methods

=head2 validate, new

See above, checks value, returns 1 if ok, undef and executes error callbacks.  If invalid params given, will excute error callbacks.

=head2 test_type

 Finfo::Validate->test_type($validation_type)
   or die;

Tests if the $validation_type is a valid type.  Returns undef on failure.

=head1 Validation Types

=head2 defined
 
Checks if value is defined.

=head2 string

Checks if value is defined and not a reference

=head2 int, integer

Checks if value is defined and an integer

=head2 pos_int, positive_integer

Checks if value is defined, integer and greater than 0

=head2 non_neg_int, non_negative_integer

Checks if value is defined, an integer and greater than or equal to 0
 
=head2 int_gt, integer_gt

Checks if value is defined, an integer and greater than the first options

=head2 int_gte, integer_gte

Checks if value is defined, an integer and greater than or equal to the first options

=head2 int_lt, integer_lt

Checks if value is defined, an integer and less than the first options

=head2 int_lte, integer_lte

Checks if value is defined, an integer and less than or equal to the first options

=head2 int_between, integer_between

Checks if value is defined, an integer, greater than or equal to the first options
 and less than or equal to the second options

=head2 float

Checks if value is defined and an float ( an integer is a float ) 

=head2 positive_float

Checks if value is defined, an float number and greater than 0
 
=head2 non_negative_float

Checks if value is defined, an float number and less than than 0

=head2 float_gt

Checks if value is defined, an float and greater than the first options

=head2 float_gte

Checks if value is defined, an float and greater than or equal to the first options

=head2 float_lt

Checks if value is defined, an float and less than the first options

=head2 float_lte

Checks if value is defined, an float and less than or equal to the first options

=head2 float_between

Checks if value is defined, an float, greater than or equal to the first options
 and less than or equal to the second options

=head2 aryref

Checks if value is defined and a reference to an array
 
=head2 non_empty_aryref

Checks if value is defined, a reference to an array and that is has values, if ref_type is defined,
 will check each value in the aryref against this type
 
=head2 hashref

Checks if value is defined and a reference to an hash, if ref_type is defined,
  will check each value in the hashref against this type
 
=head2 non_empty_hashref

Checks if value is defined, a reference to an hash and that is has values
 
=head2 file

Checks if the file is defined and not a directory

=head2 input_file

Checks if the file is defined, exists, has size and is readable
 
=head2 output_file

Checks if the file is defined, does not exist and that the directory is writable
  
=head2 rw_file

Checks if the file is defined, not a directory, and that the directory is readable/writable

=head2 path

Checks if path is defined and a directory
 
=head2 input_path

Checks if path is defined, a directory and is readable

=head2 output_path

Checks if path is defined, a directory and is writable

=head2 io_path

Checks if a path is defined, a directory, readable and writable
 
=head2 dir

Checks if directory is defined and a directory
 
=head2 input_dir

Checks if directory is defined, a directory and is readable

=head2 output_dir

Checks if directory is defined, a directory and is writable

=head2 io_dir

Checks if a directory is defined, a directory, readable and writable

=head2 executable

Checks if a file is defiend, readable and executable

=head2 object

Checks if a value is defined and that it is a blessed reference

=head2 inherits_from

Checks if a value is defined an objec and that the object inherits from 
 the list of options

=head2 inherits_from_aryref

Checks if a value is defined, an aryref of objects and checks that each object inherits from 
 the list of options

=head2 in_list

Checks if a value is defined and in the list of options

=head2 in_list_aryref

Checks if an aryref of values is defined and in the list of options

=head2 regex

Checks if a value is defined and in the list of options,
 using regex search for each option

=head2 regex_aryref

Checks if an aryref of values is defined and in the list of options,
  using regex search for each option

=head2 custom
 
Use a code reference as a custom validation type, by setting the incoming
parameter 'custom_validation' to a code reference.  The paramters that 
are passed in are the 'attr' and 'value.'  Return true for a successful
validation.

 Ex:

 Finfo::Validate->validate
 (
   attr => 'sum of x^2 + y^2',
   value => 4,
   type => 'custom',
   custom_validation => sub{ my ($attr, $value) = @_; return 1 if $value == 4 },
  )
   or die;
 
=head2 code
 
Checks is a value is defined and a reference to code
 
=head2 is_validation_type

Checks is a value is defined and is a validation method itself

=head2 y_or_n

Checks is a value is defined and starts with a Yy or Nn

=head2 project_name

Checks is a value is defined and is a GSC project name.  It returns a hashref of 
 the project info.

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details.

=head1 Author(s)

Eddie Belter <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
