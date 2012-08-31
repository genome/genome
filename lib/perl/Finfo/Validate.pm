#!/usr/bin/env genome-perl


package Finfo::Validate;

use strict;
use warnings;
no warnings 'reserved';

use base qw(Class::Accessor);

use Carp;
use Data::Dumper;
use File::Basename;
use Finfo::ClassUtils 'class';
use Finfo::Msg;
use Finfo::OldValidate;
use Regexp::Common;
use Scalar::Util 'blessed';

require Term::ANSIColor;

__PACKAGE__->mk_accessors
(qw/
    attr
    value
    ds
    is_a
    attr_opts
    empty_ok
    caller
    msg
    cb
    obj
    /);

my @data_structures = (qw/ scalar aryref hashref /);
my @isas = 
(qw/
    defined not_blank string boolean is_isa
    int real object in_list regex y_or_n
    code project_name
    file file_r file_w file_rw exe
    dir dir_r dir_w dir_rw 
    /);

sub new { validate(@_); }

sub validate
{
    my ($class, %p) = @_;

    my @caller = caller;
    if ( $p{type} )
    {
        no strict 'refs';
        warn_msg("API has changed, please update $caller[0]") unless ${ $class . '::' . $caller[0] };
        ${ $class . '::' . $caller[0] } = 1;
        return Finfo::OldValidate->validate(%p, caller => \@caller);
    }
    
    my $original_isa = delete $p{isa} || $isas[0];
    $p{ds} ||= $data_structures[0];
    $p{empty_ok} ||= 0;

    my $self = bless \%p, __PACKAGE__;

    $self->caller(\@caller) unless $self->caller;

    $self->_communicate_error("No attr given to validate")
        and return unless defined $self->attr;

    my ($isa, @attr_opts) = $self->_is_isa("Isa for " . $self->attr, $original_isa);
    $self->is_a($isa);
    $self->attr_opts(\@attr_opts) if @attr_opts;
    
    $self->_is_ds("Data structure for " . $self->attr, $self->ds);

    my $ds_method = '_validate_' . $self->ds;

    return $self->$ds_method;
}

sub is_ds
{
    my ($class, %p) = @_;

    my $self = bless \%p, __PACKAGE__;

    $self->caller([ caller ]) unless $self->caller;

    $self->_communicate_error("No attr given to validate")
        and return unless defined $self->attr;

    return $self->_is_ds($self->attr, $self->value);
}

sub _is_ds
{
    my ($self, $attr, $ds) = @_;

    return unless $self->_defined($attr, $ds);

    $self->_communicate_error("$attr is an invalid data structure ($ds)")
        and return unless grep { $ds eq $_ } @data_structures;
    
    return 1;
}

sub is_isa
{
    my ($class, %p) = @_;

    my $self = bless \%p, __PACKAGE__;

    $self->caller([ caller ]) unless $self->caller;

    $self->_communicate_error("No attribute given to validate")
        and return unless defined $self->attr;

    return $self->_is_isa($self->attr, $self->value);
}

sub _is_isa
{
    my ($self, $attr, $value) = @_;

    return unless $self->_defined($attr, $value);

    my $isa_ref = ref($value);

    my ($isa, @opts);
    if ( $isa_ref )
    {
        if ( $isa_ref eq 'ARRAY')
        {
            ($isa, @opts) = @$value;
        }
        else
        {
            $self->_commicate_error("Invalid isa: " . Dumper($value));
            return;
        }
    }
    else
    {
        ($isa, @opts) = split(/\s+/, $value);
    }
    
    $self->_communicate_error("Invalid isa ($isa) for attribute ($attr)")
        and return unless grep { $isa eq $_ } @isas;
    
    return ( wantarray ) ? ($isa, @opts) : 1;
}

sub _validate_scalar
{
    my $self = shift;

    my $attr = $self->attr;
    my $value = $self->value;
    my $isa = $self->is_a;

    my $value_ref = ref($value);
    if ( not $isa eq 'defined' and grep { $value_ref eq $_ } (qw/ HASH ARRAY /) )
    {
        $self->_communicate_error
        (
            "Data structure ($value_ref) indicated when trying to validate scalar for attribute ($attr)"
        );
        return;
    }

    return $self->_check_is_a($attr, $value);
}

sub _validate_aryref
{  
    my $self = shift;

    my $attr = $self->attr;
    my $value = $self->value;

    return unless $self->_aryref($attr, $value);

    unless ( $self->empty_ok )
    {
        $self->_communicate_error("$attr is an empty array reference")
            and return unless @$value;
        $self->_check_is_a("Value in $attr ", $value);
    }

    return 1;
}

sub _validate_hashref
{
    my $self = shift;

    my $attr = $self->attr;
    my $value = $self->value;
    
    return unless $self->_hashref($attr, $value);

    unless ( $self->empty_ok )
    {
        $self->_communicate_error("$attr is an empty hash reference")
            and return unless %$value;
        foreach my $key ( keys %$value )
        {
            $self->_check_is_a("Key ($key) in $attr", $value->{$key});
        }
    }

    return 1;
}

sub _check_is_a
{
    my ($self, $attr, $value) = @_;
    
    foreach my $value ( ( ref($value) eq 'ARRAY' ) ? @$value : $value )
    {
        my $method = '_' . $self->is_a;
        return unless $self->$method($attr, $value, $self->attr_opts);
    }

    return 1;
}

sub _communicate_error
{
    my ($self, $error) = @_;

    unless ( defined $error )
    {
        fatal_msg("No error defined to communicate", [ caller ]);
    }

    my $cb = $self->cb;
    if ( $cb )
    {
        if ( ref($cb) eq 'CODE' )
        {
            $cb->($error, { 'caller' => $self->caller });
            return 1;
        }
        else
        { 
            fatal_msg("Callback is not a code ref", { 'caller' => $self->caller });
        }
    }
    
    my $msg_level = $self->msg;
    my $obj = $self->obj;
    if ( $msg_level or $obj )
    {
        my $caller = $self->caller;
        my $pkg_to_call_msg_method = ( $obj ) 
        ? $obj
        : $caller->[0];

        $msg_level = ( $msg_level ) ? $msg_level : 'fatal';
        unless ( grep { $msg_level eq $_ } Finfo::Msg->levels )
        {
            fatal_msg("Invalid msg level ($msg_level)", $caller);
        }

        my $msg_method = $msg_level . '_msg';
        $pkg_to_call_msg_method->$msg_method($error, { 'caller' => $caller });

        return 1;
    }

    fatal_msg($error, $self->caller);
}

sub warn_msg
{
    my ($msg) = @_;

    my $msg_obj = Finfo::Msg->new
    (
        msg => $msg,
        level => 'warn',
    );

    print STDERR $msg_obj->string,"\n";

    return 1;
}

sub fatal_msg
{ 
    my ($msg, $caller) = @_;

    my $msg_obj = Finfo::Msg->new
    (
        msg => $msg,
        level => 'fatal',
        detail_level => 'more',
        'caller' => $caller,
    );

    print STDERR $msg_obj->string,"\n";

    croak;
}

sub _aryref
{
    my ($self, $attr, $value) = @_;

    return unless $self->_defined($attr, $value);

    $self->_communicate_error("$attr is not an array reference")
        and return unless ref($value) eq 'ARRAY';
    
    return 1;
}

sub _hashref
{
    my ($self, $attr, $value) = @_;

    return unless $self->_defined($attr, $value);

    $self->_communicate_error("$attr is not an hash reference")
        and return unless ref($value) eq 'HASH';

    return 1;
}

sub _defined
{
    my ($self, $attr, $value) = @_;

    $self->_communicate_error("$attr is not defined")
        and return unless defined $value;
    
    return 1;
}

sub _boolean
{
    my ($self, $attr, $value) = @_;

    return unless $self->_defined($attr, $value);

    $self->_communicate_error("$attr must be a boolean value ($value), 0 or 1)")
        and return unless grep { $value eq $_ } (qw/ 0 1 /);
    
    return 1;
}

sub _string
{
    my ($self, $attr, $value, $opts) = @_;

    return unless $self->_defined($attr, $value);

    $self->_communicate_error("$attr is not a string")
        and return if ref($value); 

    return unless $self->_validate_numeric_range('int', "length of $attr", length($value), $opts);

    return 1;
}

sub _int
{
    my ($self, $attr, $value, $opts) = @_;

    return unless $self->_defined($attr, $value);

    $self->_communicate_error("$attr ($value) needs to be an integer")
        and return unless $value =~ /^$RE{num}{int}$/;

    return unless $self->_validate_numeric_range('int', $attr, $value, $opts);

    return 1;
}

sub _real
{
    my ($self, $attr, $value, $opts) = @_;

    return unless $self->_defined($attr, $value);

    $self->_communicate_error("$attr ($value) must be a real number")
        and return unless $value =~ /^$RE{num}{real}$/;

    return unless $self->_validate_numeric_range('real', $attr, $value, $opts);

    return 1;
}

sub _validate_numeric_range
{
    my ($self, $type, $attr, $value, $opts) = @_;

    return 1 unless $opts;

    my ($function, $num1, $num2) = @$opts;

    $self->_communicate_error("No function to check numeric range for attr ($attr)")
        and return unless $function;
    
    my %functions_methods = 
    (
        'pos' => '_pos',
        'neg' => '_neg',
        'non_neg' => '_non_neg',
        'gt' => '_gt',
        '>' => '_gt',
        'gte' => '_gte',
        '>=' => '_gte',
        'lt' => '_lt',
        '<' => '_lt',
        'lte' => '_lte',
        '<=' => '_lte',
        '==' => '_num_eq',
        '!=' => '_num_ne',
        'between' => '_between',
        '<=>' => '_between',
    );
    
    my $function_method = $functions_methods{$function};

    $self->_communicate_error("Unkown function ($function) for checking numeric range for attr ($attr)")
        and return unless $function_method;

    ($num1, $num2) = (defined $num2 and $num2 < $num1 )
    ? ($num2, $num1)
    : ($num1, $num2);

    return $self->$function_method($type, $attr, $value, $num1, $num2);
}

sub _num_eq
{
    my ($self, $type, $attr, $value, $good_val) = @_;

    $self->_communicate_error("$attr ($value) must be equal to $good_val")
        and return unless $value == $good_val;

    return 1;
}

sub _num_ne
{
    my ($self, $type, $attr, $value, $bad_val) = @_;

    $self->_communicate_error("$attr ($value) must not be equal to $bad_val")
        and return unless $value != $bad_val;

    return 1;
}

sub _pos
{
    my ($self, $type, $attr, $value) = @_;

    $self->_communicate_error
    (
        sprintf('%s (%s) needs to be a positive %s', $attr, $value, $type) 
    )
        and return unless $value > 0;
    
    return 1;
}

sub _non_neg
{
    my ($self, $type, $attr, $value) = @_;

    $self->_communicate_error
    (
        sprintf('%s (%s) needs to be a non negative %s', $attr, $value, $type) 
    )
        and return unless $value >= 0;
    
    return 1;
}

sub _neg
{
    my ($self, $type, $attr, $value) = @_;

    $self->_communicate_error
    (
        sprintf('%s (%s) needs to be a negative %s', $attr, $value, $type) 
    )
        and return unless $value < 0;
    
    return 1;
}

sub _gt
{
    my ($self, $type, $attr, $value, $min) = @_;

    $self->_communicate_error("Minimum for $attr needs to be a $type")
        and return unless defined $min and $min =~ /$RE{num}{$type}/;

    $self->_communicate_error("$attr ($value) needs to be greater than $min")
        and return unless $value > $min;

    return 1;
}

sub _gte
{
    my ($self, $type, $attr, $value, $min) = @_;

    $self->_communicate_error("Minimum for $attr needs to be a $type")
        and return unless defined $min and $min =~ /$RE{num}{$type}/;

    $self->_communicate_error("$attr ($value) needs to be greater than or equal to $min")
        and return unless $value >= $min;

    return 1;
}

sub _lt
{
    my ($self, $type, $attr, $value, $max) = @_;

    $self->_communicate_error("Maximum for $attr needs to be a $type")
        and return unless defined $max and $max =~ /$RE{num}{$type}/;

    $self->_communicate_error
    (
        sprintf('%s (%s) needs to be less than %s', $attr, $value, $max) 
    )
        and return unless $value < $max;

    return 1;
}

sub _lte
{
    my ($self, $type, $attr, $value, $max) = @_;

    $self->_communicate_error("Maximum for $attr needs to be a $type")
        and return unless defined $max and $max =~ /$RE{num}{$type}/;
    
    $self->_communicate_error
    (
        sprintf('%s (%s) needs to be less than or equal to %s', $attr, $value, $max) 
    )
        and return unless $value <= $max;

    return 1;
}

sub _between
{
    my ($self, $type, $attr, $value, $min, $max) = @_;

    $self->_communicate_error("Minimum for $attr needs to be a $type")
        and return unless defined $min and $min =~ /$RE{num}{$type}/;

    $self->_communicate_error("Maximum for $attr needs to be a $type")
        and return unless defined $max and $max =~ /$RE{num}{$type}/;
    

    $self->_communicate_error
    (
        sprintf('%s (%s) needs to be between %s and %s', $attr, $value, $min, $max) 
    )
        and return unless $value >= $min and $value <= $max;

    return 1;
}

# FILE
sub _file
{
    my ($self, $attr, $value) = @_;

    return unless $self->_defined($attr, $value);

    $self->_communicate_error
    (
        sprintf("%s (%s) is a directory (should be a file)", $attr, $value)
    )
        and return if -d $value;

    $self->_communicate_error
    (
        sprintf("%s (%s) does not exist", $attr, $value)
    )
        and return unless -e $value;
 
    return 1;
}

sub _file_r
{
    my ($self, $attr, $value) = @_;

    return unless $self->_file($attr, $value);

    $self->_communicate_error
    (
        sprintf("%s (%s) is empty or does not exist", $attr, $value)
    )
        and return unless -s $value;
    
    $self->_communicate_error
    (
        sprintf("%s (%s) is not readable", $attr, $value)
    )
        and return unless -r $value;
    
    return 1;
}

sub _file_w
{
    my ($self, $attr, $value) = @_;
    
    return unless $self->_defined($attr, $value);

    $self->_communicate_error( sprintf("%s (%s) exists", $attr, $value) )
        and return if -e $value;
    
    my ($file, $dir) = fileparse($value);
    
    $self->_communicate_error("Directory ($dir) for file ($file) is not writable")
        and return unless -w $dir;
 
    return 1;
}

sub _file_rw
{
    my ($self, $attr, $value) = @_;

    return unless $self->_defined($attr, $value);

    $self->_communicate_error
    (
        sprintf("%s (%s) is a directory (should be a file)", $attr, $value)
    )
        and return if -d $value;

    my ($file, $dir) = fileparse($value);
    
    $self->_communicate_error("Directory ($dir) for file ($file) is not readable")
        and return unless -r $dir;

    $self->_communicate_error("Directory ($dir) for file ($file) is not writable")
        and return unless -w $dir;
 
    return 1;

    
    return 1;
}

sub _exe
{
    my ($self, $attr, $exe) = @_;

    return unless $self->_file_r($attr, $exe);
    
    $self->_communicate_error( sprintf("File (%s) is not executable", $attr) )
        and return unless -x $exe;
    
    return 1;
}

# DIR
sub _dir
{
    my ($self, $attr, $value) = @_;

    $self->_communicate_error
    (
        sprintf('%s (%s) does not exist', $attr, $value)
    )
        and return unless -e $value;
    
    $self->_communicate_error
    (
        sprintf('%s (%s) is not a directory', $attr, $value)
    )
        and return unless -d $value;
    
    return 1;
}

sub _dir_r
{
    my ($self, $attr, $value) = @_;

    return unless $self->_dir($attr, $value);

    $self->_communicate_error
    (
        sprintf('Can\'t read from %s (%s)', $attr, $value)
    )
        and return unless -r $value;

    return 1;
}

sub _dir_w
{
    my ($self, $attr, $value) = @_;

    return unless $self->_dir($attr, $value);

    $self->_communicate_error
    (
        sprintf('Can\'t write to %s (%s)', $attr, $value)
    )
        and return unless -w $value;

    return 1;
}

sub _dir_rw
{
    my ($self, $attr, $value) = @_;
    
    return unless $self->_dir_r($attr, $value);

    return unless $self->_dir_w($attr, $value);

    return 1;
}

sub _in_list
{
    my ($self, $attr, $value, $opts) = @_;

    return unless $self->_defined($attr, $value);

    $self->_communicate_error("Need list of valid values in isa for attribute ($attr)")
        and return unless $opts;
    
    $self->_communicate_error
    (
        sprintf
        (
            "%s (%s) is not in the list of valid values (%s)", 
            $attr,
            $value,
            join(', ', @$opts)
        )
    )
        and return unless grep { $value eq $_ } @$opts;

    return 1;
}

sub _object
{
    my ($self, $attr, $obj, $classes) = @_;

    return unless $self->_defined($attr, $obj);

    my $class = blessed($obj);

    $self->_communicate_error("$attr is not a blessed reference")
        and return unless defined $class;

    if ( $classes )
    {
        $self->_communicate_error
        (
            "$class does not inherit from " . join(', ', @$classes)
        )
            and return unless grep { UNIVERSAL::isa($obj, $_) } @$classes;
    }
    
    return $class;
}

sub _regex
{
    my ($self, $attr, $value, $opts) = @_;

    $self->_communicate_error("Need list of valid regexs in isa for attribute ($attr)")
        and return unless $opts;

    return unless $self->_defined($attr, $value);

    my @patterns = (( ref($opts) ) ? @$opts : $opts );
    
    $self->_communicate_error
    (
        sprintf
        (
            "%s (%s) does not match the given patterns (%s)",
            $attr,
            $value,
            join(', ', @patterns),
        )
    )
        and return unless grep { $value =~ /$_/i } @patterns;

    return 1;
}

sub _code
{
    my ($self, $attr, $value) = @_;

    return unless $self->_defined($attr, $value);

    $self->_communicate_error($attr . " is not a code ref")
        and return unless ref($value) eq 'CODE';
    
    return 1;
}

sub _y_or_n
{
    my ($self, $attr, $value) = @_;

    return unless $self->_defined($attr, $value);

    $self->_communicate_error("$attr needs to be [y]es or [n]o")
        and return unless $value =~ /^[yn]/i;
    
    return 1;
}

sub _not_blank
{
    my ($self, $attr, $value) = @_;

    return unless $self->_defined($attr, $value);

    $self->_communicate_error("$attr is blank")
        and return if $value eq '';
    
    return 1;
}

sub _project_name
{
    my ($self, $attr, $name) = @_;

    return unless $self->_defined($attr, $name);

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

    my $sth = $dbh->prepare("select * from projects where name = '$name'");
    $sth->execute;
    my $aryref = $sth->fetchall_arrayref;
    $dbh->disconnect;
    
    $self->_communicate_error( sprintf("Invalid %s (%s)", $attr, $name) )
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

Package with generic validation methods.  Evaluates a value for an attribute's isa.  On failure, a message will be deployed in this order:

=over

=item callback (cb), if it is defined

=item message method (msg), if defined will call the message method on the object (obj, if defiend) or the calling class.

=item create a Finfo::Msg with level of fatal, print it, then croak.
 
=back

Note that the value will be checked to see if it is defined.

=head1 Usage

 # Check if an age is valid.  Without a cb or msg param, this will
 # confess on failure.
 print "What is your age? ";
 my $age = <STDIN>;
 chomp $age;

 Finfo::Validate->validate
 (
    attr => 'age',
    value => $age,
    isa => 'int pos',
 );
 
=head1 Methods

=head2 validate

The main method, evaluates a value.  Returns 1 if ok or undef and executes error callbacks.

=head2 is_isa

 my ($isa, @opts) = Finfo::Validate->is_isa # or call as scalar
 (
    attr => 'Isa for attribute age',
    value => $isa,
    msg => 'fatal', # and obj or cb
 );

Tests if the value is a valid isa.  Returns the 'isa and valid options on success (or true if called in scalar context).Executes callbacks and returns undef on failure.  Note tha the 'isa' and 'ds' parameters are not valid here.

=head2 is_ds

 Finfo::Validate->is_ds
 (
    attr => 'data structure for attribute age',
    value => $ds,
    msg => 'fatal', # and obj or cb
 );

Tests if the value is a valid data structure.  Returns true on success.  Executes callbacks and returns undef on failure.  Note tha the 'isa' and 'ds' parameters are not valid here.

=head1 Definitions of Parameters

=head2 attr (Attribute, required)

The name of the value.  Used for message communication.

=head2 value (Value, required)

The value to check.

=head2 ds (Data Structures, optional)

Indicate a data structure by setting the 'ds' key to the incoming validation parameters.  The default assumes that scalar is the data structure.

=over 

=item scalar

=item aryref

=item hashref

=back

=head2 ISA (isa)

The isa parameter indicates what the attribute is.  This is a string that is space separated parameters.  The first is the 'isa', the rest are considered the attributes options. 

=head2 Attribute Options (included with 'isa')

The attribute options are generic options that are specific to the isa.  For example, the 'inlist' isa needs a list of valid type that the value can be.

=head2 empty_ok (aryref/hashref can be empty, optional)

This applies to data structures aryref and hashref.  If set to true, the aryref/hashref will not be checked if it has values.  The default is false.

=head2 cb (Callback)

A custom callback created by the caller.  Executed with the error msg as the first param on failure

=head2 msg (Message Level, optional)

The msg method to call, if failure occurs (should be one of Finfo::Msg::levels)

=head2 obj (Object/Class, optional)

The object/class to call the above mesage method on failrue.  If this is set and the msg is not, an attempt to call 'fatl_msg' on the object/class will be made.  If the object cannot do this method, a fatal msg will be called from Finfo::Validate to indicate this.
 
=head1 ISAs

=head2 General

=over

=item I<defined>
 
Checks if value is defined.  This is done automatically, so there isn't a need to set this isa.


=item I<code> Checks is a value is defined and a reference to code
 
=item I<is_isa> Checks is a value is defined and is an 'isa' validation itself

=item I<y_or_n> Checks is a value is defined and starts with an Y, y, N, or n

=item I<project_name> Checks is a value is defined and is a GSC project name.  It returns a hashref of 
 the project info.

=back

=head2 General with Valid Options

Valid options can be specified with these isas.  To do that, just add the options to the isa, separating each by a space.  After the value is checked to see if it is defiend, The value will be checked against the valid options.  

Ex: 
'object Finfo::Std',

This would check if the value is defined, a blessed reference and if it inherits from Finfo::Std

=over 

=item I<object> Checks if the value is a blessed reference and , if given, that it inherits from B<one> of the classes in the options.

=item I<string> Checks if the value is not a reference, and if given, checks the length of the string by adding numeric range parameters (see below)
 
=item I<in_list> Checks if the is in the list of options

=item I<regex> Checks if a value is matches B<one> of the patterns given in the options

=back 

=head2 Number ISAs

=over 

=item I<int>

Checks if value is defined and an integer.  Can also check the numeric range by adding additional parameters to the isa.  See the 'Numeric Range Paramters' below for more info.

=item I<real>

Checks if value is defined and a real number (an integer is a real number).  Can also check the numeric range by adding additional parameters to the isa.  See the 'Numeric Range Paramters' below for more info.

=item B<Numeric Range Parameters>

=over

=item pos

greater than 0

=item non_neg

greater than or equal to 0

=item 'gt #', '> #'

greater than the number (#) after the 'gt' or '>'

=item 'gte #', '>= #'

greater than or equal to the number (#) after the 'gt' or '>'

=item 'lt #', '> #'

greater than the number (#) after the 'lt' or '>'

=item 'lte #', '>= #'

greater than or equal to the number (#) after the 'lt' or '>'

=item 'between #1 #2', '<> #1 #2'

greater than or equal to the first number(#1) and less than or equal to the second number (#2)

=back

=back

=head2 Files and Directories

=over

=item I<file> Checks if the file is defined and not a directory

=item I<file_r> Checks if the file is defined, exists, has size and is readable
 
=item I<file_w> Checks if the file is defined, does not exist and that the directory is writable
  
=item I<file_rw> Checks if the file is defined, not a directory, and that the directory is readable/writable

=item I<dir> Checks if directory is defined and a directory
 
=item I<dir_r> Checks if directory is defined, a directory and is readable

=item I<dir_w> Checks if directory is defined, a directory and is writable

=item I<dir_rw> Checks if a directory is defined, a directory, and is readable and writable

=item I<exe> Checks if a file is defined, readable and executable

=back

=head1 Example(s) 
 
 # Get the value from a gui widget and validate it.  This will execute
 # the callback (cb) if the 'answer' is not in the list of valid
 # options in the isa
 my $value = $entry->get_text;

 Finfo::Validate->validtate
 (
    attr => 'answer',
    value => $value,
    isa => 'in_list yes yeah yo yep', # could use type y_or_n for this
    cb => sub{ my $msg = shift; $self->gui_popup($msg) },
 );

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
