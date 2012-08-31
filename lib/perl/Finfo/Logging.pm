package Finfo::Logging;

use strict;
use warnings;

use Carp;
use Data::Dumper;
use Finfo::ClassUtils 'class';
use Finfo::Msg;
use Finfo::Validate;
use Log::Log4perl;
use Log::Log4perl::Level;

# Init the base screen logger
my $conf = qq/
log4perl.rootLogger=DEBUG, screen
log4perl.appender.screen.Threshold=INFO
log4perl.appender.screen=Log::Dispatch::Screen
log4perl.appender.screen.stderr=1
log4perl.appender.screen.layout=PatternLayout
log4perl.appender.screen.layout.ConversionPattern=%m%n
/;

Log::Log4perl->init( \$conf );

my @levels = Finfo::Msg->levels;

my @exported_subs = 
(qw/
    debug_msg short_debug_msg 
    info_msg short_info_msg
    warn_msg short_warn_msg 
    error_msg short_error_msg
    fatal_msg short_fatal_msg
    /);

my 
(
    %error_msg, %error_msg_obj, 
    %warn_msg, %warn_msg_obj,
    %info_msg, %info_msg_obj,
    %debug_msg, %debug_msg_obj,
    %fatal_msg, %fatal_msg_obj,
    %msg_detail, %private_methods,
);

sub import 
{
    my $self = shift;

    my $caller = ( defined $_[0] and ref($_[0]) eq 'HASH' )
    ? shift->{class}
    : caller;

    my @subs_to_export = ( @_ )
    ? split(/\s+/, $_[0])
    : @exported_subs;
 
    no strict 'refs';

    for my $sub ( @subs_to_export ) 
    {
        __PACKAGE__->fatal_msg("Invalid sub ($sub) to export to $caller")
        unless grep { $sub eq $_ } @exported_subs;

        next if $caller->can($sub);

        *{ $caller . '::' . $sub } = \&{$sub}
    }

    return 1;
}

# Logger functions
sub _log_msg
{
    my $msg = shift;

    Finfo::Validate->validate
    (
        attr => 'message object to log',
        value => $msg,
        isa => 'object Finfo::Msg',
        cb => sub{ print Dumper(\@_); __PACKAGE__->fatal_msg(@_); },
    );

    unless ( Log::Log4perl->initialized )
    {
        my $not_init_msg = Finfo::Msg->new
        (
            msg => "Can't log messages, logger not initialized",
            level => 'fatal',
            detail_level => 'normal',
            'caller' => [ caller ],
        );
        print STDERR $not_init_msg->string if $not_init_msg;
        print STDERR $msg->string;
        croak;
    }

    return 1 unless _check_log_level($msg);

    my $level = $msg->level;

    return Log::Log4perl->get_logger( $msg->class_of_obj )->$level($msg->string);
}

sub set_screen_logger_level
{
    my $level = shift;

    _is_msg_level_valid($level);

    my $priority = Log::Log4perl::Level::to_priority( uc $level );
    
    return Log::Log4perl->appenders->{screen}->{level} = $priority;
}

sub _check_log_level
{
    my $msg_obj = shift;

    __PACKAGE__->fatal_msg("Need msg obj to check log level") unless $msg_obj;
    
    my $method = 'is_' . lc $msg_obj->level;

    return Log::Log4perl->get_logger( $msg_obj->class_of_obj )->$method;
}

sub add_appender
{
    my (%p) = @_;

    my $class = delete $p{class}
        or __PACKAGE__->fatal_msg("No class to add logging appender");

    my $logger = Log::Log4perl->get_logger($class);

    my $level = delete $p{level} || 'INFO';
    $logger->level($level);

    if ( exists $p{type} )
    {
        my $appender = _create_appender($p{type}, $p{params});
        $appender->layout( Log::Log4perl::Layout::PatternLayout->new('%m%n') );
        
        $logger->add_appender($appender);
    }
    elsif ( my $appender = delete $p{appender} )
    {
        $logger->add_appender($appender);
    }
    else
    {
        __PACKAGE__->fatal_msg("No logging appenders to add");
    }

    return 1;
}

sub _create_appender
{
    my ($type, $params) = @_;
   
    my %app_types_classes = 
    (
        buffer => 'Log::Log4perl::Appender::Buffer',
        dbi => 'Log::Log4perl::Appender::DBI',
        file => 'Log::Log4perl::Appender::File',
        limit => 'Log::Log4perl::Appender::Limit',
        rrds => 'Log::Log4perl::Appender::RRDs',
        screen => 'Log::Log4perl::Appender::Screen',
        screen_colored_levels => 'Log::Log4perl::Appender::ScreenColoredLevels',
        'socket' => 'Log::Log4perl::Appender::Socket',
        string => 'Log::Log4perl::Appender::String',
        synchronized => 'Log::Log4perl::Appender::Synchronized',
        test_array_buffer => 'Log::Log4perl::Appender::TestArrayBuffer',
        test_buffer => 'Log::Log4perl::Appender::TestBuffer',
        test_file_creeper => 'Log::Log4perl::Appender::TestFileCreeper',
    );

    my $class = $app_types_classes{$type};

    __PACKAGE__->fatal_msg("Unsupported appender type ($type)") unless $class;

    Finfo::Validate->validate
    (
        attr => 'appender params',
        value => $params,
        ds => 'hashref',
        cb => sub{ __PACKAGE__->fatal_msg(@_); },
    );
 
    return Log::Log4perl::Appender->new
    (
        $class,
        %$params,
    );
}

# Msg functions
sub _is_msg_level_valid
{
    my $level = shift;

    Finfo::Validate->validate
    (
        attr => 'message level',
        value => $level,
        isa => [ 'in_list',  Finfo::Msg->levels ],
        cb => sub{ __PACKAGE__->fatal_msg(@_); },
    );
}

sub get_msg_detail_level
{
    my $class = shift;
    
    Finfo::Validate->validate
    (
        attr => 'class to set detail level',
        value => $class,
        cb => sub{ $class->fatal_msg(@_); },
    );
    
    return $msg_detail{$class}; 
}

sub set_msg_detail_level
{
    my ($class, $detail_level) = @_;
    
    Finfo::Validate->validate
    (
        attr => 'class to set detail level',
        value => $class,
        cb => sub{ $class->fatal_msg(@_); },
    );

    Finfo::Validate->validate
    (
        attr => 'msg detail to set',
        value => $detail_level,
        isa => [ 'in_list', Finfo::Msg->detail_levels ],
        cb => sub{ $class->fatal_msg(@_); },
    );
    
    $msg_detail{$class} = $detail_level;

    return $msg_detail{$class};
}

# Msgs
sub debug_msg
{
    my ($self, $msg) = @_;

    return _msg
    (
        class => class($self),
        msg => $msg, 
        'caller' => [ caller ],
        level => 'debug',
    );
}

sub short_debug_msg
{
    my $self = shift;
    
    $self->fatal_msg('Cannot set debug msg string') if @_;
    
    return _msg_string
    (
        class => class($self),
        level => 'debug',
    );
}

sub info_msg
{
    my ($self, $msg) = @_;

    return _msg
    (
        class => class($self),
        msg => $msg, 
        'caller' => [ caller ],
        level => 'info',
    );
}

sub short_info_msg
{
    my $self = shift;
    
    $self->fatal_msg('Cannot set info msg string') if @_;
    
    return _msg_string
    (
        class => class($self),
        level => 'info',
    );
}

sub warn_msg
{
    my ($self, $msg) = @_;

     return _msg
    (
        class => class($self),
        msg => $msg, 
        'caller' => [ caller ],
        level => 'warn',
    );
}

sub short_warn_msg
{
    my $self = shift;
    
    $self->fatal_msg('Cannot set warn msg string') if @_;
    
    return _msg_string
    (
        class => class($self),
        level => 'warn',
    );
}

sub error_msg
{
    my ($self, $msg) = @_;

    return _msg
    (
        class => class($self),
        msg => $msg, 
        'caller' => [ caller ],
        level => 'error',
    );
}

sub short_error_msg
{
    my $self = shift;
    
    $self->fatal_msg('Cannot set error msg string') if @_;
    
    return _msg_string
    (
        class => class($self),
        level => 'error',
    );
}

sub fatal_msg
{
    my ($self, $msg, $p) = @_;

    $msg = 'Something really bad happened, but no message was provided' unless $msg;
    my $class = delete $p->{class} || class($self);
    my $caller = delete $p->{'caller'} || [ caller ];
    # TODO get error codes from Finfo::ErrorCode??
    if ( my $err_code = delete $p->{err_code})
    {
        $msg .= sprintf(' (error code %s)', $err_code)
    }
    
    my $msg_obj = Finfo::Msg->new
    (
        msg => $msg,
        level => 'fatal',
        detail_level => 'more',#get_msg_detail_level($class),
        'caller' => [ $class, $caller->[1], $caller->[2] ],
    )
        or croak "Can't create Finfo::Msg";

    my $msg_var = _msg_var('fatal');

    _log_msg($msg_obj);

    no strict 'refs';
    # Set msg on class
    ${$msg_var}{$class} = $msg_obj;

    croak "$msg\n";
}

sub short_fatal_msg
{
    my $self = shift;
    
    $self->fatal_msg('Cannot set error msg string') if @_;
    
    return _msg_string
    (
        class => class($self),
        level => 'fatal',
    );
}

sub _msg_var
{
    my $level = shift;

    croak "No level\n" unless $level;

    return sprintf('%s_msg', $level);
}

sub _msg_string
{
    my (%p) = @_;

    my $var = _msg_var($p{level});

    my $msg_obj;
    no strict 'refs';
    $msg_obj = ${ $var }{ $p{class} };

    use strict 'refs';

    return unless $msg_obj;

    return $msg_obj->msg;
}

sub _msg
{
    my (%p) = @_;

    _is_msg_level_valid($p{level});
    
    __PACKAGE__->fatal_msg("No class for logging msg")
        and return unless $p{class};

    my $msg_var = _msg_var($p{level});

    if ( defined $p{msg} )
    {
        my $msg_obj = Finfo::Msg->new
        (
            msg => $p{msg},
            level => $p{level},
            detail_level => get_msg_detail_level($p{class}),
            'caller' => [ $p{class}, $p{'caller'}->[1], $p{'caller'}->[2] ],
        )
            or croak "Can't create Finfo::Msg";

        _log_msg($msg_obj);

        no strict 'refs';
        # Set msg on class
        ${$msg_var}{ $p{class} } = $msg_obj;
    }

    no strict 'refs';
    return ${$msg_var}{ $p{class} };
}

# Cleanup
sub destroy_msgs
{
    my $self = shift;

    $self->warn_msg
    (
        "Messages are not being stored by object, so there is no reason to remove the messages"
    );
    
    return 1;
}

1;

=pod

=head1 Name

Finfo::Logging
 
=head1 Synopsis

Add messaging methods and logging to a class.  Uses Log::Log4perl to init a 'base'
screen logger and add Log::Log4perl appenders that write to files, sockects, etc.

=head1 Usage

B<Create a class, 'use' Finfo::Logging>

 pakage MyClass;
 
 use strict;
 use warnings;

 use Finfo::Logging;

 sub new
 {
    return bless {}, shift;
 }

B<later...>

 use MyClass;

 my $obj = MyClass->new();

 $obj->info_msg("Created MyClass");

 $obj->warn_msg("Look out!!");
 
 unless ( $obj->can('foo') )
 {
    $obj->fatal_msg("MyClass can't foo");
 }

B<Can be used to add logging to a script...use 'main' as the class when creating appenders, etc.  DESTROY method not needed, only 'class' msgs are stored>

 #! /usr/bin/perl

 use strict;
 use warnings;

 use Finfo::Logging;

 # use arrow method on 'main' package
 main->info_msg("Starting...");

 
B<Exported methods should be called as an object method (with arrow) only>

=head1 Methods that are exported to your class

=head2 fatal_msg

$obj->fatal_msg($msg, %p);

=over

=item I<Synopsis>   Logs the fatal error and croaks.  

=item I<Params>     msg (string), params (hash with opt keys error_code - adds error code to msg)

=item I<Returns>    croaks, but will returns true if die signal is caught

=back

=head2 error_msg warn_msg info_msg debug_msg

my $msg_obj = $obj->error_msg("Oops...I did it again");

=over

=item I<Synopsis>   Sets the msg for the class for msg level, then logs the msg 

=item I<Params>     message (string)

=item I<Returns>    Finfo::Msg (object) 

=back

=head2 short_error_msg short_warn_msg short_info_msg short_debug_msg

my $short_error_msg = $self->short_error_msg;

=over

=item I<Synopsis>   Gets the message string for the msg level

=item I<Params>     none

=item I<Returns>    msg (string) 

=back

=head1 Methods not exported (call from Finfo::Logging, as class methods)

=head2 set_screen_logger_level

 Finfo::Logging::set_screen_logger_level($level)
    or die;

=over

=item I<Synopsis>   Sets the 'base' screen logger to a new level threshold

=item I<Params>     level (string, Finfo::Msg->levels)

=item I<Returns>    boolean - true on success

=back

=head2 get_msg_detail_level

 Finfo::Logging::set_screen_logger_level($level)
    or die;

=over

=item I<Synopsis>   Get the msg detail (see Finfo::Msg) for the class

=item I<Params>     none 

=item I<Returns>    msg detail (string, Finfo::Msg::detil_levels)

=back

=head2 set_msg_detail_level

 Finfo::Logging::set_screen_logger_level($level)
    or die;

=over

=item I<Synopsis>   Set the msg detail (see Finfo::Msg) for the class

=item I<Params>     msg detail (string, Finfo::Msg::detil_levels)

=item I<Returns>    msg detail (string, Finfo::Msg::detil_levels)

=back

=head2 add_appender

 Finfo::Logging::add_appender(%p)
    or die;

=over

=item I<Synopsis>   Creates and adds an appender (Log::Log4perl::Appender classes) with a layout ( for a class Log::Log4perl::PatternLayout) for a class.  This inits the logger for the class.  Fimiliarity with the Log::Log4perl::Appenders and the Log::Log4perl::Layout::PatternLayout classes is required, but easliy obtained.  See the man pages for these classes. 

=item I<Params> params (hash)

=over

=item B<class>   req, the class to add the appender to

=item B<type>    req, the type of appender

=item B<params>  req, the specific params for the appender (see appender classes)

=item B<level>   opt, the log level (default = 'INFO')

=item B<pattern> opt, the pattern for message printing (default = '%m%n')

=back

=item I<Returns> boolean - true on success

=back

=head1 See Also

I<Finfo::Std>, I<Finfo::Msg>, I<Finfo::Validate>, I<Log::Log4perl>, I<Class::Std>

=head1 Disclaimer

Copyright (C) 2006-7 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

B<Adam Dukes> I<adukes@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
