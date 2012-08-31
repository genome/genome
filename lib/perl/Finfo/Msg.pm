package Finfo::Msg;

use strict;
use warnings;

use base 'Class::Accessor';

use Carp;
use Data::Dumper;
use Date::Format;
use Devel::StackTrace;
use Term::ANSIColor;

Finfo::Msg->mk_ro_accessors(qw/ class_of_obj msg level detail_level file_name line /);

sub new
{
    my ($class, %p) = @_;

    my $caller = ( exists $p{'caller'} ) ? delete $p{'caller'} : [ caller ];

    unless ( $p{msg} )
    {
        my $msg_obj = __PACKAGE__->new
        (
            msg => "No message sent to create object",
            level => 'fatal',
            deatil_level => 'more',
            'caller' => $caller,
        );
        print STDERR $msg_obj->string;
        croak;
    }
    
    @p{qw/ class_of_obj file_name line /} = @$caller;

    foreach my $attr (qw/ level detail_level /)
    {
        my $default_method = 'default_' . $attr;
        $p{$attr} = $class->$default_method if not defined $p{$attr};

        my $options_method =  $attr . 's';

        unless ( grep { $p{$attr} eq $_ } __PACKAGE__->$options_method )
        {
            my $msg_obj = __PACKAGE__->new
            (
                msg => "Invalid $attr ($p{$attr})",
                'caller' => $caller,
                level => 'fatal',
                detail_level => 'more',
            );
            print STDERR $msg_obj->string;
            croak;
        }
    }

    return bless \%p, $class;
}

sub string
{
    my $self = shift;

    my $string_method = sprintf('_%s_string', $self->detail_level);

    return $self->$string_method;
}

sub _normal_string
{ 
    my $self = shift;

    return sprintf
    (               
        '%s %s -> %s', 
        Term::ANSIColor::colored
        (
            sprintf(' %-6s', uc $self->level), color_params_for_level($self->level)
        ),
        $self->class_of_obj,
        $self->msg,
    );
}

sub _more_string
{ 
    my $self = shift;

    return sprintf
    (               
        "%s Msg:\t%s\n\tClass:\t%s\n\tFile:\t%s\n\tLine:\t%s\n\tDate:\t%s",
        #"%s %s -> %s\n\t at line %d in %s (%s)", 
        Term::ANSIColor::colored
        (
            sprintf(' %-6s', uc $self->level), color_params_for_level($self->level)
        ),
        $self->msg,
        $self->class_of_obj,
        $self->file_name,
        $self->line,
        Date::Format::time2str('%Y%b%d %X', time),
    );
}

sub _stack_string
{
    my $self = shift;

    my $stack_trace = Devel::StackTrace->new
    (
        ignore_class => 'Finfo::Msg',
    );

    return sprintf
    (               
        "%s Msg:\t%s\n\tClass:\t%s\n\tFile:\t%s\n\tLine:\t%s\n\tDate:\t%s\n%s",
        #"%s %s -> %s\n\t at line %d in %s (%s)", 
        Term::ANSIColor::colored
        (
            sprintf(' %-6s', uc $self->level), color_params_for_level($self->level)
        ),
        $self->msg,
        $self->class_of_obj,
        $self->file_name,
        $self->line,
        Date::Format::time2str('%Y%b%d %X', time),
        $stack_trace->as_string,
    );
}

sub levels
{
    return (qw/ debug info warn error fatal /);
}

sub detail_levels
{
    return (qw/ normal more stack /);
}

sub default_detail_level
{
    return (detail_levels())[0];
}

sub color_params_for_level
{
    my ($level) = @_;
    
    my $color_params = 
    {
        fatal => 'bold yellow on_red',
        error => 'bold white on_red',
        'warn' => 'bold white on_yellow',
        info => 'bold white on_green',
        debug => 'bold white on_blue',
    };

    return $color_params->{$level};
}

1;

=pod

=head1 Name

Finfo::Msg
 
=head1 Synopsis

A simple message class.  Used by the Finfo::Logging classes to communicate messages
and store mesage data.

=head1 Usage

 my $caller = [ caller ];
 my $msg = Finfo::Msg->new
 (
    class_of_obj => $obj->class, # req
    msg => 'got foo barred', # req
    level => 'error', # req
    line => $caller->[1], # req
    file_name => $caller->[2], #req
    msg_detail => 'more', # opt, default is 'normal'
 )
    or croak;

 print $msg->string,"\n";
 
=head1 Class Methods

=head2 levels

Get a list of the valid msg levels

=head2 detail_levels

Get a list of the valid deatil levels

=head1 Object Methods

=head2 string

Get an output string listing of the msg attributes

=head2 class_of_obj
 
Get the class of the object where the msg occurred

=head2 msg 

Get the msg
 
=head2 level 

Get the level of the msg

=head2 file_name
 
Get the actual file (module) name where the msg occurred

=head2 line 
 
Get the line of the file name where the msg occurred
 
=head1 Notes

=head2 Message Levels

Message levels pertain to the importance of the msg.  These are inspired by Log::Log4perl's msgs levels.

=head2 Message Deatil Level

The message detail level determines how much infomation is added to the output generated from the 'string' method.

=over

=item I<Normal> - Displays the level, class and msg

=item I<More> - Displays the level, msg, class, file_name, line and date in an easy to read format

=item I<Stack> - Displays the same data as 'more', but then adds a full stack trace.

=back

=head1 See Also

=over

=item Finfo::Logging

=back

=head1 Disclaimer

Copyright (C) 2006-7 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details.

=head1 Author(s)

B<Eddie Belter> <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
