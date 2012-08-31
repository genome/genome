package Finfo::CommandLineOptions;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;
use File::Basename;
use Getopt::Long;
use PP::LSF;
use Term::ANSIColor;
use Text::Wrap;

require Finfo::ClassUtils;

my %classes :name(classes:o) :ds(aryref) :default([]) :empty_ok(1);
my %add_q :name(add_q:o) :default(0);
my %q_opts :name(q_opts:o) :ds(hashref) :default({ q => 'long', u => $ENV{USER} });
my %no_clo_not_ok :name(no_clo_not_ok:o) :default(0);
my %additional_opts :name(additional_opts:o) :ds(aryref);
my %header_msg :name(header_msg:o);
my %add_to_header :name(add_to_header:o);
my %footer_msg :name(footer_msg:o) :default(None);
my %script_basename :name(_script_basename:p) :isa(string);
my %usage :name(_usage:p) :ds(hashref) :empty_ok(1) :default({});
my %req_clo :name(_req_clo:p) :ds(hashref) :empty_ok(1) :default({});
my %clo :name(_clo:p) :ds(aryref) :empty_ok(1) :default([]);
my %opts :name(_opts:p) :ds(hashref) :empty_ok(1) :default({});

sub START
{
    my $self = shift;

    $self->_script_basename( basename($0) );

    unless ( $self->header_msg )
    {
        $self->header_msg("Usage for " . $self->_script_basename);
    }

    if ( defined $self->add_to_header )
    {
        $self->header_msg( $self->header_msg . $self->add_to_header );
    }
    
    foreach my $class ( @{ $self->classes } )
    {
        my @clos = $self->_clos_for_class($class);
        $self->fatal_msg("No command line options for class ($class)") unless @clos;
        foreach my $clo ( @clos )
        {
            $self->_add_clo(%$clo);
        }
    }

    $self->_add_clo
    (
        clo => 'h|help|usage',
        cat => 'Help',
        desc => 'This help/usage message'
    )
        or return;

    if ( $self->add_q )
    {
        $self->_add_clo
        (
            clo => 'q|queue',
            cat => 'Optional',
            desc => 'Run in the LSF Queue'
        )
            or return;
    }

    $self->fatal_msg
    (
        sprintf("No command line options from classes (%s)", join(', ', @{ $self->classes }))
    )
        and return unless $self->_clo;

    if ( $self->additional_opts )
    {
        foreach my $clo_params ( @{ $self->additional_opts } )
        {
            $clo_params->{class} = 'ADDITIONAL_OPTS' unless exists $clo_params->{class};

            $self->_add_clo(%$clo_params)
                or return;
        }
    }

    return 1;
}

sub _clos_for_class : PRIVATE
{
    my ($self, $class) = @_;

    my @clos;
    foreach my $attr_method (qw/ required_attributes optional_attributes /)
    {
        $self->fatal_msg
        (
            "Can't get attributes from class ($class) using '$attr_method'.  Prehaps this class wasn't use'd or is misnamed."
        ) unless ( $class->can($attr_method) );

        foreach my $attr ( $class->$attr_method )
        {
            my $ds = $class->attributes_attribute($attr, 'ds');
            my $isa = $class->attributes_attribute($attr, 'isa');
            my ($is_a, @opts) = Finfo::Validate->is_isa
            (
                attr => $attr,
                value => $isa,
                msg => 'fatal',
            );

            if ( $is_a eq 'object' )
            {
                my $helper_class = $opts[0];
                next unless @opts == 1 and $helper_class->can('required_attributes')
                    and $helper_class->can('optional_attributes');
                push @clos, $self->_clos_for_class($helper_class);
                next;
            }

            my $desc = $class->attributes_attribute($attr, 'desc');
            next unless $desc;

            my $clo = $class->attributes_attribute($attr, 'clo') || $attr;
            $clo =~ s/\_/\-/g;
            $clo =~ s/\=.+//g;

            if ( $is_a =~ /^int/ )
            {
                $clo .= '=i';
            }
            elsif ( $is_a =~ /^real/ )
            {
                $clo .= '=f';
            }
            elsif ( $is_a ne 'boolean' )
            {
                $clo .= '=s';
            }

            if ( $ds eq 'aryref' )
            {
                $clo .= '@{,}';
            }

            push @clos, 
            {
                class => $class,
                attr => $attr,
                ds => $ds,
                isa => $isa,
                clo => $clo,
                cat => ( $class->attributes_attribute($attr, 'attr_type') eq 'r' )
                ? 'Required'
                : 'Optional',
                desc => $desc,
            };
        }
    }

    return @clos;
}

sub _add_clo : PRIVATE
{
    my ($self, %p) = @_;

    $self->fatal_msg("No command line option to add")
        and return unless $p{clo};
    
    push @{ $self->_clo }, $p{clo};

    my $usage_opt = $p{clo};
    $usage_opt =~ s/=.*$//;
    
    my $cat = ucfirst( lc ( $p{cat} || 'optional' ) );
    return unless Finfo::Validate->validate
    (
        attr => 'usage category',
        value => $cat,
        isa => [ 'in_list', $self->usage_categories ],
        err_cb => $self,
    );

    $self->_usage->{ $cat }->{$usage_opt} = $p{desc};

    foreach my $opt ( split(/\|/, $usage_opt) )
    {
        # af=s -> af
        # h|help -> h and help

        push @{ $self->_opts->{$opt} },
        {
            class => $p{class} || 'na',
            attr => $p{attr} || 'na',
        };

        next unless lc($p{cat}) eq 'required';

        $self->_add_req_clo($opt, $p{attr})
            or return;
    }

    return 1;
}

sub _add_req_clo : PRIVATE
{
    my ($self, $opt, $attr) = @_;

    $self->fatal_msg("Need cl opt and attr to add to req cl opts")
        and return unless $opt and $attr;
    
    my $req_clos = $self->_req_clo;
    push @{ $req_clos->{ $attr } }, $opt;

    return $self->_req_clo($req_clos);
}

sub _convert_command_line_options_to_params_for_classes : PRIVATE
{
    my ($self, $opts) = @_;

    my $params = {};
    foreach my $class ( @{ $self->classes } )
    {
        $params->{$class} = {};
    }

    foreach my $opt ( keys %$opts )
    {
        $self->fatal_msg("Option ($opt) not found")
            and return unless exists $self->_opts->{$opt};

        foreach my $class_ref ( @{ $self->_opts->{$opt} } )
        {
            $params->{ $class_ref->{class} }->{ $class_ref->{attr} } = $opts->{$opt};
        }
    }

    #print Dumper($opts); print Dumper($params);
    
    return $params;
}

sub _execute_script_in_queue : PRIVATE
{
    my ($self, $opts) = @_;

    Finfo::Validate->validate
    (
        attr => 'script to run in the queue',
        value => $0,
        type => 'exe',
        msg => 'fatal',
    );
    
    my $cmd = sprintf
    (
        '%s %s', 
        $0,
        ( $opts ) ?  join(' ', map { "--$_ $opts->{$_}" } keys %$opts ) : '',
    );

    $self->info_msg("Command: $cmd"); #return;

    my $job = PP::LSF->run # replace w/ finfo::pp!!!
    (
        command => $cmd,
        %{ $self->q_opts },
    );

    unless ( $job )
    {
        $self->fatal_msg("Queue submission failed.  Error, if any: $!");
    }
    else
    {
        $self->info_msg("Successfully submitted to the queue");
    }

    return;
}

# Public Methods
sub script_basename
{
    return shift->_script_basename;
}

sub usage_categories
{ 
    return (qw/ Required Optional Help Other Dev-Only /);
}

sub get_objects
{
    my $self = shift;

    my $opts = $self->get_options;

    return;
}

sub get_options
{
    my $self = shift;
    
    local $SIG{__WARN__} = sub 
    {
        my $msg = shift;
        print $self->usage;
        chomp $msg; 
        $self->fatal_msg($msg); 
        exit 1;
    };

    my %opts;
    GetOptions(\%opts, $self->command_line_options)
        or return;

    print $self->usage
        and exit 1 if not %opts and $self->no_clo_not_ok;

    foreach my $key (qw/ h help usage /)
    {
        print $self->usage
            and exit 0 if exists $opts{$key};
    }

    my $req_clos = $self->_req_clo;
    my @undef_req_opts;
    foreach my $req_attr ( keys %$req_clos )
    {
        foreach my $req_clo ( @{ $req_clos->{$req_attr} } )
        {
            push @undef_req_opts, $req_clo unless exists $opts{$req_clo};
        }
    }

    if ( @undef_req_opts )
    {
        print $self->usage;
        $self->fatal_msg
        (
            sprintf
            (
                'Command line param%s (%s) %s required',
                ( @undef_req_opts > 1 ) ? 's' : '',
                join(', ', @undef_req_opts),
                ( @undef_req_opts > 1 ) ? 'are' : 'is',
            )
        );
    }
    
    if ( $self->add_q and exists $opts{q} )
    {
        delete $opts{q};
        $self->_execute_script_in_queue(\%opts);
        exit 0;
    }

    #print Dumper(\%opts);
    
    return $self->_convert_command_line_options_to_params_for_classes(\%opts);
}

sub command_line_options
{
    return @{ shift->_clo };
}

sub usage
{
    my $self = shift;

    $Text::Wrap::columns = 100;

    my $banner;
    for (1..length($self->header_msg) + 2) { $banner .= ' '; }
    #my $usage = "\n" . colored($banner, 'bold white on_blue') . "\n";
    my $usage = "\n" . colored($self->header_msg, 'bold') . "\n";
    #$usage .= colored(' '.$self->header_msg.' ', 'bold white on_blue') . "\n";
    #$usage .= colored($banner, 'bold white on_blue') . "\n";
    #my $usage = "\n" . colored($self->header_msg, 'bold white on_green') . "\n";

    print color 'reset';

    foreach my $cat ( $self->usage_categories )
    {
        next unless exists $self->_usage->{$cat};
        $usage .= " " . colored($cat, 'underline') . "\n";
        #$usage .= colored(" $cat", 'bold') . "\n";
        foreach my $opt ( sort { $a cmp $b } keys %{ $self->_usage->{$cat} } )
        {
            $usage .= Term::ANSIColor::colored(sprintf('  %-20s', '--' . $opt), 'bold')."\n";
            #$usage .= Term::ANSIColor::colored(sprintf('  %-20s', '--' . $opt), 'bold blue')."\n";
            $usage .= Text::Wrap::wrap
            (
                "     ", # 1st line indent, 
                "      ", # all other lines indent, 
                $self->_usage->{$cat}->{$opt}, # desc
            );
            $usage .= "\n";
        }
    }

    if ( $self->footer_msg )
    {
        $usage .= " " . colored("Notes", 'underline') . "\n";
        $usage .= Text::Wrap::wrap("  ", "   ", $self->footer_msg) . "\n\n";
    }

    return  $usage;
}

1;

=pod

=head1 Name

Finfo::CommandLineOptions

=head1 Synopsis

This class combines the usefulness of the Finfo::Std and Finfo::Object attributes and Getopt::Long to get the command line options and process them for your classes.

=head1 Usage

B<In your class using Finfo::Std> 

 package Album::Reader;

 use Finfo::Std;
 
 my %input_file :name(input_file:r) :type(input_file) :clo("input-file=s") :desc("File to get albums");

B<or in your class using base Finfo::Object> 

 package Album::Reader;

 use base 'Finfo::Object';

 sub _attrs
 {
    return
    {
        'input_file:r' => 
        {
            type => 'input_file',
            clo => 'input-file=s',
            desc => 'File to get albums',
        },
    };
 }

B<later...>

# Note: clo = command line options

 my $clo = Finfo::CommandLineOptions->new
 (
    classes => [qw/ AblumReader /], 
    header_msg => 'Display you alnums!',
    additional_opts => # An aryref of hashrefs of additional clo, with these keys:
    {
        clo => 'backup=s', # Getopt::Long option
        class => 'FILES', # class(key) to set the attr in opts hashref , default (ADDITIONAL_OPTS)
        attr => 'backup', # attribute to set
        cat => 'Optional', # usage category - Required Optional 
        desc => 'File to back up album file', # usage description
    },
 )
    or die;

 # get_options - returns a hashref - keys = classes, values = hashref
 # of params to create the object. ex:
 # $opts->{AlbumReader} = { input_file => './myalbums.txt' };
 # Also, will reexecute the script in the LSF queue if the add_q param
 # is true and the q is indicated on the command line
 my $opts = $clo->get_options
     or return;
 
 # create the object(s) using the params stored in $opts
 my $reader = $class->new
 (
     %{ $opts->{'AlbumReader'} },
 )
     or die;
 
 my $album = $reader->next;
 
 # etc...

B<Required Params> (defaults)

=over

=item I<classes>    aryref of the Finfo::Std (or Finfo::Object) classes to use (req)

=back

B<Optional Params> (defaults)

=over

=item I<add_q> - boolean, add an LSF queue option (0)

=item I<q_opts - hashref of the LSF queue options ({ u => $ENV{USER}, q => 'long' })

=item I<header_msg> - string msg at top oif usage ("Usage for $0")

=item I<footer_msg> - string msg at bottom of usage (none)

=item I<additional_opts> - aryref of hashrefs of additional clo, with these keys:

=over

=item I<clo> - Getopt::Long option (req)

=item I<attr> - attribute to set (req)

=item I<class> - key in the opts hashref to set attribute (opt, ADDITIONAL_OPTS)

=item I<cat> - category to put the opt in the usage (opt, Optional)

=item I<desc> - usage description (opt, but recomended)

=back

=back

=head1 Methods

=head2 script_basename

Returns the basename of the script

=head2 usage_categories

An array of the usage categories.

=head2 get_options

Returns a hashref of keys = classes, values = hashref of params to create 
the object.  If the add_q param was used and the q option was specified
on the command line, the script will be re-executed and sent to the LSF queue.

=head2 command_line_options

The Getopt::Long params

=head2 usage

Return the usage (string) for the classes and additional options.

=head1 See Also

I<Finfo::Std>, I<Finfo::Object>, I<PP::LSF>, I<Getopt::Long>

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
