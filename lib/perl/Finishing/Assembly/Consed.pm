package Finishing::Assembly::Consed;

use warnings;
no warnings 'reserved';
use strict;

use Finfo::Std;

require Cwd;
require Data::Dumper;
require File::Basename;
require File::Copy;
use IO::File;
use IO::Handle;
require Sys::Hostname;
use XML::Simple;

my $conf;
# Load conf file
{
    my $conf_file = '/gsc/scripts/share/consed/consed-conf.xml';
    my $fh = IO::File->new('<' . $conf_file);
    __PACKAGE__->fatal_msg("Can't open conf file ($conf_file): $!")
        and return unless $fh;

    my $xs = XML::Simple->new
    (
        rootname => 'configuration',
        KeyAttr => 
        {
            dir => 'name',
            basename => 'name', 
            pkg => 'name',
            pkg_attr => 'name', 
            arg => 'name',
        },
        ForceArray => [qw/ dir basename pkg pkg_attr arg /],
        ContentKey => '-content',
        VarAttr => 'name',
    );

    $conf = $xs->XMLin($fh);
}

# Attrs
my %cwd :name(cwd:o) 
    :default(0)
    :clo('cwd')
    :desc('Change into the base directory of the acefile before running consed');
my %version :name(version:o)
    :isa([ 'in_list', __PACKAGE__->valid_versions ])
    :default(__PACKAGE__->default_version)
    :clo(sprintf('version=s', join('|', __PACKAGE__::valid_versions)))
    :desc(sprintf('Consed version to run (default - %s): %s', __PACKAGE__->default_version, join(', ', __PACKAGE__::valid_versions)));
            
# Privs            
my %_lock_file :name(_lock_file:p)
    :default('../CONSED_BUSY0');
my %_pwd :name(_pwd:p);
my %_cmd :name(_cmd:p);
my %_project :name(_project:p);
my %_login :name(_login:p);

# Args for consed
my @args = (qw/ ace add_new_reads auto_edit auto_finish do_exp new_af nav nophd socket /);
my %ace :name(ace:o) 
    :isa(file_r) 
    :clo('ace=s') 
    :desc('Ace file to open');
my %addNewReads :name(add_new_reads:o) 
    :isa(file_r) 
    :clo('addNewReads=s')
    :desc('Adds new reads to acefile and saves to next highest version');
my %auto_edit :name(auto_edit:o)
    :clo(autoEdit)
    :desc('Auto edits an acefile and saves to next higher version');
my %auto_finish :name(auto_finish:o) 
    :clo(autofinish) 
    :desc('Runs autofinish and saves acefile to next hightest version');
my %auto_report :name(auto_report:o)
    :clo(auto_report)
    :desc('Run consed auto report');
my %add_singletons :name(add_singletons:o)
    :isa(file_r)
    :clo('addSingletons=s')
    :desc('');
my %do_exp :name(do_exp:o) 
    :clo(doExperiments) 
    :desc('Add experiments when running autofinish');
my %new_af  :name(new_af:o) 
    :isa(file_w) 
    :clo('newAceFilename=s') 
    :desc('Write new acefile to this file name, instead next hightest version');
my %nav :name(nav:o) 
    :isa(file_r)
    :clo('nav=s')
    :desc('Open consed and this nav file');
my %nophd :name(nophd:o)
    :clo('nophd')
    :desc('Open consed without using phd files');
my %socket :name(socket:o)
    :isa('int <> 1024 64000')
    :clo('socket=s')
    :desc('Open consed listening to port');


sub valid_versions
{
    my $self = shift;
    $self->fatal_msg("Configuration not initialized, please 'use' this class to set config") unless $conf;
    return sort keys %{ $conf->{pkg} };
}

sub default_version
{
    my $self = shift;
    $self->fatal_msg("Configuration not initialized, please 'use' this class to set config") unless $conf;
    return $conf->{default}->{pkg};
}

sub execute
{
    my $self = shift;

    $self->fatal_msg("Configuration not initialized, please 'use' this class to set config") unless $conf;

    $self->_check_edit_dir
        or return;

    $self->_init_lock
        or return;

    $self->_init_env_vars
        or return;

    my $command = $self->_build_command
        or return;

    #print Data::Dumper::Dumper($command); return 1;
    
    my $rv = system($command);

    $self->_handle_consed_return_value($rv)
        or return;

    return 1;
}

sub DEMOLISH
{
    my $self = shift;

    if ( $self->cwd )
    {
        chdir $self->_pwd;
    }

    $self->_remove_lock_file;

    return 1;
}

# edit_dir
sub _check_edit_dir : PRIVATE
{
    my $self = shift;

   if ( $self->cwd )
    {
        $self->fatal_msg("Need acefile to use attribute cwd") unless $self->ace;

        my ($ace, $dir) = File::Basename::fileparse( $self->ace );
        $self->fatal_msg( sprintf('Can\'t get the directory for acefile (%s): ', $self->ace, $!) ) unless $dir;

        $self->fatal_msg( sprintf('Acefile (%s) is not in an edit_dir', $self->ace) ) unless $dir =~ m#edit_dir/?$#;
        
        $self->info_msg($dir);
        
        chdir $dir
            or $self->fatal_msg("Can't chdir into $dir\: $!");
    }

    my $pwd = Cwd::getcwd();
    $self->fatal_msg("Can't get current working directory") unless $pwd;

    $self->_pwd($pwd);
 
    my @dirs = split(m{/}, $pwd);
    unless ( defined $dirs[-1] and $dirs[-1] eq "edit_dir") 
    {
        $self->fatal_msg("Consed must be run in the edit_dir");
    }

    #print Data::Dumper::Dumper(\@dirs);
    
    $self->_project( $dirs[-2] );

    return 1;
}

# Foreground/promt for responses
sub _foreground_and_prompt : PRIVATE
{
    my $self = shift;

    return 
    ( 
        $conf->{pkg}->{ $self->version }->{foreground} 
            and $conf->{pkg}->{ $self->version }->{prompt} 
    );
}

# Lock
sub _use_lock : PRIVATE
{
    my $self = shift;

    return $conf->{ $self->version }->{lock};
}

sub _init_lock : PRIVATE
{
    my $self = shift;

    my $login = getpwuid($<);
    $self->_login($login);

    if ( $self->_use_lock ) 
    {
        if ( -f $self->_lock_file ) 
        {
            my $owner = getpwuid( (stat( $self->_lock_file) )[4] );
            if ($login ne $owner) 
            {
                # "\n-------------------------------\n\n",
                $self->fatal_msg
                (
                    "Sorry, this project is locked by $owner.  Contact that user or finsupport for help.",
                );
            }
            else 
            {
                $self->warn_msg
                (
                    "WARNING! You have another assembly open in this project\n",
                );
            }
        }
        else 
        {
            my $fh = IO::File->new('>' . $self->_lock_file);
            unless ( $fh )
            {
                $self->fatal_msg( sprintf("Failed to create lock file (%s):\n%s", $self->_lock_file, $!) );
                return;
            }
            $fh->print("$login\: $$\n");
            $fh->close;
        }
    }

    return 1;
}

sub _remove_lock_file : PRIVATE
{
    my $self = shift;

    if ( $self->_use_lock )
    {
        unlink( $self->_lock_file ) if -f $self->_lock_file;

        $self->warn_msg
        (
            sprintf("Failed to remove lock file (%s):\n", $self->_lock_file, $!)
        ) 
        if -f $self->_lock_file;
    }
 
    return 1;
}

# Env Vars
sub _init_env_vars :PRIVATE
{
    my $self = shift;

    $ENV{CONSED_HOME} ||= $conf->{dir}->{path};
    $ENV{CONSED_PARAMETERS} ||= sprintf
    (
        '%s/%s', 
        $conf->{dir}->{share}, $conf->{pkg}->{$self->version}->{rc}
    );

    return 1;
}

sub _build_command : PRIVATE
{
    my $self = shift;

    if ( $self->_foreground_and_prompt ) 
    {
        $self->warn_msg
        (
            "Consed will run in foreground so if it crashes, you can enter information about it. Ok? (y/n) [y]"
        );
        my $ans = <STDIN>;
        chomp($ans);
        if ($ans =~ m/^[Nn]/) 
        {
            $self->info_msg("Running in BACKground");
            $conf->{foreground} = 0;
        }
        else
        {
            $self->info_msg("Running in FOREground");
        }
    }

    unless ( $conf->{pkg}->{ $self->version }->{foreground} ) 
    {
        my $pid = fork();
        if ( $pid ) 
        {
            # parent just exits
            exit 0;
        }
        elsif ( defined $pid ) 
        {
            # child falls through
        }
        else 
        {
            $self->fatal_msg("Failed to run consed in background: $!");
            return;
        }
    }

    # construct consed command
    my $exe =  $conf->{pkg}->{ $self->version }->{exe};
    #my $exe = sprintf('%s/%s', $conf->{dir}->{path}, $conf->{pkg}->{ $self->version }->{exe} );

    my $link_count = 0;
    while (-l $exe) 
    {
        my $target = readlink($exe);
        unless ( $target ) 
        {
            $self->warn_msg("Failed to read link: $exe");
            last;
        }
        # update exe
        $exe = ($target =~ m{^/}) 
        ? $target
        : sprintf('%s/%s', $conf->{dir}->{path}, $target);

        # don't get caught in a self-referential link
        last if ++$link_count > 10;
    }

    return unless Finfo::Validate->validate
    (
        attr => 'consed executible for version ' . $self->version,
        value => $exe,
        type => 'exe',
        err_cb => $self,
    );

    my $cmd = $exe;
    foreach my $arg ( @args )
    {
        my $opt = $self->attributes_attribute($arg, 'clo');
        my $value = $self->$arg;
        next unless $value;
        $value = '' unless  $opt =~ s/=.*$//;

        $cmd .= sprintf(' -%s %s', $opt, $value);
    }

    if ( $conf->{core} ) 
    {
        $cmd = "ulimit -c 6297600; $cmd";
    }

    return $self->_cmd($cmd);
}

sub _handle_consed_return_value : PRIVATE
{
    my ($self, $rv) = @_;

    if ( $rv != 0 )
    {
        if ( $rv == -1 )
        {
            $self->fatal_msg( sprintf("Failed to execute consed (%s):\n", $self->_cmd, $!) );
            $self->_return_value(1);
            return 1;
        }

        $rv >>= 8; # exit value from consed seems to be 8 bits off
        my $es = $rv >> 8;
        #print Data::Dumper::Dumper([$rv, $es]);

        # check for signal and core
        # TODO this doesn't work properly for exiting with a cancel...
        if ( $rv & 127 ) 
        {
            my $signal = $rv & 127;
            my $core = $rv & 128;
            $self->error_msg
            (
                sprintf
                (
                    "Consed (%s) did not exit normally - exit status(%d), signal (%d)",
                    $self->_cmd,
                    $es,
                    $signal,
                )
            );

            # save cores
            if ( $core and $conf->{pkg}->{ $self->version }->{core} ) 
            {
                # save some information
                my $uname = qx(uname -a);
                chomp $uname;
                my $desc = sprintf
                (
                    "User:    %s\nProject: %s\nPath: %s\nCommand: %s\nSignal:  $%s\nExit:    %s\nUname:   %s\n\n",
                    $self->_login,
                    $self->_project,
                    $self->_pwd,
                    $self->_cmd,
                    $signal,
                    $es,
                    $uname,
                );

                # Collect information from user if running in foreground
                if ( $self->_foreground_and_prompt ) 
                {
                    $self->info_msg("Consed dumped a core");
                    $self->info_msg("Please provide description of what you were doing when consed crashed");
                    $self->info_msg("Type Ctrl-D on empty line when finished.\n");
                    
                    $desc .= join('', "User Description:\n", STDIN->getlines);
                }

                # Create desc file
                my $base = sprintf
                (
                    '/gsc/var/lib/consed/%s',
                    join
                    (
                        '-',
                        $self->version,
                        $conf->{pkg}->{ $self->version }->{exe}, 
                        Sys::Hostname::hostname, 
                        $self->_login, 
                        $$
                    ),
                );

                while (-e "$base.desc") { $base .= 'x' }
                my $desc_file = "$base.desc";
                my $desc_fh = IO::File->new(">$desc_file");
                if ( $desc_fh )
                {
                    $desc_fh->print($desc);
                    $desc_fh->close;
                }
                else 
                {
                    $self->warn_msg("Failed to create description file ($desc_file): $!");
                }

                # save core file
                if ( -f 'core' ) 
                {
                    $self->info_msg("Copying core file...");
                    my $core_file = "$base.core";
                    unless ( File::Copy::move('core', $core_file) )
                    {
                        $self->warn_msg("Failed to save core file: $!") 
                    }
                    $self->info_msg("Done");
                }
                else 
                {
                    $self->warn_msg("Core file does not exist: $!");
                }
            }

            return 1;
        }
        else 
        {
            $self->fatal_msg
            (
                sprintf('Consed (%s) exited with non-zero status (%s)', $self->_cmd, $es)
            );

            $self->_exit_status(-$es);
            
            return 1;
        }
    }

    return 1;
}

1;

=pod

=head1 Name

Finishing::Assembly::Consed::Run

=head1 Synopsis

Manages and runs the many versions of consed.

=head1 Usage

 my $consed = Finishing::Assembly::Consed::Run->new(%params)
    or die;
 my $exit_status = $consed->execute;
 exit($exit_status);

=head2 Params
 
=over

=item I<conf_file>

=item I<ace> No description available
=item I<addNewReads> Adds new reads to acefile and saves to next highest version
=item I<addSingletons> No description available
=item I<autoEdit> Auto edits an acefile and saves to next higher version
=item I<auto_report> Run consed auto report
=item I<autofinish> Runs autofinish and saves acefile to next hightest version
=item I<consed|consed_old|consed_david|consedmp|consed_davidmp Consed version to run: consed, consed_old, consed_david, consedmp, consed_davidmp> Version of consed to run, default: cs
=item I<doExperiments> Add experiments when running autofinish
=item I<nav> Open consed and this nav file
=item I<newAceFilename> Write new acefile to this file name, instead next hightest version
=item I<nophd> Open consed without using phd files
=item I<socket> Open consed listening to port

=back

=head1 Methods

=head2 execute

=over 

=item I<Synopsis>   Executes consed and returns the exit value.

=item I<Params>     none

=item I<Returns>    Consed exit code, so 0 is successful.

=back

=head1 See Also

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABI
LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Consed.pm $
#$Id: Consed.pm 29594 2007-10-29 18:31:27Z ebelter $
