package Genome::Model::Tools::Consed::ExternalNavigation;

use strict;
use warnings;

use Genome;

use Data::Dumper;
require File::Basename;
require IO::Socket;

our $initialized = 0;

sub require_requires {
    return if $initialized;
    eval {
        require Gtk2;
        Gtk2->init;
        require Gtk2::Ext::Dialogs;
        require Gtk2::Ext::EntryCrate;
        require Gtk2::Ext::PackingFactory;
        require Gtk2::Ext::Utils;
        require POE::Session;
        POE::Session->import;
        require POE::Kernel;
        POE::Kernel->import( { loop => 'Glib' });
    };
    die $@ if $@;
    $initialized = 1;
}

class Genome::Model::Tools::Consed::ExternalNavigation { 
    is => 'Command',
    has => [
    acenav => { 
        is => 'Text',
        doc => 'Multi ace navigation file, use \'nav2acenav\' to convert a regular nav file to a acenav file',
    },
    ],
    has_optional => [
    break => {
        is => 'Integer',
        doc => 'Break each navigation into sections of (x) bases each',
    },
    ]
};

# _port :default(1024);

sub execute {
    require_requires() unless $initialized;
    my $self = shift;

    $self->{_preferences} = {
        consed => 'cs',
        interval => 5, 
        use_phds => 'no', 
        warn_when_opening_acefile => 'no',
    };
    
    POE::Session->create (
        inline_states => {
            _start => sub{ shift; $self->_ui_start(@_) },
            ev_goto => sub{ shift; $self->_ui_goto(@_) },
            ev_prev => sub{ shift; $self->_ui_prev(@_) },
            ev_next => sub{ shift; $self->_ui_next(@_) },
            ev_run => sub{ shift; $self->_ui_run(@_) },
            ev_stop => sub{ shift; $self->_ui_stop(@_) },
        }
    )
        or return;

    POE::Kernel->run();

    return 1;
}

sub gtk2_utils {
    return Gtk2::Ext::Utils->instance;
}

sub gtk2_dialogs {
    return Gtk2::Ext::Dialogs->instance;
}

sub factory {
    return Gtk2::Ext::PackingFactory->instance;
}

sub _ui_start {
    my ($self, $session, $kernel) = @_[ &OBJECT, &SESSION, &KERNEL ];

    my $factory = $self->factory;

    my $window = $factory->create_window(
        title => 'Consed External Navigator',
        h => 700,
        v => 350,
        border_width => 15,
    )
        or return;

    $kernel->signal_ui_destroy($window);

    my $vbox = $factory->add_box(parent => $window, type => 'v');

    my $menu = $factory->add_menu(
        parent => $vbox,
        expand => 0,
        fill => 0,
        menu_tree => [
        File => {
            item_type => '<Branch>',
            children => [
            Quit => {
                item_type => '<StockItem>',
                extra_data => 'gtk-quit',
                callback => sub{ $window->destroy; $self->gtk2_utils->gtk2_quit; exit 0 },
                accelerator => '<ctrl>Q',
            },
            ],
        },
        Options => {
            item_type => '<Branch>',
            children => [
            'Reopen Ace' => {
                callback => sub{ $self->_reopen_ace },
                item_type => '<StockItem>',
                extra_data => 'gtk-refresh',
            },
            Separator => {
                item_type => '<Separator>'
            },
            'Preferences' => {
                callback => sub{ $self->_change_preferences_select },
                item_type => '<StockItem>',
                extra_data => 'gtk-preferences',
            },
            ],
        },
        Help => {
            item_type => '<LastBranch>',
            children => [
            About => {
                callback => sub{},
                item_type => '<StockItem>',
                extra_data => 'gtk-about',
            },
            'Wiki Fin FAQ' => {
                callback => sub
                {
                    $self->gtk2_dialogs->info_dialog("Launching firefox, please wait");
                    my $pid = fork();
                    if ($pid) { 
                        # parent falls thru
                    }
                    elsif (defined($pid)) {
                        exec("firefox http://gscweb.gsc.wustl.edu/wiki/External_nav_Fin_FAQs")
                            or $self->error_msg("Could not launch firefox: $!");
                    }
                    else {
                        $self->gtk2_error_dialog("Could not fork firefox process, bummer");
                    }
                },
                item_type => '<StockItem>',
                extra_data => 'gtk-info',
            },
            ],
        },
        ],
    );

    $window->add_accel_group( $menu->{accel_group} );

    $self->{_slist_frame} = $factory->add_frame(
        parent => $vbox,
        text => "Current File: " . File::Basename::basename($self->acenav),
    );

    $self->{_slist} = $factory->add_slist(
        parent => $factory->add_sw (
            parent => $self->{_slist_frame},
            h => 700,
            v => 150,
            fill => 1,
            expand => 1,
        ),
        columns => [qw/ Dir text Ace text Contig text Position int Comment text /],
        events => {
            'row_activated' => $session->callback('ev_goto'),
        },
    );

    $self->{_slist}->get_column(0)->set_visible(0);

    $self->_load_acenav
        or return;

    my $bbox = $factory->add_bbox(
        parent => $factory->add_frame(
            parent => $vbox,
            shadow => 'in',
            expand => 0,
        ),
        type => 'h', 
        layout => 'spread',
        border_width => 8,
        homogen => 1,
    );

    my @controls = (
        'goto' => 'gtk-jump-to',
        'next' => 'gtk-go-forward',
        run => 'gtk-media-play',
        stop => 'gtk-stop',
    );
    for (my $i = 0; $i <= $#controls; $i += 2) {
        $factory->add_button(
            parent => $bbox,
            stock => $controls[$i + 1],
            events => { 'clicked' => $session->callback('ev_'.$controls[$i]) },
        );
    }

    $self->gtk2_utils->gtk2_main;

    return 1;
}

sub _is_acefile_open {
    return $_[0]->{_acefile_ports}->{$_[1].'/'.$_[2]};
}

sub _get_socket_for_acefile {
    my ($self, $dir, $ace) = @_;

    unless ( $self->_is_acefile_open($dir, $ace) ) {
        $self->{_acefile_ports}->{$dir.'/'.$ace} = $self->_next_port_number;
        $self->_open_acefile($dir, $ace);
    }

    return IO::Socket::INET->new('localhost:' .  $self->{_acefile_ports}->{$dir.'/'.$ace});
}

sub _open_acefile {
    my ($self, $dir, $ace) = @_;

    unless ( -f $dir.'/'.$ace ) { 
        $self->gtk2_dialogs->error_dialog("Acefile ($dir/$ace) does not exist!");
        return;
    }
    
    if ( $self->{_preferences}->{warn_when_opening_acefile} =~ /^y/i ) {
        return unless $self->gtk2_dialogs->question_dialog("Continue to open acefile?") eq 'yes';
    }

    unless ( chdir $dir ) {
        $self->gtk2_dialogs->error_dialog("Cannot change to acefile directory ($dir): $!");
        return;
    }
    
    system(
        sprintf(
            "%s -ace %s -socket %s %s",
            $self->{_preferences}->{consed},
            $ace,
            $self->{_acefile_ports}->{$dir.'/'.$ace},
            ( $self->{_preferences}->{use_phds} =~ /^n/ ) ? '-nophd' : '',
        )
    );

    return 1;
}

sub _reopen_ace {
    my $self = shift;

    my ($row_data) = $self->gtk2_utils->get_selected_data_from_slist($self->{_slist});

    unless ( $row_data ) {
        $self->gtk2_dialogs->info_dialog("Please select a location");
        return;
    }

    my ($dir, $ace) = ($row_data->[0], $row_data->[1]);

    if ( $self->_is_acefile_open($dir, $ace) ) {
        if ( 
            $self->gtk2_dialogs->question_dialog ("Ace ($ace) is already open, really reopen?") eq 'yes'
        ) {
            return $self->_open_acefile($dir, $ace);
        }
    }
}

sub _next_port_number {
    my $self = shift;

    my $port = $self->{_port} || 5000;

    return $self->{_port} = ++$port;
}

sub _load_acenav {
    my $self = shift;

    my $reader = Genome::Model::Tools::Consed::Navigation::Reader->new(
        input => $self->acenav,
    )
        or return;
    
    my $count = -1;
    while ( my $nav = $reader->next ) {
        $self->fatal_msg(
            sprintf(
                'Each navigation in file (%s) is required to have an acefile.  Please add acefiles to this navigation file manually or by using \'nav2acenav\'', 
                $self->acenav
            )
        ) unless $nav->{acefile};
        
        $count++;
        my $start = $nav->{start};
        while ( 1 ) {
            my ($ace, $dir) = File::Basename::fileparse( $nav->{acefile} );

            $self->gtk2_utils->add_data_to_slist(
                $self->{_slist},
                [ $dir, $ace, $nav->{contig}, $start, $nav->{comment} ],
            );
            last unless $self->break; # not breaking into segments
            $start += $self->break;
            last if $start > $nav->{stop};
        }
    }

    $self->{_max} = $count;

    return 1;
}

sub _change_preferences_select {
    my ($self) = @_;

    my @ecrates = (
        Gtk2::Ext::EntryCrate->new(
            name => 'consed', 
            label => 'Consed Version', 
            default => $self->{_preferences}->{consed},
        ),
        Gtk2::Ext::EntryCrate->new(
            name => 'interval', 
            label => 'Time Interval', 
            is => 'int >= 1',
            default => $self->{_preferences}->{interval},
        ),
        Gtk2::Ext::EntryCrate->new(
            name => 'warn_when_opening_acefile', 
            label => 'Warn when opening an acefile?', 
            is => 'y_or_n',
            default => $self->{_preferences}->{warn_when_opening_acefile},
        ),
        Gtk2::Ext::EntryCrate->new(
            name => 'use_phds', 
            label => 'Open acefile with phds?', 
            is => 'y_or_n',
            default => $self->{_preferences}->{use_phds},
        ),
    );

    my $values = $self->gtk2_dialogs->ecrate_dialog(
        title => 'Set Preferences',
        ecrates => \@ecrates,
    );

    return unless defined $values;

    return $self->{_preferences} = ($values);
}

sub _scroll {
    my ($self, $dir, $ace, $ctg, $pos) = @_;

    my $af = $dir.'/'.$ace;
    
    my $socket = $self->_get_socket_for_acefile($dir, $ace)
        or return;
    
    return print $socket "Scroll $ctg $pos\n";

    $socket->close;

    return 1;
}

sub _ui_goto {
    my $self = shift;

    my ($row_data) = $self->gtk2_utils->get_selected_data_from_slist($self->{_slist});
    
    unless ( defined $row_data ) {
        $self->gtk2_dialogs->info_dialog("Please select a location");
        return;
    }

    return $self->_scroll(@$row_data);
}

sub _ui_next {
    my ($self) = @_;

    my ($row) = $self->gtk2_utils->get_selected_indices_from_slist($self->{_slist});

    $row = -1 unless defined $row;
    
    return if $row == $self->{_max};

    $row++;

    $self->{_slist}->set_cursor( Gtk2::TreePath->new($row) );

    return _ui_goto(@_);
}

sub _ui_prev {
    my ($self) = @_;

    my ($row) = $self->gtk2_utils->get_selected_indices_from_slist($self->{_slist});

    return if not defined $row or $row == 0;

    $row--;

    $self->{_slist}->set_cursor( Gtk2::TreePath->new($row) );

    return _ui_goto(@_);
}

sub _ui_run {
    my ($self, $kernel, $heap) = @_[ 0, &KERNEL, &HEAP ];

    $kernel->yield('ev_next');

    my ($row) = $self->gtk2_utils->get_selected_indices_from_slist($self->{_slist});

    if ( $row == $self->{_max} ) {
        _ui_stop(@_);
        $self->gtk2_dialogs->info_dialog("Reached end of navigator");
    }
    else {
        $heap->{alarm_id} = $kernel->delay_set("ev_run", $self->{_preferences}->{interval});
    }

    return 1;
}

sub _ui_stop
{
    my ($self, $kernel, $heap) = @_[ 0, &KERNEL, &HEAP ];

    $kernel->alarm_remove($heap->{alarm_id});
    
    return 1;
}

1;

=pod

=head1 Name

Genome::Model::Tools::Consed::Navigation::ExternalNavigation

=head1 Synopsis

=head1 Usage

=head1 Methods

=head2 execute

 $ex_nav->execute;

=over

=item I<Synopsis>   

=item I<Params>     

=item I<Returns>   

=back

=head1 See Also

=over

=item consed

=item Genome::Model::Tools::Consed, Genome::Model::Tools::... 

=item Gtk2

=back

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$

