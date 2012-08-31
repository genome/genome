package Genome::Model::Tools::PhredPhrap::ScfToPhd;

use strict;
use warnings;
  
use Genome;

class Genome::Model::Tools::PhredPhrap::ScfToPhd {
    is => 'Command',
    has => [
        scf_file => {
            is => 'String', #file_r
            doc => 'File of SCF',
            is_optional => 0,
        }, 
        chromat_dir => { 
            is => 'String', #dir_r
            doc => 'Directory where the SCFs are located',
            is_optional => 0,
        },
        phd_file => { 
            is => 'String', #file_w
            doc => 'File to write the most recent phd for each SCF',
            is_optional => 0,
        },
        phd_dir => {
            is => 'String', #dir_rw,
            is_optional => 0,
            doc => 'Directory to put PHDs.',
        },
        rmphd => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'Remove all phds in phd_dir, then run phred on each SCF.',
        },
        recall => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'Run phred on the SCF, naming it to a new version.',
        },
        _phd_dir => { is_optional => 1, },
    ],
};

sub help_brief {
    return 'Runs phred on SCFs to create PHDs.';
}

sub execute {
    my $self = shift;

    $self->_phd_dir( Genome::Model::Tools::Consed::PhdDirectory->create(directory => $self->phd_dir) );

    $self->_remove_phds if $self->rmphd;
    
    my $scf_fh = IO::File->new('<' . $self->scf_file)
        or ($self->error_message('Could not open scf file: ' . $self->scf_file) and return);
    my $phd_fh = IO::File->new('>' . $self->phd_file)
        or ($self->error_message( sprintr('Can\'t open file (%s) for writing', $self->phd_file)) and return);

    while ( my $scf_name = $scf_fh->getline ) {
        chomp $scf_name;

        my $phd_name;
        if ( $self->recall ) {
            $phd_name = $self->_run_phred($scf_name);
        }
        else {
            $phd_name = $self->_phd_dir->latest_phd_name($scf_name);
            $phd_name = $self->_run_phred($scf_name) unless $phd_name; 
        }

        next unless $phd_name;

        $self->_check_and_add_read_type($scf_name);

        $phd_fh->print("$phd_name\n");
    }

    $scf_fh->close;
    $phd_fh->close;

    ($self->error_message("No phds found") and return) unless -s $self->phd_file;

    return 1;
}

sub _remove_phds
{
    my $self = shift;

    return 1;

    my $dir = IO::Dir->new($self->phd_dir);
    while ( defined ($_ = $dir->read) ) {
        unlink $_ if $_ =~ /\.phd\.\d+/;
    }

    return 1;
}


sub _run_phred {
    my ($self, $scf_name) = @_;

    my $scf_file = sprintf('%s/%s', $self->chromat_dir, $scf_name);
    $self->error_message("Can't find scf ($scf_file\[.gz\])")
        and return unless -s $scf_file or -s "$scf_file.gz";
    
    my $phd_name = $self->_phd_dir->next_phd_name($scf_name);
    my $phd_file = $self->_phd_dir->phd_file($phd_name);
    my $command = "phred $scf_file -zt /tmp -nocall -p $phd_file";
    system "$command";

    unless ( -s $phd_file ) {
        # Retry w/ process_nomatch
        system "$command -process_nomatch";
    }

    $self->error_message("Phred failed on $scf_file")
        and return unless -s $phd_file;

    return $phd_name;
}

sub _check_and_add_read_type { # This will only work for GSC traces
    my ($self, $scf_name) = @_;

    my $phd = $self->_phd_dir->latest_phd($scf_name);
    return if not $phd;
    return 1 if $phd->{wr} and @{$phd->{wr}};

    my $read_name = $scf_name;
    $read_name =~ s/\.gz//;
    my ($template, $ext) = split(/\./, $scf_name);
    $ext =~ s/\d+//;

    my  ($primer_type, $lib);
    if ( $template =~ /PCR/ ) {
        $primer_type = 'pcr end';
    }
    else {
        if ( $template =~ s/_\w\d+$// ) {
            $primer_type = 'walk';
        }
        elsif ( $ext =~ /[ryg]/ ) {
            $primer_type = 'univ rev'
        }
        else {
            $primer_type = 'univ fwd'
        }
        $lib = $template;
        for (1..3) { chop $lib; }
    }

    my $date = $self->_date;

    my $wr = <<WR_TAG;
WR{
template dscript $date
name: $template
lib: $lib
}

WR{
primer dscript $date
type: $primer_type
}

WR_TAG

    my $phd_file = $self->_phd_dir->latest_phd_file($scf_name);
    my $fh = eval{ Genome::Sys->open_file_for_appending($phd_file); };
    if ( not defined $fh ) {
        Carp::confess(
            $self->error_message("Cannot open phd file ($phd_file) to append: $@")
        );
    }
    $fh->print($wr);
    $fh->close;

    return 1;
}

sub _date {
    my $self = shift;

    return $self->{_date} if defined $self->{_date};

    my ($nSecond, $nMinute, $nHour, $nDayInMonth, $nMonth, $nYear, $wday, $yday, $isdst ) = localtime;

    my $date = sprintf( 
        "%02d%02d%02d:%02d%02d%02d",
        ($nYear % 100),
        $nMonth + 1,
        $nDayInMonth,
        $nHour,
        $nMinute,
        $nSecond,
    );

    return $self->{_date} = $date;
}
1;

=pod

=head1 Name

=head1 Synopsis

=head1 Methods

=head1 Disclaimer

 Copyright (C) 2006 Washington University Genome Sequencing Center

 This module is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 License for more details.

=head1 Author(s)

 Eddie Belter <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
