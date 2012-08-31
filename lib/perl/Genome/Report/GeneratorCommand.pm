package Genome::Report::GeneratorCommand;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require File::Basename;
require File::Temp;

class Genome::Report::GeneratorCommand {
    is => 'Command',
    has_optional => [
        print_xml => {
            is => 'Boolean',
            doc => 'Print the report XML to the screen (STDOUT).',
        },
        print_datasets => {
            is => 'Boolean',
            doc => 'Print the datasets as cvs to the screen (STDOUT). Default will be to print all datasets. Indicate specific datasets with the "datasets" option.',
        },
        datasets => {
            is => 'Text',
            doc => 'Datasets to print, save or email. Defaults to all datasets.  Note that for email, this option only effects the attachments.  All datasets will still appear in the tranformed report.  Datasets will be comma sparated.  Separate dataset names by commas.',
        },
        separator => {
            is => 'Text',
            default_value => ',',
            doc => 'Separator character to use when printing, saving and email attachments. For tab separator, use "tab".',
        },
        include_headers => {
            is => 'Boolean',
            default_value => 1,
            doc => 'Include dataset headers when printing, saving and email attachments.',
        },
        email => {
            is => 'Text',
            doc => 'Email the report to these recipients.  Separate by commas.',
        },
        save => {
            is => 'Text',
            doc => 'Save report to this directory.',
        },
        force_save => {
            is => 'Boolean',
            doc => 'Force save the report, if one already exists.',
        },
        # private
        _dataset_names => {
            is => 'Array',
        },
        _datasets_svs => {
            is => 'Hash',
        },
        _datasets_files => {
            is => 'Hash',
        },
        _sv_ext => {
            is => 'Text',
            default_value => 'csv',
        },
    ],
};

#< Helps >#
sub help_detail {
    return <<STRING;
    
STRING
}

#< Generate Report >#
sub _generate_report_and_execute_functions {
    my ($self, %params) = @_;

    # Set separator stuff.  Default are ',' and 'csv'
    if ( $self->separator eq 'tab' ) {
        $self->separator("\t");
        $self->_sv_ext('tsv');
    }
    elsif ( $self->separator ne ',' ) {
        $self->_sv_ext('txt');
    }

    # Generate report
    my $generator = Genome::Model::Report::Table->create(%params);
    my $report = $generator->generate_report;
    unless ( $report ) {
        $self->error_message("Can\'t generate report.");
        return;
    }

    my @functions = (qw/ print_xml print_datasets email save /);
    my @selected_functions = grep { $self->$_ } @functions;
    unless ( @selected_functions ) {
        $self->print_datasets(1);
        @selected_functions = (qw/ print_datasets /);
    }

    $self->_resolve_dataset_names($report)
        or return;
    
    for my $function ( @selected_functions ) {
        my $method = '_'.$function;
        $self->$method($report)
            or return;
    }

    return $report;
}

sub _resolve_dataset_names {
    my ($self, $report) = @_;

    my @dataset_names = $report->get_dataset_names;
    unless ( @dataset_names ) {
        $self->error_message("Indicated to use all datasets (for print_xml, print_datasets, email save), but this report does not have any.");
        return;
    }

    if ( $self->datasets ) { # chack they exist
        my @datasets_requested;
        for my $dataset_requested ( split(',', $self->datasets) ) {
            unless ( grep { $dataset_requested eq $_ } @dataset_names ) {
                $self->error_message("Dataset requested ($dataset_requested) not found in report (".$report->name.")");
                return;
            }
            push @datasets_requested, $dataset_requested;
        }
        $self->_dataset_names(\@datasets_requested);
    }
    else { # all
        $self->_dataset_names(\@dataset_names);
    }

    return 1;
}

#< Print XML >#
sub _print_xml {
    my ($self, $report) = @_;

    return print $report->xml_string;
}

#< Print Datasets (default) >#
sub _print_datasets {
    my ($self, $report) = @_;

    my $datasets_svs = $self->_datasets_to_svs($report)
        or return;

    for my $svs ( values %$datasets_svs ) {
        print $svs;
    }

    return 1;
}

sub _datasets_to_svs {
    my ($self, $report) = @_;

    return $self->_datasets_svs if $self->_datasets_svs;

    my $dataset_names = $self->_dataset_names;
    unless ( $dataset_names ) {
        $self->_datasets_svs({});
        return $self->_datasets_svs;
    }

    my %datasets_svs;
    for my $name ( @$dataset_names ) {
        my $ds = $report->get_dataset($name);
        unless ( $ds ) { # bad
            $self->error_message("Could not get dataset ($name) from report.");
            return;
        }

        my $svs = $ds->to_separated_value_string(
            separator => $self->separator,
            include_headers => $self->include_headers,
        );
        unless ( $svs ) {
            $self->error_message("Can't get separated value string from build dataset");
            return;
        }
        $datasets_svs{$name} = $svs;
    }

    return $self->_datasets_svs(\%datasets_svs);
}

#< Save >#
sub _save {
    my ($self, $report) = @_;

    my $dir = $self->save;
    unless ( Genome::Sys->validate_existing_directory($dir) ) {
        $self->error_message("Can't save report because of problem with directory ($dir). See above error.");
        return;
    }

    unless ( $report->save($dir) ) {
        $self->error_message("Can't save report to directory ($dir).  See above errror.");
        return 1;
    }

    $self->_save_datasets($report)
        or return;

    print "Saved report to ".$self->save."\n";

    return 1;
}

sub _save_datasets {
    my ($self, $report) = @_;

    return $self->_datasets_files if $self->_datasets_files;

    my $datasets_svs = $self->_datasets_to_svs($report)
        or return;

    my $dir = ( $self->save ? $self->save : File::Temp::tempdir(CLEANUP => 1) );
    my %datasets_files;
    for my $name ( keys %$datasets_svs ) {
        my $file = sprintf('%s/%s.%s', $dir, $name, $self->_sv_ext);
        my $fh = Genome::Sys->open_file_for_writing($file)
            or return;
        $fh->print($datasets_svs->{$name});
        $fh->close;
        $datasets_files{$name} = $file;
    }

    return $self->_datasets_files(\%datasets_files);
}

#< EMail >#
sub _email {
    my ($self, $report) = @_;

    my $datasets_files = $self->_save_datasets($report)
        or return;

    my @attachments;
    for my $name ( keys %$datasets_files ) {
        my $basename = File::Basename::basename($datasets_files->{$name});
        push @attachments, {
            description => $name,
            disposition => "inline; filename=\"$basename\";\r\nContent-ID: <$name>",
            file => $datasets_files->{$name},
        };
    }

    my $confirmation = Genome::Report::Email->send_report(
        report => $report,
        to => $self->email,
        xsl_files => [ $report->generator->get_xsl_file_for_html ],
        attachments => \@attachments,
    );

    unless ( $confirmation ) {
        $self->error_message("Can't email report.");
        return;
    }

    print "Sent report to ".$self->email."\n";

    return 1;
}

1;

=pod

=head1 Name

Genome::Report::GeneratorCommand

=head1 Synopsis

Base command class for report generator commands.  Use this class as a base for yours.

=head1 Usage

=head1 Public Methods

=head2 generate_report

 my $report = $generator->generate_report
    or die;

=over

=item I<Synopsis>   Generates data and creates a Genome::Report

=item I<Arguments>  none

=item I<Returns>    Genome::Report

=back

=head1 Disclaimer

Copyright (C) 2005 - 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
