package Genome::Model::Command::Report::List;

use strict;
use warnings;

use Genome;

use Carp 'confess';
use Data::Dumper 'Dumper';
require Term::ANSIColor;
require Text::Wrap;
use Genome::Model::Report;

class Genome::Model::Command::Report::List {
    is => 'Command',
    has_optional => [
        all => {
            is => 'Boolean',
            doc => 'List generic and model type report names.  This is the default if no other options are indicated.',
        },
        type_names => {
            is => 'Boolean',
            doc => 'List report names for all model types.',
        },
        generic => {
            is => 'Boolean',
            doc => 'Reprot names that are generic and can be run on any model type.',
        },
        type_name => {
            is => 'Text', 
            is_input => 1,
            doc => 'Report names for this model type only.' 
        },
        build => { 
            is => 'Genome::Model::Build', 
            id_by => 'build_id',
        },
        build_id => {
            is => 'Text',
            is_input => 1,
            doc => 'The build id of a genome model to check which reports have been run.',
        },
    ],
};

###################################################

sub help_brief {
    return 'Lists reports for model type names';
}

sub help_synopsis {
    return <<"EOS"
List all models' available reports
 genome model report list

List reports for reference alignment models
 genome model report list --type_name 'reference alignment'

List reports for a build
 genome model report list --build-id 96426120

EOS
}

sub help_detail {
    return help_brief();
}

###################################################

sub execute {
    my $self = shift;

    my @requested_listers = grep { $self->$_ } (qw/ all type_names generic type_name build_id /);
    my $method = '_list_reports_for_';
    unless ( @requested_listers ) {
        $self->all(1);
        $method .= 'all';
    }
    elsif ( @requested_listers > 1 ) {
        $self->error_message('Please indicate one listing method.');
        return;
    }
    else {
        $method .= $requested_listers[0];
    }
    
    return $self->$method;
}

#< Helpers >#
sub _print_basic_header {
    print Term::ANSIColor::colored('Report Names', 'bold')."\n";
}

sub _print_type_and_report_names {
    my ($self, $type, @report_names) = @_;

    @report_names = (qw/ none /) unless @report_names;
    
    return print Text::Wrap::wrap(
        ' ',
        '  ',
        Term::ANSIColor::colored($type, 'red'), 
        "\n",
        join("\n", @report_names),
    )."\n";
}

#< All >#
sub _list_reports_for_all {
    my $self = shift;

    $self->_print_basic_header;
    $self->_print_reports_for_generic;
    $self->_print_reports_for_type_names;

    return 1;
}

#< Generic >#
sub _list_reports_for_generic {
    my $self = shift;

    $self->_print_basic_header;
    $self->_print_reports_for_generic;

    return 1;
}

sub _print_reports_for_generic {
    $_[0]->_print_type_and_report_names(generic => Genome::Model::Report::get_generic_report_names());
}

#< Type Names >#
sub _list_reports_for_type_name {
    my $self = shift;

    my $type_name = $self->type_name;
    unless ( grep { $type_name eq $_ } Genome::Model::Command::BaseDeprecated::get_model_type_names() ) {
        $self->error_message("Invalid type name ($type_name).");
        return;
    }
    
    $self->_print_basic_header;
    $self->_print_reports_for_type_name($self->type_name);

    return 1;
}

sub _print_reports_for_type_name {
    $_[0]->_print_type_and_report_names($_[1] => Genome::Model::Report::get_report_names_for_type_name($_[1]));
}

#< Type Names >#
sub _list_reports_for_type_names {
    my $self = shift;

    $self->_print_basic_header;
    $self->_print_reports_for_type_names;

    return 1;
}

sub _print_reports_for_type_names {
    my $self = shift;

    for my $type_name ( Genome::Model::Command::BaseDeprecated::get_model_type_names() ) {
        $self->_print_type_and_report_names(
            $type_name, Genome::Model::Report::get_report_names_for_type_name($type_name)
        );
    }

    return 1;
}

#< Build's Reports >#
sub _list_reports_for_build {
    my $self = shift;
    
    my $type_name = $self->build->model->type_name;
    my @availble_report_generators = Genome::Model::Report::get_report_names_for_type_name($type_name);

    unless ( @availble_report_generators ) { # ok
        print "Model type ($type_name) does not have any reports.\n";
        return 1;
    }
    
    my %reports = map { $_ => [] } @availble_report_generators;
    for my $report ( $self->build->reports ) {
        my ($report_subclass) = $report->get_generator =~ m#::([\w\d]+)$#;
        my $report_generator = Genome::Utility::Text::camel_case_to_string($report_subclass);
        push @{$reports{$report_generator}}, $report;
    }

    #print Dumper([keys %reports]); print $self->build->resolve_reports_directory,"\n";

    $self->_print_header_for_build
        or return;

    for my $report_name ( @availble_report_generators ) {
        $self->_print_reports($report_name, $reports{$report_name})
            or return;
    }

    return 1;
}

sub _print_header_for_build {
    my $self = shift;

    return print(
        Term::ANSIColor::colored(
            sprintf('Reports for %s Build (<Id> %s)', 
                join('', map { ucfirst } split(/\s+/, $self->build->model->type_name)),
                $self->build->id
            ), 
            'bold'
        ),
        "\n",
    );
}

sub _print_reports {
    my ($self, $report_name, $reports) = @_;

    my @strings = Term::ANSIColor::colored(sprintf('%s (%s)', $report_name, scalar(@$reports)), 'red');
    
    if ( @$reports ) {
        for my $report ( @$reports ) {
            push @strings, sprintf("%s (%s)", $report->name, $report->get_date);
        }
    }
    else {
        push @strings, "None";
    }

    return print Text::Wrap::wrap(' ', '  ', join("\n", @strings))."\n";
}

1;

#$HeadURL$
#$Id$
