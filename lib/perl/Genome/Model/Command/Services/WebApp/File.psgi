#!/usr/bin/env genome-perl

use Web::Simple 'Genome::Model::Command::Services::WebApp::File';
use MIME::Base64;

package Genome::Model::Command::Services::WebApp::File;

our $loaded = 0;
sub load_modules {
    return if $loaded;
    eval {
        use above 'Genome';
        use HTML::Tags;
        use Plack::MIME;
        use Plack::Util;
        use Plack::Request;
        use Cwd;
        use HTTP::Date;
        use JSON;
        use UR::Object::View::Default::Xsl qw/type_to_url url_to_type/;
    };
    if ($@) {
        die "failed to load required modules: $@";
    }

    # search's callbacks are expensive, web server can't change anything anyway so don't waste the time
    Genome::Search->unregister_callbacks();
}

sub dispatch_request {
    sub ( POST + /view/x/druggable-gene-go + %@families~&* + *file~ ) {
        load_modules();
        my ($self, $families, $params, $file, $env) = @_;

        my @gene_names;
        @gene_names = Genome::Sys->read_file($file->path) if $file;
        push @gene_names, split(/[\r\n]/,$params->{'genes'});
        chomp @gene_names;
        @gene_names = grep{/[\w\d]/}@gene_names;
        @gene_names = map{uc $_}@gene_names;

        my $command = Genome::DruggableGene::Command::GeneNameGroup::LookupFamilies->execute(
                gene_identifiers => \@gene_names,
                allowed_families => $families,
                );
        my %params = (
                data => $command->result,
                subject => Genome::DruggableGene::GeneNameReport->define_set,
                );
        for($params->{'return'}){
            if(/html/){
                my $html = Genome::DruggableGene::GeneNameReport::Set::View::Go::Html->create(
                        %params,
                        desired_perspective => 'status',
                        xsl_root => Genome->base_dir . '/xsl',
                        xsl_path => '/static/xsl',
                        xsl_variables => {
                            rest      => '/view',
                            resources => '/view/genome/resource.html',
                        },
                    );
                return [200, ['Content-type' => "text/html"], [$html->content]];
            } elsif(/tsv/) {
                return [200, ['Content-type' => "text/tsv"], [join("\n", $command->output)]];
            } elsif(/xml/) {
                my $xml = Genome::DruggableGene::GeneNameReport::Set::View::Go::Xml->create(%params);
                return [200, ['Content-type' => "text/xml"], [$xml->content]];
            }
        }
    },
    sub ( POST + /view/x/druggable-gene + %* + *file~ ) {
        load_modules();
        my ($self, $params, $file, $env) = @_;
        my @gene_names;
        @gene_names = Genome::Sys->read_file($file->path) if $file;
        push @gene_names, split(/[\r\n]/,$params->{'genes'});
        chomp @gene_names;
        @gene_names = grep{/[\w\d]/}@gene_names;
        @gene_names = map{uc $_}@gene_names;
        my @filter;
        for($params->{'filter'}){
            push @filter, 'drug.is_withdrawn=0,drug.is_nutraceutical=0,is_potentiator=0,(is_untyped=0 or is_known_action=1)' if /default/;
            push @filter, 'drug.is_withdrawn=0,drug.is_nutraceutical=0,is_potentiator=0,is_inhibitor=1,(is_untyped=0 or is_known_action=1)' if /1/;
            push @filter, 'drug.is_withdrawn=0,drug.is_nutraceutical=0,is_potentiator=0,gene.is_kinase=1,(is_untyped=0 or is_known_action=1)' if /2/;
            push @filter, 'drug.is_withdrawn=0,drug.is_nutraceutical=0,is_potentiator=0,drug.is_antineoplastic=1,(is_untyped=0 or is_known_action=1)' if /3/;
        }

        my @sources;
        push @sources, 'TTD' if $params->{'ttd'};
        push @sources, 'DrugBank' if $params->{'db'};
        if(@sources){
            my $filter;
            $filter .= 'source_db_name';
            $filter .= '=' if @sources == 1;
            $filter .= ':' if @sources > 1;#if we have multiple sources, we need to use : with / delimited list for boolean expr syntax
            $filter .= join '/', @sources;
            push @filter, $filter;
        }

        my $filter = join ',', @filter;

        my $command = Genome::DruggableGene::Command::GeneNameGroup::LookupInteractions->execute(
            gene_identifiers => \@gene_names,
            filter => $filter,
        );
        my %params = (
            data => $command->result,
            subject => Genome::DruggableGene::GeneNameReport->define_set,
        );
        for($params->{'return'}){
            if(/html/){
                my $html = Genome::DruggableGene::GeneNameReport::Set::View::Interaction::Html->create(
                    %params,
                    desired_perspective => 'status',
                    xsl_root => Genome->base_dir . '/xsl',
                    xsl_path => '/static/xsl',
                    xsl_variables => {
                        rest      => '/view',
                        resources => '/view/genome/resource.html',
                    },
                );
                return [200, ['Content-type' => "text/html"], [$html->content]];
            } elsif(/tsv/) {
                return [200, ['Content-type' => "text/tsv"], [join("\n", $command->output)]];
            } elsif(/xml/) {
                my $xml = Genome::DruggableGene::GeneNameReport::Set::View::Interaction::Xml->create(%params);
                return [200, ['Content-type' => "text/xml"], [$xml->content]];
            }
        }
    },
#    sub ( POST + /view/x/subject-upload + *file= ) {
    sub ( POST + /view/x/subject-upload + %* + *file= )  {

        load_modules();
        my ($self, $params, $file, $env) = @_;

        my $pathname = $file->path;

        my $c;
        {   open(my $fh, $pathname);
            undef $/;
            $c = <$fh>;
            close($fh);
        }

        my $base64 = MIME::Base64::encode_base64($c);

        my $task_params_json = encode_json( {
            nomenclature => $params->{'nomenclature'},
            subclass_name => $params->{'subclass_name'},
            project_name => $params->{'project_name'},
            content => $base64 });

        my $task_params = {
            command_class => 'Genome::Subject::Command::Import',
            user_id       => Genome::Sys->username(),
            params        => $task_params_json
        };

        my $task; eval {
            $task = Genome::Task->create(%$task_params);
            UR::Context->commit();
        };

        my $code; my $body = {};
        if ($@ || !$task) {
            $code = 200; # OK (didnt work)
            $body->{'error'} = $@ || 'Couldnt create a task with params: ' . Data::Dumper::Dumper $task_params;
            return [$code, ['Content-type' => "text/plain"], [$body->{'error'}]];
        } else {
            $code = 201; # CREATED
            $body->{'id'} = $task->id();
        }

        return [301, [ 'Location' => "/view/genome/task/status.html?id=" . $task->id() ], [encode_json($body)]];
    }
};

Genome::Model::Command::Services::WebApp::File->run_if_script;
