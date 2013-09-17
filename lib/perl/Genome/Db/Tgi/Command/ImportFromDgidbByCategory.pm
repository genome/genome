package Genome::Db::Tgi::Command::ImportFromDgidbByCategory;

use strict;
use warnings;
use Genome;

my $API_PATH = "/api/v1/genes_by_category.json";

class Genome::Db::Tgi::Command::ImportFromDgidbByCategory {
    is => 'Genome::Model::Tools::Dgidb::Base',
    has => [
        categories => {
            is => 'Text',
            is_many => 1,
        },
        output_file => {
            is => 'Text',
        },
    ],
};

sub execute {
    my $self = shift;

    my $out = Genome::Sys->open_file_for_writing($self->output_file);

    my %genes;
    for my $category ($self->categories) {
        my %params = (category => $category);
        my $resp = $self->post_request(\%params);
        if ($resp->is_success) {
            $self->append_response(decode_json($resp->content), \%genes);
        }
        else {
            die $self->error_message("Something went wrong.  Does the category $category exist?\n");
        }
    }

    for my $gene (sort keys %genes) {
        $out->print("$gene\n");
    }
    $out->close;
    return 1;
}

sub get_api_path {
    return $API_PATH;
}

sub append_response {
    my $self = shift;
    my $json_ref = shift;
    my $genes_ref = shift;

    for my $gene (@$json_ref) {
        $genes_ref->{$gene}++;
    }
    return 1;
}

1;

