package Genome::Model::Tools::RegulomeDb;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::RegulomeDb {
    is => 'Command',
};

sub fetch_large_annotation {
    my $format = shift;
    my $data = shift;

    my @data = split(/\n/, $data);
    my $content;
    while (@data) {
        my @part = splice(@data, 0, 1000);
        my $output = fetch_annotation($format, join("\n", @part));
        if ($content) {
            $content = join("\n", $content, $output);
        }
        else {
            $content = $output;
        }
    }
    return $content;
}

sub fetch_annotation {
    my $format = shift;
    my $data = shift;

    my @lines = split(/\n/, $data);
    my $num_lines = grep {!($_ =~ /^#/)} @lines;

    my $mech = WWW::Mechanize->new(stack_depth => 0);

    my $range = 240;
    my $min = 10;

    sleep(rand($range) + $min);

    $mech->get("http://regulome.stanford.edu/");
    die "Couldn't get to page: ", $mech->response->status_line
        unless $mech->success;
    
    sleep(rand($range) + $min);

    $mech->submit_form(
        fields => {data => $data,},
    );

    die "Couldn't submit data: ", $mech->response->status_line
        unless $mech->success;

    sleep(rand($range) + $min);
    
    $mech->submit_form(
        fields => {format => $format},
    );

    die "Failed to download output: ", $mech->response->status_line
        unless $mech->success;

    my $output = $mech->content;
    my @output_lines = split(/\n/, $output);
    my $output_count = grep {!($_ =~ /^#/)} @output_lines;
    unless ($output_count == $num_lines) {
        die "Output did not contain the same number of lines as the input";
    }
    return $output;
}

1;

