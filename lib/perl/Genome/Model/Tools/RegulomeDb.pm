package Genome::Model::Tools::RegulomeDb;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::RegulomeDb {
    is => 'Command',
};

sub fetch_large_annotation {
    my $self = shift;
    my $format = shift;
    my $data = shift;

    my @data = split(/\n/, $data);
    my $num_lines = scalar @data;
    $self->debug_message("Fetching annotations for $num_lines");
    my $content;
    while (@data) {
        my @part = splice(@data, 0, 1000);
        my $output = $self->fetch_annotation($format, join("\n", @part));
        if ($content) {
            my @output_lines = split(/\n/, $output);
            #remove header from subsequent lines
            shift @output_lines;
            $content = join("\n", $content, @output_lines);
        }
        else {
            $content = $output;
        }
    }
    return $content;
}

sub fetch_annotation {
    my $self = shift;
    my $format = shift;
    my $data = shift;

    my @lines = split(/\n/, $data);
    my $num_lines = grep {!($_ =~ /^#/)} @lines;

    my $mech = WWW::Mechanize->new(stack_depth => 0,
                                   timeout => 600);

    my $range = 240;
    my $min = 10;
    my $max_retries = 10;
    my $retry_count = 0;

    while ($retry_count < $max_retries and !$mech->success) {
        eval{
            $mech->get("http://regulome.stanford.edu/");
            
            $mech->submit_form(
                fields => {data => $data,},
            );    

            $mech->submit_form(
                fields => {format => $format},
            );
        };
        if ($mech->success) {
            last;
        }
        else {
            $retry_count++;
            my $sleep = rand($range) + $min;
            $self->debug_message("Sleeping $sleep seconds and retrying");
            sleep($sleep);
        }
    }

    unless ($mech->success) {
        $self->error_message("Couldn't get regulomedb info after $max_retries tries: ". $mech->response->status_line);
        die $self->error_message;
    }
    
    my $output = $mech->content;
    my @output_lines = split(/\n/, $output);
    my $output_count = grep {!($_ =~ /^#/)} @output_lines;
    #unless ($output_count == $num_lines) {
    #    $self->error_message( "Output did not contain the same number of lines as the input.  Output had $output_count and input had $num_lines.  Output: \n$output");
    #    die $self->error_message;
    #}
    $self->debug_message("Successfully got $output_count lines");
    return $output;
}

1;

