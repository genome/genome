package Genome::Model::Tools::Music::CreateVisualizations::CreateMutationDiagrams;

use warnings;
use strict;
use Genome;
use IO::File;

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::CreateVisualizations::CreateMutationDiagrams {
    is => 'Command::V2',
    has_input => [
        output_dir => {is => 'Text', doc => 'Output directory path'},
        maf_file => {is => 'Text', doc => 'final maf output file'},
        annotation_format => {is => 'Text', valid_values => ['tgi', 'vep'], doc => 'annotation file format'},
    ],
};
 
sub execute {
    my $self = shift;
    my $annotation_file = $self->create_mutation_diagram_annotation_file;
    my $cmd = Genome::Model::Tools::Graph::MutationDiagram->create(
        output_directory => $self->output_dir,
        annotation_format => $self->annotation_format,
        annotation => $annotation_file,
        annotation_format => 'tgi',
    );
    my $rv = eval{$cmd->execute()};
    if($@){
        my $error = $@;
        $self->error_message('Error running ' . $cmd->command_name . ': ' . $error);
        return;
    }

    return 1;
}

sub create_mutation_diagram_annotation_file {
    my $self = shift;
    #create a file from $self->final_maf using (1-based) columns 33-53
    my $output_dir = $self->output_dir;
    my $annotation_file = join('/', $output_dir, 'final.annotated');
    my $input_fh = IO::File->new($self->maf_file, 'r');
    my $output_fh = IO::File->new($annotation_file, 'w');
    my $header = 1;
    for my $line (<$input_fh>){
        $header-- and next if $header;
        chomp $line;
        my @columns = split("\t", $line);
        $output_fh->print(join("\t", @columns[32..52]), "\n");
    }
    $input_fh->close;
    $output_fh->close;
    return $annotation_file;
}

1;
