package Genome::DruggableGene::Command::Citation::Create;

use strict;
use warnings;

use Genome;
use IO::File;

class Genome::DruggableGene::Command::Citation::Create {
    is => 'Genome::Command::Base',
    has_input => [
        source_db_name => { is => 'Text', doc => 'The name of the druggable gene source database' },
        source_db_version => { is => 'Text', doc => 'The version identifier of the druggable gene source database' },
        citation_file => { is => 'PATH', doc => 'The path to a file containing a formatted citation' },
        base_url => {is => 'URL', default => '', doc => 'Partial url used to link back to original source entries'},
    ],
};

sub help_brief {
    return <<EOS
Create a druggable gene citation
EOS
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome druggable-gene citation create --source-db-name DrugBank --source-db-version 3 --citation-file citation.txt
EOS
}

sub help_detail {
    return <<EOS 
This takes a source database name and version along with a formated text file containing the citation for that data source.

This tool will add the citation (and its formatting) to the database and associate it with the database name and version.
EOS
}

sub execute {
    my $self = shift;
    my $source_db_name = $self->source_db_name;
    my $source_db_version = $self->source_db_version;
    my $citation_file = $self->citation_file;
    my $base_url = $self->base_url;

    my $citation_text = $self->_load_citation($citation_file);
    my $citation = Genome::DruggableGene::Citation->create(source_db_name => $source_db_name, source_db_version => $source_db_version, citation => $citation_text, base_url => $base_url);
    return $citation;
}

sub _load_citation {
   my $self = shift; 
   my $citation_file_path = shift;
   my $citation_fh = IO::File->new($citation_file_path, 'r');
   my $citation = '';
   while (my $line = <$citation_fh>){
        chomp $line;
        $citation .= $line;
        $citation .= "\n";
   }
   chomp $citation;
   return $citation;
}

1;
