package Genome::Db::DbVar::Import;

use strict;
use warnings;
use Genome;
use LWP::Simple;

class Genome::Db::DbVar::Import {
    is => "Command::V2",
    has_input => [
        filename => {
            is => 'Text',
            doc => 'file to save the db to',
        },
        species => {
            is => 'Text',
            valid_values => [qw/Bos_taurus Canis_lupus Danio_rerio Drosophila_melanogaster Equus_caballus Homo_sapiens Macaca_mulatta Mus_musculus Pan_troglodytes Sorghum_bicolor Sus_scrofa/],
        },
        assembly_name => {
            is => 'Text',
            doc => 'Name of reference assembly',
        },
    ],
};

sub execute {
    my $self = shift;

    my $url = join("/", "ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data", $self->species, "by_assembly", $self->assembly_name, "gvf", $self->assembly_name.".remap.all.germline.ucsc.gvf.gz");

    my $path = $self->filename.".gz";

    my $response = getstore($url, $path);

    die($self->error_message("Unable to download gvf file at $url")) unless $response == RC_OK;

    unless (-s $path) {
        $self->warning_message("File $path is empty, continuing");
        return 1;
    }

    `gunzip $path`;

    return 1;
}

1;

