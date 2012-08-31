package Genome::DataSource::Genes;

use Genome;

class Genome::DataSource::Genes {
    is => ['UR::DataSource::FileMux', 'UR::Singleton'],
};

sub delimiter {
    return ",";
}

sub column_order {
    return [qw(
    gene_id 
    hugo_gene_name 
    strand 
    species 
    source 
    version
    )];
}

sub sort_order {
    return ['gene_id', 'source'];
}

sub skip_first_line {
    return 0;
}

sub constant_values { [qw/data_directory reference_build_id/] };
#sub required_for_get { ['gene_id', 'species', 'source','version','data_directory'] }
sub required_for_get { ['id','data_directory','reference_build_id'] }

sub file_resolver {

    my($composite_id, $data_directory) = @_;
    #my($gene_id, $species, $source, $version, $data_directory) = @_;


    my $meta = Genome::Gene->__meta__;

    my ($gene_id) = $meta->resolve_ordered_values_from_composite_id($composite_id);
    
    my $thousand = int($gene_id / 1000);
    my $path = "$data_directory/genes/genes_" . $thousand . ".csv";
    return $path;
}


1;

