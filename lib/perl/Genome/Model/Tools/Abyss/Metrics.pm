package Genome::Model::Tools::Abyss::Stats;

use Genome;

class Genome::Model::Tools::Abyss::Stats {
    is => 'Genome::Model::Tools::Assembly::Stats',
    has => [
	    assembly_directory => {
            is => 'Text',
            doc => 'Path to soap assembly',
        },
        first_tier => {
            is => 'Number',
            doc => 'First tier value',
            is_optional => 1,
        },
        second_tier => {
            is => 'Number',
            doc => 'Second tier value',
            is_optional => 1,
        },
        major_contig_length => {
            is => 'Number',
            is_optional => 1,
            default_value => 500,
            doc => 'Cutoff value for major contig length',
        },
        output_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'Stats output file',
        },
    ],
    has_optional => [
        _metrics => { is_transient => 1, },
    ],
};

sub execute {
    # FIXME
    return 1;
}

1;
