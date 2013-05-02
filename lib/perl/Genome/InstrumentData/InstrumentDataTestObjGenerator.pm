package Genome::InstrumentData::InstrumentDataTestObjGenerator;

use strict;
use warnings;
use Genome;

our $next_instrument_data_id = "-6666";
our $sample;
our $lib;

sub create_solexa_instrument_data {
    my $bam_path = shift;
    my $id = $next_instrument_data_id;
    $next_instrument_data_id++;
    unless (defined $sample) {
        $sample = Genome::Sample->__define__(
            id => -1234,
            name => 'TEST-000',
        );
    }

    unless (defined $lib) {
        $lib = Genome::Library->__define__(
            id => -1235,
            name => $sample->name.'-testlibs1',
            sample_id => $sample->id,
            library_insert_size => 180,
        );
    }

    my $instrument_data = Genome::InstrumentData::Solexa->__define__(
        id => $id,
        sequencing_platform => 'solexa',
        read_length => 101,
        subset_name => '1-AAAA',
        index_sequence => 'AAAA',
        run_name => 'XXXXXX/1-AAAAA',
        run_type => 'Paired',
        flow_cell_id => 'XXXXX',
        lane => 1,
        library => $lib,
        bam_path => $bam_path,
        clusters => 44554,
        fwd_clusters => 44554,
        rev_clusters => 44554,
        analysis_software_version => 'not_old_illumina',
    );

    return ($instrument_data);

}

1;

