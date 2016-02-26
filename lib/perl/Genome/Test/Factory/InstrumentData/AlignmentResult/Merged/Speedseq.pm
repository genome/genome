package Genome::Test::Factory::InstrumentData::AlignmentResult::Merged::Speedseq;
use base qw(Genome::Test::Factory::Base);

use Genome;

sub generate_obj {
    my $self = shift;
    my %params = @_;

    my $idata = delete $params{instrument_data};
    my @instrument_data = $idata? @$idata : ();

    my $per_lane_params = delete $params{__per_lane_params} // {};

    my $id = delete $params{id};

    my $speedseq_result = Genome::InstrumentData::AlignmentResult::Merged::Speedseq->__define__(
        %params,
        id => $id,
    );
    for my $i (0..$#instrument_data) {
        $speedseq_result->add_input(
            name => 'instrument_data-' . $i,
            value_id => $instrument_data[$i]->id,
        );
    }
    $speedseq_result->add_param(
        name => 'instrument_data_count',
        value_id=> scalar(@instrument_data),
    );
    $speedseq_result->add_param(
        name => 'instrument_data_md5',
        value_id => Genome::Sys->md5sum_data(join(':', sort(map($_->id, @instrument_data))))
    );
    $speedseq_result->lookup_hash($speedseq_result->calculate_lookup_hash());

    for my $instrument_data (@instrument_data) {
        my $per_lane_speedseq_result = Genome::InstrumentData::AlignmentResult::Speedseq->__define__(
            instrument_data => $instrument_data,
            %params,
            %$per_lane_params,
        );
        $per_lane_speedseq_result->lookup_hash($per_lane_speedseq_result->calculate_lookup_hash());
    }

    return $speedseq_result;
}

1;
