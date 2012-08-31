package Genome::InstrumentData::Command::Align::Soap;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Align::Soap {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name => { value => 'soap' },
    ],
    has => [
        version => {is=>'String', default_value=>Genome::Model::Tools::Soap::Base->default_soap_align_version}
    ],
    doc => "align instrument data using the SOAP aligner tool (see http://soap.genomics.org.cn/soapaligner.html)",
};

sub help_synopsis {
    return <<EOS
genome instrument-data align soap

genome instrument-data align soap --params='-p NUM_CORES'

genome instrument-data align soap --params='-l SEED_LENGTH'

genome instrument-data align soap --params='-m MIN_INSERT -x MAX_INSERT'
EOS
}

sub help_detail {
return <<EOS
Launch the SOAP aligner in a standard way and produce results ready for the genome modeling pipeline.

See http://soap.genomics.org.cn/soapaligner.html for more information.
EOS
}


1;

