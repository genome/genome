package Genome::VariantReporting::Framework::Test::RunResult;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Writer;

class Genome::VariantReporting::Framework::Test::RunResult {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Result',
    has_input => [
        __input___lookup => {
            is_many => 1,
        },
    ],
    has_param => [
        __planned__ => {},
    ],
    has_transient_optional => [
        __input__ => {
            is_many => 1,
        },
    ],
};


sub output_filename {
    return '__test__.vcf.gz';
}

my $info_type = {
        id => 'TST',
        number => '1',
        type => 'String',
        description => 'Some test info showing things got passed in as expected',
};

sub _run {
    my $self = shift;

    my $reader = Genome::File::Vcf::Reader->new($self->input_vcf);
    $reader->header->info_types->{$info_type->{id}} = $info_type;

    my $output_file = File::Spec->join($self->temp_staging_directory, $self->output_filename);
    my $writer = Genome::File::Vcf::Writer->new($output_file, $reader->header);

    my $some_info = sprintf("__input___lookups|%s|__planned__|%s",
        join('|', $self->__input___lookup),
        $self->__planned__,
    );
    while (my $entry = $reader->next) {
        $entry->set_info_field('TST', $some_info);
        $writer->write($entry);
    }
    $writer->close();

    return 1;
}

sub add_to_header {
    my ($header, $key, $values) = @_;

    $header->metainfo->{$key} = Genome::File::Vcf::Header::String->new(
        content => join(', ', @$values),
    );
}


1;
