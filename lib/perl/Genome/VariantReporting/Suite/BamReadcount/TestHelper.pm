package Genome::VariantReporting::Suite::BamReadcount::TestHelper;
use Exporter 'import';
@EXPORT_OK = qw(create_default_entry create_deletion_entry create_long_deletion_entry create_no_readcount_entry);

use Genome::File::Vcf::Reader;

my $test_dir = __FILE__.".d";

sub create_deletion_entry {
    my $reader = Genome::File::Vcf::Reader->new(File::Spec->join($test_dir, "deletion.vcf"));
    return $reader->next;
}

sub create_long_deletion_entry {
    my $reader = Genome::File::Vcf::Reader->new(File::Spec->join($test_dir, "long_deletion.vcf"));
    return $reader->next;
}

sub create_default_entry {
    my $reader = Genome::File::Vcf::Reader->new(File::Spec->join($test_dir, "default.vcf"));
    return $reader->next;
}

sub create_no_readcount_entry {
    my $reader = Genome::File::Vcf::Reader->new(File::Spec->join($test_dir, "no_readcount.vcf"));
    return $reader->next;
}
1;

