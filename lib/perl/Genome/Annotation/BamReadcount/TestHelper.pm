package Genome::Annotation::BamReadcount::TestHelper;
use Exporter 'import';
@EXPORT_OK = qw(bam_readcount_line create_entry);

sub create_vcf_header {
    my $header_txt = <<EOS;
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="Passed all filters">
##FILTER=<ID=BAD,Description="This entry is bad and it should feel bad">
##INFO=<ID=A,Number=1,Type=String,Description="Info field A">
##INFO=<ID=C,Number=A,Type=String,Description="Info field C (per-alt)">
##INFO=<ID=E,Number=0,Type=Flag,Description="Info field E">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=BRCT,Number=1,Type=String,Description="Bam readcount entry">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
EOS
    my @lines = split("\n", $header_txt);
    my $header = Genome::File::Vcf::Header->create(lines => \@lines);
    return $header
}

sub create_entry {
    my $brct_value = shift;
    my $brct_string = create_bam_readcount_string($brct_value);

    my @fields = (
        '1',            # CHROM
        10,             # POS
        '.',            # ID
        'A',            # REF
        '.',            # ALT
        '10.3',         # QUAL
        'PASS',         # FILTER
        'A=B;C=8,9;E',  # INFO
        'GT:DP:BRCT',     # FORMAT
        "0/1:12:$brct_string",   # FIRST_SAMPLE
    );

    my $entry_txt = join("\t", @fields);
    my $entry = Genome::File::Vcf::Entry->new(create_vcf_header(), $entry_txt);
    return $entry;
}

sub bam_readcount_line {
    return "21	10402985	G	344	Solexa-135852	{	=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	A:2:57.50:16.00:0.00:2:0:0.35:0.02:42.50:2:0.51:231.50:0.51	C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	G:155:53.61:25.61:0.00:90:65:0.51:0.01:11.10:134:0.39:236.59:0.38	T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	+A:20:52.55:0.00:0.00:12:8:0.53:0.02:25.00:18:0.41:237.70:0.39	}	Solexa-135853	{	=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	A:1:60.00:2.00:0.00:0:1:0.76:0.10:46.00:0:-nan:250.00:0.62	C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	G:186:52.12:26.33:0.29:99:87:0.50:0.01:12.29:161:0.39:238.74:0.40	T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	}";
}

sub create_bam_readcount_string {
    my $readcount_line = shift;
    return "" unless ($readcount_line);
    return Genome::File::BamReadcount::Entry::encode($readcount_line);
}
1;

