#!/usr/bin/env perl
use strict;
use warnings;
use above 'Genome';
use Test::More tests => 13;
{
#
# find executable paths
#

my $htseq_path = Genome::Sys->sw_path('htseq','0.5.3p9','htseq-count');
ok(-e $htseq_path, "got path $htseq_path for htseq-count 0.5.3p9");
ok($htseq_path =~ /0.5.3p9/, "path contains the version number");

my $tophat_path = Genome::Sys->sw_path("tophat","2.0.7");
ok(-e $tophat_path, "got path $tophat_path for tophat 2.0.7");
ok($tophat_path =~ /2.0.7/, "path contains the version number");

my $samtools_path = Genome::Sys->sw_path("samtools","0.1.18");
ok(-e $samtools_path, "got path $samtools_path for samtools 0.1.18");
ok($samtools_path =~ /0.1.18/, "path contains the version number");

# it really sucks that this was done this way
my $bamcheck_path = Genome::Sys->sw_path("samtools/bamcheck","0.13",'bam-check');
ok(-e $bamcheck_path, "got path $bamcheck_path for bamcheck");
ok($bamcheck_path =~ /0.13/, "path contains the version number");

#
# get the entire mapping
#

my %cufflinks_version_map = Genome::Sys->sw_version_path_map('cufflinks');

my %cufflinks_version_map_expected = (
    '0.7.0' => '/gsc/pkg/bio/cufflinks/cufflinks-0.7.0.Linux_x86_64/cufflinks',
    '0.8.0' => '/gsc/pkg/bio/cufflinks/cufflinks-0.8.0.Linux_x86_64/cufflinks',
    '0.8.2' => '/gsc/pkg/bio/cufflinks/cufflinks-0.8.2.Linux_x86_64/cufflinks',
    '0.8.3' => '/gsc/pkg/bio/cufflinks/cufflinks-0.8.3.Linux_x86_64/cufflinks',
    '0.9.0' => '/gsc/pkg/bio/cufflinks/cufflinks-0.9.0.Linux_x86_64/cufflinks',
    '0.9.1' => '/gsc/pkg/bio/cufflinks/cufflinks-0.9.1.Linux_x86_64/cufflinks',
    '0.9.2' => '/gsc/pkg/bio/cufflinks/cufflinks-0.9.2.Linux_x86_64/cufflinks',
    '0.9.3' => '/gsc/pkg/bio/cufflinks/cufflinks-0.9.3.Linux_x86_64/cufflinks',
    '1.0.0' => '/gsc/pkg/bio/cufflinks/cufflinks-1.0.0.Linux_x86_64/cufflinks',
    '1.0.1' => '/gsc/pkg/bio/cufflinks/cufflinks-1.0.1.Linux_x86_64/cufflinks',
    '1.0.3' => '/gsc/pkg/bio/cufflinks/cufflinks-1.0.3.Linux_x86_64/cufflinks',
    '1.1.0' => '/gsc/pkg/bio/cufflinks/cufflinks-1.1.0.Linux_x86_64/cufflinks',
    '1.2.1' => '/gsc/pkg/bio/cufflinks/cufflinks-1.2.1.Linux_x86_64/cufflinks',
    '1.3.0' => '/usr/bin/cufflinks1.3.0',
    '2.0.2' => '/usr/bin/cufflinks2.0.2'
);


# remove any newer versions to keep the test running after new cufflinks installs
my @expected_keys = keys %cufflinks_version_map_expected;
my %extra = %cufflinks_version_map_expected;
delete @extra{@expected_keys};
for my $key (keys %extra) {
    if ($key gt '2.0.2') {
        delete $cufflinks_version_map{$key};
    }    
}

is_deeply(\%cufflinks_version_map, \%cufflinks_version_map_expected, "map for cufflinks tool matches");

#
# test the ability to get a sorted version list
#

my @versions_actual = grep { $_ le '2.0.2' } Genome::Sys->sw_versions('cufflinks');
my @versions_expected = sort keys %cufflinks_version_map_expected;
is("@versions_actual", "@versions_expected", "version list has expected content and order");
}
{
my $jar_path = Genome::Sys->jar_path('GenomeAnalysisTK','2.4');
ok(-e $jar_path, "got path $jar_path for GenomeAnalysisTK version 2.4");
ok($jar_path =~ /2.4/, "path contains the version number");

my %gatk_version_map = Genome::Sys->jar_version_path_map('GenomeAnalysisTK');

my %gatk_version_map_expected = (
    '2.4' => '/usr/share/java/GenomeAnalysisTK-2.4.jar',
);

# remove any newer versions to keep the test running after new installs
my @expected_keys = keys %gatk_version_map_expected;
my %extra = %gatk_version_map_expected;
delete @extra{@expected_keys};
for my $key (keys %extra) {
    if ($key gt '2.4') {
        delete $gatk_version_map{$key};
    }    
}

is_deeply(\%gatk_version_map, \%gatk_version_map_expected, "map for gatk tool matches");

}

