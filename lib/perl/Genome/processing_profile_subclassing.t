use above 'Genome';

use Test::More;

UR::DBI->no_commit(1);

plan tests => 31;

# Some Processing Profiles are abstract and require further info to subclass properly
# Genome::ProcessingProfile, G::PP::ReferenceAlignment and G::PP::RnaSeq

my $pp = eval { Genome::ProcessingProfile->create(name => 'will fail') };
ok(! $pp, 'Could not create an un-subclassed Processing Profile');
like($@, qr/subclassify_by calculation property 'subclass_name' requires 'type_name'/, 'Exception looks correct');

$pp = eval { Genome::ProcessingProfile::ReferenceAlignment->create(name => 'will fail') };
ok(! $pp, 'Could not create an un-subclassed ReferenceAlignment PP');
# There's no exception here because G::PP::create() does a check for G::PP::Param required values.
# since sequencing_platform is missing, it does an error_message() and returns nothing

$pp = eval { Genome::ProcessingProfile->create(name => 'will fail', type_name => 'reference alignment') };
ok(! $pp, 'Could not create an un-subclassed Processing Profile with only type_name => reference alignment');

$pp = eval { Genome::ProcessingProfile->create(name => 'will fail', type_name => 'non existent') };
ok(! $pp, 'Could not create a ProcessingProfile with type_name refering to a bad subclass');
like($@, qr/Genome::ProcessingProfile::NonExistent is not a subclass of Genome::ProcessingProfile/, 'Exception looks correct');

my @pps = Genome::ProcessingProfile->is_loaded();
ok(! scalar(@pps), 'After failed ProcessingProfile create()s, no PPs exist');

$pp = Genome::ProcessingProfile->create(name => 'test alignment 1', type_name => 'simple alignment', reference_sequence_name => 'foo');
ok($pp, 'Created a Processing Profile with type_name => simple alignment');
isa_ok($pp, 'Genome::ProcessingProfile::SimpleAlignment');
isa_ok($pp, 'Genome::ProcessingProfile');

$pp = Genome::ProcessingProfile::SimpleAlignment->create(name => 'test alignment 2', reference_sequence_name => 'foo2');
ok($pp, 'Created a Genome::ProcessingProfile::SimpleAlignment');
isa_ok($pp, 'Genome::ProcessingProfile::SimpleAlignment');
isa_ok($pp, 'Genome::ProcessingProfile');



$pp = Genome::ProcessingProfile::ReferenceAlignment->create(name => 'test refalign 1',
                                                            sequencing_platform => '454',
                                                            dna_type => 'genomic dna',
                                                            read_aligner_name => 'foo');
ok($pp, 'Created an ReferenceAlignment PP');
isa_ok($pp, 'Genome::ProcessingProfile::ReferenceAlignment::454');
isa_ok($pp, 'Genome::ProcessingProfile::ReferenceAlignment');
isa_ok($pp, 'Genome::ProcessingProfile');

$pp = Genome::ProcessingProfile->create(name => 'test refalign 2',
                                        type_name => 'reference alignment',
                                        sequencing_platform => '454',
                                        dna_type => 'genomic dna',
                                        read_aligner_name => 'foo2');
ok($pp, 'Created a Processing Profile with type_name => reference alignment');
isa_ok($pp, 'Genome::ProcessingProfile::ReferenceAlignment::454');
isa_ok($pp, 'Genome::ProcessingProfile::ReferenceAlignment');
isa_ok($pp, 'Genome::ProcessingProfile');



$pp = Genome::ProcessingProfile::RnaSeq->create(name => 'test rnaseq 1',
                                                sequencing_platform => '454',
                                                dna_type => 'cdna',
                                                read_aligner_name => 'foo');
ok($pp, 'Created an RnaSeq PP');
# Removed subclassing by sequencing platform
#isa_ok($pp, 'Genome::ProcessingProfile::RnaSeq::454');
isa_ok($pp, 'Genome::ProcessingProfile::RnaSeq');
isa_ok($pp, 'Genome::ProcessingProfile');

$pp = Genome::ProcessingProfile->create(name => 'test rnaseq 2',
                                        type_name => 'rna seq',
                                        sequencing_platform => '454',
                                        dna_type => 'cdna',
                                        read_aligner_name => 'foo2');
ok($pp, 'Created an RnaSeq PP');
# Removed subclassing by sequencing platform
#isa_ok($pp, 'Genome::ProcessingProfile::RnaSeq::454');
isa_ok($pp, 'Genome::ProcessingProfile::RnaSeq');
isa_ok($pp, 'Genome::ProcessingProfile');

SKIP:{
skip('G::PP::create() tries to check params before UR::Context::create_entity has a chance to fill it in',
     4);
$pp = Genome::ProcessingProfile::RnaSeq::454->create(name => 'test rnaseq 3',
                                                     dna_type => 'cdna',
                                                     read_aligner_name => 'foo3');
ok($pp, 'Created an RnaSeq PP');
# Removed subclassing by sequencing platform
#isa_ok($pp, 'Genome::ProcessingProfile::RnaSeq::454');
isa_ok($pp, 'Genome::ProcessingProfile::RnaSeq');
isa_ok($pp, 'Genome::ProcessingProfile');
}



