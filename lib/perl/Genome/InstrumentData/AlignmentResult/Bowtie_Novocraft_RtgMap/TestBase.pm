package Genome::InstrumentData::AlignmentResult::Bowtie_Novocraft_RtgMap::TestBase;

use Test::More;

# Common code used in Bowtie.t Novocraft.t and RtgMap.t in the parent directory

sub new {
    my $class = shift;
    return bless {}, $class;
}

sub expected_shortcut_path {
    my $self = shift;
    $self->{expected_shortcut_path} ||= do {
        my $aligner_label = $self->aligner_name() . $self->aligner_version();
        $aligner_label =~ s/\./\_/g;
        "/gscmnt/sata828/info/alignment_data/$aligner_label/TEST-human/test_run_name/4_-123456";
    };
}

sub aligner_version {
    my $self = shift;
    $self->{aligner_version} ||= do {
        my $subclass_part = Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($self->aligner_name());
        my $aligner_tools_class_name = "Genome::Model::Tools::${subclass_part}";
        $aligner_tools_class_name->default_version;
    };
}

sub reference_build {
    my $self = shift;
    $self->{reference_build} ||= do {
        my $reference_model = Genome::Model::ImportedReferenceSequence->get(name => 'TEST-human');
        ok($reference_model, "got reference model");

        my $reference_build = $reference_model->build_by_version('1');
        ok($reference_build, "got reference build");
        $reference_build;
    }
}

sub samtools_version {
    my $self = shift;
    $self->{samtools_version} ||= do {
        Genome::Model::Tools::Sam->default_samtools_version
    };
}

sub picard_version {
    my $self = shift;
    $self->{picard_version} ||= do {
        Genome::Model::Tools::Picard->default_picard_version
    };
}

sub initial_fake_instrument_data_id { -123456 }

sub next_fake_instrument_data_id {
    my $self = shift;
    $self->{next_fake_id} = defined($self->{next_fake_id})
                                ? $self->{next_fake_id} - 1
                                : $self->initial_fake_instrument_data_id;
    return $self->{next_fake_id};
}

sub execute {
    my $self = shift;

    my $reference_index = Genome::Model::Build::ReferenceSequence::AlignerIndex->create(
                                aligner_name => $self->aligner_name(),
                                aligner_params => $self->aligner_params(),
                                aligner_version => $self->aligner_version(),
                                reference_build => $self->reference_build());
    ok($reference_index, "generated reference index");

    # Uncomment this to create the dataset necessary for shorcutting to work
    #$self->test_alignment(generate_shortcut_data => 1);

    $self->test_shortcutting();
    $self->test_alignment();
    $self->test_alignment(force_fragment => 1);
}

sub test_alignment {
    my $self = shift;
    my %p = @_;
    
    my $generate_shortcut = delete $p{generate_shortcut_data};

    my $instrument_data = $self->generate_fake_instrument_data();
    my $alignment = Genome::InstrumentData::AlignmentResult->create(
                                                       instrument_data_id => $instrument_data->id,
                                                       samtools_version => $self->samtools_version(),
                                                       picard_version => $self->picard_version(),
                                                       aligner_version => $self->aligner_version(),
                                                       aligner_name => $self->aligner_name(),
                                                       reference_build => $self->reference_build(),
                                                       %p,
                                                   );

    ok($alignment, "Created Alignment");
    my $dir = $alignment->output_dir;
    ok($dir, "alignments found/generated");
    ok(-d $dir, "result is a real directory");
    ok(-s $dir . "/all_sequences.bam", "result has a bam file");

    if ($generate_shortcut) {
        print "*** Using this data to generate shortcut data! ***\n";

        my $expected_shortcut_path = $self->expected_shortcut_path();
        if (-d $expected_shortcut_path) {
            die "Expected shortcut path $expected_shortcut_path already exists, don't want to step on it";
        }
        mkpath($expected_shortcut_path);

        system("rsync -a $dir/* $expected_shortcut_path");
    }

    # clear out the temp scratch/staging paths since these normally would be auto cleaned up at completion
    my $base_tempdir = Genome::Sys->base_temp_directory;
    for (glob($base_tempdir . "/*")) {
        File::Path::rmtree($_);
    }
}

sub test_shortcutting {
    my $self = shift;

    my $fake_instrument_data = $self->generate_fake_instrument_data();

    my $alignment_result_class_name = "Genome::InstrumentData::AlignmentResult::" . Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($self->aligner_name());
    eval "use $alignment_result_class_name";

    my $alignment_result = $alignment_result_class_name->__define__(
                 id => -8765432,
                 output_dir => $self->expected_shortcut_path(),
                 instrument_data_id => $fake_instrument_data->id,
                 subclass_name => $alignment_result_class_name,
                 module_version => '12345',
                 aligner_name => $self->aligner_name(),
                 aligner_version => $self->aligner_version(),
                 samtools_version => $self->samtools_version(),
                 picard_version => $self->picard_version(),
                 reference_build => $self->reference_build(),
    );
    $alignment_result->lookup_hash($alignment_result->calculate_lookup_hash);

    # Alignment Result is a subclass of Software Result. Make sure this is true here.
    isa_ok($alignment_result, 'Genome::SoftwareResult');


    #
    # Step 1: Attempt to create an alignment that's already been created 
    # ( the one we defined up at the top of the test case )
    #
    # This ought to fail to return anything, and set the error_message property to include
    # some info about why we failed.  
    ####################################################

    Genome::InstrumentData::AlignmentResult->dump_error_messages(0);  # supress error about duplicate object
    my $bad_alignment = Genome::InstrumentData::AlignmentResult->create(
                                                              instrument_data_id => $fake_instrument_data->id,
                                                              aligner_name => $self->aligner_name(),
                                                              aligner_version => $self->aligner_version(),
                                                              samtools_version => $self->samtools_version(),
                                                              picard_version => $self->picard_version(),
                                                              reference_build => $self->reference_build(),
                                                          );
    ok(!$bad_alignment, "this should have returned undef, for attempting to create an alignment that is already created!");
    Genome::InstrumentData::AlignmentResult->dump_error_messages(1);


    #
    # Step 2: Attempt to get an alignment that's already created
    #
    #################################################
    my $alignment = Genome::InstrumentData::AlignmentResult->get_with_lock(
                                                              instrument_data_id => $fake_instrument_data->id,
                                                              aligner_name => $self->aligner_name(),
                                                              aligner_version => $self->aligner_version(),
                                                              samtools_version => $self->samtools_version(),
                                                              picard_version => $self->picard_version(),
                                                              reference_build => $self->reference_build(),
                                                              );
    ok($alignment, "got an alignment object");


    # once to find old data
    my $adir = $alignment->output_dir;
    my @list = <$adir/*>;

    ok($alignment, "Created Alignment");
    my $dir = $alignment->output_dir;
    ok($dir, "alignments found/generated");
    ok(-d $dir, "result is a real directory");
    ok(-s $dir."/all_sequences.bam", "found a bam file in there");
}

sub generate_fake_instrument_data {
    my $self = shift;

    my $fastq_directory = $ENV{GENOME_TEST_INPUTS} . '/Genome-InstrumentData-Align-Maq/test_sample_name';
    my $fake_id = $self->next_fake_instrument_data_id();
    my $instrument_data = Genome::InstrumentData::Solexa->create_mock(
                                                                      id => $fake_id,
                                                                      sequencing_platform => 'solexa',
                                                                      flow_cell_id => '12345',
                                                                      lane => '1',
                                                                      seq_id => $fake_id,
                                                                      median_insert_size => '22',
                                                                      sample_name => 'test_sample_name',
                                                                      library_id => $self->library_id(),
                                                                      run_name => 'test_run_name',
                                                                      subset_name => 4,
                                                                      run_type => 'Paired End Read 2',
                                                                      gerald_directory => $fastq_directory,
                                                                  );


    # confirm there are fastq files here, and fake the fastq_filenames method to return them
    my @in_fastq_files = glob($instrument_data->gerald_directory.'/*.txt');
    $instrument_data->set_list('dump_sanger_fastq_files',@in_fastq_files);
    $instrument_data->mock('dump_trimmed_fastq_files', sub {return Genome::InstrumentData::Solexa::dump_trimmed_fastq_files($instrument_data)});

    # fake out some properties on the instrument data
    isa_ok($instrument_data,'Genome::InstrumentData::Solexa');
    $instrument_data->set_always('sample_type','dna');
    $instrument_data->set_always('sample_id','2791246676');
    $instrument_data->set_always('is_paired_end',1);
    ok($instrument_data->is_paired_end,'instrument data is paired end');
    $instrument_data->set_always('calculate_alignment_estimated_kb_usage',10000);
    $instrument_data->set_always('resolve_quality_converter','sol2sanger');
    $instrument_data->set_always('run_start_date_formatted','Fri Jul 10 00:00:00 CDT 2009');

    return $instrument_data;
}

1;



