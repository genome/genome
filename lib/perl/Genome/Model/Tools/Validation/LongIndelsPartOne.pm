package Genome::Model::Tools::Validation::LongIndelsPartOne;

use warnings;
use strict;
use Genome;
use IO::File;

class Genome::Model::Tools::Validation::LongIndelsPartOne {
    is => 'Command',
    has_input => [
        long_indel_bed_file => {
            is => 'String',
            doc => 'unsorted, unannotated 3bp indel file in BED format!!! BED format!!!',
        },
        output_dir => {
            is => 'String',
            doc => 'directory for output files',
        },
        sample_identifier => {
            is => 'String',
            doc => 'some string to use in model names, etc, such as "BRC2"',
        },
        reference_transcripts => {
            is => 'String',
            doc => 'reference transcripts plus version to be used to annotate input indel file',
        },
    ],
    has_optional_input => [
        tumor_val_model_id => {
            is => 'Number',
            doc => 'refalign model ID for the tumor sample',
        },
        normal_val_model_id => {
            is => 'Number',
            doc => 'refalign model ID for the normal sample',
        },
        somatic_validation_model_id => {
            is => 'Number',
            doc => 'somatic-validation model ID (contains both tumor and normal)',
        },
        somatic_validation_build_id => {
            is => 'Number',
            doc => 'somatic-validation build ID (better than using the "model_id" version of this option)'
        },
        reference_sequence_build_id => {
            is => 'Text',
            doc => 'Optional reference sequence path (default is to grab it from the input models)',
        },
    ],
    has => [
        _new_tumor_model => {
            is => 'Genome::Model',
            is_optional => 1,
            doc => 'The new tumor model created by execute',
        },
        _new_normal_model => {
            is => 'Genome::Model',
            is_optional => 1,
            doc => 'The new normal model created by execute',
        },
    ],
    doc => 'Begin validation of 3bp indels.',
};

sub help_detail {
    return <<EOS
    This tool performs the first 5 steps of the 3bp indel validation process outlined on this wiki page: https://gscweb.gsc.wustl.edu/wiki/Medical_Genomics/Nimblegen_Solid_Phase_Capture_Validation/Analysis#.3E3bp_Indels. It then prints out a command which the user may use to run two builds for remapping validation data to a new reference sequence containing indel contigs, and also prints a follow-up command to run to complete the final steps of the 3bp indel process. It also prints out some details for any future manual review tickets of these indels, so save the STDOUT. Remember to address the optional parameters if you are not using HG18 build36.

Requires either paired tumor/normal models or a single somatic-validation model
EOS
}

sub execute {
    my $self = shift;

    my $output_dir = $self->output_dir;
    my $sample_id = $self->sample_identifier;
    my $input_model_type;
    my $normal_bam;
    my $tumor_bam;
    my $ref_seq_fasta;
    my $ref_seq_build;
    my $ref_seq_build_id;
    my $model;
    my $tumor_model;
    my $normal_model;
    my $tumor_build;
    my $normal_build;
    my $tumor_sample;
    my $normal_sample;
    my $tumor_sample_id;
    my $normal_sample_id;
    my $tumor_subject;
    my $normal_subject;
    my @tumor_instrument_data;
    my @normal_instrument_data;
    my @tumor_instrument_data_ids;
    my @normal_instrument_data_ids;

    #make sure we have either a somatic-validation build or a pair of tumor/normal refalign builds
    if(!(defined($self->somatic_validation_model_id || $self->somatic_validation_build_id))){
        if(!(defined($self->tumor_val_model_id)) || !(defined($self->normal_val_model_id))){
            die("ERROR: must provide either\n --somatic-validation-build-id  OR\n --tumor-val-model-id and --normal-val-model-id");
        } else {
            $input_model_type = "pairedref";
        }
    } else {
        $input_model_type = "somval";
    }

    #--------------------------------------------------------------------
    #secure BAM paths and reference sequence ids from input params

    #from somatic-validation model
    if($input_model_type eq "somval"){        
        my $build;
        if($self->somatic_validation_build_id) {
            $build = Genome::Model::Build->get($self->somatic_validation_build_id);
            die "Could not find build" unless $build;
        } else {
            my $model = Genome::Model->get($self->somatic_validation_model_id) or
            die "Could not find model ($self->somatic_validation_model_id\n";

            $build = $model->last_succeeded_build or
                die "Could not find last succeeded build from somatic model $self->somatic_validation_model_id.\n";
        }

        #get refseq build id from model unless already defined
        $ref_seq_build_id = $build->reference_sequence_build->build_id;

        if(defined($self->reference_sequence_build_id)){
            $ref_seq_build = Genome::Model::Build->get($self->reference_sequence_build_id);
        } else {
            $ref_seq_build = Genome::Model::Build->get($ref_seq_build_id);
        }
        $ref_seq_fasta = $ref_seq_build->full_consensus_path('fa');

        $tumor_sample = $build->tumor_sample;
        $tumor_sample_id = $build->tumor_sample->sample_id;
        $tumor_subject = Genome::Subject->get($tumor_sample->subject_id);
        $tumor_bam = $build->tumor_bam;

        if(defined($self->normal_val_model_id)){ #separate refalign model with normal bam
            $normal_model = Genome::Model->get($self->normal_val_model_id) or
                die "Could not find normal model with id $self->normal_val_model_id.\n";
            $normal_build = $normal_model->last_succeeded_build or
                die "Could not find last succeeded build from normal model $self->normal_val_model_id.\n";
            $normal_sample_id="";
            $normal_bam = $normal_build->whole_rmdup_bam_file or die "Cannot find normal .bam.\n";
            
        } else { #normal is in somval model
            $normal_sample = $build->normal_sample;
            $normal_sample_id = $build->normal_sample->sample_id;
            $normal_subject = Genome::Subject->get($normal_sample->subject_id);
            $normal_bam = $build->normal_bam;
        }

        #also need to grab a bunch of other stuff so that we can define refalign models later

        #instrument data ids, sample ids
        my @instrument_data_ids = $build->instrument_data_ids;
        foreach my $instid (@instrument_data_ids){
            my $inst_data = Genome::InstrumentData->get($instid);
            my $library_id = $inst_data->library_id;
            my $library = Genome::Library->get($library_id);
            my $sample_id = $library->sample_id;

            if($sample_id eq $tumor_sample_id){
                push(@tumor_instrument_data,$inst_data);
                push(@tumor_instrument_data_ids,$instid);
            } elsif ($sample_id eq $normal_sample_id){
                push(@normal_instrument_data,$inst_data);
                push(@normal_instrument_data_ids,$instid);
            } else {
                die "sample id $sample_id from instrument data does not match sample id from either tumor (" . $build->tumor_sample->sample_id . ") or normal (" . $build->normal_sample->sample_id . ")\n";
            }
        }



    #-----------
    } else { #from paired refseq builds

        $normal_model = Genome::Model->get($self->normal_val_model_id) or
            die "Could not find normal model with id $self->normal_val_model_id.\n";

        $tumor_model = Genome::Model->get($self->tumor_val_model_id) or
            die "Could not find tumor model with id $self->tumor_val_model_id.\n";

        $normal_build = $normal_model->last_succeeded_build or
            die "Could not find last succeeded build from normal model $self->normal_val_model_id.\n";

        $tumor_build = $tumor_model->last_succeeded_build or
            die "Could not find last succeeded build from tumor model $self->tumor_val_model_id.\n";

        #get refseq build id from model unless already defined
        $ref_seq_build_id = $tumor_model->reference_sequence_build->build_id;

        my $ref_seq_build;
        if(defined($self->reference_sequence_build_id)){
            $ref_seq_build = Genome::Model::Build->get($self->reference_sequence_build_id);
        } else {
            $ref_seq_build = Genome::Model::Build->get($ref_seq_build_id);
        }
        $ref_seq_fasta = $ref_seq_build->full_consensus_path('fa');

        $normal_bam = $normal_build->whole_rmdup_bam_file or die "Cannot find normal .bam.\n";
        $tumor_bam = $tumor_build->whole_rmdup_bam_file or die "Cannot find tumor .bam.\n";
    }



    my $rv = Genome::Model::Tools::Validation::LongIndelsGenerateMergedAssemblies->execute(
        long_indel_bed_file => $self->long_indel_bed_file,
        output_dir => $output_dir,
        reference_transcripts => $self->reference_transcripts,
        tumor_bam => $tumor_bam,
        normal_bam => $normal_bam,
        reference_fasta => $ref_seq_fasta,
    );

    unless ($rv) {
        $self->error_message("Failed to generate contigs_fasta file");
        return;
    }

    my $contigs_file = $rv->contigs_fasta;

    print STDERR "##### $contigs_file #####\n";

    #--------------------------------------------------------------------    
    #create reference sequence using the new contigs (define new reference and track new reference build)

    #first check for duplicate model names and if they exist, suffix the name to avoid conflicts
    my $new_ref_model_name = $sample_id . "-human";
    $new_ref_model_name = checkForDupName($new_ref_model_name);
   
    print STDERR "creating model $new_ref_model_name\n";

    my $new_ref_cmd = Genome::Model::Command::Define::ImportedReferenceSequence->create(
        species_name => 'human',
        use_default_sequence_uri => '1',
        derived_from => $ref_seq_build,
        append_to => $ref_seq_build,
        version => '500bp_assembled_contigs',
        fasta_file => $contigs_file,
        prefix => $sample_id,
        model_name => $new_ref_model_name,
    );
    unless ($new_ref_cmd->execute) {
        $self->error_message('Failed to execute the definition of the new reference sequence with added contigs.');
        return;
    }

    my $new_ref_build_id = $new_ref_cmd->result_build_id;
    my $new_ref_build = Genome::Model::Build->get($new_ref_build_id);
    my $new_ref_event = $new_ref_build->the_master_event;
    my $new_ref_event_id = $new_ref_event->id;
    my $new_ref_event_class = $new_ref_event->class;
    while ($new_ref_event->event_status eq 'Running' || $new_ref_event->event_status eq 'Scheduled') {
        sleep 120;
        $new_ref_event = $new_ref_event_class->load($new_ref_event_id);
    }
    unless ($new_ref_event->event_status eq 'Succeeded') {
        $self->error_message('New reference build not successful.');
        return;
    }


    #-----------------------------
    #create models to align data to new reference

    my ($new_tumor_model, $new_normal_model);
    if($input_model_type eq "somval"){
        #my $new_pp = "dlarson bwa0.5.9 -q 5 indel contig test picard1.42";
        my $new_pp = Genome::ProcessingProfile->get("2599983");
        my $new_pp_id = "2599983";

        my $new_tumor_model_name = $sample_id . "-Tumor-3bpIndel-Validation";
        my $new_normal_model_name = $sample_id . "-Normal-3bpIndel-Validation";
        $new_tumor_model_name = checkForDupName($new_tumor_model_name);
        $new_normal_model_name = checkForDupName($new_normal_model_name);


        #new tumor model
        # FIXME check for an existing model
        my $tumor_copy = Genome::Model::Command::Define::ReferenceAlignment->create(
            reference_sequence_build => $new_ref_build,
            auto_build_alignments => 0,
            auto_assign_inst_data => 0,
            instrument_data => \@tumor_instrument_data,
            subject => $tumor_subject,
            model_name => $new_tumor_model_name,
            processing_profile => $new_pp
            );
        $tumor_copy->dump_status_messages(1);
        $tumor_copy->execute or die "tumor define failed";
        my $new_tumor_model_id = $tumor_copy->result_model_id;
        $new_tumor_model = Genome::Model->get($new_tumor_model_id);
        print STDERR "tumor model defined: $new_tumor_model_id\n";


        #new normal model
        my $new_normal_model_id;
        if(defined($self->normal_val_model_id)){ #separate refalign model with normal bam
            #new normal model
        # FIXME check for an existing model
            my $normal_copy = Genome::Model::Command::Copy->create(
                model => $normal_model,
                overrides => [
                    'name='.$new_normal_model_name,
                    'auto_build_alignments=0',
                    'processing_profile='.$new_pp_id,
                    'reference_sequence_build='.$new_ref_build_id,
                    'annotation_reference_build=',
                    'region_of_interest_set_name=',
                    'dbsnp_build=',
                ],
                );
            $normal_copy->dump_status_messages(1);
            $normal_copy->execute or die "copy failed";
            $new_normal_model = $normal_copy->_new_model;
            $new_normal_model_id = $new_normal_model->id;

        } else {
        # FIXME check for an existing model
            my $normal_copy = Genome::Model::Command::Define::ReferenceAlignment->create(
                reference_sequence_build => $new_ref_build,
                auto_build_alignments => 0,
                auto_assign_inst_data => 0,
                instrument_data => \@normal_instrument_data,
                subject => $normal_subject,
                model_name => $new_normal_model_name,
                processing_profile => $new_pp
                );
            $normal_copy->dump_status_messages(1);
            $normal_copy->execute or die "copy failed";
            $new_normal_model_id = $normal_copy->result_model_id;
            $new_normal_model = Genome::Model->get($new_normal_model_id);
        }
        print STDERR "normal model defined: $new_normal_model_id\n";


        UR::Context->commit;

        #final notices to user
        print "\n\n\n###################### IMPORTANT INFORMATION ######################\n\n";

        #alert user to run builds for these copied models
        print "To start alignments against the new reference sequence which contains the indel contigs, please run this command from genome stable:\n\ngenome model build start $new_tumor_model_id $new_normal_model_id\n\n";

        #alert user to run LongIndelsPartTwo upon the completion of the new tumor and normal contig builds
        print "Upon the successful completion of these builds, please bsub the following command (saving STDOUT) to finish the rest of the steps for 3bp indel validation:\n\n";
        print "gmt validation long-indels-part-two --normal-val-model-copy-id $new_normal_model_id --tumor-val-model-copy-id $new_tumor_model_id --output-dir $output_dir --tier-file-location <PUT YOUR TIERING FILES HERE>\n\n";

        #print details for a manual review ticket
        my $new_ref_build_fa = $new_ref_build->full_consensus_path('fa');
        print "And lastly, for your manual review tickets, you will want to include these details, along with further info printed in the STDOUT of 'gmt validation long-indels-part-two':\n\n";
        print "Original tumor validation BAM: $tumor_bam\n";
        print "Original normal validation BAM: $normal_bam\n";
        print "New reference sequence with contigs: $new_ref_build_fa\n";



    #-----------------------------
    } else {  ##not somval
        #my $new_pp = "dlarson bwa0.5.9 -q 5 indel contig test picard1.42";
        my $new_pp = Genome::ProcessingProfile->get("2599983");
        my $new_pp_id = "2599983";

        my $new_tumor_model_name = $sample_id . "-Tumor-3bpIndel-Validation";
        my $new_normal_model_name = $sample_id . "-Normal-3bpIndel-Validation";
        $new_tumor_model_name = checkForDupName($new_tumor_model_name);
        $new_normal_model_name = checkForDupName($new_normal_model_name);

        #new tumor model
        # FIXME check for an existing model
        my $tumor_copy = Genome::Model::Command::Copy->create(
            model => $tumor_model,
            overrides => [
                'name='.$new_tumor_model_name,
                'auto_build_alignments=0',
                'processing_profile='.$new_pp_id,
                'reference_sequence_build='.$new_ref_build_id,
                'annotation_reference_build=',
                'region_of_interest_set_name=',
                'dbsnp_build=',
            ],
            );
        $tumor_copy->dump_status_messages(1);
        $tumor_copy->execute or die "copy failed";
        my $new_tumor_model = $tumor_copy->_new_model;
        my $new_tumor_model_id = $new_tumor_model->id;

        #new normal model
        # FIXME check for an existing model
        my $normal_copy = Genome::Model::Command::Copy->create(
            model => $normal_model,
            overrides => [
                'name='.$new_normal_model_name,
                'auto_build_alignments=0',
                'processing_profile='.$new_pp_id,
                'reference_sequence_build='.$new_ref_build_id,
                'annotation_reference_build=',
                'region_of_interest_set_name=',
                'dbsnp_build=',
            ],
            );
        $normal_copy->dump_status_messages(1);
        $normal_copy->execute or die "copy failed";
        my $new_normal_model = $normal_copy->_new_model;
        my $new_normal_model_id = $new_normal_model->id;

        #final notices to user
        print "\n\n\n###################### IMPORTANT INFORMATION ######################\n\n";

        #alert user to run builds for these copied models
        print "To start alignments against the new reference sequence which contains the indel contigs, please run this command from genome stable:\n\ngenome model build start $new_tumor_model_id $new_normal_model_id\n\n";

        #alert user to run LongIndelsPartTwo upon the completion of the new tumor and normal contig builds
        print "Upon the successful completion of these builds, please bsub the following command (saving STDOUT) to finish the rest of the steps for 3bp indel validation:\n\n";
        if ($ref_seq_build_id eq '101947881') { #if you are using the default, build 36, for both tools
            print "gmt validation long-indels-part-two --normal-val-model-copy-id $new_normal_model_id --tumor-val-model-copy-id $new_tumor_model_id --output-dir $output_dir\n\n";
        }
        else { #if you are using some other reference build (37)
            print "gmt validation long-indels-part-two --normal-val-model-copy-id $new_normal_model_id --tumor-val-model-copy-id $new_tumor_model_id --output-dir $output_dir --tier-file-location <PUT YOUR TIERING FILES HERE>\n\n";
        }

        #print details for a manual review ticket
        my $new_ref_build_fa = $new_ref_build->full_consensus_path('fa');
        print "And lastly, for your manual review tickets, you will want to include these details, along with further info printed in the STDOUT of 'gmt validation long-indels-part-two':\n\n";
        print "Original tumor validation BAM: $tumor_bam\n";
        print "Original normal validation BAM: $normal_bam\n";
        print "New reference sequence with contigs: $new_ref_build_fa\n";

    }

    $self->_new_tumor_model($new_tumor_model);
    $self->_new_normal_model($new_normal_model);

    return 1;
}


sub checkForDupName{
    my ($new_ref_model_name) = @_;
    my @models = Genome::Model->get("name LIKE" => "$new_ref_model_name%");    
    my $max=-1;
    if(@models > 0){   
        foreach my $m (@models){
            my $mname = $m->name;
            print STDERR $mname . "\n";
            #if we have models with a suffix already, store the highest suffix
            if ($mname=~/$new_ref_model_name-(\d+)/){
                if($1 > $max){
                        $max = $1;
                }
                #else if we have an exact match for this model name
            } elsif ($mname=~/$new_ref_model_name/){
                $max = 0;
            }
        }
        if($max > -1){
            $max++;
            $new_ref_model_name = $new_ref_model_name . "-$max";
        }
        }
    return($new_ref_model_name);
}


1;

