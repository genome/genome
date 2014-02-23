package Genome::Model::Tools::Capture::GermlineModelGroupQcIterative;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ModelGroup - Build Genome Models for Germline Capture Datasets
#
#    AUTHOR:        Will Schierding
#
#    CREATED:    2/09/2011 by W.S.
#    MODIFIED:    2/09/2011 by W.S.
#
#    NOTES:
#
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;
use Cwd;
use Genome;

class Genome::Model::Tools::Capture::GermlineModelGroupQcIterative {
    is => 'Command',

    has => [
        group_id     => { is => 'Text', doc => "ID of model group" },
    ],
    has_optional => [
        summary_file            => { is => 'Text', doc => "Outputs qc summary into this file, must be run with already finished output (turns skip-if-output-present on)" },
        output_dir              => { is => 'Text', doc => "Outputs qc into directory for each sample", default => cwd() },
        whitelist_snps_file     => { is => 'Text', doc => "File of snps to limit qc to, for example the 55 ASMS snps in ROI -- 1 rs_id per line" },
        skip_if_output_present  => { is => 'Boolean', doc => "Skip Creating new qc Files if they exist" , default => ""},
        cleanup_empty_files     => { is => 'Boolean', doc => "Delete files that pileup likely failed on so they'll re-run" , default => ""},
        use_external            => { is => 'Boolean', doc => 'Use external data source rather than internal/iscan', default_value => 0 },
    ],
};

sub help_brief {
    "Operate on germline capture model groups"
}

sub help_synopsis {
    return <<EOS
Operate on capture somatic model groups
EXAMPLE:    gmt capture germline-model-group-qc --group-id XXXX --output-dir --dbsnp-build
EOS
}

sub help_detail { '' }


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
    my $self = shift;
    my @models = Genome::ModelGroup->get($self->group_id)->models;
    my $skip_if_output_present = $self->skip_if_output_present;
    my $empty_file_cleanup = $self->skip_if_output_present;
    my $summary_file;
    if ($self->summary_file) {
        $summary_file = $self->summary_file;
        $skip_if_output_present = 1;
        unless (open(ALL_MODELS,">$summary_file")) {
            die "Could not open input file '$summary_file' for reading";
        }
        print ALL_MODELS join("\t",qw(
            Dbsnp_Build
            Sample_id
            SNPsCalled
            WithGenotype
            MetMinDepth
            Reference
            RefMatch
            RefWasHet
            RefWasHom
            Variant
            VarMatch
            HomWasHet
            HetWasHom
            VarMismatch
            VarConcord
            RareHomConcord
            OverallConcord
            )) . "\n";
    }

    # Correct the reference build name to what the database recognizes
    my $reference;
    my $build_number;
    my $db_snp_build;

    my $ref = $models[0]->reference_sequence_build;
    my $build36 = Genome::Model::Build::ReferenceSequence->get(name => "NCBI-human-build36");
    my $build37 = Genome::Model::Build::ReferenceSequence->get(name => "GRCh37-lite-build37");
    if($ref->is_compatible_with($build36)) {
        $reference = 'reference';
        $build_number = 36;
        $db_snp_build = 130;
    } elsif($ref->is_compatible_with($build37)) {
        $reference = 'GRCh37';
        $build_number = 37;
        $db_snp_build = 132;
    } else {
        die $ref->name." isn't compatible with NCBI-human-build36 or GRCh37-lite-build37\n";
    }

    my %snp_limit_hash;
    if($self->whitelist_snps_file) {
        my $snp_input = new FileHandle ($self->whitelist_snps_file);
        unless($snp_input) {
            $self->error_message("Unable to open ".$self->whitelist_snps_file);
            return;
        }

        while (my $line = <$snp_input>) {
            chomp($line);
            my ($id) = split(/\t/, $line);
            $snp_limit_hash{$id}++;
        }
    }

    my %qc_iteration_hash_genotype;
    my %qc_iteration_hash_bam_file;
    foreach my $model (@models) {
        my $model_id = $model->id;
        my $subject_name = $model->subject_name || next;
        next if $subject_name =~ m/Pooled/;
        if($model->last_succeeded_build_directory) {
            my $build = $model->last_succeeded_build;
            my $build_id = $build->id;
            my $last_build_dir = $model->last_succeeded_build_directory;
            my $bam_file = $build->whole_rmdup_bam_file;

            if($self->output_dir) {
                my $qc_dir = $self->output_dir . "/$subject_name/";
                mkdir($qc_dir);
                my $genofile = "$qc_dir/$subject_name.dbsnp$db_snp_build.genotype";
                if ($self->summary_file && ! -s $genofile) {
                    warn "You specified summary file but the script thinks there are unfinished genotype files, please run this script to finish making qc files first\nReason: file $genofile does not exist as a non-zero file\n";
                    next;
                }

                if(!$self->summary_file && ( (! -e $genofile) || !$skip_if_output_present) ) {
                    my %extract_params = (
                        output => $genofile,
                        fields => [qw(chromosome position alleles id)],
                        sample => $model->subject,
                        ($self->whitelist_snps_file?(filters => ['whitelist:whitelist_snps_file='.$self->whitelist_snps_file]):()),
                    );

                    if ($build->can('dbsnp_build') && $build->dbsnp_build) {
                        $extract_params{variation_list_build} = $build->dbsnp_build;
                    } else {
                        $extract_params{variation_list_build} = Genome::Model::ImportedVariationList->dbsnp_build_for_reference($build->reference_sequence_build),
                    }

                    if ($self->use_external) {
                        $extract_params{use_external} = $self->use_external,
                    } elsif ($build->can('genotype_microarray_build') && $build->genotype_microarray_build) {
                        $extract_params{instrument_data} = $build->genotype_microarray_build->instrument_data;
                    } elsif ($build->subject->can('default_genotype_data_id') && $build->subject->default_genotype_data_id) {
                        $extract_params{instrument_data} = Genome::InstrumentData::Imported->get($build->subject->default_genotype_data_id);
                    }

                    my $extract = Genome::InstrumentData::Command::Microarray::Extract->create(%extract_params);
                    unless ($extract) {
                        $self->error_message("Failed to create Extract Microarray for sample " . $model->subject_name);
                        return;
                    }

                    $extract->dump_status_messages(1);
                    unless ($extract->_resolve_instrument_data){
                        $self->debug_message('Skipping due to no instrument data for sample ' . $model->subject->id);
                        next;
                    }
                    unless ($extract->execute()) {
                        $self->error_message("Failed to execute Extract Microarray for sample " . $model->subject_name);
                        return;
                    }
                    unless (-s $genofile) {
                        $self->error_message("Executed Extract Microarray but geno file doesn't exist for sample " . $model->subject_name);
                        return;
                    }
                }

                $qc_iteration_hash_genotype{$subject_name}{$genofile}++;
                $qc_iteration_hash_bam_file{$subject_name}{$bam_file}++;
            }
        }
    }

    my $halt_submissions = 0;
    foreach my $subject_name1 (sort keys %qc_iteration_hash_genotype) {
        foreach my $genofile (sort keys %{$qc_iteration_hash_genotype{$subject_name1}}) {
            foreach my $subject_name2 (sort keys %qc_iteration_hash_bam_file) {
                foreach my $bam_file (sort keys %{$qc_iteration_hash_bam_file{$subject_name2}}) {
#                    print "Genosample:$subject_name1\t$genofile\nBamsample:$subject_name2\t$bam_file\n";
                    my $qc_dir = $self->output_dir . "/qc_iteration_files/";
                    mkdir($qc_dir);
                    my $qcfile = "$qc_dir/Geno_$subject_name1.Bam_$subject_name2.dbsnp$db_snp_build.qc";
                    my $output_bsub = "$qc_dir/Geno_$subject_name1.Bam_$subject_name2.dbsnp$db_snp_build.out";
                    my $error_bsub = "$qc_dir/Geno_$subject_name1.Bam_$subject_name2.dbsnp$db_snp_build.err";
                    my $bsub = "bsub -N -M 4000000 -J Geno_$subject_name1.Bam_$subject_name2.dbsnp$db_snp_build.qc -o $output_bsub -e $error_bsub -R \"select[type==LINUX64 && mem>4000 && tmp>1000] rusage[mem=4000, tmp=1000]\"";
                    my $cmd = $bsub." \'"."gmt analysis lane-qc compare-snps --genotype-file $genofile --bam-file $bam_file --output-file $qcfile --sample-name Geno_$subject_name1.Bam_$subject_name2 --min-depth-het 20 --min-depth-hom 20 --flip-alleles 1 --verbose 1 --reference-build $build_number"."\'";

                    #clean up empty qc files
                    if ($skip_if_output_present && $empty_file_cleanup && -s $qcfile &&1&&1) { #&&1&&1 to make gedit show colors correctly after a -s check
                        my $qc_input = new FileHandle ($qcfile);
                        my $qc_header = <$qc_input>;
                        my $qc_line = <$qc_input>;
                        chomp($qc_line);
                        my ($sample, $sequence_snps, $geno_snps, $covered_snps, @everything_else) = split(/\t/, $qc_line);
                        if ($geno_snps == 0) {
                            print "$qc_line\n";exit;
                            if (unlink($qcfile) == 0) {
                                print "File deleted successfully.";
                            } else {
                                print "File was not deleted.";
                            }
                            if (unlink($output_bsub) == 0) {
                                print "File deleted successfully.";
                            } else {
                                print "File was not deleted.";
                            }
                            if (unlink($error_bsub) == 0) {
                                print "File deleted successfully.";
                            } else {
                                print "File was not deleted.";
                            }
                        }
                    }
                    if ($skip_if_output_present && -s $qcfile &&1&&1) { #&&1&&1 to make gedit show colors correctly after a -s check
                    }
                    elsif ($self->summary_file) {
                        warn "You specified summary file but the script thinks there are unfinished qc files, please run this script to finish making qc files first\nReason: file $qcfile does not exist as a non-zero file\n";
                    }
                    else {
                        if ($halt_submissions > 200) {
                            $halt_submissions = 0;
                            sleep (120);
                        }
                        system("$cmd");
                        $halt_submissions++;
                    }
                }
            }
        }
    }

    if ($self->summary_file) {
        foreach my $subject_name1 (sort keys %qc_iteration_hash_genotype) {
            foreach my $genofile (sort keys %{$qc_iteration_hash_genotype{$subject_name1}}) {
                foreach my $subject_name2 (sort keys %qc_iteration_hash_bam_file) {
                    foreach my $bam_file (sort keys %{$qc_iteration_hash_bam_file{$subject_name2}}) {
                        my $qc_dir = $self->output_dir . "/qc_iteration_files/";
                        my $qcfile = "$qc_dir/Geno_$subject_name1.Bam_$subject_name2.dbsnp$db_snp_build.qc";
                        my $qc_input = new FileHandle ($qcfile);
                        my $qc_header = <$qc_input>;
                        my $qc_line = <$qc_input>;
                        chomp($qc_line);
                        print ALL_MODELS "$db_snp_build\t$qc_line\n";
                    }
                }
            }
        }
        close(ALL_MODELS);
    }

    return 1;
}

sub byChrPos
{
    my ($chr_a, $pos_a) = split(/\t/, $a);
    my ($chr_b, $pos_b) = split(/\t/, $b);

    $chr_a cmp $chr_b
        or
    $pos_a <=> $pos_b;
}

1;
