package Genome::Model::Tools::Capture::GermlineModelGroupQc;
use strict;
use warnings;
use FileHandle;
use Genome;
use Cwd;

my %stats = ();
my %already_reviewed = ();
my %wildtype_sites = my %germline_sites = ();

class Genome::Model::Tools::Capture::GermlineModelGroupQc {
    is => 'Genome::Command::Base',

    has_optional => [
        models                 => { is => 'Genome::Model', is_many => 1, shell_args_position => 1, doc => 'names or group of models to work on' },
        group_id               => { is => 'Genome::ModelGroup', doc => 'names or group of models to work on -- deprecated, can now list models/groups without --group_id' },
        use_external           => { is => 'Boolean', doc => 'Use external data source rather than internal/iscan', default_value => 0 },
        use_default            => { is => 'Boolean', doc => 'Use default data source as defined on sample (if available)', default_value => 0 },
        output_dir             => { is => 'Text', doc => "Outputs qc into directory for each sample", default => cwd() },
        summary_file           => { is => 'Text', doc => "Outputs qc summary into this file, must be run with already finished output (turns skip-if-output-present on)" },
        whitelist_snps_file    => { is => 'Text', doc => "File of snps to limit qc to, for example the 55 ASMS snps in ROI -- 1 rs_id per line" },
        skip_if_output_present => { is => 'Boolean', doc => "Skip Creating new qc Files if they exist", default => "" },
    ],
};

sub help_brief { "Operate on germline capture models and groups" }
sub help_synopsis { help_brief() . "\nEXAMPLE:    gmt capture germline-model-group-qc GROUP_ID/MODEL_NAMES --output-dir --dbsnp-build" }
sub help_detail { }

sub execute {
    my $self = shift;

    my @models;
    push @models, $self->models if $self->models;
    push @models, $self->group_id->models if $self->group_id;
    my $skip_if_output_present = $self->skip_if_output_present;
    my $summary_file = $self->summary_file;
    if ($self->summary_file) {
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
    my $build_number;
    my $db_snp_build;
    my $ref_name = $models[0]->reference_sequence_build->name;
    if($ref_name eq 'NCBI-human-build36') {
        $build_number = 36;
        $db_snp_build = 130;
        unless (grep{$_->reference_sequence_build->name eq 'NCBI-human-build36'}@models) {
            $self->error_message("Not all models are on NCBI-human-build36");
            return;
        }
    } elsif($ref_name eq 'GRCh37-lite-build37') {
        $build_number = 37;
        $db_snp_build = 132;
        unless (grep{$_->reference_sequence_build->name eq 'GRCh37-lite-build37'}@models) {
            $self->error_message("Not all models are on GRCh37-lite-build37");
            return;
        }
    } else {
        die "$ref_name isn't NCBI-human-build36 or GRCh37-lite-build37\n";
    }

    my $extract_params = $self->_resolve_gm_extract_command_params;
    return if not $extract_params;

    foreach my $model (@models) {
        my $subject_name = $model->subject_name || next;
        warn "$subject_name isn't a sample\n" and next unless $model->subject->isa('Genome::Sample');
        next if $subject_name =~ /Pooled/;
        if($model->last_succeeded_build) {
            my $bam_file = $model->last_succeeded_build->whole_rmdup_bam_file;

            my $qc_dir = $self->output_dir . "/$subject_name/";
            mkdir($qc_dir);
            my $genofile = "$qc_dir/$subject_name.dbsnp$db_snp_build.genotype";
            my $qcfile = "$qc_dir/$subject_name.dbsnp$db_snp_build.qc";

            if ($self->summary_file && -s $genofile && !-s $qcfile) {
                warn "You specified a summary file but the script thinks there are unfinished qc files, please run this script to finish making qc files first\n";
                warn "Reason: file $qcfile does not exist as a non-zero file\n";
                next;
            }
            if(!$self->summary_file && ( (! -e $genofile) || !$skip_if_output_present) ) {
                my $extract = Genome::Model::GenotypeMicroarray::Command::Extract->create(
                    sample => $model->subject,
                    variation_list_build => Genome::Model::ImportedVariationList->dbsnp_build_for_reference($model->reference_sequence_build),
                    output => $genofile.':separator=TAB:headers=chromosome,position,alleles,id:print_headers=0',
                    %$extract_params
                );
                unless ($extract->resolve_source){
                    $self->warning_message(
                        'Failed to resolve genotype build or instrumetn data for '.$model->subject->name.'skipping...'
                    );
                    next;
                }
                unless ($extract) {
                    $self->error_message("Failed to create Extract Microarray for sample " . $model->subject_name);
                    return;
                }
                $extract->dump_status_messages(1);

                unless ($extract->execute()) {
                    $self->error_message("Failed to execute Extract Microarray for sample " . $model->subject_name);
                    return;
                }

                unless (-s $genofile) {
                    $self->error_message("Executed Extract Microarray but geno file doesn't exist for sample " . $model->subject_name);
                    return;
                }
            }

            my $bsub = "bsub -N -M 4000000 -J $subject_name.dbsnp$db_snp_build.qc -o $qc_dir/$subject_name.dbsnp$db_snp_build.qc.out -e $qc_dir/$subject_name.dbsnp$db_snp_build.qc.err -R \"select[type==LINUX64 && mem>4000 && tmp>1000] rusage[mem=4000, tmp=1000]\"";
            my $cmd = $bsub." \'"."gmt analysis lane-qc compare-snps --genotype-file $genofile --bam-file $bam_file --output-file $qcfile --sample-name $subject_name --min-depth-het 20 --min-depth-hom 20 --flip-alleles 1 --verbose 1 --reference-build $build_number"."\'";
            if ($self->summary_file && !-s $genofile) {
                warn "You specified summary file but the script thinks there are unfinished qc files, please run this script to finish making qc files first\nReason: file $genofile does not exist as a non-zero file\n";
                next;
            }
            if ($self->summary_file) {
                my $qc_input = new FileHandle ($qcfile);
                my $qc_header = <$qc_input>;
                my $qc_line = <$qc_input>;
                chomp($qc_line);
                print ALL_MODELS "$db_snp_build\t$qc_line\n";
            }
            elsif(-s $genofile && ((! -e $qcfile) || !$skip_if_output_present)) {
                system("$cmd");
            }
        }
    }
    if ($self->summary_file) {
        close(ALL_MODELS);
    }
    return 1;
}

sub _resolve_gm_extract_command_params {
    my $self = shift;

    my %extract_params;

    # Sample type priority
    if ( $self->use_default and $self->use_external ) {
        $self->error_message('Specifying use_default and use_external at the same time is not allowed! Please pick one or neither.');
        return;
    }

    if ( $self->use_default ) {
        $extract_params{sample_type_priority} = [ 'default' ];
    }
    elsif ( $self->use_external ) {
        $extract_params{sample_type_priority} = [ 'external' ];
    }

    # Filters
    if ( $self->whitelist_snps_file ) {
        $extract_params{filters} = [ 'whitelist:whitelist_snps_file='.$self->whitelist_snps_file ];
    }

    return \%extract_params;
}

sub byChrPos {
    my ($chr_a, $pos_a) = split(/\t/, $a);
    my ($chr_b, $pos_b) = split(/\t/, $b);

    $chr_a cmp $chr_b
        or
    $pos_a <=> $pos_b;
}
1;
