package Genome::Model::Somatic::Report::Variant;

use strict;
use warnings;

use Genome;
use Path::Class::Dir;

class Genome::Model::Somatic::Report::Variant {
    is => 'Genome::Model::Report',
};

#variant types
use constant SNP => 'SNP'; #single nucleotide polymorphism
use constant INS => 'INS'; #insertion
use constant DEL => 'DEL'; #deletion
use constant INV => 'INV'; #inversion
use constant CTX => 'CTX'; #cross-chromosomal translocation
use constant ITX => 'ITX'; #intra-chromosomal translocation

use constant VALIDATION_TYPE => 'Official';

sub description {
    my $self = shift();

    return 'listing of variants and structural variants';
}

sub _add_to_report_xml {
    my $self = shift();

    my $model = $self->model;
    my $build = $self->build;
    
    my $doc = $self->_xml;

    my $sample = $model->tumor_model->subject;
    undef $sample unless $sample->isa('Genome::Sample');
    
    my $individual = $model->subject;

    my $individual_node = $doc->createElement('individual');
    
    my $individual_name = $model->subject_name || '';

    if($individual) {
        $individual_name = $individual->common_name || $individual_name;
        
        $individual_node->addChild( $doc->createAttribute('id', $individual->id) );
        $individual_node->addChild( $doc->createAttribute('name', $individual_name) );
        $individual_node->addChild( $doc->createAttribute('gender', $individual->gender || '') );
    } else {
        $individual_node->addChild( $doc->createAttribute('id', 0) );
        $individual_node->addChild( $doc->createAttribute('name', $individual_name) );
        $individual_node->addChild( $doc->createAttribute('gender', '') );
    }
    
    #TODO Replace this with a real datasource once one is available.
#    if($individual_name =~ m/OVC[1-4]/) {
#        my $clinical_data_node = $individual_node->addChild( $doc->createElement('clinical-data') );
#        
#        my ($diagnosis_age, $diagnosis_year, $days_survived, $alive, $treatment, $outcome, $amp);
#        
#        if($individual_name eq 'OVC1') {
#            $diagnosis_age = '56';
#            $diagnosis_year = '2003';
#            $alive = 1;
#            $treatment = 'Chemotherapy';
#            $outcome = 'Complete response';
#            $amp = 'No CCNE1/MYC amp';
#        } elsif ($individual_name eq 'OVC2') {
#            $diagnosis_age = '63';
#            $days_survived = '1204';
#            $alive = 0;
#            $treatment = 'Chemotherapy, Radiation therapy, Hormone therapy';
#            $outcome = 'Progressive disease';
#            $amp = 'CCNE1/MYC amp';
#        } elsif ($individual_name eq 'OVC3') {
#            $diagnosis_age = '53';
#            $days_survived = '223';
#            $alive = 0;
#            $treatment = 'Chemotherapy';
#            $outcome = 'Progressive disease';
#            $amp = 'MYC amp/RB1 del';
#        } elsif ($individual_name eq 'OVC4') {
#            $diagnosis_age = '50';
#            $diagnosis_year = '2001';
#            $days_survived = '1646';
#            $alive = 0;
#            $treatment = 'Chemotherapy';
#            $amp = 'No CCNE1/MYC amp';
#        } else {
#            #Something unexpected happening
#        }
#        
#        $clinical_data_node->addChild( $doc->createAttribute('diagnosis-age', $diagnosis_age) )
#            if defined $diagnosis_age;
#        $clinical_data_node->addChild( $doc->createAttribute('diagnosis-year', $diagnosis_year) )
#            if defined $diagnosis_year;
#        $clinical_data_node->addChild( $doc->createAttribute('days-survived', $days_survived) )
#            if defined $days_survived;
#        $clinical_data_node->addChild( $doc->createAttribute('alive', $alive) )
#            if defined $alive;
#        $clinical_data_node->addChild( $doc->createAttribute('treatment', $treatment) )
#            if defined $treatment;
#        $clinical_data_node->addChild( $doc->createAttribute('outcome', $outcome) )
#            if defined $outcome;
#        $clinical_data_node->addChild( $doc->createAttribute('amp', $amp) )
#            if defined $amp;
#    }
    

    my $circos_node = $individual_node->addChild( $doc->createElement('circos-images') );

    my $variants_node = $individual_node->addChild( $doc->createElement('variants') );

    my $snps_node = $variants_node->addChild( $doc->createElement('snps') );
    my $insertions_node = $variants_node->addChild( $doc->createElement('insertions') );
    my $deletions_node = $variants_node->addChild( $doc->createElement('deletions') );

    my $structural_variants_node = $individual_node->addChild( $doc->createElement('structural-variants') );

    my $sv_insertions_node = $structural_variants_node->addChild( $doc->createElement('insertions') );
    my $sv_deletions_node = $structural_variants_node->addChild( $doc->createElement('deletions') );
    my $sv_inversions_node = $structural_variants_node->addChild( $doc->createElement('inversions') );
    my $sv_translocations_node = $structural_variants_node->addChild( $doc->createElement('translocations'));

    my $samples_node = $individual_node->addChild( $doc->createElement('samples') );

    my $sample_node = $samples_node->addChild( $doc->createElement('sample') );
    $sample_node->addChild( $doc->createAttribute('id', $sample->id) )
        if $sample;
    $sample_node->addChild( $doc->createAttribute('name', $sample->name) )
        if $sample;

    my $models_node = $sample_node->addChild( $doc->createElement('models') );

    my $model_node = $models_node->addChild( $doc->createElement('model') );
    $model_node->addChild( $doc->createAttribute('id', $model->id) );
    $model_node->addChild( $doc->createAttribute('name', $model->name) );
    $model_node->addChild( $doc->createAttribute('default', ($model->is_default ? '1' : '0')) );

    my $build_node = $model_node->addChild( $doc->createElement('build') );
    $build_node->addChild( $doc->createAttribute('id', $build->id) );

    #See if this is the build with the circos images
    my $data_dir = Path::Class::Dir->new($build->data_directory);
    my $circos_large = $data_dir->file('circos_graph.png.3000x3000.png');
    my $circos_small = $data_dir->file('circos_graph.png.920x920.png');
    my $circos_server = 'http://gscweb.gsc.wustl.edu/';
    if(-e $circos_large and -e $circos_small) {
        $circos_node->addChild( $doc->createAttribute('large', $circos_server . $circos_large) );
        $circos_node->addChild( $doc->createAttribute('small', $circos_server . $circos_small) );
    }
    
    my $tumor_model = $model->tumor_model;
    my $normal_model = $model->normal_model;
        
    my $tumor_build = $tumor_model->last_succeeded_build;
    my $normal_build = $normal_model->last_succeeded_build;
        
    if($tumor_build) {
        $model_node->addChild( $doc->createAttribute('tumor-haploid-coverage', $tumor_build->get_metric('haploid_coverage')) );
    }
    if($normal_build) {
        $model_node->addChild( $doc->createAttribute('normal-haploid-coverage', $normal_build->get_metric('haploid_coverage')) );
    }
    
#    if($individual_name =~ m/OVC[1-4]/) {
#        my ($purity_estimate, $ploidy_estimate);
#        
#        if($individual_name eq 'OVC1') {
#            $purity_estimate = '0.7';
#            #$ploidy_estimate = '';
#        } elsif ($individual_name eq 'OVC2') {
#            $purity_estimate = '0.75';
#            $ploidy_estimate = '1.8-2.5';
#        } elsif ($individual_name eq 'OVC3') {
#            $purity_estimate = '0.6-0.8';
#            $ploidy_estimate = '1.8-2.5';
#        } elsif ($individual_name eq 'OVC4') {
#            $purity_estimate = '0.8-0.9';
#            $ploidy_estimate = '1.8-2.5';
#        } else {
#            #Something unexpected happening
#        }
#        
#        $model_node->addChild( $doc->createAttribute('purity-estimate', $purity_estimate) )
#            if defined $purity_estimate;
#        $model_node->addChild( $doc->createAttribute('ploidy-estimate', $ploidy_estimate) )
#            if defined $ploidy_estimate;
#    }


## Variants

    my @build_variants = Genome::Model::BuildVariant->get( build_id => $build->id );
    
    unless(@build_variants) {
        #Some very early somatic models associated the variants with the tumor build instead
        @build_variants = Genome::Model::BuildVariant->get( build_id => $tumor_build->id );
    }
    
    my $snp_validated_count = 0;
    my $ins_validated_count = 0;
    my $del_validated_count = 0;

    for my $build_variant (@build_variants) {
        my $variant = $build_variant->variant;
        my $variant_type = $self->determine_variant_type($variant);

        my $trv_type = $variant->trv_type;

        #We're only interested in Tier 1 variants
        next unless
            grep($_ eq $trv_type, qw(silent splice_site_del splice_site_ins in_frame_del frame_shift_del rna frame_shift_ins in_frame_ins missense nonsense nonstop splice_site) );


        my $validation = Genome::Model::VariantValidation->get( variant_id => $variant->variant_id, model_id => $build_variant->build->model_id, validation_type => VALIDATION_TYPE  );
        unless ($validation) {
            $self->error_message("Data problems encountered -- Could not get validation information for variant id: " . $variant->variant_id . ". Exiting.");
            die;
        }
        
        my $validation_status = $validation->validation_result;
        
        $validation_status =~ s/\n//g; #Remove extraneous newlines froms status

        my $variant_node;
        

        if($variant_type eq SNP) {
            $variant_node = $snps_node->addChild( $doc->createElement('snp') );
            $snp_validated_count++ if $self->is_validated_status($validation_status);
        } elsif($variant_type eq INS) {
            $variant_node = $insertions_node->addChild( $doc->createElement('insertion') );
            $ins_validated_count++ if $self->is_validated_status($validation_status);
        } elsif($variant_type eq DEL) {
            $variant_node = $deletions_node->addChild( $doc->createElement('deletion') );
            $del_validated_count++ if $self->is_validated_status($validation_status);
        } else {
            die("Unexpected variant type: " . $variant_type);
        }

        $variant_node->addChild( $doc->createAttribute('variant-allele', $variant->variant_allele) )
            if $variant_type ne DEL;
        $variant_node->addChild( $doc->createAttribute('reference-allele', $variant->reference_allele) )
            if $variant_type ne INS;
        $variant_node->addChild( $doc->createAttribute('chromosome', $variant->chromosome) );
        $variant_node->addChild( $doc->createAttribute('start', $variant->start_pos) );
        $variant_node->addChild( $doc->createAttribute('stop', $variant->stop_pos) );
        $variant_node->addChild( $doc->createAttribute('gene', $variant->gene_name) )
            if defined $variant->gene_name;
        $variant_node->addChild( $doc->createAttribute('trv-type', $trv_type) );
        $variant_node->addChild( $doc->createAttribute('amino-acid-change', $variant->amino_acid_change) )
            if $variant->amino_acid_change and $variant->amino_acid_change ne '-' and $variant->amino_acid_change !~ /NULL/i;
        $variant_node->addChild( $doc->createAttribute('id', $variant->id) );
        $variant_node->addChild( $doc->createAttribute('from-build-id', $build_variant->build_id) );
        $variant_node->addChild( $doc->createAttribute('validation-status', $validation_status) );
    }
    
    $snps_node->addChild( $doc->createAttribute('validated-count', $snp_validated_count) );
    $insertions_node->addChild( $doc->createAttribute('validated-count', $ins_validated_count) );
    $deletions_node->addChild( $doc->createAttribute('validated-count', $del_validated_count) );


## Structural Variants

    my @build_variant_svs = Genome::Model::BuildSV->get( build_id => $build->id );
    unless(@build_variant_svs) {
        #Some very early somatic models associated the structural variants with the tumor build instead
        @build_variant_svs = Genome::Model::BuildSV->get( build_id => $tumor_build->id );
    }
    
    my $sv_ins_validated_count = 0;
    my $sv_del_validated_count = 0;
    my $sv_translocation_validated_count = 0;
    my $sv_inversion_validated_count;

    for my $build_variant_sv (@build_variant_svs) {
        my $structural_variant = Genome::Model::SV->get(variant_id => $build_variant_sv->variant_id);

        my $variant_type = $structural_variant->event_type;
        my $validation_status = $build_variant_sv->somatic_status; #TODO After real validation table is created, stop referring to temporary column.

        #TODO further temporary measure for consistency with variants
        if($validation_status =~ m/somatic/i) {
            $validation_status = 'S';
        }

        my $sv_node;

        if($variant_type eq INS) {
            $sv_node = $sv_insertions_node->addChild( $doc->createElement('insertion') );
            $sv_ins_validated_count++ if $self->is_validated_status($validation_status);
        } elsif ($variant_type eq DEL) {
            $sv_node = $sv_deletions_node->addChild( $doc->createElement('deletion') );
            $sv_del_validated_count++ if $self->is_validated_status($validation_status);
        } elsif ($variant_type eq CTX or $variant_type eq ITX) {
            $sv_node = $sv_translocations_node->addChild( $doc->createElement('translocation') );
            $sv_translocation_validated_count++ if $self->is_validated_status($validation_status);
        } elsif ($variant_type eq INV) {
            $sv_node = $sv_inversions_node->addChild( $doc->createElement('inversion') );
            $sv_inversion_validated_count++ if $self->is_validated_status($validation_status);
        } else {
            die("Unexpected SV event type: " . $variant_type);
        }

        $sv_node->addChild( $doc->createAttribute('id', $structural_variant->id) );
        $sv_node->addChild( $doc->createAttribute('size', $structural_variant->sv_size) )
            if $variant_type ne CTX; #Value meaningless
        $sv_node->addChild( $doc->createAttribute('from-build-id', $build_variant_sv->build_id) );
        $sv_node->addChild( $doc->createAttribute('validation-status', $validation_status) );

        my $source_node = $sv_node->addChild( $doc->createElement('start') );
        $source_node->addChild( $doc->createAttribute('chromosome', $structural_variant->chromosome_1) );
        $source_node->addChild( $doc->createAttribute('position', $structural_variant->pos_predicted_1) );

        my $destination_node = $sv_node->addChild( $doc->createElement('stop') );
        $destination_node->addChild( $doc->createAttribute('chromosome', $structural_variant->chromosome_2) );
        $destination_node->addChild( $doc->createAttribute('position', $structural_variant->pos_predicted_2) );
    }

    $sv_insertions_node->addChild( $doc->createAttribute('validated-count', $sv_ins_validated_count) );
    $sv_deletions_node->addChild( $doc->createAttribute('validated-count', $sv_del_validated_count) );
    $sv_inversions_node->addChild( $doc->createAttribute('validated-count', $sv_inversion_validated_count) );
    $sv_translocations_node->addChild( $doc->createAttribute('validated-count', $sv_translocation_validated_count) );


    $self->_main_node->addChild($individual_node);
    
    return 1;
}

sub determine_variant_type {
    my ($self, $variant) = @_;

    if(!defined $variant->reference_allele or $variant->reference_allele eq '0') {
        return INS;
    } elsif(!defined $variant->variant_allele or $variant->variant_allele eq '0') {
        return DEL;
    } elsif($variant->start_pos eq $variant->stop_pos) {
        return SNP;
    } else {
        die ("Unable to determine type of variant #" . $variant->variant_id);
    }
}

sub is_validated_status {
    my ($self, $status) = @_;
    
    return unless $status;
    
    return 0 if $status eq 'P' or $status eq 'predicted';
    
    return 1;
}

1;
