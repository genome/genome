package Genome::Model::Tools::Somatic::UploadVariants;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;

class Genome::Model::Tools::Somatic::UploadVariants{
    is => 'Command',
    has => [
    variant_file => {
        is  => 'String',
        is_input => '1',
        doc => 'The file of somatic pipeline results to be uploaded. This will usually be a high confidence tier 1 or 2 snp file, or a tier 1 indel file from the somatic pipeline.',
    },
    annotation_file => {
        is  => 'String',
        is_input => '1',
        doc => 'The file containing the annotation of all of the variants from the corresponding variant file. This will usually be the annotation output for snps or indels from the somatic pipeline.',
    },
    output_file => {
        is  => 'String',
        is_input => '1',
        is_output => '1',
        doc => 'The output file containing all of the annotation lines from the annotation_file for each of the variants from the variant_file',
    },
    build_id => {
        is => 'Number',
        is_input => '1',
        doc => 'The build id that should be linked to the variant. Enter a number <= 0 to skip uploading (test mode).',
    },
    _skip => {
        is => 'Boolean',
        default => '0',
        is_input => 1,
        is_optional => 1,
        doc => "If set to true... this will do nothing! Fairly useless, except this is necessary for workflow.",
    },
    ],
};

sub help_brief {
    "Adds results from the somatic pipeline to the database tables for tracking all known tier 1 and 2 variants.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
    gmt somatic upload-results --variant-file high_confidence_file.out --annotation-file annotation_file.out --output-file upload.out
EOS
}

sub help_detail {                           
    return <<EOS 
Adds results from the somatic pipeline to the database tables for tracking all known tier 1 and 2 variants. 
EOS
}

sub execute {
    my $self = shift;
    
    if ($self->_skip) {
        $self->status_message("Skipping execution: Skip flag set");
        return 1;
    }
    
    my @accession_list = @{$self->accession_list};
    my $accession_lookup = $self->make_hash_lookup(@accession_list);
      
    my $variant_fh = IO::File->new($self->variant_file);
    unless ($variant_fh) {
        $self->error_message("Could not open variant file: " . $self->variant_file . " for reading. $!");
        die;
    }

    my $annotation_fh = IO::File->new($self->annotation_file);
    unless ($annotation_fh) {
        $self->error_message("Could not open annotation file: " . $self->annotation_file . " for reading. $!");
        die;
    }

    my $ofh = IO::File->new($self->output_file, "w");
    unless($ofh) {
        $self->error_message("Unable to open " . $self->output_file . " for writing. $!");
        die;
    }

    # Fill the variant hash to later lookup and grab ALL annotation lines for each variant... help... i've turned into dlarson (just kidding dave, I <3 you)
    # FIXME this is pretty hacky going through both files and using the hash and everything... it would be nice to clean this up but for now lets get it working
    my %annotation;
    while(my $line = $annotation_fh->getline) {
        chomp $line;
        my ($chr, $start, $stop, $reference, $variant, $variation_type, $gene, $transcript, $species, $transcript_source, $transcript_version, $strand, $transcript_status, $type, $aa_string) = split "\t", $line;
        push ( @{$annotation{$chr}{$start}{$stop}{$reference}{$variant}}, $line);
    }
    
    # Go through each line in the variant file and get each annotation line that matches from the annotation file
    # For each line, print it to the output file and upload it to the database
    while (my $line = $variant_fh->getline) {
        my ($chr, $start, $stop, $reference, $variant) = split "\t", $line;
        
        # Get each possible variant from IUB code unless the variant is an indel
        my @variant_alleles;
        if($reference eq "0" || $variant eq "0") {
            @variant_alleles = ($variant);
        }else {
            @variant_alleles = Genome::Info::IUB->variant_alleles_for_iub($reference, $variant);
        }

        for my $variant_allele (@variant_alleles) {
            # There should be annotation for each variant line, or something went wrong in the pipeline
            unless (defined $annotation{$chr}{$start}{$stop}{$reference}{$variant_allele}) {
                $self->error_message("Could not find annotation for variant: $chr $start $stop $reference $variant");
                die;
            }

            # This should hold the entire annotation line from the transcript annotation file
            for my $annotation (@{$annotation{$chr}{$start}{$stop}{$reference}{$variant_allele}}) {
                $ofh->print("$annotation\n");

                # If the build id is <= 0 we are in test mode and are not uploading
                next if ($self->build_id <= 0);
                
                my ($chr, $start, $stop, $reference, $variant, $variation_type, $gene, $transcript, $species, $transcript_source, $transcript_version, $strand, $transcript_status, $trv_type, $c_position, $amino_acid_change, $ucsc_cons, $domain) = split("\t", $annotation);
                if (length($amino_acid_change) > 255) {
                    $amino_acid_change = substr($amino_acid_change,0,240)."...truncated";
                }
                if (length($ucsc_cons) > 255) {
                    $ucsc_cons = substr($ucsc_cons,0,240)."...truncated";
                }
                my $accession_domains = $self->make_domain_into_ids($domain, $accession_lookup);
                my $new_variant;
                my $variant_already_exists = Genome::Model::Variant->get(
                    chromosome      => $chr,
                    start_pos       => $start,
                    stop_pos        => $stop,
                    reference_allele=>$reference,
                    variant_allele  =>$variant
                );
                if($variant_already_exists) {
                    $new_variant = $variant_already_exists;
                }
                else{
                    $new_variant = Genome::Model::Variant->create(
                        chromosome         => $chr,
                        start_pos          => $start,
                        stop_pos           => $stop,
                        reference_allele   => $reference,
                        variant_allele     => $variant,
                        gene_name          => $gene,
                        transcript_name    => $transcript,
                        transcript_source  => $transcript_source,
                        transcript_version => $transcript_version,
                        strand             => $strand,
                        transcript_status  => $transcript_status,
                        trv_type           => $trv_type,
                        c_position         => $c_position,
                        amino_acid_change  => $amino_acid_change,
                        ucsc_cons          => $ucsc_cons,
                        domain             => $accession_domains,
                    );
                }
                unless($new_variant) {
                    $self->error_message("UR was unable to instantiate a variant object for this line:");
                    $self->error_message($annotation);
                    die;
                }

                my $new_build_variant = Genome::Model::BuildVariant->get_or_create(
                    variant => $new_variant,
                    build_id => $self->build_id,
                );
                
                my $build = Genome::Model::Build->get($self->build_id);
                unless (defined $build) {
                    $self->error_message("Could not get a build for build id " . $self->build_id . ". Please use a valid build id.");
                    die;
                }
                
                my $model = $build->model;
                my $new_variant_validation = Genome::Model::VariantValidation->get_or_create(
                    variant=>$new_variant,
                    validation_type=>'Official',
                    validation_result=>'P',
                    model_id =>$model->id,
                );
                unless($new_build_variant && $new_variant_validation) {
                    $self->error_message("Unable to create Build-Variant Link OR VariantValidation Status");
                    $self->error_message("Problem line: $annotation");
                    die;
                }
            }
        }
    }

    return 1;

}

sub make_domain_into_ids {
    my ($self, $domain, $accession_list) = @_;
    my $new_list;
    my @domains = split /,/, $domain;
    for my $key (@domains) {
        my $new_domain = $accession_list->{$key};
        if (!defined($new_domain)) {
            my $new_key = "HMMPfam_".$key;
            $new_domain = $accession_list->{$new_key};
        }
        if (!defined($new_domain)) {
            my $new_key = "superfamily_".$key;
            $new_domain = $accession_list->{$new_key};
        }
        $new_list .= "$new_domain ";
    }
    
    return $new_list;
} 

sub make_hash_lookup {
    my $self=shift;
    my @accession_list = @_;
    my %hash_lookup;
    for my $list_item (@accession_list) {
        my ($access, $human_readable) = split /\t/, $list_item;
        $hash_lookup{$human_readable}=$access;
    }
    $hash_lookup{"-"}="-";
    $hash_lookup{"NULL"}="NULL";
    return \%hash_lookup;
}

sub accession_list {
    my $self=shift;
    my $accession_list = <<EOF;
SSF46585	superfamily_HR1 repeat
PF00168	HMMPfam_C2
PF00388	HMMPfam_PI-PLC-X
PF00017	HMMPfam_SH2
PF00018	HMMPfam_SH3_1
SSF50044	superfamily_SH3-domain
PF00387	HMMPfam_PI-PLC-Y
PF00169	HMMPfam_PH
SSF49562	superfamily_C2 domain (Calcium/lipid-binding domain CaLB)
SSF47473	superfamily_EF-hand
SSF50729	superfamily_PH domain-like
SSF51695	superfamily_PLC-like phosphodiesterases
SSF55550	superfamily_SH2 domain
PF00047	HMMPfam_ig
SSF48726	superfamily_Immunoglobulin
PF03066	HMMPfam_Nucleoplasmin
SSF69203	superfamily_Nucleoplasmin-like core domain
PF04908	HMMPfam_SH3BGR
SSF52833	superfamily_Thioredoxin-like
SSF53137	superfamily_Translational machinery components
PF07686	HMMPfam_V-set
SSF47454	superfamily_A DNA-binding domain in eukaryotic transcription factors
PF00170	HMMPfam_bZIP_1
PF01769	HMMPfam_MgtE
SSF52151	superfamily_FabD/lysophospholipase-like
PF02996	HMMPfam_Prefoldin
SSF46579	superfamily_Prefoldin
PF00076	HMMPfam_RRM_1
PF00658	HMMPfam_PABP
SSF54928	superfamily_RNA-binding domain RBD
PF00454	HMMPfam_PI3_PI4_kinase
PF02260	HMMPfam_FATC
PF07714	HMMPfam_Pkinase_Tyr
PF09027	HMMPfam_GTPase_binding
PF00794	HMMPfam_PI3K_rbd
PF00613	HMMPfam_PI3Ka
HMMPfam_PF00787	HMMPfam_PX
PF00792	HMMPfam_PI3K_C2
PF00621	HMMPfam_RhoGEF
PF00307	HMMPfam_CH
PF00130	HMMPfam_C1_1
PF00786	HMMPfam_PBD
PF00069	HMMPfam_Pkinase
PF00433	HMMPfam_Pkinase_C
PF00780	HMMPfam_CNH
PF08826	HMMPfam_DMPK_coil
PF00488	HMMPfam_MutS_V
PF01624	HMMPfam_MutS_I
PF05192	HMMPfam_MutS_III
PF05188	HMMPfam_MutS_II
PF05190	HMMPfam_MutS_IV
SSF57196	superfamily_EGF/Laminin
PF00226	HMMPfam_DnaJ
PF00320	HMMPfam_GATA
PF00349	HMMPfam_Hexokinase_1
PF03727	HMMPfam_Hexokinase_2
PF06390	HMMPfam_NESP55
PF00788	HMMPfam_RA
PF08947	HMMPfam_BPS
PF00618	HMMPfam_RasGEF_N
PF00617	HMMPfam_RasGEF
PF00855	HMMPfam_PWWP
PF02020	HMMPfam_W2
PF02854	HMMPfam_MIF4G
PF02847	HMMPfam_MA3
PF00071	HMMPfam_Ras
PF07710	HMMPfam_P53_tetramer
PF07647	HMMPfam_SAM_2
PF00870	HMMPfam_P53
PF07653	HMMPfam_SH3_2
PF02828	HMMPfam_L27
PF07679	HMMPfam_I-set
PF02985	HMMPfam_HEAT
PF08064	HMMPfam_UME
PF08311	HMMPfam_Mad3_BUB1_I
PF02192	HMMPfam_PI3K_p85B
PF00754	HMMPfam_F5_F8_type_C
PF00041	HMMPfam_fn3
PF00010	HMMPfam_HLH
SSF47459	superfamily_HLH helix-loop-helix DNA-binding domain
PF01056	HMMPfam_Myc_N
PF02344	HMMPfam_Myc-LZ
PF00612	HMMPfam_IQ
PF00063	HMMPfam_Myosin_head
SSF56112	superfamily_Protein kinase-like (PK-like)
SSF52540	superfamily_P-loop containing nucleoside triphosphate hydrolases
PF00373	HMMPfam_Band_41
PF00641	HMMPfam_zf-RanBP
PF02201	HMMPfam_SWIB
SSF47592	superfamily_SWIB/MDM2 domain
SSF57850	superfamily_RING/U-box
SSF90209	superfamily_Znf265 first zinc-finger domain
PF05053	HMMPfam_Menin
SSF49265	superfamily_Fibronectin type III
SSF52821	superfamily_Rhodanese/Cell cycle control phosphatase
PF00008	HMMPfam_EGF
SSF49899	superfamily_Concanavalin A-like lectins/glucanases
PF06816	HMMPfam_NOD
PF07684	HMMPfam_NODP
PF07645	HMMPfam_EGF_CA
PF07974	HMMPfam_EGF_2
PF00755	HMMPfam_Carn_acyltransf
PF00066	HMMPfam_Notch
PF00023	HMMPfam_Ank
SSF48097	superfamily_Regulator of G-protein signaling RGS
PF00096	HMMPfam_zf-C2H2
SSF109640	superfamily_KRAB domain (Kruppel-associated box Pfam 01352)
SSF57667	superfamily_C2H2 and C2HC zinc fingers
PF00001	HMMPfam_7tm_1
SSF81321	superfamily_SSF81321
PF05624	HMMPfam_LSR
PF02038	HMMPfam_ATP1G1_PLM_MAT8
PF07988	HMMPfam_Wos2
PF00249	HMMPfam_Myb_DNA-binding
PF09316	HMMPfam_Cmyb_C
PF00595	HMMPfam_PDZ
PF08926	HMMPfam_DUF1908
PF02185	HMMPfam_HR1
PF00640	HMMPfam_PID
PF00400	HMMPfam_WD40
SSF50978	superfamily_WD40 repeat-like
SSF48371	superfamily_ARM repeat
SSF57889	superfamily_Cysteine-rich domain
SSF90193	superfamily_Notch domain
PF08447	HMMPfam_PAS_3
PF00989	HMMPfam_PAS
SSF55785	superfamily_PYP-like sensor domain (PAS domain)
SSF81321	superfamily_Family A G protein-coupled receptor-like
PF07885	HMMPfam_Ion_trans_2
SSF81324	superfamily_Voltage-gated potassium channels
SSF50156	superfamily_PDZ domain-like
PF01825	HMMPfam_GPS
PF00002	HMMPfam_7tm_2
PF00751	HMMPfam_DM
SSF82927	superfamily_Cysteine-rich DNA binding domain (DM domain)
PF08917	HMMPfam_ecTbetaR2
PF01414	HMMPfam_DSL
PF07657	HMMPfam_MNNL
PF01352	HMMPfam_KRAB
PF00090	HMMPfam_TSP_1
PF01421	HMMPfam_Reprolysin
PF01562	HMMPfam_Pep_M12B_propep
PF05986	HMMPfam_ADAM_spacer1
PF08686	HMMPfam_PLAC
SSF55486	superfamily_Metalloproteases ("zincins") catalytic domain
SSF57059	superfamily_omega toxin-like
SSF82895	superfamily_TSP-1 type 1 repeat
PF00058	HMMPfam_Ldl_recept_b
PF02757	HMMPfam_YLP
PF00757	HMMPfam_Furin-like
PF01392	HMMPfam_Fz
SSF46565	superfamily_Chaperone J-domain
SSF52799	superfamily_(Phosphotyrosine protein) phosphatases II
SSF57716	superfamily_Glucocorticoid receptor-like (DNA-binding domain)
PF00167	HMMPfam_FGF
SSF50353	superfamily_Cytokine
PF01833	HMMPfam_TIG
PF01347	HMMPfam_Vitellogenin_N
PF06448	HMMPfam_DUF1081
PF09172	HMMPfam_DUF1943
PF00054	HMMPfam_Laminin_G_1
PF02210	HMMPfam_Laminin_G_2
PF02065	HMMPfam_Melibiase
PF00200	HMMPfam_Disintegrin
PF08516	HMMPfam_ADAM_CR
PF00102	HMMPfam_Y_phosphatase
PF00629	HMMPfam_MAM
PF00059	HMMPfam_Lectin_C
PF01477	HMMPfam_PLAT
PF02010	HMMPfam_REJ
PF08016	HMMPfam_PKD_channel
PF00084	HMMPfam_Sushi
PF00178	HMMPfam_Ets
PF02198	HMMPfam_SAM_PNT
PF00611	HMMPfam_FCH
PF00206	HMMPfam_Lyase_1
PF00250	HMMPfam_Fork_head
PF07716	HMMPfam_bZIP_2
PF03798	HMMPfam_LAG1
SSF46689	superfamily_Homeodomain-like
PF00646	HMMPfam_F-box
PF04300	HMMPfam_FBA
SSF49785	superfamily_Galactose-binding domain-like
SSF81383	superfamily_F-box domain
SSF48403	superfamily_Ankyrin repeat
PF00620	HMMPfam_RhoGAP
PF01436	HMMPfam_NHL
PF00057	HMMPfam_Ldl_recept_a
PF01030	HMMPfam_Recep_L_domain
SSF57184	superfamily_Growth factor receptor domain
SSF52058	superfamily_L domain-like
SSF47769	superfamily_SAM/Pointed domain
PF01404	HMMPfam_Ephrin_lbd
PF00536	HMMPfam_SAM_1
SSF46785	superfamily_"Winged helix" DNA-binding domain
SSF47668	superfamily_N-terminal domain of cbl (N-cbl)
PF02984	HMMPfam_Cyclin_C
PF00134	HMMPfam_Cyclin_N
SSF47954	superfamily_Cyclin-like
PF06617	HMMPfam_M-inducer_phosp
PF00581	HMMPfam_Rhodanese
PF02234	HMMPfam_CDI
PF00498	HMMPfam_FHA
SSF49879	superfamily_SMAD/FHA domain
SSF57184	superfamily_Grow_fac_recept
SSF56112	superfamily_Kinase_like
SSF52058	superfamily_SSF52058
SSF81296	superfamily_E set domains
SSF101912	superfamily_Sema domain
SSF103575	superfamily_Plexin repeat
PF00779	HMMPfam_BTK
PF00531	HMMPfam_Death
PF06733	HMMPfam_DEAD_2
PF06777	HMMPfam_DUF1227
PF02174	HMMPfam_IRS
SSF56672	superfamily_DNA/RNA polymerases
PF01391	HMMPfam_Collagen
SSF56496	superfamily_Fibrinogen C-terminal domain-like
PF00092	HMMPfam_VWA
SSF53300	superfamily_vWA-like
PF00112	HMMPfam_Peptidase_C1
PF08246	HMMPfam_Inhibitor_I29
SSF54001	superfamily_Cysteine proteinases
PF01302	HMMPfam_CAP_GLY
PF00965	HMMPfam_TIMP
SSF50242	superfamily_TIMP-like
PF00643	HMMPfam_zf-B_box
PF00439	HMMPfam_Bromodomain
PF00097	HMMPfam_zf-C3HC4
PF00628	HMMPfam_PHD
SSF57903	superfamily_FYVE/PHD zinc finger
SSF47370	superfamily_Bromodomain
SSF57845	superfamily_B-box zinc-binding domain
PF03836	HMMPfam_RasGAP_C
PF00397	HMMPfam_WW
PF00616	HMMPfam_RasGAP
SSF48350	superfamily_GTPase activation domain GAP
SSF47576	superfamily_Calponin-homology domain CH-domain
PF00377	HMMPfam_Prion
PF03991	HMMPfam_Prion_octapep
SSF54897	superfamily_Protease propeptides/inhibitors
SSF52743	superfamily_Subtilisin-like
PF04597	HMMPfam_Ribophorin_I
SSF53271	superfamily_PRTase-like
PF00609	HMMPfam_DAGK_acc
PF00781	HMMPfam_DAGK_cat
PF00036	HMMPfam_efhand
SSF52518	superfamily_Thiamin diphosphate-binding fold (THDP-binding)
PF00533	HMMPfam_BRCT
PF00452	HMMPfam_Bcl-2
PF00651	HMMPfam_BTB
PF01064	HMMPfam_Activin_recp
PF08515	HMMPfam_TGF_beta_GS
PF02196	HMMPfam_RBD
PF00634	HMMPfam_BRCA2
PF09103	HMMPfam_BRCA-2_OB1
PF09104	HMMPfam_BRCA-2_OB3
PF09121	HMMPfam_Tower
PF09169	HMMPfam_BRCA-2_helical
PF00653	HMMPfam_BIR
SSF57924	superfamily_Inhibitor of apoptosis (IAP) repeat
PF00111	HMMPfam_Fer2
SSF54292	superfamily_2Fe-2S ferredoxin-like
SSF54117	superfamily_Interleukin 8-like chemokines
PF00683	HMMPfam_TB
PF08563	HMMPfam_P53_TAD
PF01857	HMMPfam_RB_B
PF01858	HMMPfam_RB_A
PF08934	HMMPfam_Rb_C
PF00051	HMMPfam_Kringle
PF00089	HMMPfam_Trypsin
PF00024	HMMPfam_PAN_1
PF01585	HMMPfam_G-patch
PF00564	HMMPfam_PB1
PF00630	HMMPfam_Filamin
SSF54236	superfamily_Ubiquitin-like
SSF64268	superfamily_PX domain
PF02518	HMMPfam_HATPase_c
SSF55874	superfamily_ATPase domain of HSP90 chaperone/DNA topoisomerase II/histidine kinase
SSF69012	superfamily_alpha-ketoacid dehydrogenase kinase N-terminal domain
PF00046	HMMPfam_Homeobox
SSF63570	superfamily_PABC (PABP) domain
PF01929	HMMPfam_Ribosomal_L14e
PF08912	HMMPfam_Rho_Binding
PF00560	HMMPfam_LRR_1
PF07723	HMMPfam_LRR_2
PF08477	HMMPfam_Miro
PF00053	HMMPfam_Laminin_EGF
PF01630	HMMPfam_Glyco_hydro_56
PF00240	HMMPfam_ubiquitin
PF00688	HMMPfam_TGFb_propeptide
PF00019	HMMPfam_TGF_beta
PF06218	HMMPfam_NPR2
PF00341	HMMPfam_PDGF
PF03128	HMMPfam_CXCXC
PF02460	HMMPfam_Patched
PF00627	HMMPfam_UBA
PF02865	HMMPfam_STAT_int
PF01017	HMMPfam_STAT_alpha
PF04388	HMMPfam_Hamartin
PF02145	HMMPfam_Rap_GAP
PF03542	HMMPfam_Tuberin
SSF48678	superfamily_Moesin tail domain
PF00769	HMMPfam_ERM
SSF47031	superfamily_Second domain of FERM
SSF53150	superfamily_DNA repair protein MutS domain II
SSF55271	superfamily_DNA repair protein MutS domain I
SSF63748	superfamily_Tudor/PWWP/MBT
PF01403	HMMPfam_Sema
PF01437	HMMPfam_PSI
PF02178	HMMPfam_AT_hook
PF00856	HMMPfam_SET
PF02008	HMMPfam_zf-CXXC
PF05964	HMMPfam_FYRN
PF05965	HMMPfam_FYRC
SSF57424	superfamily_LDL receptor-like module
PF00622	HMMPfam_SPRY
PF08513	HMMPfam_LisH
PF01753	HMMPfam_zf-MYND
PF05175	HMMPfam_MTS
PF01576	HMMPfam_Myosin_tail_1
PF02736	HMMPfam_Myosin_N
SSF50084	superfamily_Myosin S1 fragment N-terminal domain
SSF57501	superfamily_Cystine-knot cytokines
PF00505	HMMPfam_HMG_box
SSF47095	superfamily_HMG-box
PF01823	HMMPfam_MACPF
PF07648	HMMPfam_Kazal_2
PF06311	HMMPfam_NumbF
SSF109640	superfamily_SSF109640
SSF57667	superfamily_SSF57667
SSF46934	superfamily_UBA-like
SSF101278	superfamily_N-terminal domain of adenylylcyclase associated protein CAP
SSF52047	superfamily_RNI-like
PF06017	HMMPfam_Myosin_TH1
PF02204	HMMPfam_VPS9
SSF109993	superfamily_VPS9 domain (Pfam 02204)
SSF57283	superfamily_PMP inhibitors
PF00503	HMMPfam_G-alpha
PF02437	HMMPfam_Ski_Sno
PF08782	HMMPfam_c-SKI_SMAD_bind
PF02012	HMMPfam_BNR
PF02014	HMMPfam_Reeler
SSF53335	superfamily_S-adenosyl-L-methionine-dependent methyltransferases
SSF50952	superfamily_Soluble quinoprotein glucose dehydrogenase
PF00027	HMMPfam_cNMP_binding
PF00520	HMMPfam_Ion_trans
PF08412	HMMPfam_Ion_trans_N
PF00005	HMMPfam_ABC_tran
PF05466	HMMPfam_BASP1
PF00405	HMMPfam_Transferrin
PF03166	HMMPfam_MH2
PF03165	HMMPfam_MH1
SSF55234	superfamily_Cyanase C-terminal domain
PF06920	HMMPfam_Ded_cyto
PF00515	HMMPfam_TPR_1
PF07719	HMMPfam_TPR_2
PF09311	HMMPfam_Rab5-bind
SSF48452	superfamily_TPR-like
PF00501	HMMPfam_AMP-binding
SSF56801	superfamily_Acetyl-CoA synthetase-like
PF00022	HMMPfam_Actin
SSF53067	superfamily_Actin-like ATPase domain
PF01390	HMMPfam_SEA
PF04503	HMMPfam_SSDP
PF08561	HMMPfam_Ribosomal_L37
PF02864	HMMPfam_STAT_bind
PF00028	HMMPfam_Cadherin
SSF57586	superfamily_TNF receptor-like
PF00561	HMMPfam_Abhydrolase_1
PF06441	HMMPfam_EHN
SSF53474	superfamily_alpha/beta-Hydrolases
SSF52777	superfamily_CoA-dependent acyltransferases
PF03188	HMMPfam_Cytochrom_B561
PF02319	HMMPfam_E2F_TDP
SSF54160	superfamily_Chromo domain-like
PF00079	HMMPfam_Serpin
SSF52087	superfamily_CRAL/TRIO domain
SSF49265	superfamily_FN_III-like
PF09132	HMMPfam_BmKX
SSF52087	superfamily_CRAL_TRIO_C
SSF48350	superfamily_Rho_GAP
SSF48371	superfamily_SSF48371
SSF74924	superfamily_Cap-Gly domain
PF01119	HMMPfam_DNA_mis_repair
PF08676	HMMPfam_MutL_C
SSF54211	superfamily_Ribosomal protein S5 domain 2-like
PF09247	HMMPfam_TBP-binding
SSF47055	superfamily_TAF(II)230 TBP-binding fragment
PF04437	HMMPfam_RINT1_TIP1
SSF101898	superfamily_NHL repeat
PF01049	HMMPfam_Cadherin_C
SSF49313	superfamily_Cadherin-like
PF00435	HMMPfam_Spectrin
SSF57581	superfamily_TB module/8-cys domain
PF01088	HMMPfam_Peptidase_C12
PF00514	HMMPfam_Arm
PF00070	HMMPfam_Pyr_redox
PF00355	HMMPfam_Rieske
PF07992	HMMPfam_Pyr_redox_2
PF00864	HMMPfam_P2X_receptor
PF08919	HMMPfam_F_actin_bind
PF00615	HMMPfam_RGS
SSF54277	superfamily_CAD  PB1 domains
PF00014	HMMPfam_Kunitz_BPTI
SSF56399	superfamily_ADP-ribosylation
SSF50494	superfamily_Trypsin-like serine proteases
SSF57440	superfamily_Kringle-like
SSF57414	superfamily_Hairpin loop containing domain-like
SSF63501	superfamily_Frizzled cysteine-rich domain
SSF82199	superfamily_SET domain
SSF48334	superfamily_DNA repair protein MutS domain III
SSF63825	superfamily_YWTD domain
SSF50969	superfamily_YVTN repeat-like/Quinoprotein amine dehydrogenase
PF02262	HMMPfam_Cbl_N
PF02761	HMMPfam_Cbl_N2
PF02762	HMMPfam_Cbl_N3
PF08758	HMMPfam_Cadherin_pro
PF04836	HMMPfam_IFRD_C
PF05004	HMMPfam_IFRD
PF00583	HMMPfam_Acetyltransf_1
PF01582	HMMPfam_TIR
PF04434	HMMPfam_SWIM
PF00194	HMMPfam_Carb_anhydrase
PF01847	HMMPfam_VHL
PF01834	HMMPfam_XRCC1_N
PF03957	HMMPfam_Jun
PF01153	HMMPfam_Glypican
PF01593	HMMPfam_Amino_oxidase
SSF51905	superfamily_FAD/NAD(P)-binding domain
SSF54373	superfamily_FAD-linked reductases C-terminal domain
SSF47986	superfamily_DEATH domain
SSF48557	superfamily_L-aspartase-like
SSF47895	superfamily_Transducin (alpha subunit) insertion domain
PF00021	HMMPfam_UPAR_LY6
SSF57302	superfamily_Snake toxin-like
PF00067	HMMPfam_p450
SSF48264	superfamily_Cytochrome P450
PF07707	HMMPfam_BACK
PF08005	HMMPfam_PHR
SSF54695	superfamily_POZ domain
PF01569	HMMPfam_PAP2
SSF48317	superfamily_Acid phosphatase/Vanadium-dependent haloperoxidase
SSF51206	superfamily_cAMP-binding domain-like
PF01734	HMMPfam_Patatin
SSF82671	superfamily_SEA domain
PF01534	HMMPfam_Frizzled
PF05001	HMMPfam_RNA_pol_Rpb1_R
PF00193	HMMPfam_Xlink
SSF51445	superfamily_(Trans)glycosidases
SSF55729	superfamily_Acyl-CoA N-acyltransferases (Nat)
SSF48726	superfamily_SSF48726
SSF100895	superfamily_Kazal-type serine protease inhibitors
SSF46966	superfamily_Spectrin repeat
PF00207	HMMPfam_A2M
SSF48239	superfamily_Terpenoid cyclases/Protein prenyltransferases
PF07678	HMMPfam_A2M_comp
SSF56854	superfamily_Bcl-2 inhibitors of programmed cell death
PF02187	HMMPfam_GAS2
PF06001	HMMPfam_DUF902
PF09030	HMMPfam_Creb_binding
SSF47040	superfamily_Kix domain of CBP (creb binding protein)
SSF57933	superfamily_TAZ domain
PF08781	HMMPfam_DP
PF02019	HMMPfam_WIF
PF00554	HMMPfam_RHD
PF00625	HMMPfam_Guanylate_kin
PF08685	HMMPfam_GON
PF00619	HMMPfam_CARD
PF00656	HMMPfam_Peptidase_C14
PF00004	HMMPfam_AAA
PF05889	HMMPfam_SLA_LP_auto_ag
PF02932	HMMPfam_Neur_chan_memb
PF02931	HMMPfam_Neur_chan_LBD
PF00176	HMMPfam_SNF2_N
PF00271	HMMPfam_Helicase_C
PF07529	HMMPfam_HSA
PF07533	HMMPfam_TCH
PF08880	HMMPfam_QLQ
SSF82866	superfamily_Multidrug efflux transporter AcrB transmembrane domain
PF03623	HMMPfam_Focal_AT
SSF68993	superfamily_FAT domain of focal adhesion kinase
PF00318	HMMPfam_Ribosomal_S2
PF02218	HMMPfam_HS1_rep
PF01926	HMMPfam_MMR_HSR1
PF07001	HMMPfam_BAT2_N
PF02205	HMMPfam_WH2
SSF54236	superfamily_SSF54236
SSF49417	superfamily_p53-like transcription factors
SSF48092	superfamily_Transcription factor STAT-4 N-domain
PF01853	HMMPfam_MOZ_SAS
PF01344	HMMPfam_Kelch_1
SSF50965	superfamily_Galactose oxidase central domain
PF07646	HMMPfam_Kelch_2
PF00415	HMMPfam_RCC1
PF00039	HMMPfam_fn1
SSF50985	superfamily_RCC1/BLIP-II
PF02259	HMMPfam_FAT
PF08163	HMMPfam_NUC194
PF05923	HMMPfam_APC_crr
PF05924	HMMPfam_SAMP
PF05937	HMMPfam_EB1_binding
PF05956	HMMPfam_APC_basic
PF05972	HMMPfam_APC_15aa
SSF48366	superfamily_Ras GEF
SSF47220	superfamily_alpha-catenin/vinculin
SSF101447	superfamily_Formin homology 2 domain (FH2 domain)
SSF47655	superfamily_STAT
SSF50044	superfamily_SH3
SSF46626	superfamily_Cytochrome_c
PF05044	HMMPfam_Prox1
SSF46689	superfamily_Homeodomain_like
SSF57552	superfamily_Blood coagulation inhibitor (disintegrin)
PF00188	HMMPfam_SCP
SSF55797	superfamily_PR-1-like
PF01284	HMMPfam_MARVEL
PF04050	HMMPfam_Upf2
SSF54098	superfamily_Prion-like
PF03167	HMMPfam_UDG
SSF52141	superfamily_DNA glycosylase
SSF52200	superfamily_Toll/Interleukin receptor TIR domain
PF02165	HMMPfam_WT1
PF01490	HMMPfam_Aa_trans
PF00040	HMMPfam_fn2
SSF53850	superfamily_Periplasmic binding protein-like II
PF00324	HMMPfam_AA_permease
PF02743	HMMPfam_Cache_1
PF08399	HMMPfam_VWA_N
PF08473	HMMPfam_VGCC_alpha2
PF07392	HMMPfam_P19Arf_N
PF06211	HMMPfam_BAMBI
SSF48340	superfamily_Interferon-induced guanylate-binding protein 1 (GBP1) C-terminal domain
SSF50022	superfamily_ISP domain
SSF55424	superfamily_FAD/NAD-linked reductases dimerisation (C-terminal) domain
PF08205	HMMPfam_C2-set_2
PF01504	HMMPfam_PIP5K
SSF56104	superfamily_SAICAR synthase-like
PF04851	HMMPfam_ResIII
PF00337	HMMPfam_Gal-bind_lectin
PF01603	HMMPfam_B56
PF01294	HMMPfam_Ribosomal_L13e
PF00149	HMMPfam_Metallophos
PF02994	HMMPfam_Transposase_22
PF02257	HMMPfam_RFX_DNA_binding
PF00209	HMMPfam_SNF
PF00956	HMMPfam_NAP
PF05826	HMMPfam_Phospholip_A2_2
PF01434	HMMPfam_Peptidase_M41
PF06480	HMMPfam_FtsH_ext
PF06461	HMMPfam_DUF1086
SSF57603	superfamily_Fibronectin type I module
SSF57095	superfamily_Scorpion toxin-like
SSF46955	superfamily_Putative DNA-binding domain
SSF50249	superfamily_Nucleic acid-binding proteins
SSF81872	superfamily_BRCA2 helical domain
SSF81878	superfamily_BRCA2 tower domain
PF06512	HMMPfam_Na_trans_assoc
PF00335	HMMPfam_Tetraspannin
SSF48652	superfamily_Tetraspanin
PF00003	HMMPfam_7tm_3
PF00020	HMMPfam_TNFR_c6
PF00385	HMMPfam_Chromo
PF06465	HMMPfam_DUF1087
PF02194	HMMPfam_PXA
PF08628	HMMPfam_Nexin_C
PF00880	HMMPfam_Nebulin
SSF46565	superfamily_DnaJ_N
PF02338	HMMPfam_OTU
SSF54001	superfamily_SSF54001
PF00060	HMMPfam_Lig_chan
PF01094	HMMPfam_ANF_receptor
SSF53822	superfamily_Periplasmic binding protein-like I
PF03953	HMMPfam_Tubulin_C
SSF52490	superfamily_SSF52490
SSF55307	superfamily_SSF55307
PF00091	HMMPfam_Tubulin
SSF52490	superfamily_Tubulin nucleotide-binding domain-like
SSF55307	superfamily_Tubulin C-terminal domain-like
PF00782	HMMPfam_DSPc
SSF53756	superfamily_UDP-Glycosyltransferase/glycogen phosphorylase
PF00221	HMMPfam_PAL
SSF48508	superfamily_Nuclear receptor ligand-binding domain
PF04617	HMMPfam_Hox9_act
PF05296	HMMPfam_TAS2R
PF00357	HMMPfam_Integrin_alpha
PF01839	HMMPfam_FG-GAP
PF08441	HMMPfam_Integrin_alpha2
SSF69179	superfamily_Integrin domains
SSF69318	superfamily_Integrin alpha N-terminal domain
PF05587	HMMPfam_Anth_Ig
SSF53300	superfamily_SSF53300
PF00191	HMMPfam_Annexin
SSF47874	superfamily_Annexin
PF01412	HMMPfam_ArfGap
SSF57863	superfamily_Pyk2-associated protein beta ARF-GAP domain
SSF50729	superfamily_SSF50729
PF02138	HMMPfam_Beach
SSF81837	superfamily_BEACH domain
SSF55073	superfamily_Adenylyl and guanylyl cyclase catalytic domain
SSF75399	superfamily_Plakin repeat
PF00681	HMMPfam_Plectin
SSF51069	superfamily_Carbonic anhydrase
PF00572	HMMPfam_Ribosomal_L13
SSF52161	superfamily_Ribosomal_L13
PF04547	HMMPfam_DUF590
SSF48065	superfamily_SSF48065
PF00858	HMMPfam_ASC
PF03247	HMMPfam_Prothymosin
PF03171	HMMPfam_2OG-FeII_Oxy
PF03587	HMMPfam_Nep1
PF00431	HMMPfam_CUB
SSF47719	superfamily_p53 tetramerization domain
SSF51011	superfamily_Glycosyl hydrolase domain
SSF111347	superfamily_Rap/Ran-GAP (Pfam 02145)
PF08161	HMMPfam_NUC173
PF00701	HMMPfam_DHDPS
SSF51569	superfamily_Aldolase
SSF56399	superfamily_SSF56399
PF00571	HMMPfam_CBS
PF01595	HMMPfam_DUF21
SSF54631	superfamily_CBS-domain
PF01991	HMMPfam_vATP-synt_E
PF07531	HMMPfam_TAFH
PF08788	HMMPfam_NHR2
SSF54768	superfamily_dsRNA-binding domain-like
PF00183	HMMPfam_HSP90
PF01299	HMMPfam_Lamp
PF00530	HMMPfam_SRCR
SSF56487	superfamily_SRCR-like
PF00896	HMMPfam_Mtap_PNP
SSF49599	superfamily_TRAF domain-like
PF01501	HMMPfam_Glyco_transf_8
SSF53448	superfamily_Nucleotide-diphospho-sugar transferases
PF00153	HMMPfam_Mito_carr
SSF103506	superfamily_SSF103506
PF00748	HMMPfam_Calpain_inhib
SSF52047	superfamily_SSF52047
PF00225	HMMPfam_Kinesin
PF00110	HMMPfam_wnt
SSF48065	superfamily_DBL homology domain (DH-domain)
SSF47598	superfamily_Ribbon-helix-helix
SSF81301	superfamily_Nucleotidyltransferase
SSF81631	superfamily_PAP/OAS1 substrate-binding domain
PF01909	HMMPfam_NTP_transf_2
PF01805	HMMPfam_Surp
SSF90123	superfamily_Multidrug resistance ABC transporter MsbA N-terminal domain
PF02453	HMMPfam_Reticulon
PF02197	HMMPfam_RIIa
SSF88633	superfamily_Positive stranded ssRNA viruses
SSF46561	superfamily_Ribosomal protein L29 (L29p)
SSF56574	superfamily_Serpins
PF01835	HMMPfam_A2M_N
PF00362	HMMPfam_Integrin_beta
PF07965	HMMPfam_Integrin_B_tail
PF08725	HMMPfam_Integrin_b_cyt
SSF69687	superfamily_Integrin beta tail domain
PF01462	HMMPfam_LRRNT
SSF100895	superfamily_SSF100895
SSF103473	superfamily_SSF103473
SSF68906	superfamily_SAP domain
PF08416	HMMPfam_PTB
PF03137	HMMPfam_OATP
SSF103473	superfamily_MFS general substrate transporter
PF01529	HMMPfam_zf-DHHC
PF03810	HMMPfam_IBN_N
PF05391	HMMPfam_Lsm_interact
PF00365	HMMPfam_PFK
SSF53784	superfamily_Phosphofructokinase
SSF47762	superfamily_PAH2 domain
PF07525	HMMPfam_SOCS_box
SSF52113	superfamily_BRCT domain
PF00664	HMMPfam_ABC_membrane
PF02919	HMMPfam_Topoisom_I_N
SSF57903	superfamily_FYVE_PHD_ZnF
PF08074	HMMPfam_CHDCT2
PF08073	HMMPfam_CHDNT
SSF52540	superfamily_SSF52540
SSF54160	superfamily_SSF54160
PF08771	HMMPfam_Rapamycin_bind
SSF47212	superfamily_FKBP12-rapamycin-binding domain of FKBP-rapamycin-associated protein (FRAP)
PF00147	HMMPfam_Fibrinogen_C
SSF63763	superfamily_SAND domain-like
SSF56366	superfamily_SMAD MH1 domain
SSF57262	superfamily_Leech antihemostatic proteins
PF00038	HMMPfam_Filament
PF07539	HMMPfam_DRIM
PF00610	HMMPfam_DEP
PF00009	HMMPfam_GTP_EFTU
PF03144	HMMPfam_GTP_EFTU_D2
SSF50447	superfamily_Translation proteins
SSF50465	superfamily_EF-Tu/eEF-1alpha/eIF2-gamma C-terminal domain
PF09173	HMMPfam_eIF2_C
SSF48097	superfamily_SSF48097
SSF48479	superfamily_SSF48479
SSF64268	superfamily_SSF64268
PF00157	HMMPfam_Pou
SSF47413	superfamily_Lambda_like_DNA
SSF49854	superfamily_Spermadhesin CUB domain
PF02793	HMMPfam_HRM
PF02180	HMMPfam_BH4
SSF48439	superfamily_Prenyl_trans
SSF74653	superfamily_SSF74653
SSF54928	superfamily_SSF54928
SSF46785	superfamily_SSF46785
PF03607	HMMPfam_DCX
SSF89837	superfamily_SSF89837
PF00348	HMMPfam_polyprenyl_synt
SSF48576	superfamily_Terpenoid synthases
SSF48403	superfamily_ANK
PF03828	HMMPfam_PAP_assoc
SSF81301	superfamily_SSF81301
SSF81631	superfamily_SSF81631
SSF49854	superfamily_CUB
PF07479	HMMPfam_NAD_Gly3P_dh_C
SSF48179	superfamily_6DGDH_C_like
PF01210	HMMPfam_NAD_Gly3P_dh_N
SSF51735	superfamily_SSF51735
SSF47413	superfamily_lambda repressor-like DNA-binding domains
PF00104	HMMPfam_Hormone_recep
PF00105	HMMPfam_zf-C4
PF05631	HMMPfam_DUF791
PF01463	HMMPfam_LRRCT
PF02214	HMMPfam_K_tetra
SSF56436	superfamily_C-type lectin-like
SSF50978	superfamily_WD40_like
SSF52075	superfamily_Outer arm dynein light chain 1
PF00676	HMMPfam_E1_dh
PF02779	HMMPfam_Transket_pyr
PF02466	HMMPfam_Tim17
PF04734	HMMPfam_Ceramidase_alk
SSF50494	superfamily_Pept_Ser_Cys
PF00061	HMMPfam_Lipocalin
SSF50814	superfamily_Calycin
PF00083	HMMPfam_Sugar_tr
PF08443	HMMPfam_RimK
SSF56059	superfamily_Glutathione synthetase ATP-binding domain-like
PF07677	HMMPfam_A2M_recep
SSF49410	superfamily_AM_receptor_bind
SSF69203	superfamily_Nucleoplasmin
SSF56436	superfamily_SSF56436
PF00013	HMMPfam_KH_1
PF01145	HMMPfam_Band_7
SSF50156	superfamily_PDZ
PF01369	HMMPfam_Sec7
SSF48425	superfamily_Sec7
PF07177	HMMPfam_Neuralized
PF01336	HMMPfam_tRNA_anti
PF09170	HMMPfam_DUF1879
SSF50249	superfamily_Nucleic_acid_OB
PF01039	HMMPfam_Carboxyl_trans
PF00364	HMMPfam_Biotin_lipoyl
PF02786	HMMPfam_CPSase_L_D2
PF00289	HMMPfam_CPSase_L_chain
PF02785	HMMPfam_Biotin_carb_C
SSF51230	superfamily_Single hybrid motif
SSF52266	superfamily_Esterase_SGNH_hydro-type
SSF109905	superfamily_Surp module (SWAP domain Pfam 01805)
PF00777	HMMPfam_Glyco_transf_29
PF00801	HMMPfam_PKD
SSF49299	superfamily_PKD domain
PF01496	HMMPfam_V_ATPase_I
PF02902	HMMPfam_Peptidase_C48
PF03020	HMMPfam_LEM
SSF63451	superfamily_LEM domain
PF07502	HMMPfam_MANEC
SSF49373	superfamily_Invasin/intimin cell-adhesion fragments
SSF49410	superfamily_Alpha-macroglobulin receptor domain
PF07703	HMMPfam_A2M_N_2
SSF56487	superfamily_Srcr_receptor
PF01749	HMMPfam_IBB
PF08599	HMMPfam_Nbs1_C
PF00175	HMMPfam_NAD_binding_1
PF00667	HMMPfam_FAD_binding_1
PF02898	HMMPfam_NO_synthase
PF00258	HMMPfam_Flavodoxin_1
SSF52218	superfamily_Flavoproteins
SSF52343	superfamily_Ferredoxin reductase-like C-terminal NADP-linked domain
PF05716	HMMPfam_AKAP_110
PF06741	HMMPfam_Ataxin-2_N
PF07145	HMMPfam_PAM2
SSF50182	superfamily_Sm-like ribonucleoproteins
PF00270	HMMPfam_DEAD
PF01822	HMMPfam_WSC
PF01363	HMMPfam_FYVE
SSF55753	superfamily_Actin depolymerizing proteins
PF00652	HMMPfam_Ricin_B_lectin
SSF50370	superfamily_Ricin B-like lectins
PF04666	HMMPfam_Glyco_transf_54
PF07766	HMMPfam_LETM1
PF06246	HMMPfam_Isy1
SSF57535	superfamily_Complement control module/SCR domain
SSF46596	superfamily_Eukaryotic DNA topoisomerase I dispensable insert domain
SSF56349	superfamily_DNA breaking-rejoining enzymes
PF01028	HMMPfam_Topoisom_I
SSF56741	superfamily_Eukaryotic DNA topoisomerase I N-terminal DNA-binding fragment
PF00082	HMMPfam_Peptidase_S8
PF01483	HMMPfam_P_proprotein
SSF52313	superfamily_Ribosomal protein S2
PF00098	HMMPfam_zf-CCHC
SSF56784	superfamily_HAD-like
SSF81653	superfamily_Calcium ATPase transduction domain A
SSF81660	superfamily_Metal cation-transporting ATPase ATP-binding domain N
SSF81665	superfamily_Calcium ATPase transmembrane domain M
PF00006	HMMPfam_ATP-synt_ab
PF00306	HMMPfam_ATP-synt_ab_C
SSF47917	superfamily_C-terminal domain of alpha and beta subunits of F1 ATP synthase
PF04714	HMMPfam_BCL_N
PF00903	HMMPfam_Glyoxalase
SSF54593	superfamily_Glyoxalase/Bleomycin resistance protein/Dihydroxybiphenyl dioxygenase
PF00288	HMMPfam_GHMP_kinases_N
PF08544	HMMPfam_GHMP_kinases_C
SSF55060	superfamily_GHMP Kinase C-terminal domain
PF07654	HMMPfam_C1-set
PF06464	HMMPfam_DMAP_binding
PF00293	HMMPfam_NUDIX
SSF55811	superfamily_NUDIX_hydrolase
PF00248	HMMPfam_Aldo_ket_red
SSF51430	superfamily_NAD(P)-linked oxidoreductase
SSF51430	superfamily_Aldo/ket_red
PF00996	HMMPfam_GDI
PF02820	HMMPfam_MBT
SSF47769	superfamily_SAM_homology
SSF63748	superfamily_SSF63748
PF09036	HMMPfam_Bcr-Abl_Oligo
SSF69036	superfamily_Bcr-Abl oncoprotein oligomerization domain
SSF63491	superfamily_BAG domain
PF01410	HMMPfam_COLFI
SSF52129	superfamily_Caspase-like
PF00093	HMMPfam_VWC
PF02135	HMMPfam_zf-TAZ
PF00569	HMMPfam_ZZ
PF02172	HMMPfam_KIX
PF00246	HMMPfam_Peptidase_M14
SSF49464	superfamily_Carboxypeptidase regulatory domain
SSF53187	superfamily_Zn-dependent exopeptidases
PF00285	HMMPfam_Citrate_synt
SSF48256	superfamily_Citrate synthase
SSF49785	superfamily_Gal_bind_like
PF00888	HMMPfam_Cullin
SSF74788	superfamily_Cullin repeat
SSF75632	superfamily_Cullin homology domain
PF09282	HMMPfam_Mago-bind
SSF101931	superfamily_Pym (Within the bgcn gene intron protein WIBG) N-terminal domain
SSF49313	superfamily_Cadherin
PF00062	HMMPfam_Lys
SSF53955	superfamily_Lysozyme-like
PF05781	HMMPfam_MRVI1
SSF54791	superfamily_Eukaryotic type KH-domain (KH-domain type I)
PF00025	HMMPfam_Arf
PF01044	HMMPfam_Vinculin
PF02759	HMMPfam_RUN
PF00282	HMMPfam_Pyridoxal_deC
SSF53383	superfamily_PyrdxlP-dep_Trfase_major
PF00791	HMMPfam_ZU5
SSF51246	superfamily_Rudiment single hybrid motif
PF08326	HMMPfam_ACC_central
SSF52096	superfamily_ClpP/crotonase
SSF52440	superfamily_PreATP-grasp domain
PF00441	HMMPfam_Acyl-CoA_dh_1
PF02770	HMMPfam_Acyl-CoA_dh_M
PF02771	HMMPfam_Acyl-CoA_dh_N
SSF47203	superfamily_Acyl-CoA dehydrogenase C-terminal domain-like
SSF56645	superfamily_Acyl-CoA dehydrogenase NM domain-like
PF00211	HMMPfam_Guanylate_cyc
PF06327	HMMPfam_DUF1053
PF05240	HMMPfam_APOBEC_C
PF08210	HMMPfam_APOBEC_N
SSF53927	superfamily_Cytidine deaminase-like
PF00230	HMMPfam_MIP
SSF81338	superfamily_Aquaporin-like
SSF110942	superfamily_HSP90 C-terminal domain (C-terminal part of Pfam 00183)
PF01085	HMMPfam_HH_signal
PF01079	HMMPfam_Hint
SSF55166	superfamily_Hedgehog/DD-peptidase
SSF51294	superfamily_Hedgehog/intein (Hint) domain
PF08824	HMMPfam_Serine_rich
SSF53167	superfamily_Purine and uridine phosphorylases
PF02209	HMMPfam_VHP
SSF47050	superfamily_VHP Villin headpiece domain
PF00626	HMMPfam_Gelsolin
PF04916	HMMPfam_Laminin_A
PF02798	HMMPfam_GST_N
SSF47616	superfamily_GST_C_like
SSF52833	superfamily_Thiordxn-like_fd
PF01553	HMMPfam_Acyltransferase
SSF69593	superfamily_Glycerol-3-phosphate (1)-acyltransferase
PF00100	HMMPfam_Zona_pellucida
PF08347	HMMPfam_CTNNB1_binding
PF07522	HMMPfam_DRMBL
SSF56281	superfamily_Metallo-hydrolase/oxidoreductase
PF00151	HMMPfam_Lipase
SSF49723	superfamily_Lipase/lipooxygenase domain (PLAT/LH2 domain)
PF04791	HMMPfam_LMBR1
PF00481	HMMPfam_PP2C
SSF81606	superfamily_SSF81606
PF07065	HMMPfam_D123
PF00050	HMMPfam_Kazal_1
PF09289	HMMPfam_FOLN
PF02809	HMMPfam_UIM
PF02920	HMMPfam_Integrase_DNA
PF00233	HMMPfam_PDEase_I
SSF109604	superfamily_HD-domain/PDEase-like
SSF53098	superfamily_Ribonuclease H-like
PF00929	HMMPfam_Exonuc_X-T
PF04821	HMMPfam_TIMELESS
PF05029	HMMPfam_TIMELESS_C
PF01105	HMMPfam_EMP24_GP25L
SSF101576	superfamily_Supernatant protein factor (SPF) C-terminal domain
PF01092	HMMPfam_Ribosomal_S6e
PF03568	HMMPfam_Peptidase_C50
SSF82754	superfamily_C-terminal gelsolin-like domain of Sec23/24
PF04184	HMMPfam_ST7
PF07691	HMMPfam_PA14
SSF48439	superfamily_Protein prenylyltransferase
PF00094	HMMPfam_VWD
PF01826	HMMPfam_TIL
SSF57567	superfamily_Serine proterase inhibitors
PF08742	HMMPfam_DUF1787
PF00143	HMMPfam_Interferon
SSF47266	superfamily_4-helical cytokines
PF03522	HMMPfam_KCl_Cotrans_1
PF01846	HMMPfam_FF
SSF81698	superfamily_FF domain
PF02874	HMMPfam_ATP-synt_ab_N
SSF50615	superfamily_N-terminal domain of alpha and beta subunits of F1 ATP synthase
PF08983	HMMPfam_DUF1856
PF01063	HMMPfam_Aminotran_4
SSF56752	superfamily_D-aminoacid aminotransferase-like PLP-dependent enzymes
PF08763	HMMPfam_Ca_chan_IQ
PF03690	HMMPfam_UPF0160
SSF74784	superfamily_Translin
SSF48239	superfamily_Terp_cyc_toroid
PF02141	HMMPfam_DENN
PF03455	HMMPfam_dDENN
PF07524	HMMPfam_Bromo_TP
SSF47113	superfamily_Histone-fold
PF00566	HMMPfam_TBC
SSF47923	superfamily_Ypt/Rab-GAP domain of gyp1p
PF00344	HMMPfam_SecY
SSF103491	superfamily_SecY
SSF109925	superfamily_Lissencephaly-1 protein (Lis-1 PAF-AH alpha) N-terminal domain
SSF47031	superfamily_SSF47031
PF08433	HMMPfam_KTI12
PF01592	HMMPfam_NifU_N
SSF82649	superfamily_SufE/NifU
PF00644	HMMPfam_PARP
PF06473	HMMPfam_FGF-BP1
SSF47923	superfamily_RabGAP_TBC
SSF56801	superfamily_SSF56801
PF00412	HMMPfam_LIM
PF00106	HMMPfam_adh_short
PF03372	HMMPfam_Exo_endo_phos
SSF56219	superfamily_DNase I-like
PF00313	HMMPfam_CSD
PF02825	HMMPfam_WWE
PF04950	HMMPfam_DUF663
PF08142	HMMPfam_AARP2CN
PF05679	HMMPfam_CHGN
PF00086	HMMPfam_Thyroglobulin_1
SSF47473	superfamily_SSF47473
SSF57610	superfamily_SSF57610
PF02158	HMMPfam_Neuregulin
PF02115	HMMPfam_Rho_GDI
PF01129	HMMPfam_ART
PF00690	HMMPfam_Cation_ATPase_N
PF00702	HMMPfam_Hydrolase
PF00689	HMMPfam_Cation_ATPase_C
PF00122	HMMPfam_E1-E2_ATPase
PF01494	HMMPfam_FAD_binding_3
PF04083	HMMPfam_Abhydro_lipase
SSF50998	superfamily_Quinoprotein alcohol dehydrogenase-like
PF03062	HMMPfam_MBOAT
PF09324	HMMPfam_DUF1981
PF00241	HMMPfam_Cofilin_ADF
PF08623	HMMPfam_TIP120
PF01062	HMMPfam_Bestrophin
PF01429	HMMPfam_MBD
PF02791	HMMPfam_DDT
SSF54171	superfamily_DNA-binding domain
PF01740	HMMPfam_STAS
SSF52091	superfamily_Anti-sigma factor antagonist SpoIIaa
PF00916	HMMPfam_Sulfate_transp
PF00535	HMMPfam_Glycos_transf_2
PF00113	HMMPfam_Enolase_C
PF03952	HMMPfam_Enolase_N
SSF51604	superfamily_SSF51604
SSF54826	superfamily_SSF54826
SSF50923	superfamily_Hemopexin-like domain
PF00413	HMMPfam_Peptidase_M10
PF01471	HMMPfam_PG_binding_1
SSF47090	superfamily_PGBD-like
PF00056	HMMPfam_Ldh_1_N
PF02866	HMMPfam_Ldh_1_C
SSF56327	superfamily_LDH C-terminal domain-like
SSF51735	superfamily_NAD(P)-binding Rossmann-fold domains
PF00632	HMMPfam_HECT
SSF56204	superfamily_Hect E3 ligase catalytic domain
SSF49562	superfamily_C2_CaLB
PF05383	HMMPfam_La
PF08409	HMMPfam_DUF1736
PF02121	HMMPfam_IP_trans
PF02862	HMMPfam_DDHD
PF08235	HMMPfam_LNS2
SSF49468	superfamily_VHL
SSF53659	superfamily_Isocitrate/Isopropylmalate dehydrogenase-like
SSF54695	superfamily_SSF54695
PF06312	HMMPfam_Neurexophilin
PF01424	HMMPfam_R3H
SSF82708	superfamily_R3H domain
SSF51045	superfamily_WW domain
SSF56512	superfamily_Nitric oxide (NO) synthase oxygenase domain
SSF63380	superfamily_Riboflavin synthase domain-like
PF00292	HMMPfam_PAX
SSF54999	superfamily_Ribosomal protein S10
PF01556	HMMPfam_DnaJ_C
SSF49493	superfamily_HSP40/DnaJ peptide-binding domain
PF03099	HMMPfam_BPL_LipA_LipB
SSF46934	superfamily_UBA_like
PF00889	HMMPfam_EF_TS
SSF54713	superfamily_EF_TS
SSF47802	superfamily_DNA polymerase beta N-terminal domain-like
SSF101288	superfamily_L27 domain
SSF57196	superfamily_SSF57196
PF02179	HMMPfam_BAG
SSF57850	superfamily_SSF57850
PF01222	HMMPfam_ERG4_ERG24
PF00953	HMMPfam_Glycos_transf_4
PF03637	HMMPfam_Mob1_phocein
SSF101152	superfamily_Mob1/phocein
PF00254	HMMPfam_FKBP_C
SSF54534	superfamily_FKBP-like
PF03456	HMMPfam_uDENN
PF05033	HMMPfam_Pre-SET
SSF56281	superfamily_SSF56281
PF02535	HMMPfam_Zip
SSF57424	superfamily_SSF57424
PF00685	HMMPfam_Sulfotransfer_1
SSF56784	superfamily_SSF56784
SSF46589	superfamily_tRNA-binding arm
PF01979	HMMPfam_Amidohydro_1
SSF51338	superfamily_Metalo_hydrolase
SSF51556	superfamily_SSF51556
SSF56219	superfamily_SSF56219
PF00443	HMMPfam_UCH
SSF56300	superfamily_SSF56300
PF08916	HMMPfam_Phe_ZIP
SSF109805	superfamily_SSF109805
SSF55550	superfamily_SSF55550
PF02318	HMMPfam_RPH3A_effector
SSF69125	superfamily_Nuclear receptor coactivator interlocking domain
PF06010	HMMPfam_DUF906
PF01157	HMMPfam_Ribosomal_L21e
SSF50104	superfamily_Transl_SH3_like
PF05739	HMMPfam_SNARE
PF00804	HMMPfam_Syntaxin
SSF47661	superfamily_t-snare proteins
PF00044	HMMPfam_Gp_dh_N
PF02800	HMMPfam_Gp_dh_C
PF08518	HMMPfam_GIT_SHD
PF00118	HMMPfam_Cpn60_TCP1
SSF48592	superfamily_GroEL equatorial domain-like
SSF52029	superfamily_GroEL apical domain-like
SSF54849	superfamily_GroEL-intermediate domain like
PF00464	HMMPfam_SHMT
SSF53383	superfamily_PLP-dependent transferases
PF00854	HMMPfam_PTR2
PF00534	HMMPfam_Glycos_transf_1
SSF53756	superfamily_SSF53756
SSF48425	superfamily_Sec7 domain
PF01237	HMMPfam_Oxysterol_BP
PF00386	HMMPfam_C1q
SSF49842	superfamily_TNF-like
PF01633	HMMPfam_Choline_kinase
PF04408	HMMPfam_HA2
PF07717	HMMPfam_DUF1605
PF04046	HMMPfam_PSP
SSF57756	superfamily_Retrovirus zinc finger-like domains
PF05395	HMMPfam_DARPP-32
PF07139	HMMPfam_DUF1387
PF01896	HMMPfam_DNA_primase_S
SSF56747	superfamily_Prim-pol domain
PF07970	HMMPfam_DUF1692
PF00351	HMMPfam_Biopterin_H
PF00045	HMMPfam_Hemopexin
PF02376	HMMPfam_CUT
PF02148	HMMPfam_zf-UBP
PF01900	HMMPfam_RNase_P_Rpp14
PF01031	HMMPfam_Dynamin_M
SSF55961	superfamily_Bet v1-like
PF04812	HMMPfam_HNF-1B_C
PF00474	HMMPfam_SSF
PF05625	HMMPfam_PAXNEB
PF00594	HMMPfam_Gla
SSF57630	superfamily_SSF57630
PF05794	HMMPfam_Tcp11
PF05127	HMMPfam_DUF699
PF08351	HMMPfam_DUF1726
PF00198	HMMPfam_2-oxoacid_dh
PF01204	HMMPfam_Trehalase
SSF48208	superfamily_Six-hairpin glycosidases
PF08165	HMMPfam_FerA
PF08150	HMMPfam_FerB
PF08151	HMMPfam_FerI
PF01661	HMMPfam_A1pp
PF00125	HMMPfam_Histone
SSF52949	superfamily_Macro domain-like
SSF55455	superfamily_SRF-like
PF02463	HMMPfam_SMC_N
PF06470	HMMPfam_SMC_hinge
SSF75553	superfamily_Smc hinge domain
PF00108	HMMPfam_Thiolase_N
PF02803	HMMPfam_Thiolase_C
SSF53901	superfamily_Thiolase-like
PF00328	HMMPfam_Acid_phosphat_A
SSF53254	superfamily_Phosphoglycerate mutase-like
PF08726	HMMPfam_efhand_Ca_insen
PF02023	HMMPfam_SCAN
PF02191	HMMPfam_OLF
PF00035	HMMPfam_dsrm
PF09091	HMMPfam_CenpB-DNA-bind
PF01248	HMMPfam_Ribosomal_L7Ae
SSF55315	superfamily_SSF55315
SSF52075	superfamily_SSF52075
SSF50939	superfamily_Sialidases (neuraminidases)
PF03009	HMMPfam_GDPD
PF03982	HMMPfam_DAGAT
PF05875	HMMPfam_aPHC
PF01762	HMMPfam_Galactosyl_T
PF02567	HMMPfam_PhzC-PhzF
SSF54506	superfamily_Diaminopimelate epimerase-like
PF05193	HMMPfam_Peptidase_M16_C
SSF63411	superfamily_LuxS/MPP-like metallohydrolase
PF03671	HMMPfam_UPF0185
PF03160	HMMPfam_Calx-beta
PF04880	HMMPfam_NUDE_C
PF03114	HMMPfam_BAR
PF08773	HMMPfam_CathepsinC_exc
SSF75001	superfamily_Dipeptidyl peptidase I (cathepsin C) exclusion domain
PF02845	HMMPfam_CUE
PF00031	HMMPfam_Cystatin
SSF54403	superfamily_Cystatin/monellin
SSF81585	superfamily_DNA polymerase beta-like second domain
PF00378	HMMPfam_ECH
SSF49899	superfamily_ConA_like_lec_gl
PF01733	HMMPfam_Nucleoside_tran
SSF55811	superfamily_Nudix
PF06668	HMMPfam_ITI_HC_C
PF08487	HMMPfam_VIT
PF03493	HMMPfam_BK_channel_a
PF06665	HMMPfam_DUF1172
SSF110296	superfamily_Oligoxyloglucan reducing end-specific cellobiohydrolase
PF00586	HMMPfam_AIRS
PF02769	HMMPfam_AIRS_C
SSF109736	superfamily_FGAM synthase PurL linker domain
SSF52317	superfamily_Class I glutamine amidotransferase-like
SSF55326	superfamily_PurM N-terminal domain-like
SSF56042	superfamily_PurM C-terminal domain-like
SSF82697	superfamily_PurS-like
SSF48305	superfamily_Class II MHC-associated invariant chain ectoplasmic trimerization domain
SSF49503	superfamily_Cupredoxin
PF07731	HMMPfam_Cu-oxidase_2
PF07732	HMMPfam_Cu-oxidase_3
PF02170	HMMPfam_PAZ
PF02171	HMMPfam_Piwi
SSF101690	superfamily_PAZ domain
PF04889	HMMPfam_Cwf_Cwc_15
SSF51182	superfamily_RmlC-like cupins
PF07847	HMMPfam_DUF1637
PF08081	HMMPfam_RBM1CTR
PF04091	HMMPfam_Sec15
PF06858	HMMPfam_NOG1
PF08155	HMMPfam_NOGCT
PF04130	HMMPfam_Spc97_Spc98
SSF101059	superfamily_B-form DNA mimic Ocr
PF07690	HMMPfam_MFS_1
PF04810	HMMPfam_zf-Sec23_Sec24
PF04811	HMMPfam_Sec23_trunk
PF02225	HMMPfam_PA
PF04253	HMMPfam_TFR_dimer
PF04389	HMMPfam_Peptidase_M28
SSF47672	superfamily_Transferrin receptor ectodomain C-terminal domain
SSF52025	superfamily_Transferrin receptor ectodomain apical domain
PF00007	HMMPfam_Cys_knot
PF00852	HMMPfam_Glyco_transf_10
PF08355	HMMPfam_EF_assoc_1
PF08356	HMMPfam_EF_assoc_2
PF04103	HMMPfam_CD20
PF06333	HMMPfam_TRAP_240kDa
PF04813	HMMPfam_HNF-1A_C
PF04814	HMMPfam_HNF-1_N
SSF100957	superfamily_Dimerization cofactor of HNF-1 alpha
PF00719	HMMPfam_Pyrophosphatase
SSF50324	superfamily_Inorganic pyrophosphatase
PF02891	HMMPfam_zf-MIZ
PF07593	HMMPfam_UnbV_ASPIC
PF06468	HMMPfam_Spond_N
PF09014	HMMPfam_Sushi_2
PF00339	HMMPfam_Arrestin_N
PF02752	HMMPfam_Arrestin_C
PF02990	HMMPfam_EMP70
SSF55781	superfamily_GAF domain-like
PF04326	HMMPfam_AAA_4
PF08154	HMMPfam_NLE
PF04193	HMMPfam_PQ-loop
PF09305	HMMPfam_TACI-CRD2
PF00567	HMMPfam_TUDOR
PF05008	HMMPfam_V-SNARE
PF07546	HMMPfam_EMI
SSF57863	superfamily_ArfGAP
PF04992	HMMPfam_RNA_pol_Rpb1_6
PF04997	HMMPfam_RNA_pol_Rpb1_1
PF04998	HMMPfam_RNA_pol_Rpb1_5
PF05000	HMMPfam_RNA_pol_Rpb1_4
SSF64484	superfamily_beta and beta-prime subunits of DNA dependent RNA-polymerase
PF01161	HMMPfam_PBP
SSF49777	superfamily_PEBP-like
PF01266	HMMPfam_DAO
PF01571	HMMPfam_GCV_T
PF08669	HMMPfam_GCV_T_C
SSF101790	superfamily_Aminomethyltransferase beta-barrel domain
SSF103025	superfamily_Aminomethyltransferase folate-binding domain
PF08318	HMMPfam_COG4
PF08397	HMMPfam_IMD
SSF49354	superfamily_PapD-like
PF01728	HMMPfam_FtsJ
SSF53335	superfamily_SSF53335
PF05822	HMMPfam_UMPH-1
PF03803	HMMPfam_Scramblase
PF07780	HMMPfam_Spb1_C
PF03098	HMMPfam_An_peroxidase
SSF48113	superfamily_Heme-dependent peroxidases
PF03184	HMMPfam_DDE
SSF56204	superfamily_SSF56204
PF01138	HMMPfam_RNase_PH
PF03725	HMMPfam_RNase_PH_C
SSF55666	superfamily_3_ExoRNase
SSF54211	superfamily_SSF54211
SSF50405	superfamily_Actin-crosslinking proteins
PF08598	HMMPfam_Sds3
SSF46458	superfamily_Globin-like
PF03535	HMMPfam_Paxillin
PF07856	HMMPfam_DUF1650
PF07851	HMMPfam_TMPIT
SSF56821	superfamily_Prismane protein-like
PF01842	HMMPfam_ACT
SSF55021	superfamily_SSF55021
SSF56534	superfamily_SSF56534
PF02991	HMMPfam_MAP1_LC3
PF00467	HMMPfam_KOW
SSF51045	superfamily_WW_Rsp5_WWP
PF00350	HMMPfam_Dynamin_N
PF02212	HMMPfam_GED
PF02817	HMMPfam_E3_binding
SSF51230	superfamily_Hybrid_motif
SSF47005	superfamily_SSF47005
SSF52777	superfamily_SSF52777
PF00155	HMMPfam_Aminotran_1_2
PF02655	HMMPfam_ATP-grasp_3
PF08953	HMMPfam_DUF1899
PF08954	HMMPfam_DUF1900
PF00907	HMMPfam_T-box
PF05217	HMMPfam_STOP
SSF63501	superfamily_SSF63501
PF00171	HMMPfam_Aldedh
SSF53720	superfamily_ALDH-like
PF00962	HMMPfam_A_deaminase
SSF51556	superfamily_Metallo-dependent hydrolases
PF02137	HMMPfam_A_deamin
PF01442	HMMPfam_Apolipoprotein
SSF47162	superfamily_Apolipoprotein
PF04430	HMMPfam_DUF498
SSF64076	superfamily_Hypothetical protein MT938 (MTH938)
PF00152	HMMPfam_tRNA-synt_2
SSF55681	superfamily_Class II aaRS and biotin synthetases
PF05593	HMMPfam_RHS_repeat
PF06484	HMMPfam_Ten_N
SSF63570	superfamily_SSF63570
PF06337	HMMPfam_DUF1055
SSF55205	superfamily_EPT/RTPC-like
SSF64356	superfamily_SNARE-like
PF03133	HMMPfam_TTL
PF04157	HMMPfam_EAP30
PF01530	HMMPfam_zf-C2HC
PF00822	HMMPfam_PMP22_Claudin
PF00174	HMMPfam_Oxidored_molyb
SSF56524	superfamily_Sulfite oxidase middle catalytic domain
PF00173	HMMPfam_Cyt-b5
SSF55856	superfamily_Cytochrome b5-like heme/steroid binding domain
PF03404	HMMPfam_Mo-co_dimer
PF07803	HMMPfam_GSG-1
PF06374	HMMPfam_NDUF_C2
PF09279	HMMPfam_efhand_like
PF03719	HMMPfam_Ribosomal_S5_C
PF00333	HMMPfam_Ribosomal_S5
SSF54768	superfamily_SSF54768
PF05270	HMMPfam_AbfB
PF02870	HMMPfam_Methyltransf_1N
PF03736	HMMPfam_EPTP
PF00438	HMMPfam_S-AdoMet_synt_N
PF02772	HMMPfam_S-AdoMet_synt_M
PF02773	HMMPfam_S-AdoMet_synt_C
PF02493	HMMPfam_MORN
SSF82185	superfamily_Histone H3 K4-specific methyltransferase SET7/9 N-terminal domain
PF01035	HMMPfam_DNA_binding_1
SSF46767	superfamily_Methylated DNA-protein cysteine methyltransferase C-terminal domain
SSF53155	superfamily_Methylated DNA-protein cysteine methyltransferase domain
PF08065	HMMPfam_K167R
PF02516	HMMPfam_STT3
PF01416	HMMPfam_PseudoU_synth_1
SSF53822	superfamily_SSF53822
PF04815	HMMPfam_Sec23_helical
PF08033	HMMPfam_Sec23_BS
SSF81811	superfamily_Helical domain of Sec23/24
SSF81995	superfamily_beta-sandwich domain of Sec23/24
SSF82919	superfamily_Zn-finger domain of Sec23/24
PF01590	HMMPfam_GAF
PF03028	HMMPfam_Dynein_heavy
PF05148	HMMPfam_Methyltransf_8
SSF81383	superfamily_SSF81383
PF02906	HMMPfam_Fe_hyd_lg_C
SSF53920	superfamily_Fe-only hydrogenase
PF06060	HMMPfam_Mesothelin
PF06762	HMMPfam_DUF1222
PF05399	HMMPfam_EVI2A
SSF47240	superfamily_Ferritin-like
PF00975	HMMPfam_Thioesterase
SSF111418	superfamily_SSF111418
PF04258	HMMPfam_Peptidase_A22B
SSF52025	superfamily_SSF52025
PF00418	HMMPfam_Tubulin-binding
PF01500	HMMPfam_Keratin_B2
PF07562	HMMPfam_NCD3G
PF00043	HMMPfam_GST_C
SSF47616	superfamily_Glutathione S-transferase (GST) C-terminal domain
PF03941	HMMPfam_INCENP_ARK-bind
PF07941	HMMPfam_K_channel_TID
PF02469	HMMPfam_Fasciclin
PF04863	HMMPfam_EGF_alliinase
SSF50911	superfamily_Mannose 6-phosphate receptor domain
PF04564	HMMPfam_U-box
PF08606	HMMPfam_Prp19
PF05224	HMMPfam_NDT80_PhoG
SSF55856	superfamily_Cyt_B5
PF00487	HMMPfam_FA_desaturase
PF03954	HMMPfam_Lectin_N
PF04952	HMMPfam_AstE_AspA
PF01019	HMMPfam_G_glu_transpept
PF06367	HMMPfam_Drf_FH3
PF06371	HMMPfam_Drf_GBD
PF00892	HMMPfam_DUF6
SSF103481	superfamily_SSF103481
SSF55287	superfamily_RPB5-like RNA polymerase subunit
PF00229	HMMPfam_TNF
SSF49842	superfamily_TNF_like
PF04191	HMMPfam_PEMT
SSF55347	superfamily_Glyceraldehyde-3-phosphate dehydrogenase-like C-terminal domain
PF00884	HMMPfam_Sulfatase
SSF53649	superfamily_Alkaline phosphatase-like
PF05693	HMMPfam_Glycogen_syn
PF01198	HMMPfam_Ribosomal_L31e
SSF54575	superfamily_Ribosomal_L31e
SSF48464	superfamily_ENTH/VHS domain
PF00790	HMMPfam_VHS
SSF48464	superfamily_ENTH_VHS
SSF53955	superfamily_SSF53955
PF00048	HMMPfam_IL8
SSF54117	superfamily_SSF54117
PF01852	HMMPfam_START
SSF55961	superfamily_SSF55961
PF06292	HMMPfam_DUF1041
PF07354	HMMPfam_Sp38
SSF57362	superfamily_BPTI-like
PF00095	HMMPfam_WAP
SSF50630	superfamily_Pept_Aspartic
PF03859	HMMPfam_CG-1
PF03256	HMMPfam_APC10
PF04901	HMMPfam_RAMP
PF01417	HMMPfam_ENTH
SSF51101	superfamily_Mannose-binding lectins
SSF47384	superfamily_Homodimeric domain of signal transducing histidine kinase
SSF52129	superfamily_SSF52129
SSF52313	superfamily_Ribosomal_S2
PF04116	HMMPfam_FA_hydroxylase
PF08914	HMMPfam_Myb_DNA-bind_2
PF00642	HMMPfam_zf-CCCH
SSF90229	superfamily_CCCH zinc finger
PF03328	HMMPfam_HpcH_HpaI
SSF51621	superfamily_Phosphoenolpyruvate/pyruvate domain
PF02434	HMMPfam_Fringe
PF07773	HMMPfam_DUF1619
PF08385	HMMPfam_DHC_N1
PF08393	HMMPfam_DHC_N2
PF01759	HMMPfam_NTR
PF00055	HMMPfam_Laminin_N
SSF56601	superfamily_beta-lactamase/transpeptidase-like
PF01849	HMMPfam_NAC
PF00243	HMMPfam_NGF
SSF51126	superfamily_Pectin lyase-like
SSF101912	superfamily_Sema
PF08337	HMMPfam_Plexin_cytopl
SSF103575	superfamily_SSF103575
SSF51338	superfamily_Composite domain of metallo-dependent hydrolases
SSF47576	superfamily_SSF47576
PF05761	HMMPfam_5_nucleotid
PF00462	HMMPfam_Glutaredoxin
PF02852	HMMPfam_Pyr_redox_dim
PF01167	HMMPfam_Tub
SSF54518	superfamily_Transcriptional factor tubby C-terminal domain
SSF90213	superfamily_NZF domain
PF01608	HMMPfam_I_LWEQ
PF07651	HMMPfam_ANTH
SSF109885	superfamily_I/LWEQ domain (Pfam 01608)
PF04360	HMMPfam_Serglycin
PF08152	HMMPfam_GUCT
SSF103491	superfamily_Preprotein translocase SecY subunit
PF01186	HMMPfam_Lysyl_oxidase
PF02383	HMMPfam_Syja_N
SSF49417	superfamily_P53_like_DNA_bnd
PF00711	HMMPfam_Defensin_beta
PF04706	HMMPfam_Dickkopf_N
PF04676	HMMPfam_CwfJ_C_2
PF04677	HMMPfam_CwfJ_C_1
PF08241	HMMPfam_Methyltransf_11
PF09004	HMMPfam_DUF1891
SSF57630	superfamily_GLA-domain
PF02629	HMMPfam_CoA_binding
PF00549	HMMPfam_Ligase_CoA
SSF52210	superfamily_Succinyl-CoA synthetase domains
PF01756	HMMPfam_ACOX
PF00300	HMMPfam_PGAM
SSF53254	superfamily_SSF53254
PF02893	HMMPfam_GRAM
PF04387	HMMPfam_PTPLA
PF03630	HMMPfam_Fumble
PF01202	HMMPfam_SKI
PF00291	HMMPfam_PALP
SSF53686	superfamily_Tryptophan synthase beta subunit-like PLP-dependent enzymes
SSF55021	superfamily_ACT-like
PF04438	HMMPfam_zf-HIT
PF01764	HMMPfam_Lipase_3
PF00214	HMMPfam_Calc_CGRP_IAPP
SSF50960	superfamily_TolB C-terminal domain
PF02628	HMMPfam_COX15-CtaA
PF00389	HMMPfam_2-Hacid_dh
PF02826	HMMPfam_2-Hacid_dh_C
SSF52283	superfamily_Formate/glycerate dehydrogenase catalytic domain-like
SSF103637	superfamily_SSF103637
SSF55729	superfamily_SSF55729
PF04055	HMMPfam_Radical_SAM
PF06969	HMMPfam_HemN_C
SSF102114	superfamily_SSF102114
PF03190	HMMPfam_DUF255
SSF48113	superfamily_Peroxidase_super
PF00447	HMMPfam_HSF_DNA-bind
PF00354	HMMPfam_Pentaxin
PF02359	HMMPfam_CDC48_N
PF02933	HMMPfam_CDC48_2
SSF50692	superfamily_ADC-like
SSF54585	superfamily_Cdc48 domain 2-like
PF00489	HMMPfam_IL6
PF06407	HMMPfam_BDV_P40
SSF101399	superfamily_P40 nucleoprotein
PF03931	HMMPfam_Skp1_POZ
PF07798	HMMPfam_DUF1640
SSF50475	superfamily_FMN-binding split barrel
PF01243	HMMPfam_Pyridox_oxidase
PF00850	HMMPfam_Hist_deacetyl
SSF52768	superfamily_Arginase/deacetylase
PF00789	HMMPfam_UBX
PF07998	HMMPfam_DUF1695
PF02758	HMMPfam_PAAD_DAPIN
PF05729	HMMPfam_NACHT
PF00887	HMMPfam_ACBP
SSF47027	superfamily_Acyl-CoA binding protein
PF06553	HMMPfam_BNIP3
PF00774	HMMPfam_Ca_channel_B
SSF110221	superfamily_AbfB domain (Pfam 05270)
PF08718	HMMPfam_GLTP
SSF110004	superfamily_Glycolipid_transfer_prot
SSF103506	superfamily_Mitoch_carrier
PF00277	HMMPfam_SAA
SSF56327	superfamily_Lactate_DH/Glyco_hydro_4_C
PF08367	HMMPfam_M16C_assoc
PF05350	HMMPfam_GSK-3_bind
SSF55973	superfamily_S-adenosylmethionine synthetase
PF02351	HMMPfam_GDNF
SSF110035	superfamily_GDNF receptor-like (Pfam 02351)
SSF55120	superfamily_SSF55120
PF00849	HMMPfam_PseudoU_synth_2
SSF56854	superfamily_SSF56854
PF07004	HMMPfam_DUF1309
PF06916	HMMPfam_DUF1279
PF03632	HMMPfam_Glyco_hydro_65m
SSF48208	superfamily_Glyco_trans_6hp
PF04505	HMMPfam_CD225
PF00578	HMMPfam_AhpC-TSA
PF05879	HMMPfam_RHD3
PF05721	HMMPfam_PhyH
PF04376	HMMPfam_ATE_N
PF04377	HMMPfam_ATE_C
PF03826	HMMPfam_OAR
PF02840	HMMPfam_Prp18
PF08799	HMMPfam_PRP4
SSF47938	superfamily_Functional domain of the splicing factor Prp18
PF03148	HMMPfam_Tektin
SSF109854	superfamily_YfiT-like putative metal-dependent hydrolases
PF07162	HMMPfam_B9
PF00575	HMMPfam_S1
PF00550	HMMPfam_PP-binding
SSF47336	superfamily_ACP-like
SSF50129	superfamily_GroES-like
PF00107	HMMPfam_ADH_zinc_N
PF08242	HMMPfam_Methyltransf_12
PF00109	HMMPfam_ketoacyl-synt
PF02801	HMMPfam_Ketoacyl-synt_C
PF00698	HMMPfam_Acyl_transf_1
SSF55048	superfamily_Probable ACP-binding domain of malonyl-CoA ACP transacylase
SSF51971	superfamily_Nucleotide-binding domain
PF03577	HMMPfam_Peptidase_C69
PF01217	HMMPfam_Clat_adaptor_s
PF00918	HMMPfam_Gastrin
PF01055	HMMPfam_Glyco_hydro_31
PF00088	HMMPfam_Trefoil
SSF57492	superfamily_Trefoil
PF01007	HMMPfam_IRK
PF03520	HMMPfam_KCNQ_channel
PF00605	HMMPfam_IRF
SSF47370	superfamily_SSF47370
SSF46988	superfamily_Tubulin chaperone cofactor A
PF00932	HMMPfam_IF_tail
SSF74853	superfamily_Lamin A/C globular tail domain
PF01365	HMMPfam_RYDR_ITPR
PF02815	HMMPfam_MIR
PF08454	HMMPfam_RIH_assoc
PF08709	HMMPfam_Ins145_P3_rec
SSF100909	superfamily_IP3 receptor type 1 binding core domain 2
SSF82109	superfamily_MIR domain (Pfam 02815)
PF06428	HMMPfam_Sec2p
PF04970	HMMPfam_NC
PF01448	HMMPfam_ELM2
SSF52949	superfamily_SSF52949
PF02181	HMMPfam_FH2
PF02189	HMMPfam_ITAM
PF05517	HMMPfam_p25-alpha
PF00808	HMMPfam_CBFD_NFYB_HMF
PF01244	HMMPfam_Peptidase_M19
PF01207	HMMPfam_Dus
SSF51395	superfamily_SSF51395
PF00179	HMMPfam_UQ_con
SSF54495	superfamily_UBC-like
PF00654	HMMPfam_Voltage_CLC
SSF54631	superfamily_SSF54631
SSF81340	superfamily_SSF81340
PF01467	HMMPfam_CTP_transf_2
SSF52374	superfamily_Nucleotidylyl transferase
PF03283	HMMPfam_PAE
PF08240	HMMPfam_ADH_N
PF03127	HMMPfam_GAT
PF02883	HMMPfam_Alpha_adaptinC2
SSF49348	superfamily_Clathrin adaptor appendage domain
SSF89009	superfamily_GAT-like domain
SSF57256	superfamily_Elafin-like
PF04080	HMMPfam_Per1
PF03031	HMMPfam_NIF
PF04906	HMMPfam_Tweety
PF05020	HMMPfam_zf-NPL4
PF05021	HMMPfam_NPL4
PF00728	HMMPfam_Glyco_hydro_20
PF03881	HMMPfam_Fructosamin_kin
SSF50965	superfamily_Gal_oxid_central
SSF81464	superfamily_Cytochrome c oxidase subunit II-like transmembrane region
PF00735	HMMPfam_Septin
SSF51395	superfamily_FMN-linked oxidoreductases
PF04598	HMMPfam_Gasdermin
PF02485	HMMPfam_Branch
PF01121	HMMPfam_CoaE
PF06237	HMMPfam_DUF1011
PF08583	HMMPfam_UPF0287
PF05238	HMMPfam_CHL4
PF00778	HMMPfam_DIX
PF08833	HMMPfam_Axin_b-cat_bind
PF00370	HMMPfam_FGGY_N
PF01256	HMMPfam_Carb_kinase
SSF53613	superfamily_Ribokinase-like
PF09079	HMMPfam_Cdc6_C
SSF82153	superfamily_FAS1 domain
PF01066	HMMPfam_CDP-OH_P_transf
PF08391	HMMPfam_Ly49
PF01923	HMMPfam_Cob_adeno_trans
PF08366	HMMPfam_LLGL
SSF47954	superfamily_Cyclin_like
PF00135	HMMPfam_COesterase
SSF53474	superfamily_SSF53474
PF00131	HMMPfam_Metallothio
SSF57868	superfamily_Metallothionein
SSF57868	superfamily_SSF57868
PF04097	HMMPfam_NIC
PF08403	HMMPfam_AA_permease_N
PF08703	HMMPfam_PLC-beta_C
PF09333	HMMPfam_ATG_C
PF03399	HMMPfam_SAC3_GANP
SSF103506	superfamily_Mitochondrial carrier
PF03147	HMMPfam_FDX-ACB
SSF54991	superfamily_Anticodon-binding domain of PheRS
PF03367	HMMPfam_zf-ZPR1
PF00305	HMMPfam_Lipoxygenase
SSF48484	superfamily_Lipoxigenase
PF03662	HMMPfam_Glyco_hydro_79n
PF06602	HMMPfam_Myotub-related
SSF89837	superfamily_Doublecortin (DC)
PF08344	HMMPfam_TRP_2
PF04136	HMMPfam_Sec34
PF00232	HMMPfam_Glyco_hydro_1
PF04790	HMMPfam_Sarcoglycan_1
SSF49348	superfamily_Clath_adapt
PF05292	HMMPfam_MCD
PF03992	HMMPfam_ABM
SSF54909	superfamily_Dimer_A_B_barrel
SSF47729	superfamily_IHF-like DNA-binding proteins
PF08022	HMMPfam_FAD_binding_8
PF08030	HMMPfam_NAD_binding_6
PF01794	HMMPfam_Ferric_reduct
PF04968	HMMPfam_CHORD
PF04969	HMMPfam_CS
SSF49764	superfamily_HSP20_chap
SSF56712	superfamily_Prokaryotic type I DNA topoisomerase
PF01751	HMMPfam_Toprim
PF06839	HMMPfam_zf-GRF
PF01131	HMMPfam_Topoisom_bac
PF01396	HMMPfam_zf-C4_Topoisom
SSF57783	superfamily_Zinc beta-ribbon
PF04218	HMMPfam_CENP-B_N
PF08332	HMMPfam_CaMKII_AD
SSF54427	superfamily_NTF2-like
PF05010	HMMPfam_TACC
PF02208	HMMPfam_Sorb
PF01423	HMMPfam_LSM
PF03820	HMMPfam_Mtc
PF05832	HMMPfam_DUF846
PF04579	HMMPfam_Keratin_matx
PF00202	HMMPfam_Aminotran_3
PF09006	HMMPfam_Surfac_D-trimer
PF07814	HMMPfam_WAPL
PF02037	HMMPfam_SAP
PF03572	HMMPfam_Peptidase_S41
PF01399	HMMPfam_PCI
PF00704	HMMPfam_Glyco_hydro_18
SSF51445	superfamily_SSF51445
PF00614	HMMPfam_PLDc
SSF56024	superfamily_Phospholipase D/nuclease
PF04433	HMMPfam_SWIRM
SSF56300	superfamily_Metallo-dependent phosphatases
PF07894	HMMPfam_DUF1669
SSF56024	superfamily_SSF56024
PF01554	HMMPfam_MatE
PF00156	HMMPfam_Pribosyltran
PF04221	HMMPfam_RelB
PF04732	HMMPfam_Filament_head
PF03762	HMMPfam_VOMI
SSF51092	superfamily_VOMI
PF00219	HMMPfam_IGFBP
SSF57610	superfamily_Thyroglobulin type-1 domain
SSF57802	superfamily_Rubredoxin-like
PF00177	HMMPfam_Ribosomal_S7
SSF47973	superfamily_Ribosomal protein S7
PF03189	HMMPfam_DUF270
PF06941	HMMPfam_NT5C
PF01073	HMMPfam_3Beta_HSD
SSF52402	superfamily_Adenine nucleotide alpha hydrolases-like
PF08839	HMMPfam_CDT1
PF01636	HMMPfam_APH
PF00588	HMMPfam_SpoU_methylase
PF08032	HMMPfam_SpoU_sub_bind
SSF75217	superfamily_SSF75217
PF04057	HMMPfam_Rep-A_N
PF08646	HMMPfam_Rep_fac-A_C
PF01127	HMMPfam_Sdh_cyt
SSF81343	superfamily_Fumarate reductase respiratory complex transmembrane subunits
PF01866	HMMPfam_Diphthamide_syn
PF07973	HMMPfam_tRNA_SAD
SSF55186	superfamily_SSF55186
PF07292	HMMPfam_NID
PF07334	HMMPfam_IFP_35_N
PF07815	HMMPfam_Abi_HHR
SSF63712	superfamily_Nicotinic receptor ligand binding domain-like
SSF90112	superfamily_Neurotransmitter-gated ion-channel transmembrane pore
PF05881	HMMPfam_CNPase
SSF55144	superfamily_LigT-like
PF01040	HMMPfam_UbiA
PF08167	HMMPfam_NUC201
PF08166	HMMPfam_NUC202
SSF54076	superfamily_RNase A-like
PF00074	HMMPfam_RnaseA
PF04960	HMMPfam_Glutaminase
PF04675	HMMPfam_DNA_ligase_A_N
PF04679	HMMPfam_DNA_ligase_A_C
PF01068	HMMPfam_DNA_ligase_A_M
SSF56091	superfamily_DNA ligase/mRNA capping enzyme catalytic domain
SSF103107	superfamily_Hypothetical protein c14orf129 hspc210
SSF88723	superfamily_PIN domain-like
PF01280	HMMPfam_Ribosomal_L19e
SSF48140	superfamily_Ribosomal protein L19 (L19e)
PF04931	HMMPfam_DNA_pol_phi
PF02176	HMMPfam_zf-TRAF
PF03441	HMMPfam_FAD_binding_7
PF00875	HMMPfam_DNA_photolyase
SSF52425	superfamily_Cryptochrome/photolyase N-terminal domain
SSF48173	superfamily_Cryptochrome/photolyase FAD-binding domain
PF00472	HMMPfam_RF-1
PF03462	HMMPfam_PCRF
SSF75620	superfamily_Release factor (Pfam 00472)
PF04089	HMMPfam_BRICHOS
SSF54736	superfamily_ClpS-like
PF00077	HMMPfam_RVP
SSF50630	superfamily_Acid proteases
PF01965	HMMPfam_DJ-1_PfpI
PF01268	HMMPfam_FTHFS
PF03096	HMMPfam_Ndr
SSF81822	superfamily_RuBisCo LSMT C-terminal substrate-binding domain
PF01883	HMMPfam_DUF59
PF00648	HMMPfam_Peptidase_C2
PF01067	HMMPfam_Calpain_III
SSF49758	superfamily_Calpain large subunit middle domain (domain III)
SSF47323	superfamily_Anticodon-binding domain of a subclass of class I aminoacyl-tRNA synthetases
PF01406	HMMPfam_tRNA-synt_1e
PF00199	HMMPfam_Catalase
SSF56634	superfamily_Heme-dependent catalases
PF06657	HMMPfam_DUF1167
SSF50199	superfamily_Staphylococcal nuclease
PF01401	HMMPfam_Peptidase_M2
PF07728	HMMPfam_AAA_5
PF00103	HMMPfam_Hormone_1
SSF47266	superfamily_4_helix_cytokine
PF00631	HMMPfam_G-gamma
SSF48670	superfamily_Transducin (heterotrimeric G protein) gamma chain
PF08168	HMMPfam_NUC205
SSF53649	superfamily_SSF53649
PF04494	HMMPfam_TFIID_90kDa
PF04379	HMMPfam_DUF525
SSF110069	superfamily_ApaG-like (Pfam 04379)
PF00939	HMMPfam_Na_sulph_symp
SSF54814	superfamily_Prokaryotic type KH domain (KH-domain type II)
PF04280	HMMPfam_Tim44
PF00227	HMMPfam_Proteasome
SSF56235	superfamily_N-terminal nucleophile aminohydrolases (Ntn hydrolases)
PF04849	HMMPfam_HAP1_N
PF00159	HMMPfam_Hormone_3
SSF63712	superfamily_Neur_chan_LBD
SSF90112	superfamily_SSF90112
PF07810	HMMPfam_TMC
PF06079	HMMPfam_Apyrase
SSF101887	superfamily_Apyrase
PF07742	HMMPfam_BTG
PF08375	HMMPfam_Rpn3_C
PF05805	HMMPfam_L6_membrane
PF08164	HMMPfam_TRAUB
PF01694	HMMPfam_Rhomboid
PF07699	HMMPfam_GCC2_GCC3
PF04922	HMMPfam_DIE2_ALG10
PF07002	HMMPfam_Copine
PF04124	HMMPfam_Dor1
PF01327	HMMPfam_Pep_deformylase
SSF56420	superfamily_Fmet_deformylase
PF08772	HMMPfam_NOB1_Zn_bind
PF06419	HMMPfam_COG6
PF02578	HMMPfam_DUF152
PF00623	HMMPfam_RNA_pol_Rpb1_2
PF04983	HMMPfam_RNA_pol_Rpb1_3
PF06959	HMMPfam_RecQ5
SSF90229	superfamily_SSF90229
SSF88659	superfamily_Sigma3 and sigma4 domains of RNA polymerase sigma factors
SSF109993	superfamily_SSF109993
PF01096	HMMPfam_TFIIS_C
PF02150	HMMPfam_RNA_POL_M_15KD
SSF57783	superfamily_SSF57783
PF00085	HMMPfam_Thioredoxin
PF00955	HMMPfam_HCO3_cotransp
PF07565	HMMPfam_Band_3_cyto
SSF55804	superfamily_SSF55804
PF04588	HMMPfam_HIG_1_N
PF00679	HMMPfam_EFG_C
PF03764	HMMPfam_EFG_IV
SSF50447	superfamily_Translat_factor
SSF54980	superfamily_EFG_III_V
PF01233	HMMPfam_NMT
PF02799	HMMPfam_NMT_C
PF00265	HMMPfam_TK
SSF89028	superfamily_Cobalamin adenosyltransferase-like
PF02263	HMMPfam_GBP
PF03359	HMMPfam_GKAP
PF04731	HMMPfam_Caudal_act
PF01432	HMMPfam_Peptidase_M3
PF08266	HMMPfam_Cadherin_2
PF08374	HMMPfam_Protocadherin
PF00650	HMMPfam_CRAL_TRIO
PF04707	HMMPfam_MSF1
PF03765	HMMPfam_CRAL_TRIO_N
SSF46938	superfamily_CRAL/TRIO N-terminal domain
PF00829	HMMPfam_Ribosomal_L21p
SSF69572	superfamily_Activating enzymes of the ubiquitin-like proteins
SSF111347	superfamily_SSF111347
PF05041	HMMPfam_Pecanex_C
PF08939	HMMPfam_DUF1917
SSF48452	superfamily_SSF48452
SSF55785	superfamily_SSF55785
PF03571	HMMPfam_Peptidase_M49
SSF57027	superfamily_Plant inhibitors of proteinases and amylases
PF04139	HMMPfam_Rad9
SSF55979	superfamily_DNA clamp
PF02149	HMMPfam_KA1
SSF103243	superfamily_Kinase associated domain 1 KA1
PF05060	HMMPfam_MGAT2
PF04509	HMMPfam_CheC
SSF103111	superfamily_SSF103111
PF03109	HMMPfam_ABC1
PF08736	HMMPfam_FA
PF04696	HMMPfam_Pinin_SDK_memA
PF04697	HMMPfam_Pinin_SDK_N
SSF52743	superfamily_Pept_S8_S53
PF01758	HMMPfam_SBF
PF00160	HMMPfam_Pro_isomerase
SSF50891	superfamily_CSA_PPIase
SSF53613	superfamily_SSF53613
SSF47323	superfamily_tRNAsyn_1a_bind
SSF52374	superfamily_SSF52374
PF08534	HMMPfam_Redoxin
PF02161	HMMPfam_Prog_receptor
PF02824	HMMPfam_TGS
PF02377	HMMPfam_Dishevelled
SSF51604	superfamily_Enolase C-terminal domain-like
SSF54826	superfamily_Enolase N-terminal domain-like
PF03815	HMMPfam_LCCL
SSF69848	superfamily_LCCL domain
PF06730	HMMPfam_DUF1208
PF04128	HMMPfam_Psf2
PF01812	HMMPfam_5-FTHF_cyc-lig
SSF100950	superfamily_NagB/RpiA/CoA transferase-like
PF02373	HMMPfam_JmjC
SSF88713	superfamily_Glycoside hydrolase/deacetylase
SSF51197	superfamily_Clavaminate synthase-like
SSF55073	superfamily_SSF55073
PF02937	HMMPfam_COX6C
SSF81415	superfamily_SSF81415
SSF54518	superfamily_SSF54518
PF01398	HMMPfam_Mov34
PF00493	HMMPfam_MCM
PF06297	HMMPfam_PET
SSF57716	superfamily_SSF57716
SSF90123	superfamily_SSF90123
PF07469	HMMPfam_DUF1518
PF08815	HMMPfam_Nuc_rec_co-act
PF08832	HMMPfam_SRC-1
PF01433	HMMPfam_Peptidase_M1
SSF53187	superfamily_SSF53187
PF04893	HMMPfam_Yip1
PF01033	HMMPfam_Somatomedin_B
SSF90188	superfamily_Somatomedin B domain
PF07575	HMMPfam_Nucleopor_Nup85
SSF57889	superfamily_SSF57889
PF01151	HMMPfam_ELO
SSF63737	superfamily_Leukotriene A4 hydrolase N-terminal domain
PF02775	HMMPfam_TPP_enzyme_C
PF00205	HMMPfam_TPP_enzyme_M
PF02776	HMMPfam_TPP_enzyme_N
SSF52467	superfamily_SSF52467
SSF52518	superfamily_SSF52518
SSF54897	superfamily_Prot_inh_propept
PF01586	HMMPfam_Basic
SSF47459	superfamily_HLH_basic
PF00749	HMMPfam_tRNA-synt_1c
PF03950	HMMPfam_tRNA-synt_1c_C
PF04557	HMMPfam_tRNA_synt_1c_R2
PF04558	HMMPfam_tRNA_synt_1c_R1
SSF50715	superfamily_Ribosomal_L25rel
SSF50104	superfamily_Translation proteins SH3-like domain
PF08420	HMMPfam_zf-C4_C
PF01166	HMMPfam_TSC22
SSF50129	superfamily_GroES_like
PF02580	HMMPfam_Tyr_Deacylase
PF06608	HMMPfam_DUF1143
PF02244	HMMPfam_Propep_M14
PF05198	HMMPfam_IF3_N
SSF54364	superfamily_Translation initiation factor IF3 N-terminal domain
SSF55200	superfamily_Translation initiation factor IF3 C-terminal domain
PF00812	HMMPfam_Ephrin
SSF49503	superfamily_Cupredoxins
PF00752	HMMPfam_XPG_N
PF00867	HMMPfam_XPG_I
PF00368	HMMPfam_HMG-CoA_red
SSF55035	superfamily_NAD-binding domain of HMG-CoA reductase
SSF56542	superfamily_Substrate-binding domain of HMG-CoA reductase
PF08286	HMMPfam_Spc24
SSF64438	superfamily_CNF1/YfiH-like putative cysteine hydrolases
PF08442	HMMPfam_ATP-grasp_2
PF05241	HMMPfam_EBP
PF07700	HMMPfam_HNOB
PF05063	HMMPfam_MT-A70
PF04000	HMMPfam_Sas10_Utp3
PF02696	HMMPfam_UPF0061
PF02668	HMMPfam_TauD
PF08610	HMMPfam_Pex16
PF06046	HMMPfam_Sec6
PF09325	HMMPfam_Vps5
PF07286	HMMPfam_DUF1445
PF08766	HMMPfam_DEK_C
SSF109715	superfamily_DEK C-terminal domain
PF05210	HMMPfam_Sprouty
PF06487	HMMPfam_SAP18
PF07202	HMMPfam_Tcp10_C
PF08558	HMMPfam_TRF
SSF63600	superfamily_Telomeric repeat binding factor (TRF) dimerisation domain
PF02099	HMMPfam_Josephin
PF01699	HMMPfam_Na_Ca_ex
PF00992	HMMPfam_Troponin
PF00917	HMMPfam_MATH
PF00928	HMMPfam_Adap_comp_sub
SSF49447	superfamily_Second domain of Mu2 adaptin subunit (ap50) of ap2 adaptor
PF01283	HMMPfam_Ribosomal_S26e
PF06462	HMMPfam_Hyd_WA
PF06031	HMMPfam_SERTA
PF08424	HMMPfam_DUF1740
SSF46915	superfamily_Polynucleotide phosphorylase/guanosine pentaphosphate synthase (PNPase/GPSI) domain 3
PF08075	HMMPfam_NOPS
SSF81442	superfamily_Cytochrome c oxidase subunit I-like
PF05168	HMMPfam_HEPN
SSF81593	superfamily_Nucleotidyltransferase substrate binding subunit/domain
PF00710	HMMPfam_Asparaginase
SSF53774	superfamily_Asp/Glutamnse
PF01074	HMMPfam_Glyco_hydro_38
SSF74650	superfamily_Gal_mut_like
SSF88713	superfamily_Glyco_hydro/deAcase_b/a-brl
PF07748	HMMPfam_Glyco_hydro_38C
PF09261	HMMPfam_Alpha-mann_mid
SSF88688	superfamily_SSF88688
PF03604	HMMPfam_DNA_RNApol_7kD
SSF63393	superfamily_RNA polymerase subunits
PF04144	HMMPfam_SCAMP
PF08510	HMMPfam_PIG-P
PF00072	HMMPfam_Response_reg
SSF52172	superfamily_CheY-like
PF08629	HMMPfam_PDE8
SSF49447	superfamily_AP50
SSF64356	superfamily_Longin_like
SSF55486	superfamily_SSF55486
SSF63737	superfamily_SSF63737
PF06951	HMMPfam_PLA2G12
SSF48619	superfamily_Phospholipase A2 PLA2
PF00795	HMMPfam_CN_hydrolase
SSF56317	superfamily_Carbon-nitrogen hydrolase
PF00278	HMMPfam_Orn_DAP_Arg_deC
SSF50923	superfamily_Hemopexin
PF02574	HMMPfam_S-methyl_trans
SSF82282	superfamily_S_methyl_trans
PF08514	HMMPfam_STAG
PF00244	HMMPfam_14-3-3
SSF48445	superfamily_14-3-3 protein
PF01652	HMMPfam_IF4E
SSF55418	superfamily_TIF_eIF_4E
PF00957	HMMPfam_Synaptobrevin
SSF81324	superfamily_SSF81324
SSF53271	superfamily_SSF53271
PF00853	HMMPfam_Runt
PF08504	HMMPfam_RunxI
PF00310	HMMPfam_GATase_2
SSF53448	superfamily_SSF53448
PF00976	HMMPfam_ACTH_domain
PF08035	HMMPfam_Op_neuropeptide
PF08384	HMMPfam_NPP
PF01409	HMMPfam_tRNA-synt_2d
PF08279	HMMPfam_HTH_11
PF00262	HMMPfam_Calreticulin
PF02237	HMMPfam_BPL_C
SSF69500	superfamily_DTD-like (Pfam 02580)
PF00352	HMMPfam_TBP
SSF55945	superfamily_TATA-box binding protein-like
PF06119	HMMPfam_NIDO
PF03381	HMMPfam_CDC50
SSF54495	superfamily_SSF54495
SSF53720	superfamily_SSF53720
PF01034	HMMPfam_Syndecan
SSF50809	superfamily_XRCC4 N-terminal domain
PF03644	HMMPfam_Glyco_hydro_85
SSF57048	superfamily_Gurmarin-like
PF08647	HMMPfam_BRE1
PF04493	HMMPfam_Endonuclease_5
SSF56770	superfamily_Nickel-iron hydrogenase small subunit
PF03528	HMMPfam_Rabaptin
PF07888	HMMPfam_CALCOCO1
SSF47807	superfamily_5' to 3' exonuclease C-terminal subdomain
PF03155	HMMPfam_Alg6_Alg8
PF02375	HMMPfam_JmjN
PF02540	HMMPfam_NAD_synthase
PF02878	HMMPfam_PGM_PMM_I
PF02879	HMMPfam_PGM_PMM_II
SSF53738	superfamily_Phosphoglucomutase first 3 domains
SSF55315	superfamily_L30e-like
SSF81296	superfamily_Ig_E-set
SSF101898	superfamily_SSF101898
PF07701	HMMPfam_HNOBA
SSF111126	superfamily_H-NOX domain
PF04615	HMMPfam_Utp14
PF04990	HMMPfam_RNA_pol_Rpb1_7
PF06472	HMMPfam_ABC_membrane_2
PF04698	HMMPfam_MOBP
SSF111418	superfamily_Hormone receptor domain (HRM Pfam 02793)
PF01150	HMMPfam_GDA1_CD39
PF07975	HMMPfam_C1_4
PF04056	HMMPfam_Ssl1
PF08466	HMMPfam_IRK_N
PF06080	HMMPfam_DUF938
PF06747	HMMPfam_CHCH
PF03901	HMMPfam_Glyco_transf_22
PF02892	HMMPfam_zf-BED
PF05699	HMMPfam_hATC
PF02630	HMMPfam_SCO1-SenC
PF05914	HMMPfam_RIB43A
PF00806	HMMPfam_PUF
PF08338	HMMPfam_DUF1731
PF02136	HMMPfam_NTF2
PF03943	HMMPfam_TAP_C
PF09162	HMMPfam_Tap-RNA_bind
PF09128	HMMPfam_RGS-like
PF01388	HMMPfam_ARID
PF08169	HMMPfam_RBB1NT
SSF46774	superfamily_ARID-like
PF02732	HMMPfam_ERCC4
SSF52980	superfamily_Restriction endonuclease-like
PF00012	HMMPfam_HSP70
SSF100920	superfamily_Heat shock protein 70kD (HSP70) peptide-binding domain
SSF100934	superfamily_Heat shock protein 70kD (HSP70) C-terminal subdomain
PF05386	HMMPfam_TEP1_N
PF05731	HMMPfam_TROVE
PF05652	HMMPfam_DcpS
SSF102860	superfamily_mRNA decapping enzyme DcpS N-terminal domain
SSF54197	superfamily_HIT-like
SSF46966	superfamily_Spectrin
PF00821	HMMPfam_PEPCK
SSF53795	superfamily_PEP carboxykinase-like
SSF68923	superfamily_PEP carboxykinase N-terminal domain
PF02366	HMMPfam_PMT
PF08232	HMMPfam_Striatin
PF03832	HMMPfam_PkinA_anch
PF08557	HMMPfam_Lipid_DES
PF03357	HMMPfam_Snf7
PF05089	HMMPfam_NAGLU
SSF56574	superfamily_Prot_inh_serpin
PF00287	HMMPfam_Na_K-ATPase
PF00403	HMMPfam_HMA
SSF55008	superfamily_HMA heavy metal-associated domain
PF04710	HMMPfam_Pellino
SSF81653	superfamily_SSF81653
SSF81660	superfamily_SSF81660
SSF81665	superfamily_SSF81665
SSF47587	superfamily_Domain of poly(ADP-ribose) polymerase
PF01011	HMMPfam_PQQ
PF06479	HMMPfam_Ribonuc_2-5A
SSF63411	superfamily_Metalloenz_metal-bd
SSF53098	superfamily_RNaseH_fold
PF00029	HMMPfam_Connexin
PF00568	HMMPfam_WH1
SSF50370	superfamily_RicinB_like
PF04824	HMMPfam_Rad21_Rec8
PF04825	HMMPfam_Rad21_Rec8_N
PF00725	HMMPfam_3HCDH
PF02737	HMMPfam_3HCDH_N
SSF48179	superfamily_6-phosphogluconate dehydrogenase C-terminal domain-like
SSF57302	superfamily_SSF57302
SSF48225	superfamily_Seven-hairpin glycosidases
PF02969	HMMPfam_TAF
PF07571	HMMPfam_DUF1546
PF03999	HMMPfam_MAP65_ASE1
SSF47090	superfamily_PGBD_like
PF00557	HMMPfam_Peptidase_M24
SSF55920	superfamily_SSF55920
PF02784	HMMPfam_Orn_Arg_deC_N
SSF50621	superfamily_Alanine racemase C-terminal domain-like
SSF51419	superfamily_PLP-binding barrel
PF05110	HMMPfam_AF-4
PF02202	HMMPfam_Tachykinin
PF01351	HMMPfam_RNase_HII
PF00268	HMMPfam_Ribonuc_red_sm
PF00733	HMMPfam_Asn_synthase
SSF52402	superfamily_SSF52402
SSF56235	superfamily_SSF56235
PF00145	HMMPfam_DNA_methylase
SSF55681	superfamily_SSF55681
PF00570	HMMPfam_HRDC
SSF47819	superfamily_HRDC_like
PF08072	HMMPfam_BDHCT
SSF47986	superfamily_DEATH_like
PF00080	HMMPfam_Sod_Cu
SSF49329	superfamily_CuZn superoxide dismutase-like
PF03451	HMMPfam_HELP
PF00763	HMMPfam_THF_DHG_CYH
PF02882	HMMPfam_THF_DHG_CYH_C
SSF53223	superfamily_Aminoacid dehydrogenase-like N-terminal domain
PF06632	HMMPfam_XRCC4
PF00330	HMMPfam_Aconitase
PF03747	HMMPfam_ADP_ribosyl_GH
SSF101478	superfamily_ADP-ribosylglycohydrolase
PF03097	HMMPfam_BRO1
SSF52799	superfamily_SSF52799
SSF48334	superfamily_SSF48334
SSF53150	superfamily_SSF53150
SSF55271	superfamily_SSF55271
PF09280	HMMPfam_XPC-binding
SSF101238	superfamily_SSF101238
PF02271	HMMPfam_UCR_14kD
SSF81524	superfamily_14 kDa protein of cytochrome bc1 complex (Ubiquinol-cytochrome c reductase)
PF01279	HMMPfam_Parathyroid
SSF49879	superfamily_SMAD_FHA
SSF50242	superfamily_TIMP_like
PF00635	HMMPfam_Motile_Sperm
PF02039	HMMPfam_Adrenomedullin
SSF48366	superfamily_Ras_GEF
PF06189	HMMPfam_5-nucleotidase
PF08625	HMMPfam_Utp13
PF07959	HMMPfam_Fucokinase
PF00180	HMMPfam_Iso_dh
SSF53659	superfamily_SSF53659
SSF47672	superfamily_SSF47672
PF00882	HMMPfam_Zn_dep_PLPC
PF03223	HMMPfam_V-ATPase_C
SSF82199	superfamily_SSF82199
PF07779	HMMPfam_Cas1p
PF02109	HMMPfam_DAD
PF06621	HMMPfam_SIM_C
PF01101	HMMPfam_HMG14_17
SSF54189	superfamily_L23_L15e_core
PF00276	HMMPfam_Ribosomal_L23
SSF50952	superfamily_Quino_gluc_DH
PF02877	HMMPfam_PARP_reg
PF05406	HMMPfam_WGR
PF02731	HMMPfam_SKIP_SNW
SSF47364	superfamily_Domain of the SRP/SRP receptor G-proteins
PF01650	HMMPfam_Peptidase_C13
PF08704	HMMPfam_GCD14
PF02475	HMMPfam_Met_10
PF08952	HMMPfam_DUF1866
PF05997	HMMPfam_Nop52
PF06875	HMMPfam_PRF
PF08423	HMMPfam_Rad51
PF01612	HMMPfam_3_5_exonuc
SSF57552	superfamily_Disintegrin
PF05670	HMMPfam_DUF814
PF00458	HMMPfam_WHEP-TRS
PF00579	HMMPfam_tRNA-synt_1b
SSF47060	superfamily_S15/NS1 RNA-binding domain
PF06583	HMMPfam_Neogenin_C
SSF56534	superfamily_Aromatic aminoacid monoxygenases catalytic and oligomerization domains
PF00264	HMMPfam_Tyrosinase
SSF48056	superfamily_Di-copper centre-containing domain
SSF56672	superfamily_SSF56672
PF04775	HMMPfam_Bile_Hydr_Trans
PF08840	HMMPfam_BAAT_C
PF01271	HMMPfam_Granin
PF00217	HMMPfam_ATP-gua_Ptrans
PF02807	HMMPfam_ATP-gua_PtransN
SSF48034	superfamily_Guanido kinase N-terminal domain
SSF55931	superfamily_Glutamine synthetase/guanido kinase
SSF53167	superfamily_SSF53167
PF05458	HMMPfam_Siva
PF03782	HMMPfam_AMOP
PF05186	HMMPfam_Dpy-30
PF08371	HMMPfam_PLD_envelope
PF00709	HMMPfam_Adenylsucc_synt
SSF55957	superfamily_Phosphoglucomutase C-terminal domain
PF08243	HMMPfam_SPT2
SSF52788	superfamily_Phosphotyrosine protein phosphatases I
PF04722	HMMPfam_Ssu72
PF01379	HMMPfam_Porphobil_deam
PF03900	HMMPfam_Porphobil_deamC
SSF54782	superfamily_Porphobilinogen deaminase (hydroxymethylbilane synthase) C-terminal domain
PF06881	HMMPfam_Elongin_A
PF08711	HMMPfam_TFIIS
SSF47676	superfamily_Conserved domain common to transcription factors TFIIS elongin A CRSP70
PF07474	HMMPfam_G2F
SSF54511	superfamily_GFP-like
PF00675	HMMPfam_Peptidase_M16
PF01747	HMMPfam_ATP-sulfurylase
PF02544	HMMPfam_Steroid_dh
PF03016	HMMPfam_Exostosin
PF09258	HMMPfam_EXTL2
PF05676	HMMPfam_NDUF_B7
SSF51206	superfamily_cNMP_binding
PF01148	HMMPfam_CTP_transf_1
PF07527	HMMPfam_Hairy_orange
PF07528	HMMPfam_DZF
SSF63825	superfamily_SSF63825
PF07496	HMMPfam_zf-CW
PF02412	HMMPfam_TSP_3
PF05735	HMMPfam_TSP_C
SSF103647	superfamily_SSF103647
PF02816	HMMPfam_Alpha_kinase
PF00390	HMMPfam_malic
PF03949	HMMPfam_Malic_M
SSF101576	superfamily_SSF101576
PF02838	HMMPfam_Glyco_hydro_20b
SSF55545	superfamily_beta-N-acetylhexosaminidase-like domain
PF02046	HMMPfam_COX6A
SSF81411	superfamily_COX6A
SSF53850	superfamily_SSF53850
SSF53223	superfamily_SSF53223
SSF50939	superfamily_Sialidase
PF01413	HMMPfam_C4
PF01099	HMMPfam_Uteroglobin
SSF48201	superfamily_Uteroglobin-like
PF01619	HMMPfam_Pro_dh
SSF51730	superfamily_FAD-linked oxidoreductase
PF02230	HMMPfam_Abhydrolase_2
PF05153	HMMPfam_DUF706
SSF53732	superfamily_Aconitase iron-sulfur domain
PF00334	HMMPfam_NDK
SSF54919	superfamily_Nucleoside diphosphate kinases
PF08674	HMMPfam_AChE_tetra
SSF57501	superfamily_SSF57501
PF07985	HMMPfam_SRR1
PF04511	HMMPfam_DER1
PF01602	HMMPfam_Adaptin_N
PF02251	HMMPfam_PA28_alpha
PF02252	HMMPfam_PA28_beta
SSF47216	superfamily_Proteasome activator reg(alpha)
SSF56366	superfamily_MAD_MH1
PF01918	HMMPfam_Alba
SSF82704	superfamily_AlbA-like
PF00758	HMMPfam_EPO_TPO
PF04212	HMMPfam_MIT
PF00026	HMMPfam_Asp
SSF47862	superfamily_Saposin
SSF50911	superfamily_Man6php_recept
PF07915	HMMPfam_PRKCSH
PF03034	HMMPfam_PSS
PF00375	HMMPfam_SDF
PF07574	HMMPfam_SMC_Nse1
PF08746	HMMPfam_zf-RING-like
PF07039	HMMPfam_DUF1325
PF01735	HMMPfam_PLA2_B
SSF48034	superfamily_ATP-gua_Ptrans
SSF55931	superfamily_SSF55931
PF07809	HMMPfam_RTP801_C
PF06567	HMMPfam_Neural_ProG_Cyt
PF06566	HMMPfam_Chon_Sulph_att
PF03061	HMMPfam_4HBT
SSF54637	superfamily_SSF54637
PF01125	HMMPfam_G10
PF01535	HMMPfam_PPR
PF07263	HMMPfam_DMP1
PF01387	HMMPfam_Synuclein
PF02733	HMMPfam_Dak1
PF02734	HMMPfam_Dak2
SSF101473	superfamily_Citrobacter dihydroxyacetone kinase extra ATP-binding domain
SSF82549	superfamily_DAK1/DegV-like
PF06373	HMMPfam_CART
SSF64546	superfamily_Satiety factor CART (cocaine and amphetamine regulated transcript)
PF01992	HMMPfam_vATP-synt_AC39
SSF103486	superfamily_ATPase_V0/A0_c/d
PF00414	HMMPfam_MAP1B_neuraxin
PF02867	HMMPfam_Ribonuc_red_lgC
PF03477	HMMPfam_ATP-cone
SSF48168	superfamily_R1 subunit of ribonucleotide reductase N-terminal domain
PF00317	HMMPfam_Ribonuc_red_lgN
SSF51998	superfamily_PFL-like glycyl radical enzymes
PF04047	HMMPfam_PWP2
SSF50993	superfamily_Prolyl oligopeptidase N-terminal domain
PF02140	HMMPfam_Gal_Lectin
PF06345	HMMPfam_Drf_DAD
SSF101447	superfamily_SSF101447
PF00773	HMMPfam_RNB
PF06427	HMMPfam_UDP-g_GGTase
PF00762	HMMPfam_Ferrochelatase
PF09236	HMMPfam_AHSP
SSF109751	superfamily_Alpha-hemoglobin stabilizing protein AHSP
PF03164	HMMPfam_DUF254
PF04072	HMMPfam_LCM
PF00753	HMMPfam_Lactamase_B
PF06083	HMMPfam_IL17
PF03343	HMMPfam_SART-1
PF00970	HMMPfam_FAD_binding_6
PF02060	HMMPfam_ISK_Channel
PF01108	HMMPfam_Tissue_fac
PF04818	HMMPfam_DUF618
PF08451	HMMPfam_A_deaminase_N
PF05185	HMMPfam_PRMT5
PF00817	HMMPfam_IMS
SSF100879	superfamily_Lesion bypass DNA polymerase (Y-family) little finger domain
SSF48508	superfamily_Str_ncl_receptor
PF05432	HMMPfam_BSP_II
PF00476	HMMPfam_DNA_pol_A
PF09067	HMMPfam_EpoR_lig-bind
PF01285	HMMPfam_TEA
SSF49764	superfamily_HSP20-like chaperones
SSF53901	superfamily_SSF53901
PF00030	HMMPfam_Crystall
SSF49695	superfamily_gamma-Crystallin-like
SSF51069	superfamily_Euk_COanhd
PF03921	HMMPfam_ICAM_N
PF09038	HMMPfam_53-BP1_Tudor
SSF56059	superfamily_SSF56059
PF01591	HMMPfam_6PF2K
PF04729	HMMPfam_Anti-silence
SSF101546	superfamily_Anti-silence
PF00319	HMMPfam_SRF-TF
PF07834	HMMPfam_RanGAP1_C
SSF69099	superfamily_Ran-GTPase activating protein 1 (RanGAP1) C-terminal domain
PF00835	HMMPfam_SNAP-25
PF00297	HMMPfam_Ribosomal_L3
PF03145	HMMPfam_Sina
PF01583	HMMPfam_APS_kinase
SSF88697	superfamily_PUA domain-like
SSF55257	superfamily_RBP11-like subunits of RNA polymerase
PF01193	HMMPfam_RNA_pol_L
PF00129	HMMPfam_MHC_I
SSF75412	superfamily_Hypothetical protein MTH1880
PF05287	HMMPfam_PMG
PF03801	HMMPfam_Ndc80_HEC
PF01663	HMMPfam_Phosphodiest
PF04987	HMMPfam_PigN
PF02714	HMMPfam_DUF221
PF08327	HMMPfam_AHSA1
PF09229	HMMPfam_Aha1_N
SSF53800	superfamily_Chelatase
PF00343	HMMPfam_Phosphorylase
PF01411	HMMPfam_tRNA-synt_2c
PF02272	HMMPfam_DHHA1
SSF101353	superfamily_Putative anticodon-binding domain of alanyl-tRNA synthetase (AlaRS)
PF00260	HMMPfam_Protamine_P1
PF04054	HMMPfam_Not1
PF02057	HMMPfam_Glyco_hydro_59
PF01227	HMMPfam_GTP_cyclohydroI
SSF55620	superfamily_Tetrahydrobiopterin biosynthesis enzymes-like
PF08430	HMMPfam_Fork_head_N
PF06482	HMMPfam_Endostatin
PF01770	HMMPfam_Folate_carrier
PF08180	HMMPfam_BAGE
SSF81606	superfamily_Protein serine/threonine phosphatase 2C catalytic domain
PF04749	HMMPfam_PLAC8
PF01731	HMMPfam_Arylesterase
SSF63829	superfamily_SSF63829
PF02044	HMMPfam_Bombesin
PF03416	HMMPfam_Peptidase_C54
PF07884	HMMPfam_VKOR
SSF63829	superfamily_Calcium-dependent phosphotriesterase
PF07160	HMMPfam_DUF1395
PF01928	HMMPfam_CYTH
PF00303	HMMPfam_Thymidylat_synt
SSF55831	superfamily_Thymidylate synthase/dCMP hydroxymethylase
SSF47391	superfamily_SSF47391
PF01223	HMMPfam_Endonuclease_NS
SSF54060	superfamily_His-Me finger endonucleases
PF04678	HMMPfam_DUF607
PF05160	HMMPfam_DSS1_SEM1
SSF54980	superfamily_EF-G C-terminal domain-like
PF03178	HMMPfam_CPSF_A
PF01729	HMMPfam_QRPTase_C
PF02749	HMMPfam_QRPTase_N
SSF54675	superfamily_Nicotinate/Quinolinate PRTase C-terminal domain-like
PF08312	HMMPfam_cwf21
PF00123	HMMPfam_Hormone_2
PF04988	HMMPfam_AKAP95
PF05831	HMMPfam_GAGE
PF06990	HMMPfam_Gal-3-0_sulfotr
SSF48264	superfamily_Cytochrome_P450
PF00836	HMMPfam_Stathmin
SSF101494	superfamily_Stathmin
PF03561	HMMPfam_Allantoicase
PF04589	HMMPfam_RFX1_trans_act
PF02935	HMMPfam_COX7C
SSF81427	superfamily_Mitochondrial cytochrome c oxidase subunit VIIc (aka VIIIa)
SSF50324	superfamily_Pyrophosphatase
PF08512	HMMPfam_Rtt106
PF08644	HMMPfam_SPT16
SSF55920	superfamily_Creatinase/aminopeptidase
PF00042	HMMPfam_Globin
PF01765	HMMPfam_RRF
SSF55194	superfamily_Ribosome_recyc_fac
PF06765	HMMPfam_HS6ST
PF02836	HMMPfam_Glyco_hydro_2_C
PF02585	HMMPfam_PIG-L
SSF102588	superfamily_SSF102588
SSF50814	superfamily_Lipocalins
PF05572	HMMPfam_Peptidase_M43
PF01061	HMMPfam_ABC2_membrane
PF00133	HMMPfam_tRNA-synt_1
SSF50677	superfamily_ValRS_IleRS_edit
PF08264	HMMPfam_Anticodon_1
PF09334	HMMPfam_tRNA-synt_1g
PF03234	HMMPfam_CDC37_N
PF08564	HMMPfam_CDC37_C
PF08565	HMMPfam_CDC37_M
SSF101391	superfamily_Hsp90 co-chaperone CDC37
SSF47781	superfamily_RuvA domain 2-like
PF08190	HMMPfam_PIH1
PF04042	HMMPfam_DNA_pol_E_B
PF06409	HMMPfam_NPIP
PF05738	HMMPfam_Cna_B
SSF49452	superfamily_Starch-binding domain-like
PF01194	HMMPfam_RNA_pol_N
SSF46924	superfamily_RNA polymerase subunit RPB10
PF01122	HMMPfam_Cobalamin_bind
PF07200	HMMPfam_Mod_r
PF03055	HMMPfam_RPE65
PF08687	HMMPfam_ASD2
SSF56808	superfamily_Ribosomal protein L1
PF02450	HMMPfam_LACT
PF01773	HMMPfam_Nucleos_tra2_N
PF07670	HMMPfam_Gate
PF07662	HMMPfam_Nucleos_tra2_C
PF08282	HMMPfam_Hydrolase_3
PF08969	HMMPfam_DUF1873
SSF52821	superfamily_SSF52821
PF00011	HMMPfam_HSP20
PF00525	HMMPfam_Crystallin
PF00692	HMMPfam_dUTPase
SSF51283	superfamily_SSF51283
PF03253	HMMPfam_UT
SSF46955	superfamily_Putativ_DNA_bind
SSF63763	superfamily_SAND_like
PF00596	HMMPfam_Aldolase_II
SSF53639	superfamily_AraD-like aldolase/epimerase
PF00052	HMMPfam_Laminin_B
PF05349	HMMPfam_GATA-N
PF00459	HMMPfam_Inositol_P
SSF56655	superfamily_Carbohydrate phosphatase
SSF55186	superfamily_ThrRS/AlaRS common domain
PF05039	HMMPfam_Agouti
SSF57055	superfamily_Agouti-related protein
PF03556	HMMPfam_DUF298
PF06701	HMMPfam_MIB_HERC2
SSF82927	superfamily_DM_DNA_bd
PF03474	HMMPfam_DMA
PF08144	HMMPfam_CPL
PF07064	HMMPfam_RIC1
SSF50969	superfamily_Amine_DH_B_like
SSF50182	superfamily_Sm_like_riboprot
PF08444	HMMPfam_Gly_acyl_tr_C
PF06021	HMMPfam_Gly_acyl_tr_N
PF03414	HMMPfam_Glyco_transf_6
PF07572	HMMPfam_BCNT
PF00266	HMMPfam_Aminotran_5
PF03473	HMMPfam_MOSC
PF03476	HMMPfam_MOSC_N
SSF50800	superfamily_PK beta-barrel domain-like
PF01342	HMMPfam_SAND
PF03172	HMMPfam_Sp100
PF00213	HMMPfam_OSCP
SSF47928	superfamily_N-terminal domain of the delta subunit of the F1F0-ATP synthase
PF05024	HMMPfam_Gpi1
PF09068	HMMPfam_efhand_1
PF09069	HMMPfam_efhand_2
PF01158	HMMPfam_Ribosomal_L36e
SSF103642	superfamily_Sec-C motif
PF00478	HMMPfam_IMPDH
SSF51412	superfamily_SSF51412
PF06870	HMMPfam_RNA_pol_I_A49
SSF51126	superfamily_Pectin_lyas_like
PF09049	HMMPfam_SNN_transmemb
PF09050	HMMPfam_SNN_linker
PF09051	HMMPfam_SNN_cytoplasm
PF05327	HMMPfam_RRN3
PF08238	HMMPfam_Sel1
SSF81901	superfamily_HCP-like
SSF82895	superfamily_SSF82895
PF01454	HMMPfam_MAGE
PF00865	HMMPfam_Osteopontin
PF05460	HMMPfam_ORC6_1
PF02492	HMMPfam_cobW
PF07683	HMMPfam_CobW_C
SSF90002	superfamily_SSF90002
SSF81338	superfamily_MIP
PF05672	HMMPfam_MAP7
PF07458	HMMPfam_SPAN-X
SSF46992	superfamily_Ribosomal protein S20
SSF55103	superfamily_FAD-linked oxidases C-terminal domain
SSF90002	superfamily_Hypothetical protein YjiA C-terminal domain
PF08672	HMMPfam_APC2
PF07324	HMMPfam_DGCR6
PF05091	HMMPfam_eIF-3_zeta
PF01000	HMMPfam_RNA_pol_A_bac
SSF56553	superfamily_Insert subdomain of RNA polymerase alpha subunit
SSF52091	superfamily_STAS
PF08688	HMMPfam_ASD1
PF02166	HMMPfam_Androgen_recep
SSF54452	superfamily_MHC_I/II-like_Ag-recog
PF01546	HMMPfam_Peptidase_M20
PF07687	HMMPfam_M20_dimer
SSF55031	superfamily_SSF55031
SSF48163	superfamily_tRNA-synt_bind
PF07534	HMMPfam_TLD
PF03232	HMMPfam_COQ7
PF07989	HMMPfam_Microtub_assoc
PF01428	HMMPfam_zf-AN1
PF01754	HMMPfam_zf-A20
PF06008	HMMPfam_Laminin_I
PF06009	HMMPfam_Laminin_II
PF03388	HMMPfam_Lectin_leg-like
PF01400	HMMPfam_Astacin
PF09045	HMMPfam_L27_2
SSF101288	superfamily_SSF101288
PF02936	HMMPfam_COX4
SSF81406	superfamily_Mitochondrial cytochrome c oxidase subunit IV
PF02423	HMMPfam_OCD_Mu_crystall
PF05038	HMMPfam_Cytochrom_B558a
PF01180	HMMPfam_DHO_dh
SSF55666	superfamily_Ribonuclease PH domain 2-like
PF01597	HMMPfam_GCV_H
PF00839	HMMPfam_Cys_rich_FGFR
PF03635	HMMPfam_Vps35
PF00208	HMMPfam_ELFV_dehydrog
PF02812	HMMPfam_ELFV_dehydrog_N
SSF56496	superfamily_SSF56496
PF02182	HMMPfam_YDG_SRA
SSF48317	superfamily_AcPase_VanPerase
SSF82109	superfamily_SSF82109
PF00485	HMMPfam_PRK
PF05648	HMMPfam_PEX11
SSF51004	superfamily_C-terminal (heme d1) domain of cytochrome cd1-nitrite reductase
PF02017	HMMPfam_CIDE-N
PF00490	HMMPfam_ALAD
PF00406	HMMPfam_ADK
PF00764	HMMPfam_Arginosuc_synth
SSF69864	superfamily_Argininosuccinate synthetase C-terminal domain
PF06012	HMMPfam_DUF908
PF06025	HMMPfam_DUF913
SSF110783	superfamily_Hypothetical protein MTH677
PF04826	HMMPfam_DUF634
PF06220	HMMPfam_zf-U1
PF07850	HMMPfam_Renin_r
PF00999	HMMPfam_Na_H_Exchanger
PF00210	HMMPfam_Ferritin
PF01425	HMMPfam_Amidase
SSF75304	superfamily_Amidase_sig_enz
SSF51695	superfamily_SSF51695
PF03732	HMMPfam_Retrotrans_gag
PF00838	HMMPfam_TCTP
SSF51316	superfamily_Mss4-like
PF04487	HMMPfam_CITED
PF07984	HMMPfam_DUF1693
SSF47928	superfamily_ATPsynt_OSCP
PF04427	HMMPfam_Brix
PF01485	HMMPfam_IBR
PF00075	HMMPfam_RnaseH
SSF55658	superfamily_L9_N_like
PF00429	HMMPfam_TLV_coat
PF06003	HMMPfam_SMN
PF02173	HMMPfam_pKID
PF04145	HMMPfam_Ctr
PF03920	HMMPfam_TLE_N
PF04777	HMMPfam_Evr1_Alr
SSF69000	superfamily_FAD-dependent thiol oxidase
PF00316	HMMPfam_FBPase
PF05682	HMMPfam_PHK_AB
SSF74650	superfamily_Galactose mutarotase-like
SSF88688	superfamily_Families 57/38 glycoside transferase middle domain
PF02926	HMMPfam_THUMP
PF02494	HMMPfam_HYR
SSF90193	superfamily_SSF90193
PF06809	HMMPfam_NPDC1
PF01704	HMMPfam_UDPGP
PF03821	HMMPfam_Mtp
PF04382	HMMPfam_SAB
PF05902	HMMPfam_4_1_CTD
SSF89807	superfamily_Flavin-binding protein dodecin
SSF55874	superfamily_ATP_bd_ATPase
SSF110942	superfamily_SSF110942
SSF103637	superfamily_CCHHC domain
PF06535	HMMPfam_RGM_N
SSF53182	superfamily_SSF53182
PF01843	HMMPfam_DIL
SSF49870	superfamily_Osmotin thaumatin-like protein
PF04727	HMMPfam_ELMO_CED12
PF00684	HMMPfam_DnaJ_CXXCXGXG
SSF51679	superfamily_Bacterial luciferase-like
SSF48695	superfamily_Multiheme cytochromes
PF07035	HMMPfam_Mic1
PF00861	HMMPfam_Ribosomal_L18p
SSF53137	superfamily_SSF53137
PF08774	HMMPfam_VRR_NUC
PF06417	HMMPfam_DUF1077
PF00049	HMMPfam_Insulin
SSF56994	superfamily_Insulin-like
SSF109805	superfamily_Phenylalanine zipper
PF01163	HMMPfam_RIO1
SSF47391	superfamily_Dimerization-anchoring domain of cAMP-dependent type II PK regulatory subunit
PF03131	HMMPfam_bZIP_Maf
PF08383	HMMPfam_Maf_N
PF08333	HMMPfam_DUF1725
PF05640	HMMPfam_DUF798
PF07778	HMMPfam_Mis6
PF04118	HMMPfam_Dopey_N
SSF103481	superfamily_Multidrug resistance efflux transporter EmrE
PF01370	HMMPfam_Epimerase
SSF46946	superfamily_S13-like H2TH domain
PF06831	HMMPfam_H2TH
PF02229	HMMPfam_PC4
SSF54447	superfamily_ssDNA_bind_regul
PF09238	HMMPfam_IL4Ra_N
PF08783	HMMPfam_DWNN
PF02690	HMMPfam_Na_Pi_cotrans
SSF57938	superfamily_DnaJ/Hsp40 cysteine-rich domain
PF01687	HMMPfam_Flavokinase
SSF82114	superfamily_Riboflavin kinase-like
PF06413	HMMPfam_Neugrin
PF06534	HMMPfam_RGM_C
PF07221	HMMPfam_GlcNAc_2-epim
PF02678	HMMPfam_Pirin
PF05726	HMMPfam_Pirin_C
PF02190	HMMPfam_LON
PF07259	HMMPfam_ProSAAS
PF01479	HMMPfam_S4
PF08071	HMMPfam_RS4NT
PF00900	HMMPfam_Ribosomal_S4e
PF01426	HMMPfam_BAH
SSF55174	superfamily_SSF55174
PF09262	HMMPfam_PEX-1N
PF09263	HMMPfam_PEX-2N
PF03143	HMMPfam_GTP_EFTU_D3
SSF50465	superfamily_Elong_init_C
PF01135	HMMPfam_PCMT
PF06309	HMMPfam_Torsin
PF02513	HMMPfam_Spin-Ssty
PF06828	HMMPfam_Fukutin-related
PF01491	HMMPfam_Frataxin_Cyay
SSF55387	superfamily_Frataxin-like
PF00879	HMMPfam_Defensin_propep
PF00323	HMMPfam_Defensin_1
SSF57392	superfamily_Defensin-like
PF01238	HMMPfam_PMI_typeI
PF05577	HMMPfam_Peptidase_S28
SSF101238	superfamily_XPC-binding domain
SSF51306	superfamily_LexA/Signal peptidase
SSF47819	superfamily_HRDC-like
PF04685	HMMPfam_DUF608
PF07052	HMMPfam_Hep_59
SSF81273	superfamily_H-NS histone-like proteins
PF04762	HMMPfam_IKI3
SSF75011	superfamily_3-carboxy-ciscis-mucoante lactonizing enzyme
PF06278	HMMPfam_DUF1032
PF03501	HMMPfam_S10_plectin
PF05324	HMMPfam_Sperm_Ag_HE2
PF03521	HMMPfam_Kv2channel
PF03151	HMMPfam_TPT
SSF49384	superfamily_Carbohydrate-binding domain
PF01139	HMMPfam_UPF0027
SSF49493	superfamily_HSP40_DnaJ_pep
SSF57938	superfamily_SSF57938
PF03446	HMMPfam_NAD_binding_2
PF08389	HMMPfam_Xpo1
PF05967	HMMPfam_DUF887
PF00995	HMMPfam_Sec1
SSF56815	superfamily_Sec1/munc18-like (SM) proteins
SSF47055	superfamily_TAF_II_230
PF08585	HMMPfam_DUF1767
PF01875	HMMPfam_UPF0103
PF03939	HMMPfam_Ribosomal_L23eN
PF09257	HMMPfam_BCMA-Tall_bind
PF04146	HMMPfam_YTH
PF00298	HMMPfam_Ribosomal_L11
SSF55031	superfamily_Bacterial exopeptidase dimerisation domain
PF06818	HMMPfam_Fez1
SSF81382	superfamily_Skp1 dimerisation domain-like
SSF81271	superfamily_TGS-like
PF03946	HMMPfam_Ribosomal_L11_N
SSF46906	superfamily_Ribosomal protein L11 C-terminal domain
SSF54747	superfamily_Ribosomal protein L11 N-terminal domain
SSF63600	superfamily_Telo_rept_bnd_D
PF01205	HMMPfam_UPF0029
PF05773	HMMPfam_RWD
PF01532	HMMPfam_Glyco_hydro_47
PF07297	HMMPfam_DPM2
PF08913	HMMPfam_VBS
PF09141	HMMPfam_Talin_middle
SSF109880	superfamily_A middle domain of Talin 1
PF00638	HMMPfam_Ran_BP1
PF00366	HMMPfam_Ribosomal_S17
PF05477	HMMPfam_SURF2
PF01082	HMMPfam_Cu2_monooxygen
SSF49742	superfamily_PHM/PNGase F
PF03351	HMMPfam_DOMON
PF03712	HMMPfam_Cu2_monoox_C
PF08474	HMMPfam_MYT1
PF09174	HMMPfam_Maf1
PF00826	HMMPfam_Ribosomal_L10e
PF00639	HMMPfam_Rotamase
PF02078	HMMPfam_Synapsin_N
PF02750	HMMPfam_Synapsin_C
PF00891	HMMPfam_Methyltransf_2
PF05361	HMMPfam_PP1_inhibitor
SSF81790	superfamily_Myosin phosphatase inhibitor 17kDa protein CPI-17
PF06682	HMMPfam_DUF1183
PF04133	HMMPfam_Vps55
PF03600	HMMPfam_CitMHS
PF06529	HMMPfam_Vert_IL3-reg_TF
PF02598	HMMPfam_DUF171
PF02545	HMMPfam_Maf
SSF52972	superfamily_Maf/Ham1
PF06209	HMMPfam_COBRA1
PF05019	HMMPfam_Coq4
SSF75217	superfamily_alpha/beta knot
PF05742	HMMPfam_DUF833
PF00694	HMMPfam_Aconitase_C
PF02221	HMMPfam_E1_DerP2_DerF2
SSF52016	superfamily_Aconitase
PF01257	HMMPfam_Complex1_24kDa
SSF47895	superfamily_Transducn_insert
PF00935	HMMPfam_Ribosomal_L44
SSF57829	superfamily_Ribosomal_zn-bd
SSF82615	superfamily_Polo-box domain
PF08699	HMMPfam_DUF1785
SSF103365	superfamily_UPF0027
PF01273	HMMPfam_LBP_BPI_CETP
PF02886	HMMPfam_LBP_BPI_CETP_C
SSF55394	superfamily_SSF55394
PF05461	HMMPfam_ApoL
SSF54762	superfamily_Signal recognition particle alu RNA binding heterodimer SRP9/14
PF06105	HMMPfam_Aph-1
PF03036	HMMPfam_Perilipin
PF05250	HMMPfam_UPF0193
SSF101262	superfamily_Methenyltetrahydrofolate cyclohydrolase-like
SSF68906	superfamily_SSF68906
SSF48592	superfamily_GroEL-ATPase
SSF52029	superfamily_SSF52029
SSF54849	superfamily_SSF54849
PF00163	HMMPfam_Ribosomal_S4
PF03126	HMMPfam_Plus-3
PF08737	HMMPfam_Rgp1
SSF55174	superfamily_Alpha-L RNA-binding motif
PF00573	HMMPfam_Ribosomal_L4
SSF52166	superfamily_Ribosomal protein L4
PF02026	HMMPfam_RyR
PF04757	HMMPfam_Pex2_Pex12
PF02347	HMMPfam_GDC-P
PF07817	HMMPfam_GLE1
PF01465	HMMPfam_GRIP
SSF101283	superfamily_GRIP domain
PF04923	HMMPfam_Ninjurin
SSF50346	superfamily_PRC-barrel domain
PF08683	HMMPfam_DUF1781
PF02188	HMMPfam_GoLoco
PF05180	HMMPfam_zf-DNL
SSF47175	superfamily_Cytochromes
PF07142	HMMPfam_DUF1388
PF04874	HMMPfam_Mak16
PF09047	HMMPfam_MEF2_binding
PF04508	HMMPfam_Pox_A_type_inc
PF07718	HMMPfam_Coatamer_beta_C
SSF47240	superfamily_Ferritin/RR_like
SSF55711	superfamily_Subdomain of clathrin and coatomer appendage domain
PF09066	HMMPfam_B2-adapt-app_C
PF06393	HMMPfam_BID
PF03073	HMMPfam_TspO_MBR
PF01596	HMMPfam_Methyltransf_3
PF00591	HMMPfam_Glycos_transf_3
PF02885	HMMPfam_Glycos_trans_3N
PF01246	HMMPfam_Ribosomal_L24e
SSF101173	superfamily_Docking domain B of the erythromycin polyketide synthase (DEBS)
PF05653	HMMPfam_DUF803
PF01149	HMMPfam_Fapy_DNA_glyco
PF09292	HMMPfam_Neil1-DNA_bind
SSF81624	superfamily_N-terminal domain of MutM-like DNA repair proteins
PF03700	HMMPfam_Sorting_nexin
PF05485	HMMPfam_THAP
PF02666	HMMPfam_PS_Dcarbxylase
PF07962	HMMPfam_Swi3
PF04587	HMMPfam_ADP_PFK_GK
PF08547	HMMPfam_CIA30
PF09336	HMMPfam_Vps4_C
PF08976	HMMPfam_DUF1880
PF01103	HMMPfam_Bac_surface_Ag
PF08697	HMMPfam_TFP11
PF04499	HMMPfam_SAPS
PF00587	HMMPfam_tRNA-synt_2b
PF03129	HMMPfam_HGTP_anticodon
SSF52954	superfamily_Anticodon_bd
SSF47060	superfamily_S15/NS1_bind
SSF48445	superfamily_14-3-3
SSF81340	superfamily_Clc chloride channel
PF06121	HMMPfam_DUF959
PF01419	HMMPfam_Jacalin
PF06221	HMMPfam_zf-C2HC5
PF00868	HMMPfam_Transglut_N
PF01841	HMMPfam_Transglut_core
PF00927	HMMPfam_Transglut_C
SSF49309	superfamily_Transglutaminase two C-terminal domains
PF08778	HMMPfam_HIF-1a_CTAD
PF01275	HMMPfam_Myelin_PLP
PF05760	HMMPfam_IER
PF06814	HMMPfam_Lung_7-TM_R
PF05236	HMMPfam_TAF4
PF05071	HMMPfam_NDUFA12
PF01564	HMMPfam_Spermine_synth
PF04658	HMMPfam_TAFII55_N
SSF57492	superfamily_P_trefoil
PF06423	HMMPfam_GWT1
PF02167	HMMPfam_Cytochrom_C1
SSF46626	superfamily_Cytochrome c
PF00383	HMMPfam_dCMP_cyt_deam_1
SSF53927	superfamily_SSF53927
PF06459	HMMPfam_RR_TM4-6
PF00231	HMMPfam_ATP-synt
SSF52943	superfamily_ATP synthase (F1-ATPase) gamma subunit
PF03770	HMMPfam_IPK
SSF56104	superfamily_SSF56104
PF05118	HMMPfam_Asp_Arg_Hydrox
SSF51197	superfamily_SSF51197
PF02274	HMMPfam_Amidinotransf
SSF55909	superfamily_Pentein
SSF52210	superfamily_CoA_ligase
PF04158	HMMPfam_Sof1
SSF54686	superfamily_SSF54686
PF00410	HMMPfam_Ribosomal_S8
SSF56047	superfamily_Ribosomal_S8
PF01545	HMMPfam_Cation_efflux
PF02036	HMMPfam_SCP2
SSF55718	superfamily_Sterol carrier protein SCP
SSF52954	superfamily_Anticodon-binding domain of Class II aaRS
PF01137	HMMPfam_RTC
SSF55205	superfamily_RNA3apos_cycl/enolpyr_transf_A/B
PF05189	HMMPfam_RTC_insert
PF04049	HMMPfam_APC8
SSF89069	superfamily_N-terminal cytoplasmic domain of anti-sigmaE factor RseA
SSF52467	superfamily_DHS-like NAD/FAD-binding domain
SSF55347	superfamily_SSF55347
SSF52113	superfamily_SSF52113
SSF56091	superfamily_SSF56091
SSF48056	superfamily_Di-copper_centre
SSF46589	superfamily_tRNA_binding_arm
PF02403	HMMPfam_Seryl_tRNA_N
PF00326	HMMPfam_Peptidase_S9
PF03938	HMMPfam_OmpH
PF06346	HMMPfam_Drf_FH1
PF01702	HMMPfam_TGT
SSF51713	superfamily_tRNA_ribo_trans
PF00705	HMMPfam_PCNA_N
PF01470	HMMPfam_Peptidase_C15
SSF48201	superfamily_SSF48201
SSF81406	superfamily_COX4
SSF63451	superfamily_LEM_like
PF01641	HMMPfam_SelR
SSF51316	superfamily_Mss4_like
PF03726	HMMPfam_PNPase
SSF46915	superfamily_3_ExoRNase
SSF54791	superfamily_SSF54791
PF01669	HMMPfam_Myelin_MBP
PF05281	HMMPfam_Secretogranin_V
PF02155	HMMPfam_GCR
PF03045	HMMPfam_DAN
PF04201	HMMPfam_TPD52
PF01980	HMMPfam_UPF0066
SSF57288	superfamily_Midkine
PF01323	HMMPfam_DSBA
PF00714	HMMPfam_IFN-gamma
PF00494	HMMPfam_SQS_PSY
PF07978	HMMPfam_NIPSNAP
PF08491	HMMPfam_SE
PF07159	HMMPfam_DUF1394
PF00068	HMMPfam_Phospholip_A2_1
PF05057	HMMPfam_DUF676
PF04121	HMMPfam_Nup84_Nup100
PF03531	HMMPfam_SSrecog
SSF52218	superfamily_SSF52218
SSF52343	superfamily_SSF52343
SSF56512	superfamily_SSF56512
SSF63380	superfamily_SSF63380
PF06052	HMMPfam_3-HAO
PF02163	HMMPfam_Peptidase_M50
PF03645	HMMPfam_Tctex-1
PF03962	HMMPfam_Mnd1
PF02921	HMMPfam_UCR_TM
SSF56568	superfamily_AlphaBeta_subunt
PF09165	HMMPfam_Ubiq-Cytc-red_N
SSF50022	superfamily_Rieske_dom
SSF81502	superfamily_SSF81502
PF01091	HMMPfam_PTN_MK_C
PF05196	HMMPfam_PTN_MK_N
SSF57288	superfamily_SSF57288
SSF52200	superfamily_TIR
SSF47912	superfamily_WASP_C
PF07831	HMMPfam_PYNP_C
SSF47648	superfamily_Nucleoside phosphorylase/phosphoribosyltransferase N-terminal domain
SSF52418	superfamily_Nucleoside phosphorylase/phosphoribosyltransferase catalytic domain
SSF54680	superfamily_Pyrimidine nucleoside phosphorylase C-terminal domain
PF01821	HMMPfam_ANATO
SSF49309	superfamily_Transglut_C
PF02059	HMMPfam_IL3
PF02275	HMMPfam_CBAH
PF07258	HMMPfam_HCaRG
PF00144	HMMPfam_Beta-lactamase
PF08127	HMMPfam_Propeptide_C1
SSF49464	superfamily_CarboxypepD_reg
SSF100920	superfamily_SSF100920
SSF100934	superfamily_SSF100934
SSF53067	superfamily_SSF53067
PF02822	HMMPfam_Antistasin
PF07792	HMMPfam_DUF1630
PF01218	HMMPfam_Coprogen_oxidas
SSF102886	superfamily_Coprogen_oxidas
PF03226	HMMPfam_Yippee
PF06807	HMMPfam_Clp1
PF00810	HMMPfam_ER_lumen_recept
PF09040	HMMPfam_H-K_ATPase_N
PF02897	HMMPfam_Peptidase_S9_N
SSF50993	superfamily_SSF50993
PF02747	HMMPfam_PCNA_C
PF02297	HMMPfam_COX6B
SSF47694	superfamily_SSF47694
PF09184	HMMPfam_PPP4R2
PF08620	HMMPfam_RPAP1_C
PF08621	HMMPfam_RPAP1_N
PF06662	HMMPfam_C5-epim_C
PF08499	HMMPfam_PDEase_I_N
SSF109604	superfamily_SSF109604
PF02345	HMMPfam_TIL_assoc
PF01112	HMMPfam_Asparaginase_2
PF02404	HMMPfam_SCF
PF04573	HMMPfam_SPC22
PF00184	HMMPfam_Hormone_5
PF00220	HMMPfam_Hormone_4
SSF49606	superfamily_Neurhyp_horm
PF00456	HMMPfam_Transketolase_N
PF02780	HMMPfam_Transketolase_C
PF04111	HMMPfam_APG6
PF01431	HMMPfam_Peptidase_M13
PF05649	HMMPfam_Peptidase_M13_N
PF03876	HMMPfam_RNA_pol_Rpb7_N
SSF56988	superfamily_Anthrax protective antigen
PF05756	HMMPfam_S-antigen
PF01170	HMMPfam_UPF0020
SSF52922	superfamily_Transketo_C_like
PF01265	HMMPfam_Cyto_heme_lyase
PF01234	HMMPfam_NNMT_PNMT_TEMT
PF01027	HMMPfam_UPF0005
PF04819	HMMPfam_DUF716
SSF54277	superfamily_SSF54277
PF09110	HMMPfam_HAND
PF09111	HMMPfam_SLIDE
SSF101224	superfamily_SSF101224
PF06077	HMMPfam_LR8
PF03511	HMMPfam_Fanconi_A
PF00294	HMMPfam_PfkB
PF02390	HMMPfam_Methyltransf_4
SSF51621	superfamily_SSF51621
PF08007	HMMPfam_Cupin_4
PF08336	HMMPfam_P4Ha_N
PF01231	HMMPfam_IDO
PF00797	HMMPfam_Acetyltransf_2
PF09034	HMMPfam_TRADD_N
SSF55044	superfamily_SSF55044
PF09296	HMMPfam_NUDIX-like
PF09297	HMMPfam_zf-NADH-PPase
PF03823	HMMPfam_Neurokinin_B
PF04733	HMMPfam_Coatomer_E
PF02372	HMMPfam_IL15
PF04549	HMMPfam_CD47
PF08204	HMMPfam_V-set_CD47
PF08321	HMMPfam_PPP5
PF04095	HMMPfam_NAPRTase
PF03807	HMMPfam_F420_oxidored
PF02441	HMMPfam_Flavoprotein
PF05279	HMMPfam_Asp-B-Hydro_N
PF05132	HMMPfam_RNA_pol_Rpc4
PF06650	HMMPfam_DUF1162
PF02290	HMMPfam_SRP14
SSF56524	superfamily_Oxidored_molyb
SSF81837	superfamily_SSF81837
SSF52096	superfamily_SSF52096
PF00201	HMMPfam_UDPGT
PF03523	HMMPfam_Macscav_rec
PF01968	HMMPfam_Hydantoinase_A
PF02538	HMMPfam_Hydantoinase_B
PF05378	HMMPfam_Hydant_A_N
PF04063	HMMPfam_DUF383
PF04064	HMMPfam_DUF384
SSF57262	superfamily_Antihaemostatic
PF04840	HMMPfam_Vps16_C
PF04841	HMMPfam_Vps16_N
PF00687	HMMPfam_Ribosomal_L1
SSF56808	superfamily_SSF56808
PF06979	HMMPfam_DUF1301
PF01476	HMMPfam_LysM
SSF69848	superfamily_SSF69848
PF03366	HMMPfam_YEATS
PF03104	HMMPfam_DNA_pol_B_exo
PF00136	HMMPfam_DNA_pol_B
PF07830	HMMPfam_PP2C_C
SSF81601	superfamily_Protein serine/threonine phosphatase 2C C-terminal domain
PF01599	HMMPfam_Ribosomal_S27
SSF81790	superfamily_SSF81790
PF02709	HMMPfam_Galactosyl_T_2
PF02207	HMMPfam_zf-UBR
PF01182	HMMPfam_Glucosamine_iso
SSF100950	superfamily_SSF100950
PF06699	HMMPfam_PIG-F
SSF47917	superfamily_ATPase_a/b_C
SSF50615	superfamily_ATPase_a/b_N
PF05915	HMMPfam_DUF872
PF04621	HMMPfam_ETS_PEA3_N
SSF88798	superfamily_N-terminal heterodimerisation domain of RBP4 (RpoE)
PF07225	HMMPfam_NDUF_B4
PF04209	HMMPfam_HgmA
PF01179	HMMPfam_Cu_amine_oxid
SSF49998	superfamily_Amine oxidase catalytic domain
PF02727	HMMPfam_Cu_amine_oxidN2
PF02728	HMMPfam_Cu_amine_oxidN3
SSF54416	superfamily_Amine oxidase N-terminal region
PF00708	HMMPfam_Acylphosphatase
SSF54975	superfamily_Acylphosphatase
SSF49998	superfamily_CuNH_oxidase
SSF54416	superfamily_SSF54416
PF08145	HMMPfam_BOP1NT
PF00732	HMMPfam_GMC_oxred_N
PF05199	HMMPfam_GMC_oxred_C
PF08357	HMMPfam_SEFIR
PF01176	HMMPfam_eIF-1a
PF06546	HMMPfam_Vert_HS_TF
SSF56655	superfamily_SSF56655
SSF50677	superfamily_ValRS/IleRS/LeuRS editing domain
SSF53686	superfamily_SSF53686
PF08288	HMMPfam_PIGA
SSF101904	superfamily_DNA gyrase A C-terminal domain
PF01725	HMMPfam_Ham1p_like
PF09026	HMMPfam_Cenp-B_dimeris
SSF101160	superfamily_SSF101160
PF03185	HMMPfam_CaKB
SSF50353	superfamily_Cytok_IL1_like
PF05191	HMMPfam_ADK_lid
PF06446	HMMPfam_Hepcidin
SSF57756	superfamily_SSF57756
SSF103486	superfamily_V-type ATP synthase subunit C
PF00715	HMMPfam_IL2
PF07303	HMMPfam_Occludin_ELL
PF04048	HMMPfam_Sec8_exocyst
SSF52507	superfamily_Homo-oligomeric flavin-containing Cys decarboxylases HFCD
PF05154	HMMPfam_TM2
PF06789	HMMPfam_UPF0258
PF05712	HMMPfam_MRG
SSF101353	superfamily_SSF101353
PF08059	HMMPfam_SEP
SSF102848	superfamily_SSF102848
PF00261	HMMPfam_Tropomyosin
PF05686	HMMPfam_DUF821
PF06743	HMMPfam_FAST_1
PF08368	HMMPfam_FAST_2
PF08373	HMMPfam_RAP
SSF54106	superfamily_LysM domain
PF03271	HMMPfam_EB1
PF01109	HMMPfam_GM_CSF
PF00657	HMMPfam_Lipase_GDSL
PF08996	HMMPfam_zf-DNA_Pol
SSF90234	superfamily_SSF90234
PF00931	HMMPfam_NB-ARC
PF01070	HMMPfam_FMN_dh
SSF109775	superfamily_Mannose-6-phosphate receptor binding protein 1 (Tip47) C-terminal domain
SSF48256	superfamily_SSF48256
PF07421	HMMPfam_Pro-NT_NN
PF00930	HMMPfam_DPPIV_N
SSF82171	superfamily_Dipeptidyl peptidase IV/CD26 N-terminal domain
SSF90188	superfamily_SSF90188
PF00446	HMMPfam_GnRH
PF02928	HMMPfam_zf-C5HC2
PF08429	HMMPfam_PLU-1
PF01146	HMMPfam_Caveolin
PF00132	HMMPfam_Hexapep
PF02922	HMMPfam_Isoamylase_N
PF00128	HMMPfam_Alpha-amylase
PF02806	HMMPfam_Alpha-amylase_C
PF00117	HMMPfam_GATase
PF06418	HMMPfam_CTP_synth_N
PF00432	HMMPfam_Prenyltrans
PF02525	HMMPfam_Flavodoxin_2
SSF74853	superfamily_SSF74853
PF08701	HMMPfam_GN3L_Grn1
PF01102	HMMPfam_Glycophorin_A
PF00274	HMMPfam_Glycolytic
SSF51569	superfamily_SSF51569
PF03370	HMMPfam_CBM_21
SSF101790	superfamily_SSF101790
SSF103025	superfamily_SSF103025
PF03909	HMMPfam_BSD
PF08567	HMMPfam_TFIIH_BTF_p62_N
SSF57774	superfamily_SSF57774
PF01263	HMMPfam_Aldose_epim
PF02238	HMMPfam_COX7a
SSF81419	superfamily_Mitochondrial cytochrome c oxidase subunit VIIa
PF04263	HMMPfam_TPK_catalytic
PF04265	HMMPfam_TPK_B1_binding
SSF63862	superfamily_SSF63862
SSF63999	superfamily_SSF63999
PF08648	HMMPfam_DUF1777
SSF101233	superfamily_PWI domain
PF05527	HMMPfam_DUF758
PF01575	HMMPfam_MaoC_dehydratas
SSF55718	superfamily_SCP2
PF03435	HMMPfam_Saccharop_dh
PF01262	HMMPfam_AlaDh_PNT_C
PF05222	HMMPfam_AlaDh_PNT_N
PF06469	HMMPfam_DUF1088
PF00645	HMMPfam_zf-PARP
PF01160	HMMPfam_Opiods_neuropep
PF07706	HMMPfam_TAT_ubiq
PF02765	HMMPfam_Telo_bind
SSF111069	superfamily_Hypothetical protein yfbM
PF01291	HMMPfam_LIF_OSM
SSF52922	superfamily_TK C-terminal domain-like
PF02025	HMMPfam_IL5
SSF52317	superfamily_SSF52317
SSF52266	superfamily_SGNH hydrolase
PF00402	HMMPfam_Calponin
PF02146	HMMPfam_SIR2
PF00483	HMMPfam_NTP_transferase
SSF51161	superfamily_Trimeric LpxA-like enzymes
PF06663	HMMPfam_DUF1170
PF07106	HMMPfam_TBPIP
SSF47162	superfamily_SSF47162
PF03360	HMMPfam_Glyco_transf_43
PF03153	HMMPfam_TFIIA
SSF47396	superfamily_TFIIA_helical
SSF50784	superfamily_TFIIA_betabarrel
PF05051	HMMPfam_COX17
PF02536	HMMPfam_mTERF
PF03870	HMMPfam_RNA_pol_Rpb8
PF04065	HMMPfam_Not3
PF04153	HMMPfam_NOT2_3_5
PF04152	HMMPfam_Mre11_DNA_bind
SSF101233	superfamily_SSF101233
SSF49599	superfamily_Traf_like
SSF53784	superfamily_SSF53784
PF02939	HMMPfam_UcrQ
SSF81508	superfamily_UcrQ
PF01598	HMMPfam_Sterol_desat
PF00329	HMMPfam_Complex1_30kDa
SSF51182	superfamily_RmlC_like_cupin
SSF63887	superfamily_Calret_calnex_P
PF06637	HMMPfam_PV-1
PF01384	HMMPfam_PHO4
SSF50475	superfamily_FMN_binding
PF06886	HMMPfam_TPX2
PF09041	HMMPfam_Aurora-A_bind
SSF57362	superfamily_Prot_inh_Kunz-m
SSF47203	superfamily_AcylCoADH_C_like
SSF56645	superfamily_AcylCoA_dehyd_NM
PF02301	HMMPfam_HORMA
SSF56019	superfamily_The spindle assembly checkpoint protein mad2
PF02312	HMMPfam_CBF_beta
SSF50723	superfamily_CBF_beta
PF05454	HMMPfam_DAG1
SSF111006	superfamily_Dystroglycan domain 2
PF05715	HMMPfam_zf-piccolo
PF00137	HMMPfam_ATP-synt_C
SSF81333	superfamily_F1F0 ATP synthase subunit C
PF06631	HMMPfam_DUF1154
PF02285	HMMPfam_COX8
SSF81431	superfamily_Mitochondrial cytochrome c oxidase subunit VIIIb (aka IX)
SSF52440	superfamily_SSF52440
SSF50998	superfamily_Quin_alc_DH_like
PF03487	HMMPfam_IL13
PF06212	HMMPfam_GRIM-19
PF04114	HMMPfam_Gaa1
PF04506	HMMPfam_Rft-1
SSF55060	superfamily_SSF55060
PF02029	HMMPfam_Caldesmon
PF02913	HMMPfam_FAD-oxidase_C
PF01565	HMMPfam_FAD_binding_4
SSF56176	superfamily_FAD-binding domain
PF03134	HMMPfam_TB2_DP1_HVA22
PF04423	HMMPfam_Rad50_zn_hook
PF06625	HMMPfam_DUF1151
SSF52335	superfamily_SSF52335
PF03227	HMMPfam_GILT
PF05995	HMMPfam_CDO_I
SSF69179	superfamily_SSF69179
SSF47807	superfamily_5_3_exo_C
SSF88723	superfamily_SSF88723
PF05047	HMMPfam_L51_S25_CI-B8
SSF48484	superfamily_Lipoxygenase
SSF49723	superfamily_SSF49723
SSF48670	superfamily_SSF48670
SSF56634	superfamily_Catalase_N
SSF54452	superfamily_MHC antigen-recognition domain
PF08061	HMMPfam_P68HR
PF03735	HMMPfam_ENT
PF07859	HMMPfam_Abhydrolase_3
PF04856	HMMPfam_Securin
PF03089	HMMPfam_RAG2
PF06058	HMMPfam_DCP1
PF03261	HMMPfam_CDK5_activator
PF04190	HMMPfam_DUF410
PF01380	HMMPfam_SIS
SSF53697	superfamily_SIS domain
PF03850	HMMPfam_Tfb4
PF08645	HMMPfam_PNK3P
PF01531	HMMPfam_Glyco_transf_11
PF00037	HMMPfam_Fer4
SSF46548	superfamily_alpha-helical ferredoxin
SSF53671	superfamily_Asp/Orn_carbamoyltranf
PF00185	HMMPfam_OTCace
PF02729	HMMPfam_OTCace_N
PF00988	HMMPfam_CPSase_sm_chain
SSF52021	superfamily_CP_synthsmall
PF02787	HMMPfam_CPSase_L_D3
PF02142	HMMPfam_MGS
SSF48108	superfamily_SSF48108
PF07452	HMMPfam_CHRD
PF01974	HMMPfam_tRNA_int_endo
SSF53032	superfamily_tRNA_int_endo_C
PF01851	HMMPfam_PC_rep
PF00340	HMMPfam_IL1
PF02394	HMMPfam_IL1_propep
SSF57392	superfamily_SSF57392
PF05706	HMMPfam_CDKN3
PF02024	HMMPfam_Leptin
PF08418	HMMPfam_Pol_alpha_B_N
PF01712	HMMPfam_dNK
PF05743	HMMPfam_UEV
PF05347	HMMPfam_Complex1_LYR
PF09202	HMMPfam_Rio2_N
SSF82185	superfamily_SSF82185
PF01124	HMMPfam_MAPEG
SSF57277	superfamily_Granulin repeat
PF01015	HMMPfam_Ribosomal_S3Ae
PF08604	HMMPfam_Nup153
SSF109880	superfamily_Talin_cent
SSF109885	superfamily_SSF109885
SSF52151	superfamily_SSF52151
PF01287	HMMPfam_eIF-5a
SSF81923	superfamily_Double Clp-N motif
PF07744	HMMPfam_SPOC
SSF100939	superfamily_SPOC domain-like
SSF57586	superfamily_SSF57586
PF05195	HMMPfam_AMP_N
SSF53092	superfamily_SSF53092
PF00727	HMMPfam_IL4
SSF51717	superfamily_Dihydropteroate synthetase-like
SSF69318	superfamily_SSF69318
SSF51905	superfamily_SSF51905
SSF54373	superfamily_SSF54373
SSF75553	superfamily_SSF75553
PF01648	HMMPfam_ACPS
SSF56214	superfamily_4-PPT_transf
SSF111126	superfamily_SSF111126
PF09007	HMMPfam_EBP50_C-term
PF05282	HMMPfam_AAR2
PF01086	HMMPfam_Clathrin_lg_ch
PF01023	HMMPfam_S_100
PF05438	HMMPfam_TRH
SSF89028	superfamily_SSF89028
PF00255	HMMPfam_GSHPx
PF03006	HMMPfam_HlyIII
PF04081	HMMPfam_DNA_pol_delta_4
SSF55979	superfamily_SSF55979
PF01408	HMMPfam_GFO_IDH_MocA
PF02894	HMMPfam_GFO_IDH_MocA_C
PF06068	HMMPfam_TIP49
PF01422	HMMPfam_zf-NF-X1
PF01126	HMMPfam_Heme_oxygenase
SSF48613	superfamily_Heme oxygenase-like
SSF82671	superfamily_SSF82671
PF04979	HMMPfam_IPP-2
SSF48678	superfamily_Moesin
PF01048	HMMPfam_PNP_UDP_1
PF05493	HMMPfam_ATP_synt_H
PF08577	HMMPfam_PI31_Prot_Reg
SSF69687	superfamily_Integrin_bsu_tail
PF03567	HMMPfam_Sulfotransfer_2
SSF109715	superfamily_SSF109715
PF00832	HMMPfam_Ribosomal_L39
SSF48662	superfamily_Ribosomal_L39
SSF54919	superfamily_NDK
PF07986	HMMPfam_TBCC
PF05918	HMMPfam_API5
PF01990	HMMPfam_ATP-synt_F
PF03485	HMMPfam_Arg_tRNA_synt_N
PF05746	HMMPfam_DALR_1
PF00750	HMMPfam_tRNA-synt_1d
SSF55190	superfamily_Arginyl-tRNA synthetase (ArgRS) N-terminal 'additional' domain
PF02782	HMMPfam_FGGY_C
PF00215	HMMPfam_OMPdecase
SSF51366	superfamily_Ribulose-phoshate binding barrel
PF02296	HMMPfam_Alpha_adaptin_C
SSF47170	superfamily_Aspartate receptor ligand-binding domain
PF06365	HMMPfam_CD34_antigen
PF01267	HMMPfam_F-actin_cap_A
SSF90096	superfamily_Subunits of heterodimeric actin filament capping protein Capz
PF02888	HMMPfam_CaMBD
PF03530	HMMPfam_SK_channel
SSF81327	superfamily_Small-conductance potassium channel
SSF57770	superfamily_SSF57770
PF07857	HMMPfam_DUF1632
PF01315	HMMPfam_Ald_Xan_dh_C
PF09024	HMMPfam_Sak_Polo
PF01393	HMMPfam_Chromo_shadow
PF01230	HMMPfam_HIT
PF03301	HMMPfam_Trp_dioxygenase
SSF54665	superfamily_Aldxan_dh_hamm
SSF54292	superfamily_Ferredoxin
PF00941	HMMPfam_FAD_binding_5
PF01799	HMMPfam_Fer2_2
SSF47741	superfamily_2Fe-2S_bind
PF03450	HMMPfam_CO_deh_flav_C
SSF55447	superfamily_CO_deh_flav_C
PF02738	HMMPfam_Ald_Xan_dh_C2
SSF56003	superfamily_SSF56003
SSF56176	superfamily_SSF56176
PF04868	HMMPfam_PDE6_gamma
PF05366	HMMPfam_Sarcolipin
SSF48019	superfamily_Pol_clamp_load_C
PF08542	HMMPfam_Rep_fac_C
PF00480	HMMPfam_ROK
PF02350	HMMPfam_Epimerase_2
SSF81811	superfamily_SSF81811
SSF81995	superfamily_SSF81995
SSF82754	superfamily_SSF82754
SSF82919	superfamily_SSF82919
SSF54403	superfamily_SSF54403
PF00466	HMMPfam_Ribosomal_L10
PF00428	HMMPfam_Ribosomal_60s
PF00450	HMMPfam_Peptidase_S10
PF00637	HMMPfam_Clathrin
PF01394	HMMPfam_Clathrin_propel
SSF49299	superfamily_PKD
PF00958	HMMPfam_GMP_synt_C
PF06508	HMMPfam_ExsB
SSF54810	superfamily_GMP synthetase C-terminal dimerisation domain
SSF82708	superfamily_SSF82708
PF00633	HMMPfam_HHH
PF03834	HMMPfam_Rad10
SSF47781	superfamily_RuvA_2_like
PF06888	HMMPfam_Put_Phosphatase
SSF47454	superfamily_Euk_transcr_DNA
PF05046	HMMPfam_Img2
PF08831	HMMPfam_MHCassoc_trimer
PF09307	HMMPfam_MHC2-interact
PF04799	HMMPfam_Fzo_mitofusin
PF00562	HMMPfam_RNA_pol_Rpb2_6
PF00538	HMMPfam_Linker_histone
PF01192	HMMPfam_RNA_pol_Rpb6
SSF63562	superfamily_RPB6/omega subunit-like
PF08320	HMMPfam_PIG-X
PF05328	HMMPfam_CybS
SSF81343	superfamily_SSF81343
PF01902	HMMPfam_ATP_bind_4
PF04051	HMMPfam_TRAPP_Bet3
PF05158	HMMPfam_RNA_pol_Rpc34
PF07228	HMMPfam_SpoIIE
PF05557	HMMPfam_MAD
SSF57924	superfamily_SSF57924
PF06645	HMMPfam_SPC12
PF07966	HMMPfam_A1_Propeptide
PF08324	HMMPfam_PUL
PF09070	HMMPfam_PFU
PF03002	HMMPfam_Somatostatin
PF03874	HMMPfam_RNA_pol_Rpb4
PF08147	HMMPfam_DBP10CT
PF03159	HMMPfam_XRN_N
PF01840	HMMPfam_TCL1_MTCP1
SSF50904	superfamily_Oncogene products
SSF102712	superfamily_SSF102712
PF00730	HMMPfam_HhH-GPD
SSF48150	superfamily_DNA-glycosylase
PF01869	HMMPfam_BcrAD_BadFG
PF00189	HMMPfam_Ribosomal_S3_C
SSF54821	superfamily_Ribosomal protein S3 C-terminal domain
PF07650	HMMPfam_KH_2
PF09029	HMMPfam_Preseq_ALAS
SSF90209	superfamily_SSF90209
PF00682	HMMPfam_HMGL-like
PF02436	HMMPfam_PYC_OADA
SSF51246	superfamily_Rudmnt_hyb_motif
SSF89000	superfamily_SSF89000
SSF49777	superfamily_SSF49777
PF01335	HMMPfam_DED
SSF48619	superfamily_PhospholipaseA2
PF08797	HMMPfam_HIRAN
PF01209	HMMPfam_Ubie_methyltran
PF04692	HMMPfam_PDGF_N
PF06026	HMMPfam_Rib_5-P_isom_A
SSF75445	superfamily_D-ribose-5-phosphate isomerase (RpiA) lid domain
PF00394	HMMPfam_Cu-oxidase
PF09268	HMMPfam_Clathrin-link
SSF50989	superfamily_Clathrin heavy-chain terminal domain
PF03651	HMMPfam_Psf1
PF07738	HMMPfam_Sad1_UNC
PF00322	HMMPfam_Endothelin
PF06752	HMMPfam_E_Pc_C
PF01873	HMMPfam_eIF-5_eIF-2B
SSF100966	superfamily_Translation initiation factor 2 beta aIF2beta N-terminal domain
SSF75689	superfamily_Zinc-binding domain of translation initiation factor 2 beta
PF05873	HMMPfam_Mt_ATP-synt_D
PF08925	HMMPfam_DUF1907
PF05221	HMMPfam_AdoHcyase
PF00670	HMMPfam_AdoHcyase_NAD
SSF52283	superfamily_SSF52283
PF01215	HMMPfam_COX5B
SSF53633	superfamily_Aa_kinase
PF04768	HMMPfam_DUF619
PF03071	HMMPfam_GNT-I
PF07576	HMMPfam_BRAP2
PF04321	HMMPfam_RmlD_sub_bind
PF03259	HMMPfam_Robl_LC7
PF09310	HMMPfam_PD-C2-AF1
PF02245	HMMPfam_Pur_DNA_glyco
SSF50486	superfamily_FMT C-terminal domain-like
PF06728	HMMPfam_PIG-U
PF04560	HMMPfam_RNA_pol_Rpb2_7
PF04561	HMMPfam_RNA_pol_Rpb2_2
PF04563	HMMPfam_RNA_pol_Rpb2_1
PF04565	HMMPfam_RNA_pol_Rpb2_3
PF04567	HMMPfam_RNA_pol_Rpb2_5
PF06883	HMMPfam_RNA_pol_Rpa2_4
PF09058	HMMPfam_L27_1
SSF54184	superfamily_Penicillin-binding protein 2x (pbp-2x) c-terminal domain
PF06821	HMMPfam_DUF1234
PF07649	HMMPfam_C1_3
PF03199	HMMPfam_GSH_synthase
PF03917	HMMPfam_GSH_synth_ATP
SSF51971	superfamily_SSF51971
SSF55194	superfamily_Ribosome recycling factor RRF
PF01469	HMMPfam_Pentapeptide_2
PF05104	HMMPfam_Rib_recp_KP_reg
PF08755	HMMPfam_YccV-like
PF02199	HMMPfam_SapA
PF05184	HMMPfam_SapB_1
PF03489	HMMPfam_SapB_2
SSF47862	superfamily_Saposin_like
PF09064	HMMPfam_Tme5_EGF_like
PF02184	HMMPfam_HAT
PF03298	HMMPfam_Stanniocalcin
PF08603	HMMPfam_CAP_C
SSF69340	superfamily_CARP
PF01213	HMMPfam_CAP_N
SSF101278	superfamily_SSF101278
PF09230	HMMPfam_DFF40
PF02961	HMMPfam_BAF
SSF47798	superfamily_Barrier-to-autointegration factor BAF
PF06248	HMMPfam_Zw10
SSF54991	superfamily_Fdx_AntiC_bd
SSF57625	superfamily_Invertebrate chitin-binding proteins
PF06573	HMMPfam_Churchill
SSF47802	superfamily_DNApol_B_N_like
SSF81585	superfamily_SSF81585
PF04912	HMMPfam_Dynamitin
PF00342	HMMPfam_PGI
SSF53697	superfamily_SSF53697
PF03088	HMMPfam_Str_synth
PF01512	HMMPfam_Complex1_51K
PF02071	HMMPfam_NSF
SSF54762	superfamily_SRP9/14
PF02755	HMMPfam_RPEL
PF01175	HMMPfam_Urocanase
SSF111326	superfamily_Urocanase
PF09329	HMMPfam_zf-primase
PF09332	HMMPfam_Mcm10
PF05680	HMMPfam_ATP-synt_E
PF08450	HMMPfam_SGL
PF06907	HMMPfam_Latexin
PF03024	HMMPfam_Folate_rec
PF03039	HMMPfam_IL12
PF08292	HMMPfam_RNA_pol_Rbc25
SSF52141	superfamily_UDNA_glycsylseSF
SSF48225	superfamily_Glyco_hydro_47
SSF46548	superfamily_Helical_ferredxn
PF00479	HMMPfam_G6PD_N
PF02781	HMMPfam_G6PD_C
PF01776	HMMPfam_Ribosomal_L22e
PF04140	HMMPfam_ICMT
PF00993	HMMPfam_MHC_II_alpha
SSF55120	superfamily_Pseudouridine synthase
PF03402	HMMPfam_V1R
SSF64182	superfamily_DHH phosphoesterases
PF04005	HMMPfam_Hus1
PF04909	HMMPfam_Amidohydro_2
PF00393	HMMPfam_6PGD
PF03299	HMMPfam_TF_AP-2
PF08490	HMMPfam_DUF1744
SSF100909	superfamily_SSF100909
PF08767	HMMPfam_CRM1_C
SSF56425	superfamily_Succinate dehydrogenase/fumarate reductase flavoprotein catalytic domain
PF02947	HMMPfam_Flt3_lig
SSF46774	superfamily_ARID
PF01229	HMMPfam_Glyco_hydro_39
SSF51011	superfamily_SSF51011
SSF69500	superfamily_DTyrtRNA_deacyls
PF08776	HMMPfam_VASP_tetra
PF07047	HMMPfam_OPA3
SSF82171	superfamily_SSF82171
PF04538	HMMPfam_BEX
PF04756	HMMPfam_OST3_OST6
SSF55973	superfamily_S-AdoMet_synt
SSF89009	superfamily_SSF89009
PF01290	HMMPfam_Thymosin
PF01242	HMMPfam_PTPS
SSF55620	superfamily_SSF55620
PF04113	HMMPfam_Gpi16
PF05392	HMMPfam_COX7B
SSF81423	superfamily_SSF81423
SSF50405	superfamily_Actin_crosslink
PF06268	HMMPfam_Fascin
PF08661	HMMPfam_Rep_fac-A_3
PF09341	HMMPfam_Pcc1
PF08390	HMMPfam_TRAM1
PF08160	HMMPfam_NUC156
SSF55116	superfamily_Formiminotr
PF04961	HMMPfam_FTCD_C
SSF101262	superfamily_Cyclodeamin/cyclohydro
PF07837	HMMPfam_FTCD_N
PF02971	HMMPfam_FTCD
SSF51161	superfamily_Trimer_LpxA_like
PF00969	HMMPfam_MHC_II_beta
PF01513	HMMPfam_NAD_kinase
SSF111331	superfamily_SSF111331
PF03348	HMMPfam_Serinc
SSF52768	superfamily_SSF52768
PF02827	HMMPfam_PKI
PF03483	HMMPfam_B3_4
PF03484	HMMPfam_B5
SSF56037	superfamily_B3/B4 domain of PheRS PheT
PF04724	HMMPfam_Glyco_transf_17
SSF101152	superfamily_Mob1_phocein
PF05994	HMMPfam_FragX_IP
PF05076	HMMPfam_SUFU
SSF103359	superfamily_Suppressor of Fused N-terminal domain
SSF57933	superfamily_TAZ_finger
PF09033	HMMPfam_DFF-C
SSF81783	superfamily_C-terminal domain of DFF45/ICAD (DFF-C domain)
PF00235	HMMPfam_Profilin
SSF55770	superfamily_Profilin (actin-binding protein)
PF05422	HMMPfam_SIN1
PF02223	HMMPfam_Thymidylate_kin
PF08213	HMMPfam_DUF1713
PF07092	HMMPfam_DUF1356
PF03730	HMMPfam_Ku_C
PF03731	HMMPfam_Ku_N
PF02735	HMMPfam_Ku
PF05827	HMMPfam_ATP-synt_S1
PF00696	HMMPfam_AA_kinase
SSF53633	superfamily_Carbamate kinase-like
PF05817	HMMPfam_Ribophorin_II
PF06623	HMMPfam_MHC_I_C
PF09042	HMMPfam_Titin_Z
PF06137	HMMPfam_TFA
PF08785	HMMPfam_Ku_PK_bind
SSF101420	superfamily_C-terminal domain of Ku80
PF00396	HMMPfam_Granulin
PF07443	HMMPfam_HARP
PF02177	HMMPfam_A4_EXTRA
SSF89811	superfamily_Amyloid beta a4 protein copper binding domain (domain 2)
SSF56491	superfamily_A heparin-binding domain
SSF109843	superfamily_CAPPD an extracellular domain of amyloid beta A4 protein
SSF47794	superfamily_Rad51_N
PF04566	HMMPfam_RNA_pol_Rpb2_4
SSF64484	superfamily_SSF64484
PF08234	HMMPfam_Spindle_Spc25
PF00384	HMMPfam_Molybdopterin
PF09326	HMMPfam_DUF1982
SSF53706	superfamily_SSF53706
SSF54862	superfamily_SSF54862
PF05571	HMMPfam_DUF766
SSF54427	superfamily_SSF54427
PF04886	HMMPfam_PT
PF01115	HMMPfam_F_actin_cap_B
SSF90096	superfamily_SSF90096
PF05005	HMMPfam_Ocnus
PF01381	HMMPfam_HTH_3
PF00245	HMMPfam_Alk_phosphatase
PF03179	HMMPfam_V-ATPase_G
PF03154	HMMPfam_Atrophin-1
SSF69125	superfamily_Nuc_recept_coact
SSF47040	superfamily_KIX
PF05030	HMMPfam_SSXT
SSF57535	superfamily_SSF57535
PF00883	HMMPfam_Peptidase_M17
PF05090	HMMPfam_VKG_Carbox
SSF52021	superfamily_Carbamoyl phosphate synthetase small subunit N-terminal domain
PF00162	HMMPfam_PGK
SSF53748	superfamily_Phosphoglycerate kinase
PF09088	HMMPfam_MIF4G_like
PF09090	HMMPfam_MIF4G_like_2
SSF74788	superfamily_SSF74788
SSF75632	superfamily_SSF75632
SSF47050	superfamily_VHP
PF04664	HMMPfam_OGFr_N
PF04680	HMMPfam_OGFr_III
PF01642	HMMPfam_MM_CoA_mutase
PF02310	HMMPfam_B12-binding
SSF52242	superfamily_Cobalamin (vitamin B12)-binding domain
SSF51703	superfamily_Cobalamin (vitamin B12)-dependent enzymes
SSF53732	superfamily_Aconitase_N
SSF52016	superfamily_Aconitase/3IPM_dehydase_swvl
PF01549	HMMPfam_ShK
PF05507	HMMPfam_MAGP
PF01808	HMMPfam_AICARFT_IMPCHas
SSF52335	superfamily_Methylglyoxal synthase-like
SSF64197	superfamily_AICAR transformylase domain of bifunctional purine biosynthesis enzyme ATIC
PF05474	HMMPfam_Semenogelin
SSF57256	superfamily_WAP
PF04739	HMMPfam_AMPKBI
PF08122	HMMPfam_NDUF_B12
PF01111	HMMPfam_CKS
SSF55637	superfamily_Cell cycle regulatory proteins
PF04845	HMMPfam_PurA
PF02953	HMMPfam_zf-Tim10_DDP
PF05741	HMMPfam_zf-nanos
PF00491	HMMPfam_Arginase
PF03074	HMMPfam_GCS
PF03332	HMMPfam_PMM
SSF47668	superfamily_Adaptor_Cbl_N
PF08202	HMMPfam_Mis12_component
SSF50974	superfamily_N2O_reductase_N
PF03533	HMMPfam_SPO11_like
PF04406	HMMPfam_TP6A_N
SSF56726	superfamily_SSF56726
SSF54637	superfamily_Thioesterase/thiol ester dehydrase-isomerase
SSF54171	superfamily_SSF54171
SSF47226	superfamily_Histidine-containing phosphotransfer domain HPT domain
PF05786	HMMPfam_Barren
PF08523	HMMPfam_MBF1
PF02602	HMMPfam_HEM4
SSF69618	superfamily_HEM4_synth
PF03345	HMMPfam_DDOST_48kD
SSF81698	superfamily_FF
PF04402	HMMPfam_DUF541
PF01937	HMMPfam_DUF89
SSF111321	superfamily_SSF111321
PF05172	HMMPfam_MPPN
PF02106	HMMPfam_Fanconi_C
PF04629	HMMPfam_ICA69
PF06456	HMMPfam_Arfaptin
SSF57593	superfamily_Heparin-binding domain from vascular endothelial growth factor
SSF49758	superfamily_Peptidase_C2
PF03344	HMMPfam_Daxx
SSF48108	superfamily_Carbamoyl phosphate synthetase large subunit connection domain
PF02219	HMMPfam_MTHFR
SSF51730	superfamily_SSF51730
PF00726	HMMPfam_IL10
PF05511	HMMPfam_ATP-synt_F6
SSF111357	superfamily_Mitochondrial ATP synthase coupling factor 6 (Pfam 05511)
PF05351	HMMPfam_GMP_PDE_delta
PF00212	HMMPfam_ANP
PF04858	HMMPfam_TH1
SSF100939	superfamily_SSF100939
SSF53748	superfamily_PGK
SSF54593	superfamily_SSF54593
PF01286	HMMPfam_XPA_N
PF05181	HMMPfam_XPA_C
PF08784	HMMPfam_RPA_C
SSF55804	superfamily_Phoshotransferase/anion transport protein
SSF48150	superfamily_DNA_glycsylse
SSF55945	superfamily_TFIID_C/glycos_N
PF07934	HMMPfam_OGG_N
PF00497	HMMPfam_SBP_bac_3
PF01566	HMMPfam_Nramp
PF04704	HMMPfam_Zfx_Zfy_act
PF06090	HMMPfam_DUF941
SSF53328	superfamily_formyl_transf
SSF50486	superfamily_FMT_C_like
PF05450	HMMPfam_Nicastrin
PF03177	HMMPfam_Nucleoporin
PF02875	HMMPfam_Mur_ligase_C
SSF53244	superfamily_MurD-like peptide ligases peptide-binding domain
SSF53623	superfamily_MurD-like peptide ligases catalytic domain
PF02724	HMMPfam_CDC45
SSF53795	superfamily_SSF53795
SSF68923	superfamily_SSF68923
PF06391	HMMPfam_MAT1
SSF56019	superfamily_SSF56019
PF04101	HMMPfam_Glyco_tran_28_C
PF02789	HMMPfam_Peptidase_M17_N
PF04627	HMMPfam_ATP-synt_Eps
SSF48690	superfamily_ATP_synth_E
PF03102	HMMPfam_NeuB
PF08666	HMMPfam_SAF
SSF51269	superfamily_SSF51269
PF04440	HMMPfam_Dysbindin
PF04718	HMMPfam_ATP-synt_G
PF05301	HMMPfam_DUF738
PF05456	HMMPfam_eIF_4EBP
PF05770	HMMPfam_Ins134_P3_kin
PF05147	HMMPfam_LANC_like
PF07491	HMMPfam_PPI_Ypi1
SSF47396	superfamily_Transcription factor IIA (TFIIA) alpha-helical domain
SSF50784	superfamily_Transcription factor IIA (TFIIA) beta-barrel domain
SSF47592	superfamily_MDM2
SSF48019	superfamily_DNA polymerase III clamp loader subunits C-terminal domain
PF08519	HMMPfam_RFC1
PF00401	HMMPfam_ATP-synt_DE
PF02823	HMMPfam_ATP-synt_DE_N
SSF46604	superfamily_Epsilon subunit of F1F0-ATP synthase C-terminal domain
SSF51344	superfamily_Epsilon subunit of F1F0-ATP synthase N-terminal domain
PF04900	HMMPfam_Fcf1
SSF88633	superfamily_SSF88633
PF03403	HMMPfam_PAF-AH_p_II
PF03932	HMMPfam_CutC
SSF110395	superfamily_CutC-like (Pfam 03932)
PF05028	HMMPfam_PARG_cat
PF04098	HMMPfam_Rad52_Rad22
SSF56741	superfamily_TopoI_DNA_bd_euk
SSF46596	superfamily_Topismrse_insert
SSF56349	superfamily_DNA_brk_join_enz
PF05544	HMMPfam_Pro_racemase
PF03781	HMMPfam_DUF323
SSF54534	superfamily_SSF54534
PF00834	HMMPfam_Ribul_P_3_epim
PF07933	HMMPfam_DUF1681
PF08040	HMMPfam_NADH_oxidored
PF00584	HMMPfam_SecE
SSF103456	superfamily_Preprotein translocase SecE subunit
SSF50346	superfamily_PRCH_cytoplasmic
PF00120	HMMPfam_Gln-synt_C
PF03951	HMMPfam_Gln-synt_N
PF05821	HMMPfam_NDUF_B8
SSF50800	superfamily_PK_B_barrel_like
PF00224	HMMPfam_PK
PF02887	HMMPfam_PK_C
SSF52935	superfamily_Pyruvate_kinase
SSF48576	superfamily_Terpenoid_synth
PF02267	HMMPfam_Rib_hydrolayse
SSF56629	superfamily_Rib_hydrolayse
SSF50891	superfamily_Cyclophilin-like
SSF81423	superfamily_Mitochondrial cytochrome c oxidase subunit VIIb
SSF110296	superfamily_SSF110296
SSF69099	superfamily_RanGAP1_C
PF00878	HMMPfam_CIMR
PF07926	HMMPfam_TPR_MLP1_2
PF02889	HMMPfam_Sec63
PF03800	HMMPfam_Nuf2
PF03792	HMMPfam_PBC
PF00203	HMMPfam_Ribosomal_S19
SSF54570	superfamily_Ribosomal protein S19
PF03721	HMMPfam_UDPG_MGDP_dh_N
PF00984	HMMPfam_UDPG_MGDP_dh
PF03720	HMMPfam_UDPG_MGDP_dh_C
SSF52413	superfamily_UDP-Glc/GDP-Man_DH_C
PF09270	HMMPfam_Beta-trefoil
SSF110217	superfamily_DNA-binding protein LAG-1 (CSL)
PF09271	HMMPfam_LAG1-DNAbind
PF00940	HMMPfam_RNA_pol
PF05622	HMMPfam_HOOK
SSF63867	superfamily_SSF63867
PF05808	HMMPfam_Podoplanin
PF05337	HMMPfam_CSF-1
PF09180	HMMPfam_ProRS-C_1
SSF64586	superfamily_SSF64586
PF06827	HMMPfam_zf-FPG_IleRS
PF05602	HMMPfam_CLPTM1
SSF103612	superfamily_SBT domain
PF05790	HMMPfam_C2-set
PF09191	HMMPfam_CD4-extracel
SSF48092	superfamily_STAT
PF01588	HMMPfam_tRNA_bind
PF03694	HMMPfam_Erg28
SSF111384	superfamily_OmpH-like (Pfam 03938)
SSF50974	superfamily_Nitrous oxide reductase N-terminal domain
PF01120	HMMPfam_Alpha_L_fucos
PF00703	HMMPfam_Glyco_hydro_2
SSF49303	superfamily_beta-Galactosidase/glucuronidase domain
PF02837	HMMPfam_Glyco_hydro_2_N
PF00551	HMMPfam_Formyl_trans_N
SSF110004	superfamily_Glycolipid transfer protein GLTP
PF07140	HMMPfam_IFNGR1
SSF81431	superfamily_SSF81431
PF01527	HMMPfam_Transposase_8
SSF54368	superfamily_Gln_synt_beta
PF00904	HMMPfam_Involucrin
PF02159	HMMPfam_Oest_recep
PF02320	HMMPfam_UCR_hinge
SSF81531	superfamily_SSF81531
PF03980	HMMPfam_Nnf1
SSF56629	superfamily_ADP ribosyl cyclase-like
PF05007	HMMPfam_Mannosyl_trans
PF05179	HMMPfam_CDC73
PF04188	HMMPfam_Mannosyl_trans2
PF04831	HMMPfam_Popeye
PF09011	HMMPfam_DUF1898
PF05177	HMMPfam_RCSD
PF03650	HMMPfam_UPF0041
SSF49695	superfamily_G_crystallin_SF
PF00636	HMMPfam_Ribonuclease_3
SSF69065	superfamily_RNase III catalytic domain-like (Pfam 00636)
PF02354	HMMPfam_Latrophilin
PF06729	HMMPfam_NRIF3
PF04553	HMMPfam_Tis11B_N
PF00408	HMMPfam_PGM_PMM_IV
PF02880	HMMPfam_PGM_PMM_III
SSF53738	superfamily_SSF53738
PF07541	HMMPfam_EIF_2_alpha
SSF110993	superfamily_eIF-2-alpha C-terminal domain
PF00312	HMMPfam_Ribosomal_S15
PF08999	HMMPfam_SP_C-Propep
SSF54556	superfamily_Chitinase insertion domain
PF01607	HMMPfam_CBM_14
PF05405	HMMPfam_Mt_ATP-synt_B
PF01008	HMMPfam_IF-2B
PF09120	HMMPfam_EGF-like_subdom
PF03494	HMMPfam_Beta-APP
PF04488	HMMPfam_Gly_transf_sug
PF04572	HMMPfam_Gb3_synth
PF05483	HMMPfam_SCP-1
PF06638	HMMPfam_Strabismus
SSF57774	superfamily_Microbial and mitochondrial ADK insert "zinc finger" domain
SSF50621	superfamily_Racem_decarbox_C
SSF51419	superfamily_SSF51419
PF04272	HMMPfam_Phospholamban
SSF101546	superfamily_ASF1-like
SSF47587	superfamily_SSF47587
PF08082	HMMPfam_PRO8NT
PF08083	HMMPfam_PROCN
PF08084	HMMPfam_PROCT
PF05536	HMMPfam_Neurochondrin
PF01813	HMMPfam_ATP-synt_D
PF06140	HMMPfam_Ifi-6-16
PF01368	HMMPfam_DHH
PF02833	HMMPfam_DHHA2
SSF64182	superfamily_SSF64182
PF05764	HMMPfam_YL1
PF02070	HMMPfam_NMU
PF08641	HMMPfam_Mis14
PF02058	HMMPfam_Guanylin
SSF89890	superfamily_Guanylin
PF04127	HMMPfam_DFP
PF00923	HMMPfam_Transaldolase
PF01080	HMMPfam_Presenilin
PF03054	HMMPfam_tRNA_Me_trans
PF01521	HMMPfam_Fe-S_biosyn
SSF89360	superfamily_SSF89360
PF01093	HMMPfam_Clusterin
PF02389	HMMPfam_Cornifin
PF02760	HMMPfam_HIN
PF01216	HMMPfam_Calsequestrin
PF04614	HMMPfam_Pex19
SSF54556	superfamily_SSF54556
PF04636	HMMPfam_PA26
SSF69118	superfamily_AhpD-like
PF01536	HMMPfam_SAM_decarbox
SSF56276	superfamily_SSF56276
PF04044	HMMPfam_Nup133
PF08801	HMMPfam_Nup133_N
PF08434	HMMPfam_CLCA_N
PF09315	HMMPfam_DUF1973
PF05040	HMMPfam_HS2ST
SSF55957	superfamily_SSF55957
PF06328	HMMPfam_Lep_receptor_Ig
PF01133	HMMPfam_ER
SSF49742	superfamily_PHM_PNGase_F
SSF56317	superfamily_Ntlse/CNhydtse
PF08063	HMMPfam_PADR1
PF01480	HMMPfam_PWI
PF04709	HMMPfam_AMH_N
PF02100	HMMPfam_ODC_AZ
PF02144	HMMPfam_Rad1
PF00121	HMMPfam_TIM
SSF51351	superfamily_Triosephosphate isomerase (TIM)
PF00186	HMMPfam_DHFR_1
SSF53597	superfamily_SSF53597
PF05182	HMMPfam_Fip1
SSF48557	superfamily_L-Aspartase-like
PF08265	HMMPfam_YL1_C
PF03152	HMMPfam_UFD1
PF09166	HMMPfam_Biliv-reduc_cat
PF01191	HMMPfam_RNA_pol_Rpb5_C
SSF55287	superfamily_RNApol_RPB5_like
PF03871	HMMPfam_RNA_pol_Rpb5_N
SSF53036	superfamily_RNA_pol_Rpb5_N
PF05478	HMMPfam_Prominin
PF04934	HMMPfam_MED6
PF00743	HMMPfam_FMO-like
PF06432	HMMPfam_GPI2
PF09294	HMMPfam_Interfer-bind
PF08365	HMMPfam_IGF2_C
SSF48613	superfamily_Heme_oxygenase
PF06467	HMMPfam_zf-MYM
PF01963	HMMPfam_TraB
SSF102645	superfamily_SSF102645
PF05694	HMMPfam_SBP56
SSF51004	superfamily_Cyt_cd1_haem_C
PF00837	HMMPfam_T4_deiodinase
SSF102114	superfamily_Radical SAM enzymes
PF00994	HMMPfam_MoCF_biosynth
SSF53218	superfamily_Molybdenum cofactor biosynthesis proteins
PF01507	HMMPfam_PAPS_reduct
PF06331	HMMPfam_Tbf5
PF08694	HMMPfam_DUF1782
PF01738	HMMPfam_DLH
PF00346	HMMPfam_Complex1_49kDa
SSF56762	superfamily_Nickel-iron hydrogenase large subunit
PF04711	HMMPfam_ApoA-II
SSF53218	superfamily_MoCF_biosynth
PF03453	HMMPfam_MoeA_N
SSF63882	superfamily_MoeA_N
PF03454	HMMPfam_MoeA_C
PF00236	HMMPfam_Hormone_6
SSF55190	superfamily_SSF55190
PF00809	HMMPfam_Pterin_bind
PF02607	HMMPfam_B12-binding_2
SSF47644	superfamily_Met_synth_Cbl-bd
PF02965	HMMPfam_Met_synt_B12
SSF52242	superfamily_Cbl-bd
SSF51717	superfamily_DHP_synth_like
SSF56507	superfamily_SSF56507
PF02841	HMMPfam_GBP_C
SSF48340	superfamily_GBP
PF03508	HMMPfam_Connexin43
PF06202	HMMPfam_GDE_C
PF02063	HMMPfam_MARCKS
PF08997	HMMPfam_UCR_6-4kD
SSF81518	superfamily_Ub-CytCRdtase_cplx_6.4kD
PF07569	HMMPfam_Hira
SSF69593	superfamily_SSF69593
PF01997	HMMPfam_Translin
PF01154	HMMPfam_HMG_CoA_synt_N
PF08540	HMMPfam_HMG_CoA_synt_C
PF01144	HMMPfam_CoA_trans
PF01715	HMMPfam_IPPT
PF02089	HMMPfam_Palm_thioest
PF08660	HMMPfam_Alg14
PF06758	HMMPfam_DUF1220
SSF56209	superfamily_Nitrile hydratase alpha chain
PF03281	HMMPfam_Mab-21
PF08911	HMMPfam_NUP50
PF06087	HMMPfam_Tyr-DNA_phospho
PF05793	HMMPfam_TFIIF_alpha
SSF50916	superfamily_Rap30/74 interaction domains
PF03849	HMMPfam_Tfb2
PF07741	HMMPfam_BRF1
PF08271	HMMPfam_TFIIB_Zn_Ribbon
PF00382	HMMPfam_TFIIB
SSF89890	superfamily_Proguanylin
PF01261	HMMPfam_AP_endonuc_2
SSF51658	superfamily_SSF51658
PF09340	HMMPfam_NuA4
PF08153	HMMPfam_NGP1NT
PF02055	HMMPfam_Glyco_hydro_30
PF02938	HMMPfam_GAD
PF05499	HMMPfam_DMAP1
PF00659	HMMPfam_POLO_box
SSF47686	superfamily_Anaphylotoxins (complement system)
SSF48552	superfamily_Serum albumin-like
PF00273	HMMPfam_Serum_albumin
PF09164	HMMPfam_VitD-bind_III
PF02295	HMMPfam_z-alpha
PF02872	HMMPfam_5_nucleotid_C
SSF55816	superfamily_5'-nucleotidase (syn. UDP-sugar hydrolase) C-terminal domain
SSF55261	superfamily_SSF55261
SSF57414	superfamily_SSF57414
SSF47794	superfamily_Rad51 N-terminal domain-like
SSF55781	superfamily_SSF55781
PF05645	HMMPfam_RNA_pol_Rpc82
PF08221	HMMPfam_HTH_9
PF02114	HMMPfam_Phosducin
PF04670	HMMPfam_Gtr1_RagA
PF01239	HMMPfam_PPTA
SSF55257	superfamily_RNAP_RBP11-like
PF05719	HMMPfam_GPP34
PF00034	HMMPfam_Cytochrom_C
SSF103647	superfamily_TSP type-3 repeat
PF01187	HMMPfam_MIF
SSF55331	superfamily_Tautomerase/MIF
PF00081	HMMPfam_Sod_Fe_N
PF02777	HMMPfam_Sod_Fe_C
SSF46609	superfamily_FeMn superoxide dismutase (SOD) N-terminal domain
SSF54719	superfamily_FeMn superoxide dismutase (SOD) C-terminal domain
PF04857	HMMPfam_CAF1
PF02233	HMMPfam_PNTB
PF04178	HMMPfam_Got1
PF09060	HMMPfam_L27_N
SSF47741	superfamily_CO dehydrogenase ISP C-domain like
SSF55447	superfamily_CO dehydrogenase flavoprotein C-terminal domain-like
SSF56003	superfamily_Molybdenum cofactor-binding domain
PF08377	HMMPfam_MAP2_projctn
SSF57770	superfamily_Methionyl-tRNA synthetase (MetRS) Zn-domain
PF08123	HMMPfam_DOT1
PF01208	HMMPfam_URO-D
SSF51726	superfamily_SSF51726
PF02368	HMMPfam_Big_2
SSF49373	superfamily_Invasin_intimin
PF04275	HMMPfam_P-mevalo_kinase
SSF101160	superfamily_Dimerisation domain of CENP-B
PF01053	HMMPfam_Cys_Met_Meta_PP
PF07770	HMMPfam_SFT2
PF00445	HMMPfam_Ribonuclease_T2
SSF55895	superfamily_Ribonuclease Rh-like
SSF54975	superfamily_Acylphosphatase-like
SSF54368	superfamily_Glutamine synthetase N-terminal domain
SSF51412	superfamily_Inosine monophosphate dehydrogenase (IMPDH)
SSF54862	superfamily_4Fe-4S ferredoxins
SSF47005	superfamily_Peripheral subunit-binding domain of 2-oxo acid dehydrogenase complex
PF01172	HMMPfam_SBDS
SSF109728	superfamily_Hypothetical protein AF0491 middle domain
SSF89895	superfamily_FYSH domain
PF01451	HMMPfam_LMWPc
SSF55418	superfamily_Translation initiation factor eIF4e
SSF48690	superfamily_Epsilon subunit of mitochondrial F1F0-ATP synthase
PF00576	HMMPfam_Transthyretin
SSF49472	superfamily_Transthyretin (prealbumin)
SSF53182	superfamily_Pyrrolidone carboxyl peptidase (pyroglutamate aminopeptidase)
PF05355	HMMPfam_Apo-CII
PF08212	HMMPfam_Lipocalin_2
PF00473	HMMPfam_CRF
PF07883	HMMPfam_Cupin_2
PF01557	HMMPfam_FAA_hydrolase
SSF56529	superfamily_FAH
PF09298	HMMPfam_DUF1969
SSF63433	superfamily_Fumarylacetoacetate hydrolase FAH N-terminal domain
SSF55044	superfamily_TRADD N-terminal domain
SSF50904	superfamily_TCL1_MTCP1
SSF56276	superfamily_S-adenosylmethionine decarboxylase
PF01505	HMMPfam_Vault
PF01130	HMMPfam_CD36
SSF50715	superfamily_Ribosomal protein L25-like
SSF64586	superfamily_C-terminal domain of ProRS
PF01087	HMMPfam_GalP_UDP_transf
PF02744	HMMPfam_GalP_UDP_tr_C
PF09273	HMMPfam_Rubis-subs-bind
PF01071	HMMPfam_GARS_A
PF02843	HMMPfam_GARS_C
PF02844	HMMPfam_GARS_N
SSF53328	superfamily_Formyltransferase
SSF81411	superfamily_Mitochondrial cytochrome c oxidase subunit VIa
SSF49606	superfamily_Neurophysin II
SSF56994	superfamily_SSF56994
SSF53597	superfamily_Dihydrofolate reductases
SSF50723	superfamily_Core binding factor beta CBF
SSF50916	superfamily_TFIIF_interactn
SSF69065	superfamily_RNase_III
PF03368	HMMPfam_DUF283
SSF101690	superfamily_SSF101690
SSF52935	superfamily_PK C-terminal domain-like
PF04592	HMMPfam_SelP_N
PF04593	HMMPfam_SelP_C
PF02671	HMMPfam_PAH
PF08295	HMMPfam_HDAC_interact
SSF47694	superfamily_Cytochrome c oxidase subunit h
SSF81415	superfamily_Mitochondrial cytochrome c oxidase subunit VIc
SSF102886	superfamily_Coproporphyrinogen III oxidase
SSF82282	superfamily_Homocysteine S-methyltransferase
PF01214	HMMPfam_CK_II_beta
SSF53671	superfamily_Aspartate/ornithine carbamoyltransferase
PF03095	HMMPfam_PTPA
PF01169	HMMPfam_UPF0016
PF07500	HMMPfam_TFIIS_M
SSF46942	superfamily_Elongation factor TFIIS domain 2
PF08148	HMMPfam_DSHCT
SSF89000	superfamily_post-HMGL domain-like
PF09243	HMMPfam_Rsm22
SSF103196	superfamily_Roadblock/LC7 domain
SSF81827	superfamily_Porin chaperone SurA peptide-binding domain
SSF51713	superfamily_tRNA-guanine transglycosylase
PF05495	HMMPfam_zf-CHY
SSF102848	superfamily_NSFL1 (p97 ATPase) cofactor p47 SEP domain
SSF57798	superfamily_Casein kinase II beta subunit
PF07535	HMMPfam_zf-DBF
PF02348	HMMPfam_CTP_transf_3
SSF51283	superfamily_dUTPase-like
SSF69340	superfamily_C-terminal domain of adenylylcyclase associated protein
SSF63887	superfamily_P-domain of calnexin/calreticulin
PF04716	HMMPfam_ETC_C1_NDUFA5
SSF54665	superfamily_CO dehydrogenase molybdoprotein N-domain-like
SSF82615	superfamily_SSF82615
SSF47789	superfamily_C-terminal domain of RNA polymerase alpha subunit
SSF47912	superfamily_Wiscott-Aldrich syndrome protein WASP C-terminal domain
PF02911	HMMPfam_Formyl_trans_C
PF02515	HMMPfam_CoA_transf_3
SSF89796	superfamily_CoA-transferase family III (CaiB/BaiF)
SSF102712	superfamily_JAB1/MPN domain
SSF53706	superfamily_Formate dehydrogenase/DMSO reductase domains 1-3
SSF103111	superfamily_Activator of Hsp90 ATPase Aha1
PF05837	HMMPfam_CENP-H
PF06387	HMMPfam_Calcyon
PF07326	HMMPfam_DUF1466
PF03835	HMMPfam_Rad4
PF00814	HMMPfam_Peptidase_M22
SSF55116	superfamily_Formiminotransferase domain of formiminotransferase-cyclodeaminase.
PF05365	HMMPfam_UCR_UQCRX_QCR9
PF01713	HMMPfam_Smr
PF08590	HMMPfam_DUF1771
SSF53092	superfamily_Creatinase/prolidase N-terminal domain
PF04800	HMMPfam_ETC_C1_NDUFA4
SSF111331	superfamily_NAD kinase (Pfam 01513)
SSF74653	superfamily_TolA/TonB C-terminal domain
SSF51658	superfamily_Xylose isomerase-like
PF04801	HMMPfam_Sin_N
PF04096	HMMPfam_Nucleoporin2
SSF82215	superfamily_C-terminal autoproteolytic domain of nucleoporin nup98
SSF56214	superfamily_4'-phosphopantetheinyl transferase
PF03215	HMMPfam_Rad17
SSF102588	superfamily_LmbE-like
PF00731	HMMPfam_AIRC
SSF52255	superfamily_N5-CAIR mutase (phosphoribosylaminoimidazole carboxylase PurE)
PF01259	HMMPfam_SAICAR_synt
PF00338	HMMPfam_Ribosomal_S10
PF04855	HMMPfam_SNF5
SSF63862	superfamily_Thiamin pyrophosphokinase substrate-binding domain
SSF63999	superfamily_Thiamin pyrophosphokinase catalytic domain
PF02284	HMMPfam_COX5A
SSF48479	superfamily_Cytochrome c oxidase subunit E
PF07787	HMMPfam_DUF1625
PF05835	HMMPfam_Synaphin
PF05348	HMMPfam_UMP1
PF09031	HMMPfam_SRI_2
SSF46942	superfamily_SSF46942
PF00521	HMMPfam_DNA_topoisoIV
PF08070	HMMPfam_DTHCT
PF00204	HMMPfam_DNA_gyraseB
SSF56719	superfamily_Type II DNA topoisomerase
SSF81531	superfamily_Non-heme 11 kDa protein of cytochrome bc1 complex (Ubiquinol-cytochrome c reductase)
PF07819	HMMPfam_PGAP1
SSF55277	superfamily_GYF domain
SSF51726	superfamily_UROD/MetE-like
PF05208	HMMPfam_ALG3
PF01020	HMMPfam_Ribosomal_L40e
PF06775	HMMPfam_DUF1226
SSF63882	superfamily_MoeA N-terminal region -like
SSF63867	superfamily_MoeA C-terminal domain-like
PF08944	HMMPfam_p47_phox_C
PF03224	HMMPfam_V-ATPase_H
SSF55261	superfamily_Prokaryotic AspRS insert domain
PF07304	HMMPfam_SRA1
SSF109635	superfamily_DnaK suppressor protein DksA alpha-hairpin domain
PF07767	HMMPfam_Nop53
PF07557	HMMPfam_Shugoshin_C
SSF56529	superfamily_Fumarylacetoacetase_C-rel
SSF63433	superfamily_Fumarylacetoacetase_N
PF00890	HMMPfam_FAD_binding_2
PF02910	HMMPfam_Succ_DH_flav_C
SSF46977	superfamily_Succinate dehydrogenase/fumarate reductase flavoprotein C-terminal domain
PF01798	HMMPfam_Nop
PF08060	HMMPfam_NOSIC
SSF89124	superfamily_Nop domain
PF07994	HMMPfam_NAD_binding_5
PF01658	HMMPfam_Inos-1-P_synth
SSF47644	superfamily_Methionine synthase domain
SSF56507	superfamily_Methionine synthase activation domain-like
PF01791	HMMPfam_DeoC
PF05859	HMMPfam_Mis12
PF09302	HMMPfam_XLF
SSF54076	superfamily_RNaseA
PF01459	HMMPfam_Porin_3
PF03200	HMMPfam_Glyco_hydro_63
SSF56726	superfamily_DNA topoisomerase IV alpha subunit
PF02291	HMMPfam_TFIID-31kDa
SSF50118	superfamily_Cell growth inhibitor/plasmid maintenance toxic component
SSF48552	superfamily_Serum_albumin
SSF111321	superfamily_Hypothetical protein At2g17340
SSF51269	superfamily_AFP III-like domain (Pfam 01354)
SSF52413	superfamily_UDP-glucose/GDP-mannose dehydrogenase C-terminal domain
SSF69618	superfamily_Uroporphyrinogen III synthase (U3S HemD)
SSF55394	superfamily_Bactericidal permeability-increasing protein BPI
SSF48168	superfamily_Ribonucleo_red_N
SSF51998	superfamily_SSF51998
SSF90234	superfamily_Zinc finger domain of DNA polymerase-alpha
SSF54713	superfamily_Elongation factor Ts (EF-Ts) dimerisation domain
SSF82704	superfamily_SSF82704
PF04926	HMMPfam_PAP_RNA-bind
PF04928	HMMPfam_PAP_central
SSF55003	superfamily_PAP/Archaeal CCA-adding enzyme C-terminal domain
PF01301	HMMPfam_Glyco_hydro_35
SSF89360	superfamily_HesB-like domain
PF09272	HMMPfam_Hepsin-SRCR
PF05161	HMMPfam_MOFRL
SSF82544	superfamily_GckA/TtuD-like
SSF46938	superfamily_Sec14p_like_N
SSF48662	superfamily_Ribosomal protein L39e
SSF47762	superfamily_SSF47762
PF01755	HMMPfam_Glyco_transf_25
SSF101224	superfamily_HAND domain of the nucleosome remodeling ATPase ISWI
PF01058	HMMPfam_Oxidored_q6
SSF55658	superfamily_L9 N-domain-like
SSF111326	superfamily_Urocanase (Pfam 01175)
SSF110849	superfamily_ParB/Sulfiredoxin
PF07722	HMMPfam_Peptidase_C26
PF06773	HMMPfam_Bim_N
PF08945	HMMPfam_Bclx_interact
PF07774	HMMPfam_DUF1620
SSF48163	superfamily_An anticodon-binding domain of class I aminoacyl-tRNA synthetases
PF03959	HMMPfam_FSH1
PF07947	HMMPfam_YhhN
PF00465	HMMPfam_Fe-ADH
SSF56796	superfamily_Dehydroquinate synthase-like
SSF81872	superfamily_SSF81872
SSF81878	superfamily_SSF81878
PF09127	HMMPfam_Leuk-A4-hydro_C
SSF49363	superfamily_Purple acid phosphatase N-terminal domain
PF04515	HMMPfam_DUF580
SSF55797	superfamily_SSF55797
SSF82607	superfamily_YbaB-like
SSF46585	superfamily_PKN_effector
PF07713	HMMPfam_DUF1604
SSF48479	superfamily_Cyt_c_ox5A
PF06229	HMMPfam_FRG1
SSF53036	superfamily_Eukaryotic RPB5 N-terminal domain
PF07347	HMMPfam_CI-B14_5a
PF08572	HMMPfam_PRP3
PF01296	HMMPfam_Galanin
PF08743	HMMPfam_Nse4
PF05178	HMMPfam_Krr1
PF06239	HMMPfam_ECSIT
PF03378	HMMPfam_CAS_CSE1
PF08506	HMMPfam_Cse1
PF01510	HMMPfam_Amidase_2
SSF55846	superfamily_SSF55846
PF06375	HMMPfam_BLVR
PF01199	HMMPfam_Ribosomal_L34e
SSF81573	superfamily_F1F0 ATP synthase subunit B membrane domain
SSF82215	superfamily_SSF82215
PF03878	HMMPfam_YIF1
PF03285	HMMPfam_Paralemmin
PF02134	HMMPfam_UBACT
PF00899	HMMPfam_ThiF
SSF56935	superfamily_Porins
PF04502	HMMPfam_DUF572
PF06244	HMMPfam_DUF1014
PF01958	HMMPfam_DUF108
PF03447	HMMPfam_NAD_binding_3
PF01171	HMMPfam_ATP_bind_3
SSF53032	superfamily_tRNA-intron endonuclease catalytic domain-like
PF06350	HMMPfam_HSL_N
SSF55144	superfamily_Cyc_nuc_Pdiester
PF02096	HMMPfam_60KD_IMP
PF03398	HMMPfam_DUF292
PF09085	HMMPfam_Adhes-Ig_like
SSF56568	superfamily_Non-globular alpha beta subunits of globular proteins
SSF81502	superfamily_ISP transmembrane anchor
PF05064	HMMPfam_Nsp1_C
PF03439	HMMPfam_Supt5
PF05049	HMMPfam_IIGP
SSF109775	superfamily_SSF109775
PF00380	HMMPfam_Ribosomal_S9
SSF82679	superfamily_N-utilization substance G protein NusG N-terminal domain
PF01269	HMMPfam_Fibrillarin
PF01012	HMMPfam_ETF
EOF
    my @return = split "\n", $accession_list;
return \@return;
}
    
1;
