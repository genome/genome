package Genome::Info::SnapModelFileAbbreviations;

use strict;
use warnings;

my %model_files = (
    'A.canium.hmm'                          => 'sn_acan',
    'A_ceylanicum_1.3.ec.cg.pg.cgmkv2.hmm'  => 'sn_ac13cgmkv2',
    'A_ceylanicum_1.3.ec.cg.pg.mkv2.hmm'    => 'sn_ac13mkv2',
    'A_duodenale.v2.2.cg.mk2.hmm'           => 'sn_ad22cgmkv2',
    'A_duodenale.v2.2.mk2.hmm'              => 'sn_ad22mkv2',
    'A_duodenale.v2.2.cg.mk2a.hmm'          => 'sn_ad22cgmkv2a',
    'A_duodenale.v2.2.mk2a.hmm'             => 'sn_ad22mkv2a',
    'Acaninum.9.3.2.cegma.hmm '             => 'sn_accg932',
    'Acaninum.9.3.2.cgv2.hmm'               => 'sn_accg932v2',
    'Acaninum.9.3.2.mkv2.hmm'               => 'sn_acmk932v2',
    'B.malayi.hmm'                          => 'sn_bmal',
    'C.elegans.hmm'                         => 'sn_cele',
    'D.viviparus_cegma.hmm'                 => 'sn_dvcg',
    'D_viviparus_9.3.1.ec.cgmkv2.hmm'       => 'sn_dv931cgmkv2',
    'D_viviparus_9.3.1.ec.mkv2.hmm'         => 'sn_dv931mkv2',
    'N.americanus_cegma.hmm'                => 'sn_nacg',
    'Namericanus.4.2.1.cegma.hmm'           => 'sn_nacg421',
    'Namericanus_4.2.2.ec.cg.v4.hmm'        => 'sn_nacg422ecv4',
    'Namericanus_4.2.2.ec.mkcb.v3.hmm'      => 'sn_namkcb422ecv3',
    'Namericanus_4.2.2.ec.mkis.v3.hmm'      => 'sn_namkis422ecv3',
    'Namericanus_4.2.2.ec.mkrd.v3.hmm'      => 'sn_namkrd422ecv3',
    'S_cerevisiae_ref_cegma.hmm'            => 'sn_screfcg',
    'T_circumcincta.v14.cg.mk2.hmm'         => 'sn_tc14cgmkv2',
    'T_circumcincta.v14.mk2.hmm'            => 'sn_tc14mkv2',
    'T_suis.v1.cgmk2.hmm'                   => 'n_ts1cgmk2',
    'T_suis.v1.mk1.hmm'                     => 'sn_ts1mk1',
    'bmal.intronset.hmm'                    => 'sn_bmal30',
    'caninum_cegma.hmm'                     => 'sn_accg',
    'caninum_maker.hmm'                     => 'sn_acmk',
    'spiralis_cegma_mod.hmm'                => 'sn_tscg',
    'spiralis_cegma_mod_intron_min30.hmm'   => 'sn_tscg30',
    'O_dentatum_10.0.ec.cg.mk2.hmm'         => 'sn_od10cgmkv2',
    'O_dentatum_10.0.cg.mk2.hmm'            => 'sn_od10mkv2',
    'C_aethiops.hmm'			    => 'sn_vercgmk2',
    'F_hepatica_v1.0.allpaths.pg.cgmk2.hmm' => 'sn_fh1cgmk2',
    'F_hepatica_v1.0.allpaths.pg.mk2.hmm'   => 'sn_fh1mk2',
);


sub abbreviation_for_model_file {
    my $model_file = shift;
    return $model_files{$model_file};
}

1;

