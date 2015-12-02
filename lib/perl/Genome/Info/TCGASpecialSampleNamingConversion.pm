package Genome::Info::TCGASpecialSampleNamingConversion;

use strict;
use warnings;

# This module holds the TCGA special sample naming conversions, which
# are not handled properly on LIMS. For example, LAML samples from
# Ley's lab have caTissue not TCGA as nomenclature and they don't have
# TCGA- names as sample name or extraction_label, and their sample
# attributes' external_name often contain multiple TCGA- names from
# different sources. The only way to find correct TCGA names is to
# parse TCGA data inventory .dat files. Refer to jira issue:
# https://jira.gsc.wustl.edu/browse/APIPE-2742

my %conversion = (
    'H_GV-933124G-S.9043'   =>  'TCGA-AB-2993-03A-01D-0741-05',
    'H_KA-103342-0814958'   =>  'TCGA-AB-2978-03A-01D-0742-05',
    'H_KA-113971-0814714'   =>  'TCGA-AB-2966-03A-01D-0741-05',
    'H_KA-123172G-S.3627'   =>  'TCGA-AB-2963-03A-01D-0741-05',
    'H_KA-142074-0814712'   =>  'TCGA-AB-2965-03A-01D-0741-05',
    'H_KA-179223-0814937'   =>  'TCGA-AB-2972-03A-01D-0741-05',
    'H_KA-202127-0903575'   =>  'TCGA-AB-2998-03A-01D-0742-05',
    'H_KA-224143-0903576'   =>  'TCGA-AB-2988-03A-01D-0742-05',
    'H_KA-225373-0814956'   =>  'TCGA-AB-2977-03A-01D-0742-05',
    'H_KA-246634-0814715'   =>  'TCGA-AB-2968-03A-01D-0741-05',
    'H_KA-254137-0909964'   =>  'TCGA-AB-2986-03A-01D-0739-09',
    'H_KA-255421-0927564'   =>  'TCGA-AB-2918-11A-01W-0761-09',
    'H_KA-273919-S.18840'   =>  'TCGA-AB-3000-03A-01D-0741-05',
    'H_KA-321258-S.16585'   =>  'TCGA-AB-3001-03A-01D-0741-05',
    'H_KA-335640-0815277'   =>  'TCGA-AB-2974-03A-01D-0742-05',
    'H_KA-400220-0814727'   =>  'TCGA-AB-2970-03A-01D-0741-05',
    'H_KA-440422-0909966'   =>  'TCGA-AB-2907-03A-01D-0742-05',
    'H_KA-445045-S.18834'   =>  'TCGA-AB-2996-03A-01D-0741-05',
    'H_KA-452198-0814719'   =>  'TCGA-AB-2969-03A-01D-0739-09',
    'H_KA-455499-0816986'   =>  'TCGA-AB-2982-03A-01D-0739-09',
    'H_KA-456892-0815956'   =>  'TCGA-AB-2967-03A-01D-0742-05',
    'H_KA-501944G-S.9049'   =>  'TCGA-AB-2991-03A-01D-0741-05',
    'H_KA-529205-0903574'   =>  'TCGA-AB-2906-03A-01D-0742-05',
    'H_KA-545259-0815275'   =>  'TCGA-AB-2979-03A-01D-0742-05',
    'H_KA-548327-0903586'   =>  'TCGA-AB-2990-03A-01D-0742-05',
    'H_KA-573988-0814940'   =>  'TCGA-AB-2973-03A-01D-0742-05',
    'H_KA-617776-1007333'   =>  'TCGA-AB-2930-03A-01W-0761-09',
    'H_KA-673778G-S.15614'  =>  'TCGA-AB-3012-03A-01D-0741-05',
    'H_KA-700717-0816989'   =>  'TCGA-AB-2983-03A-01D-0742-05',
    'H_KA-702808-0909968'   =>  'TCGA-AB-2987-03A-01D-0742-05',
    'H_KA-753374-0909960'   =>  'TCGA-AB-2989-03A-01D-0742-05',
    'H_KA-775109-S.16152'   =>  'TCGA-AB-3005-03A-01D-0739-09',
    'H_KA-808642-S.18835'   =>  'TCGA-AB-3006-03A-01D-0741-05',
    'H_KA-816067-0815276'   =>  'TCGA-AB-2981-03A-01D-0742-05',
    'H_KA-817156-0814950'   =>  'TCGA-AB-2975-03A-01D-0739-09',
    'H_KA-831711-S.22465'   =>  'TCGA-AB-2964-03A-01D-0741-05',
    'H_KA-849660G-S.16586'  =>  'TCGA-AB-3008-03A-01D-0741-05',
    'H_KA-869586G-S.16427'  =>  'TCGA-AB-3009-03A-01D-0741-05',
    'H_KA-906708-0814938'   =>  'TCGA-AB-2971-03A-01D-0741-05',
    'H_KA-907786-0909970'   =>  'TCGA-AB-2985-03A-01D-0742-05',
    'H_KA-943309-S.18828'   =>  'TCGA-AB-3007-03A-01D-0741-05',
    'H_KA-991612-S.18829'   =>  'TCGA-AB-2995-03A-01D-0741-05',
);

sub tcga_naming_conversion {
    return %conversion; 
}

1;

