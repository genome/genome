package Genome::Model::Tools::CgHub::UpdateFileForImport;

use strict;
use warnings;

use Genome;

use List::MoreUtils;
use Params::Validate ':types';
use Text::CSV;

class Genome::Model::Tools::CgHub::UpdateFileForImport {
    is => 'Command::V2',
    has_input => {
        file => {
            is => 'Text',
            doc => 'Comma (.csv) or tab (.tsv) separated file of that has an "analysis_id". Separator is determined by file extension.',
        },
    },
    has_output => {
        output_file => { 
            is => 'Text',
            doc => 'Output file. Will automatically use the separator of the incoming file.',
        },
    },
    has_optional_transient => {
        _incoming_data => { is => 'ARRAY', },
        _incoming_headers => { is => 'ARRAY', },
        _incoming_sep_char => { is => 'ARRAY', },
        _metadata => { is => 'Genome::Model::Tools::CgHub::Metadata', },
    },
};

sub tcga_attribute_names {
    return (qw/
        individual.name individual.taxon individual.upn
        instdata.analysis_id
        library.name
        sample.common_name sample.disease_abbr
        sample.extraction_type sample.name sample.nomenclature
        sample.tissue_desc
        /);
}

sub help_brief {
    return 'Derive info from tcga legacy sample names';
}

sub help_detail {
    my $help = <<HELP;
This command takes a CSV/TSV file that minimally contains a column with the CG Hub 'analysis_id'. It outputs a CSV/TSV file adding the attributes needed for import. Th original data is retained, except if the column is one of the attributes listed below. In that case, the old value is overwritten. The attributes are queried directly from CG Hub and translated as necessary into TGI attributes and values. These attributes are:

HELP
    $help .= join("\n\n", tcga_attribute_names());
    return $help;
}

sub execute {
    my $self = shift;

    $self->_load_file;
    $self->_query;
    $self->_process;

    return 1;
}

sub _load_file {
    my $self = shift;

    my $file = $self->file;
    my ($dir, $basename, $ext) = File::Basename::fileparse($file, 'csv', 'tsv');
    die $self->error_message("Cannot determine type for file: %s. It needs to end with .csv or .tsv.", $file) if not $ext;
    my $sep_char = ( $ext eq 'csv' ? ',' : "\t" );
    my $parser = Text::CSV->new({
            sep_char => $sep_char,
            empty_is_undef => 1,
        });
    die $self->error_message('Failed to create Text::CSV parser!') if not $parser;
    $self->_incoming_sep_char($sep_char);

    my $fh = Genome::Sys->open_file_for_reading($file);
    my $incoming_headers = $parser->getline($fh);
    die $self->error_message('File (%s) does not have any headers! Is it empty?', $file) if not $incoming_headers;
    $parser->column_names($incoming_headers);
    $self->_incoming_headers($incoming_headers);

    if ( not List::MoreUtils::any { $_ eq 'analysis_id' } @$incoming_headers ) {
        die $self->error_message("No 'analysis_id' column found.");
    }

    my @incoming_data;
    while ( my $data = $parser->getline_hr($fh) ) {
        die $self->error_message(
            'No analysis_id found on line: %s', join($sep_char, map { $data->{$_} } @$incoming_headers)
        ) if not $data->{analysis_id};
        push @incoming_data, $data;
    }
    $self->_incoming_data(\@incoming_data);

    return 1;
}

sub _query {
    my $self = shift;

    my $metadata = Genome::Model::Tools::CgHub::Metadata->create;
    $self->_metadata($metadata);

    my @batches;
    my ($idx, $length) = (qw/ 0 0 /);
    for my $analysis_id ( map { $_->{analysis_id} } @{$self->_incoming_data} ) {
        push @{$batches[$idx]}, $analysis_id;
        $length += length($analysis_id);
        if ( $length > 800 ) {
            $idx++;
            $length = 0;
        }
    }

    for my $batch ( @batches ) {
        $self->_query_batch($batch);
    }

    return 1;
}

sub _query_batch {
    my ($self, $batch) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => ARRAYREF});

    my $query = 'analysis_id=('.join(' OR ', @$batch).')';
    my $query_cmd = Genome::Model::Tools::CgHub::Query->execute(
        query => $query,
    );
    die $self->error_message('Failed to execute query command!') if not $query_cmd->result;
    $self->_metadata->add_xml($query_cmd->metadata_xml);

    return 1;
}

sub _process {
    my $self = shift;

    my $output_fh = Genome::Sys->open_file_for_writing( $self->output_file );
    my $incoming_headers = $self->_incoming_headers;
    my @headers = sort { $a cmp $b } List::MoreUtils::uniq(@$incoming_headers, $self->tcga_attribute_names);
    my $sep_char = $self->_incoming_sep_char;
    $output_fh->print( join($sep_char, @headers)."\n" );

    for my $incoming_data ( @{$self->_incoming_data} ) {
        my $analysis_id = $incoming_data->{analysis_id};
        die $self->error_message(
            'No analysis_id found on line: %s', join($sep_char, map { $incoming_data->{$_} } @$incoming_headers)
        ) if not $analysis_id;
        my %data = ( %$incoming_data, $self->process_analysis_id($analysis_id) );
        $output_fh->print( join($sep_char, map { $data{$_} } @headers)."\n" );
    }
    $output_fh->close if $self->output_file ne '-';

    return 1;
}

sub process_analysis_id {
    my ($self, $analysis_id) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => SCALAR});

    my $legacy_sample_id = $self->_metadata->get_attribute_value($analysis_id, 'legacy_sample_id');
    if ( not $legacy_sample_id ) {
        die $self->error_message('No legacy sample id for %s!', $analysis_id);
    }

    my @tokens = split('-', $legacy_sample_id);

    my %attrs = (
        'individual.name' => join('-', @tokens[0..2]),
        'individual.taxon' => 'human',
        'individual.upn' => join('-', @tokens[1..2]),
        'instdata.analysis_id' => $analysis_id,
        'sample.name' => $legacy_sample_id,
        'sample.nomenclature' => 'TCGA',
        'sample.common_name' => $self->_resolve_common_name($legacy_sample_id),
        'sample.disease_abbr' => ( $self->_metadata->get_attribute_value($analysis_id, 'disease_abbr') || 'unknown' ),
        'sample.extraction_type' => $self->_resolve_extraction_type($legacy_sample_id),
        'sample.tissue_desc' => $self->_resolve_tissue_desc($analysis_id),
        'library.name' => $legacy_sample_id.'-extlibs',
    );


    return %attrs;
}

sub _resolve_common_name {
    my ($self, $analysis_id) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => SCALAR});

    my $sample_type = $self->_metadata->get_attribute_value($analysis_id, 'sample_type');
    if ( not $sample_type ) {
        die $self->error_message("No sample type for analysis id: $analysis_id");
    }

    if ( $sample_type >= 1 and $sample_type <= 9 ) {
        return 'tumor';
    }
    elsif ( $sample_type >= 10 and $sample_type <= 19 ) {
        return 'normal';
    }
    elsif ( $sample_type >= 20 and $sample_type <= 29 ) {
        return 'control';
    }
    elsif ( $sample_type >= 60 and $sample_type <= 69 ) {
        return 'xenograft';
    }

    return 'unknown';
}

sub _resolve_extraction_type {
    my ($self, $analysis_id) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => SCALAR});

    my $analyte_code = $self->_metadata->get_attribute_value($analysis_id, 'analyte_code');
    if ( not $analyte_code ) {
        die $self->error_message("No analyte code for analysis id: $analysis_id");
    }

    my $extraction_type = $self->analyte_codes_and_extraction_types->{$analyte_code};
    if ( not $extraction_type ) {
        die $self->error_message("No extraction type for analyte code: $analyte_code");
    }

    return $extraction_type;
}

sub _resolve_tissue_desc {
    my ($self, $analysis_id) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => SCALAR});

    my $tss_id = $self->_metadata->get_attribute_value($analysis_id, 'tss_id');
    return 'unknown' if not $tss_id;

    my $tissues_and_codes = tissues_and_codes();
    for my $tissue ( keys %$tissues_and_codes ) {
        return $tissue if List::MoreUtils::any { $tss_id eq $_ } @{$tissues_and_codes->{$tissue}};
    }

    return 'unknown';
}

sub analyte_codes_and_extraction_types {
    return {
        D => 'genomic dna',
        G => 'ipr product',
        R => 'rna',
        T => 'total rna',
        W => 'ipr product',
        X => 'ipr product',
    };
}

sub tissues_and_codes {
    return {
        'adrenal gland' => [qw/ OR OS OT OU P6 PA PF PK U6 /],
        'bile duct' => [qw/ 1X 3X 4G 5A W5 W6 W7 WD XZ YR ZD ZH ZK ZU /],
        'bladder' => [qw/ 2B 2F 4M 4Z 5N BL BT C4 CF CU DK E5 E7 FD FJ FT G2 GC GD GU GV H4 HQ J6 JK K4 KQ LC LT LU MV O4 PQ Q5 R3 S5 SY TB UY VU X5 XF YC YF ZF /],
        'blood' => [qw/ 3Y 5J AB DN G1 GP ID IJ IV JJ K2 M1 O7 PH Q4 TU Y3 /],
        'brain' => [qw/ 02 06 08 12 14 15 16 19 26 27 28 32 41 45 4W 65 74 76 81 84 87 BV CP CS DB DH DU E1 EZ F6 FG FN HK HT HW IK IX KT NV OA OX P5 QH R8 RR RY S9 TM TQ VM VV VW W9 WH WY YK /],
        'breast' => [qw/ 2R 3C 3J 4H 5L 5T A1 A2 A7 A8 AC AN AO AQ AR B6 BH C7 C8 D8 E2 E9 EW GI GM HN JL K3 LD LL LQ M2 MS OK OL PE PL Q6 S3 T4 TV UL UU V7 W8 WT XX YE YM Z7 ZZ /],
        'cervix' => [qw/ 1Y 2W 4J BI C5 DG DR DS EA EK EN EX FU GH HG HM IR IT JM JW JX K5 LE LP LS M3 MA MU MY O5 PN PW Q1 QV R2 RA T5 TD TW UC VS WL XS ZJ ZX /],
        'colon' => [qw/ 2S 3L 4N 4T 5M A6 AA AD AE AM AU AY AZ CA CK CM CY D5 DM F4 G4 GT NH Q8 QG QL R1 RU SK SS T9 VE VK WS WU Y2 /],
        'esophagus' => [qw/ 1U 2C 2H 2M IC IG JP JY JZ KA KH L5 L7 LF LN LW M9 NT NU Q9 QY R4 R6 RE RF S8 TF UG V5 VR X8 XB XP XW Z6 ZR /],
        'eye' => [qw/ 5H RZ V3 V4 VD VY WC YV YZ /],
        'head' => [qw/ 4P BA BB BE BY C9 CN CQ CR CV CX D6 DQ F7 H7 HD HL IQ IY JA KU MT MZ NM O3 P3 PV QK R7 RH RJ RS ST T2 T3 TG TN TY UF UP VJ WA /],
        'kidney' => [qw/ 1W 2K 2T 2Z 3Z 4A 5P 6C 6D A3 A4 AK AL AS AT B0 B1 B2 B3 B4 B8 B9 BP BQ CB CJ CW CZ DV DW DZ EU EV F8 F9 G6 G7 GK GL GY HE IA IZ J7 JB JC KL KM KN KO KV L8 LM MD MH MM MW N2 N3 N4 NP O9 P4 PJ PX Q2 QZ SX T7 T8 U9 UA UK UM UN UW UZ V8 V9 WM WN Y7 Y8 /],
        'liver' => [qw/ 2V 2Y 3K 4R 5C 5R 5Y BC BD BW CC DD ED EP ES FV G3 GJ HP JR K7 KR LG LV MI MR NI O8 PD PY QA RC RG SF T1 TH TZ UB UV WJ WQ WX XR YA ZP ZS ZY /],
        'lung' => [qw/ 03 05 11 17 18 21 22 33 34 35 37 38 39 3F 3H 3U 3V 40 43 44 46 48 49 4B 50 51 52 53 55 56 58 60 62 63 64 66 67 68 69 6A 6B 70 71 73 75 77 78 79 80 82 83 85 86 88 89 90 91 92 93 94 95 96 97 98 99 J1 J2 JD JE KW KX L3 L4 L9 LA LK ME MF ML MN MP MQ NB NC NJ NK NN NQ NW NX O1 O2 OB OC OZ PM PS S2 SC SH SZ T6 TS U1 U2 U7 UD UH UJ UT WG XC XT Y9 YS ZE ZN /],
        'lymph node' => [qw/ C2 FA FF FM G8 GR GS GZ H3 J3 JF JT KY M5 O6 OQ OV PB RQ VB /],
        'normal' => [qw/ 07 AV /],
        'ovary' => [qw/ 01 04 09 10 13 20 23 24 25 29 2Q 30 31 36 3P 42 47 4D 57 59 5E 5X 61 6F 72 OY QJ R9 VG WR Y4 /],
        'pancreas' => [qw/ 2D 2J 2L 2P 3A 3E 5Q 5Z F2 FB FQ FZ H6 H8 HV HZ IB JG KG L1 LB M8 NY OE ON P9 PU PZ Q3 RB RL RV S4 SD U3 US WF XD XN YB YH YY Z5 Z8 ZW /],
        'prostate' => [qw/ 2A 4L 4S C3 CH EJ FC G9 GX H9 HC HI J4 J9 JH KC KK L2 M7 MG OF QU SU TK TP UR UX V1 V2 VN VP WW X4 XA XJ XK XQ Y6 YJ YL ZG /],
        'rectum' => [qw/ 6G AF AG AH AI BM CI CL DC DT DY EF EG EH EI F5 G5 QP TA VL /],
        'skin' => [qw/ 1T 3N BF D3 D9 DA EB EE ER FR FS FW GF GN HR IH JQ K8 LH M4 NS OD QB RP S1 TE TR W3 WE WV XV YD YG YN Z2 /],
        'stomach' => [qw/ 1V 2N 3D 3M B7 BR BX CD CG D7 EQ F1 FP GW H1 HA HF HH HJ HU IN IP J5 JI KB KZ MX OH OP QW QX R5 RD RK SM SW TL V6 VA VC VQ VX YX ZA ZQ /],
        'testes' => [qw/ 2G 2X 4K S6 SB SN SO U8 VF W4 WZ X3 XE XL XY YU ZM /],
        'thymus' => [qw/ 1Z 3G 3Q 3S 3T 4V 4X 5G 5K 5U 5V X7 XH XM XU YT ZB ZC ZL ZT /],
        'thyroid' => [qw/ 2U 4C 5F BJ BZ CE DE DJ DO E3 E8 EL EM ET FE FH FK FY GE H2 IM J8 JS K9 KS L6 LJ LR LZ MC MK OJ OM QD U5 UQ /],
        'unknown' => [qw/ 3B 3R 3W 4Q 4Y 5D DX FX HB HS HY IE IF IS IU IW JN JV K1 KD KF LI LY M6 MB MJ MO N1 OG OW P7 P8 PC PR PT QC QQ QR QT RM RN RT RW RX S7 SA SE SG SI SP SQ SR SV TJ TT U4 UE VH VT VZ W2 WB WK WP X2 X6 X9 XG Y5 YW Z3 Z4 Z9 /],
        'uterus' => [qw/ 2E 4E 4F 4U 5B 5S 5W 6E A5 AJ AP AW AX B5 BG BK BS D1 DF DI E6 EC EO EY FI FL GG H5 JO JU K6 KE KJ KP LX N5 N6 N7 N8 N9 NA ND NE NF NG NL NR NZ P1 P2 PG Q7 QE QF QM QN QS SJ SL TX /],
    };
}

1;

