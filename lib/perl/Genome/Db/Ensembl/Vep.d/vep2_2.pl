#!/usr/bin/env genome-perl

=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Variant Effect Predictor - a script to predict the consequences of genomic variants

http://www.ensembl.org/info/docs/variation/vep/vep_script.html

Version 2.2

by Will McLaren (wm2@ebi.ac.uk)
=cut

use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use Bio::EnsEMBL::Variation::Utils::VEP qw(
    parse_line
    vf_to_consequences
    validate_vf
    load_dumped_adaptor_cache
    dump_adaptor_cache
    get_all_consequences
    get_slice
    build_full_cache
    read_cache_info
    get_time
    debug
    @OUTPUT_COLS
    @REG_FEAT_TYPES
);

# global vars
my $VERSION = '2.2';

# set output autoflush for progress bars
$| = 1;

# configure from command line opts
my $config = &configure(scalar @ARGV);

# run the main sub routine
&main($config);

# this is the main sub-routine - it needs the configured $config hash
sub main {
    my $config = shift;

    debug("Starting...") unless defined $config->{quiet};

    my $tr_cache = {};
    my $rf_cache = {};

    # create a hash to hold slices so we don't get the same one twice
    my %slice_cache = ();

    my @vfs;
    my ($vf_count, $total_vf_count);
    my $in_file_handle = $config->{in_file_handle};

    # initialize line number in config
    $config->{line_number} = 0;

    # read the file
    while(<$in_file_handle>) {
        chomp;

        $config->{line_number}++;

        # header line?
        next if /^\#/;

        # some lines (pileup) may actually parse out into more than one variant
        foreach my $vf(@{&parse_line($config, $_)}) {

            # validate the VF
            next unless validate_vf($config, $vf);

            # now get the slice
            if(!defined($vf->{slice})) {
                my $slice;

                # don't get slices if we're using cache
                # we can steal them from transcript objects later
                if((!defined($config->{cache}) && !defined($config->{whole_genome})) || defined($config->{check_ref}) || defined($config->{convert})) {

                    # check if we have fetched this slice already
                    if(defined $slice_cache{$vf->{chr}}) {
                        $slice = $slice_cache{$vf->{chr}};
                    }

                    # if not create a new one
                    else {

                        $slice = &get_slice($config, $vf->{chr});

                        # if failed, warn and skip this line
                        if(!defined($slice)) {
                            warn("WARNING: Could not fetch slice named ".$vf->{chr}." on line ".$config->{line_number}."\n") unless defined $config->{quiet};
                            next;
                        }

                        # store the hash
                        $slice_cache{$vf->{chr}} = $slice;
                    }
                }

                $vf->{slice} = $slice;
            }

            # make a name if one doesn't exist
            $vf->{variation_name} ||= $vf->{chr}.'_'.$vf->{start}.'_'.$vf->{allele_string};

            # jump out to convert here
            if(defined($config->{convert})) {
                &convert_vf($config, $vf);
                next;
            }

            if(defined $config->{whole_genome}) {
                push @vfs, $vf;
                $vf_count++;
                $total_vf_count++;

                if($vf_count == $config->{buffer_size}) {
                    debug("Read $vf_count variants into buffer") unless defined($config->{quiet});

                    print_line($config, $_) foreach @{get_all_consequences($config, \@vfs, $tr_cache, $rf_cache)};

                    debug("Processed $total_vf_count total variants") unless defined($config->{quiet});

                    @vfs = ();
                    $vf_count = 0;
                }
            }
            else {
                print_line($config, $_) foreach @{vf_to_consequences($config, $vf)};
                $vf_count++;
                $total_vf_count++;
                debug("Processed $vf_count variants") if $vf_count =~ /0$/ && defined($config->{verbose});
            }
        }
    }

    # if in whole-genome mode, finish off the rest of the buffer
    if(defined $config->{whole_genome} && scalar @vfs) {
        debug("Read $vf_count variants into buffer") unless defined($config->{quiet});

        print_line($config, $_) foreach @{get_all_consequences($config, \@vfs, $tr_cache, $rf_cache)};

        debug("Processed $total_vf_count total variants") unless defined($config->{quiet});
    }

    debug("Executed ", defined($Bio::EnsEMBL::DBSQL::StatementHandle::count_queries) ? $Bio::EnsEMBL::DBSQL::StatementHandle::count_queries : 'unknown number of', " SQL statements") if defined($config->{count_queries}) && !defined($config->{quiet});

    debug("Finished!") unless defined $config->{quiet};
}

# sets up configuration hash that is used throughout the script
sub configure {
    my $args = shift;

    my $config = {};

    GetOptions(
        $config,
        'help',                    # displays help message

        # input options,
        'config=s',                # config file name
        'input_file=s',            # input file name
        'format=s',                # input file format

        # DB options
        'species=s',               # species e.g. human, homo_sapiens
        'registry=s',              # registry file
        'host=s',                  # database host
        'port=s',                  # database port
        'user=s',                  # database user name
        'password=s',              # database password
        'db_version=i',            # Ensembl database version to use e.g. 62
        'genomes',                 # automatically sets DB params for e!Genomes
        'refseq',                  # use otherfeatures RefSeq DB instead of Ensembl
        #'no_disconnect',           # disables disconnect_when_inactive

        # runtime options
        'most_severe',             # only return most severe consequence
        'summary',                 # only return one line per variation with all consquence types
        'per_gene',                # only return most severe per gene
        'buffer_size=i',           # number of variations to read in before analysis
        'chunk_size=s',            # size in bases of "chunks" used in internal hash structure
        'failed=i',                # include failed variations when finding existing
        'no_whole_genome',         # disables now default whole-genome mode
        'whole_genome',            # proxy for whole genome mode - now just warns user
        'gp',                      # read coords from GP part of INFO column in VCF (probably only relevant to 1KG)
        'chr=s',                   # analyse only these chromosomes, e.g. 1-5,10,MT
        'check_ref',               # check supplied reference allele against DB
        'check_existing',          # find existing co-located variations
        'check_alleles',           # only attribute co-located if alleles are the same
        'check_frequency',         # enable frequency checking
        'freq_filter=s',           # exclude or include
        'freq_freq=f',             # frequency to filter on
        'freq_gt_lt=s',            # gt or lt (greater than or less than)
        'freq_pop=s',              # population to filter on

        # verbosity options
        'verbose',                 # print out a bit more info while running
        'quiet',                   # print nothing to STDOUT (unless using -o stdout)
        'no_progress',             # don't display progress bars

        # output options
        'output_file=s',           # output file name
        'force_overwrite',         # force overwrite of output file if already exists
        'terms=s',                 # consequence terms to use e.g. NCBI, SO
        'coding_only',             # only return results for consequences in coding regions
        'canonical',               # indicates if transcript is canonical
        'protein',                 # add e! protein ID to extra column
        'hgnc',                    # add HGNC gene ID to extra column
        'hgvs',                    # add HGVS names to extra column
        'sift=s',                  # SIFT predictions
        'polyphen=s',              # PolyPhen predictions
        'condel=s',                # Condel predictions
        'gene',                    # force gene column to be populated (disabled by default, enabled when using cache)
        'regulatory',              # enable regulatory stuff
        'convert=s',               # convert input to another format (doesn't run VEP)
        'no_intergenic',           # don't print out INTERGENIC consequences

        # cache stuff
        'cache',                   # use cache
        'write_cache',             # enables writing to the cache
        'build=s',                 # builds cache from DB from scratch; arg is either all (all top-level seqs) or a list of chrs
        'prefetch',                # prefetch exons, translation, introns, codon table etc for each transcript
        'strip',                   # strips adaptors etc from objects before caching them
        'rebuild=s',               # rebuilds cache by reading in existing then redumping - probably don't need to use this any more
        'dir=s',                   # dir where cache is found (defaults to $HOME/.vep/)
        'cache_region_size=i',     # size of region in bases for each cache file
        'no_slice_cache',          # tell API not to cache features on slice
        'standalone',              # standalone mode uses minimal set of modules installed in same dir, no DB connection
        'skip_db_check',           # don't compare DB parameters with cached
        'compress=s',              # by default we use zcat to decompress; user may want to specify gzcat or "gzip -dc"

        # debug
        'cluck',                   # these two need some mods to Bio::EnsEMBL::DBSQL::StatementHandle to work. Clucks callback trace and SQL
        'count_queries',           # counts SQL queries executed
        'admin',                   # allows me to build off public hosts
        'debug',                   # print out debug info
    );

    # print usage message if requested or no args supplied
    if(defined($config->{help}) || !$args) {
        &usage;
        exit(0);
    }

    # config file?
    if(defined $config->{config}) {

        open CONFIG, $config->{config} or die "ERROR: Could not open config file \"".$config->{config}."\"\n";

        while(<CONFIG>) {
            next if /^\#/;
            my ($key, $value) = split /\s+|\=/;
            $key =~ s/^\-//g;
            $config->{$key} = $value unless defined $config->{$key};
        }

        close CONFIG;
    }

    # can't be both quiet and verbose
    die "ERROR: Can't be both quiet and verbose!" if defined($config->{quiet}) && defined($config->{verbose});

    # check file format
    if(defined $config->{format}) {
        die "ERROR: Unrecognised input format specified \"".$config->{format}."\"\n" unless $config->{format} =~ /pileup|vcf|guess|hgvs|ensembl|id/i;
    }

    # check convert format
    if(defined $config->{convert}) {
        die "ERROR: Unrecognised output format for conversion specified \"".$config->{convert}."\"\n" unless $config->{convert} =~ /vcf|ensembl|pileup/i;
    }

    # connection settings for Ensembl Genomes
    if($config->{genomes}) {
        $config->{host} ||= 'mysql.ebi.ac.uk';
        $config->{port} ||= 4157;
    }

    # connection settings for main Ensembl
    else {
        $config->{species} ||= "homo_sapiens";
        $config->{host}    ||= 'ensembldb.ensembl.org';
        $config->{port}    ||= 5306;
    }

    # refseq or core?
    if(defined($config->{refseq})) {
        die "ERROR: SIFT, PolyPhen and Condel predictions not available fore RefSeq transcripts" if defined $config->{sift} || defined $config->{polyphen} || defined $config->{condel};

        $config->{core_type} = 'otherfeatures';
    }
    else {
        $config->{core_type} = 'core';
    }

    # output term
    if(defined $config->{terms}) {
        die "ERROR: Unrecognised consequence term type specified \"".$config->{terms}."\" - must be one of ensembl, so, ncbi\n" unless $config->{terms} =~ /ensembl|display|so|ncbi/i;
        if($config->{terms} =~ /ensembl|display/i) {
            $config->{terms} = 'display';
        }
        else {
            $config->{terms} = uc($config->{terms});
        }
    }

    # check nsSNP tools
    foreach my $tool(grep {defined $config->{lc($_)}} qw(SIFT PolyPhen Condel)) {
        die "ERROR: Unrecognised option for $tool \"", $config->{lc($tool)}, "\" - must be one of p (prediction), s (score) or b (both)\n" unless $config->{lc($tool)} =~ /^(s|p|b)/;

        die "ERROR: $tool not available for this species\n" unless $config->{species} =~ /human|homo/i;

        die "ERROR: $tool not available in standalone mode\n" if defined($config->{standalone});

        # use V2 of the Condel algorithm, possibly gives fewer false positives
        if($tool eq 'Condel' && $config->{lc($tool)} =~ /1$/) {
            $Bio::EnsEMBL::Variation::Utils::Condel::USE_V2 = 0;
        }
    }

    # force quiet if outputting to STDOUT
    if(defined($config->{output_file}) && $config->{output_file} =~ /stdout/i) {
        delete $config->{verbose} if defined($config->{verbose});
        $config->{quiet} = 1;
    }

    # summarise options if verbose
    if(defined $config->{verbose}) {
        my $header =<<INTRO;
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

version $VERSION

By Will McLaren (wm2\@ebi.ac.uk)

Configuration options:

INTRO
        print $header;

        my $max_length = (sort {$a <=> $b} map {length($_)} keys %$config)[-1];

        foreach my $key(sort keys %$config) {
            print $key.(' ' x (($max_length - length($key)) + 4)).$config->{$key}."\n";
        }

        print "\n".("-" x 20)."\n\n";
    }

    # set defaults
    $config->{user}              ||= 'anonymous';
    $config->{buffer_size}       ||= 5000;
    $config->{chunk_size}        ||= '50kb';
    $config->{output_file}       ||= "variant_effect_output.txt";
    $config->{tmpdir}            ||= '/tmp';
    $config->{format}            ||= 'guess';
    $config->{terms}             ||= 'display';
    $config->{gene}              ||= 1 unless defined($config->{whole_genome});
    $config->{cache_region_size} ||= 1000000;
    $config->{dir}               ||= join '/', ($ENV{'HOME'}, '.vep');
    $config->{compress}          ||= 'zcat';

    # frequency filtering
    if(defined($config->{check_frequency})) {
        foreach my $flag(qw(freq_freq freq_filter freq_pop freq_gt_lt)) {
            die "ERROR: To use --check_frequency you must also specify flag --$flag" unless defined $config->{$flag};
        }

        # need to set check_existing
        $config->{check_existing} = 1;
    }

    $config->{check_existing} = 1 if defined $config->{check_alleles};

    # warn users still using whole_genome flag
    if(defined($config->{whole_genome})) {
        debug("INFO: Whole-genome mode is now the default run-mode for the script. To disable it, use --no_whole_genome") unless defined($config->{quiet});
    }

    $config->{whole_genome}      = 1 unless defined $config->{no_whole_genome};
    $config->{include_failed}    = 1 unless defined $config->{include_failed};
    $config->{chunk_size}        =~ s/mb?/000000/i;
    $config->{chunk_size}        =~ s/kb?/000/i;
    $config->{cache_region_size} =~ s/mb?/000000/i;
    $config->{cache_region_size} =~ s/kb?/000/i;

    # cluck and display executed SQL?
    $Bio::EnsEMBL::DBSQL::StatementHandle::cluck = 1 if defined($config->{cluck});

    # standalone needs cache, can't use HGVS
    if(defined($config->{standalone})) {
        $config->{cache} = 1;

        die("ERROR: Cannot generate HGVS coordinates in standalone mode") if defined($config->{hgvs});
        die("ERROR: Cannot use HGVS as input in standalone mode") if $config->{format} eq 'hgvs';
        die("ERROR: Cannot use variant identifiers as input in standalone mode") if $config->{format} eq 'id';
        die("ERROR: Cannot do frequency filtering in standalone mode") if defined($config->{check_frequency});
    }

    # write_cache needs cache
    $config->{cache} = 1 if defined $config->{write_cache};

    # no_slice_cache, prefetch and whole_genome have to be on to use cache
    if(defined($config->{cache})) {
        $config->{prefetch} = 1;
        $config->{no_slice_cache} = 1;
        $config->{whole_genome} = 1;
        $config->{strip} = 1;
    }

    $config->{build} = $config->{rebuild} if defined($config->{rebuild});

    # force options for full build
    if(defined($config->{build})) {
        $config->{prefetch} = 1;
        $config->{gene} = 1;
        $config->{hgnc} = 1;
        $config->{no_slice_cache} = 1;
        $config->{cache} = 1;
        $config->{strip} = 1;
        $config->{write_cache} = 1;
    }

    # connect to databases
    $config->{reg} = &connect_to_dbs($config);

    # complete dir with species name and db_version
    $config->{dir} .= '/'.(
        join '/', (
            defined($config->{standalone}) ? $config->{species} : ($config->{reg}->get_alias($config->{species}) || $config->{species}),
            $config->{db_version} || $config->{reg}->software_version
        )
    );

    if(defined($config->{cache})) {
        # read cache info
        if(read_cache_info($config)) {
            debug("Read existing cache info") unless defined $config->{quiet};
        }
    }

    # include regulatory modules if requested
    if(defined($config->{regulatory})) {
        # do the use statements here so that users don't have to have the
        # funcgen API install to use the rest of the script
        use Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryFeatureAdaptor;
        use Bio::EnsEMBL::Funcgen::DBSQL::MotifFeatureAdaptor;
        use Bio::EnsEMBL::Funcgen::MotifFeature;
        use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
        use Bio::EnsEMBL::Funcgen::BindingMatrix;
    }

    # warn user cache directory doesn't exist
    if(!-e $config->{dir}) {

        # if using write_cache
        if(defined($config->{write_cache})) {
            debug("INFO: Cache directory ", $config->{dir}, " not found - it will be created") unless defined($config->{quiet});
        }

        # want to read cache, not found
        elsif(defined($config->{cache})) {
            die("ERROR: Cache directory ", $config->{dir}, " not found");
        }
    }

    # suppress warnings that the FeatureAdpators spit if using no_slice_cache
    Bio::EnsEMBL::Utils::Exception::verbose(1999) if defined($config->{no_slice_cache});

    # get adaptors
    if(defined($config->{cache}) && !defined($config->{write_cache})) {

        # try and load adaptors from cache
        if(!&load_dumped_adaptor_cache($config)) {
            &get_adaptors($config);
            &dump_adaptor_cache($config) if defined($config->{write_cache});
        }

        # check cached adaptors match DB params
        else {
            my $dbc = $config->{sa}->{dbc};

            my $ok = 1;

            if($dbc->{_host} ne $config->{host}) {

                # ens-livemirror, useastdb and ensembldb should all have identical DBs
                unless(
                    (
                        $dbc->{_host} eq 'ens-livemirror'
                        || $dbc->{_host} eq 'ensembldb.ensembl.org'
                        || $dbc->{_host} eq 'useastdb.ensembl.org'
                    ) && (
                        $config->{host} eq 'ens-livemirror'
                        || $config->{host} eq 'ensembldb.ensembl.org'
                        || $config->{host} eq 'useastdb.ensembl.org'
                    )
                ) {
                    $ok = 0;
                }

                # but we still need to reconnect
                debug("INFO: Defined host ", $config->{host}, " is different from cached ", $dbc->{_host}, " - reconnecting to host") unless defined($config->{quiet});

                &get_adaptors($config);
            }

            if(!$ok) {
                if(defined($config->{skip_db_check})) {
                    debug("INFO: Defined host ", $config->{host}, " is different from cached ", $dbc->{_host}) unless defined($config->{quiet});
                }
                else {
                    die "ERROR: Defined host ", $config->{host}, " is different from cached ", $dbc->{_host}, ". If you are sure this is OK, rerun with -skip_db_check flag set";
                }
            }
        }
    }
    else {
        &get_adaptors($config);
        &dump_adaptor_cache($config) if defined($config->{write_cache})
    }

    # reg adaptors (only fetches if not retrieved from cache already)
    &get_reg_adaptors($config) if defined($config->{regulatory});

    # get terminal width for progress bars
    unless(defined($config->{quiet})) {
        my $width;

        # module may not be installed
        eval {
            use Term::ReadKey;
        };

        if(!$@) {
            my ($w, $h);

            # module may be installed, but e.g.
            eval {
                ($w, $h) = GetTerminalSize();
            };

            $width = $w if defined $w;
        }

        $width ||= 60;
        $width -= 12;
        $config->{terminal_width} = $width;
    }

    # jump out to build cache if requested
    if(defined($config->{build})) {

        if($config->{host} =~ /^(ensembl|useast)db\.ensembl\.org$/ && !defined($config->{admin})) {
            die("ERROR: Cannot build cache using public database server ", $config->{host}, "\n");
        }

        # build the cache
        debug("Building cache for ".$config->{species}) unless defined($config->{quiet});
        build_full_cache($config);

        # exit script
        debug("Finished building cache") unless defined($config->{quiet});
        exit(0);
    }

    # warn user DB will be used for SIFT/PolyPhen/Condel/HGVS/frequency
    if(defined($config->{cache})) {

        # these two def depend on DB
        foreach my $param(grep {defined $config->{$_}} qw(hgvs check_frequency)) {
            debug("INFO: Database will be accessed when using --$param") unless defined($config->{quiet});
        }

        # as does using HGVS or IDs as input
        debug("INFO: Database will be accessed when using --format ", $config->{format}) if ($config->{format} eq 'id' || $config->{format} eq 'hgvs') && !defined($config->{quiet});

        # the rest may be in the cache
        foreach my $param(grep {defined $config->{$_}} qw(sift polyphen condel regulatory)) {
            next if defined($config->{'cache_'.$param});
            debug("INFO: Database will be accessed when using --$param; consider using the complete cache containing $param data (see documentation for details)") unless defined($config->{quiet});
        }
    }

    # get list of chrs if supplied
    if(defined($config->{chr})) {
        my %chrs;

        foreach my $val(split /\,/, $config->{chr}) {
            my @nnn = split /\-/, $val;

            foreach my $chr($nnn[0]..$nnn[-1]) {
                $chrs{$chr} = 1;
            }
        }

        $config->{chr} = \%chrs;
    }

    # get input file handle
    $config->{in_file_handle} = &get_in_file_handle($config);

    # configure output file
    $config->{out_file_handle} = &get_out_file_handle($config);

    return $config;
}

# connects to DBs; in standalone mode this just loads registry module
sub connect_to_dbs {
    my $config = shift;

    # get registry
    my $reg = 'Bio::EnsEMBL::Registry';

    unless(defined($config->{standalone})) {
        # load DB options from registry file if given
        if(defined($config->{registry})) {
            debug("Loading DB config from registry file ", $config->{registry}) unless defined($config->{quiet});
            $reg->load_all(
                $config->{registry},
                $config->{verbose},
                undef,
                $config->{no_slice_cache}
            );
        }

        # otherwise manually connect to DB server
        else {
            $reg->load_registry_from_db(
                -host       => $config->{host},
                -user       => $config->{user},
                -pass       => $config->{password},
                -port       => $config->{port},
                -db_version => $config->{db_version},
                -species    => $config->{species} =~ /^[a-z]+\_[a-z]+/i ? $config->{species} : undef,
                -verbose    => $config->{verbose},
                -no_cache   => $config->{no_slice_cache},
            );
        }

        eval { $reg->set_reconnect_when_lost() };

        if(defined($config->{verbose})) {
            # get a meta container adaptors to check version
            my $core_mca = $reg->get_adaptor($config->{species}, 'core', 'metacontainer');
            my $var_mca = $reg->get_adaptor($config->{species}, 'variation', 'metacontainer');

            if($core_mca && $var_mca) {
                debug(
                    "Connected to core version ", $core_mca->get_schema_version, " database ",
                    "and variation version ", $var_mca->get_schema_version, " database"
                );
            }
        }
    }

    return $reg;
}

# get adaptors from DB
sub get_adaptors {
    my $config = shift;

    die "ERROR: No registry" unless defined $config->{reg};

    $config->{vfa}   = $config->{reg}->get_adaptor($config->{species}, 'variation', 'variationfeature');
    $config->{tva}   = $config->{reg}->get_adaptor($config->{species}, 'variation', 'transcriptvariation');
    $config->{pfpma} = $config->{reg}->get_adaptor($config->{species}, 'variation', 'proteinfunctionpredictionmatrix');
    $config->{va}    = $config->{reg}->get_adaptor($config->{species}, 'variation', 'variation');

    # get fake ones for species with no var DB
    if(!defined($config->{vfa})) {
        $config->{vfa} = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake($config->{species});
        $config->{tva} = Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor->new_fake($config->{species});
    }
    else {
        $config->{vfa}->db->include_failed_variations($config->{include_failed}) if defined($config->{vfa}->db) && $config->{vfa}->db->can('include_failed_variations');
    }

    $config->{sa}  = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'slice');
    $config->{ga}  = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'gene');
    $config->{ta}  = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'transcript');
    $config->{mca} = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'metacontainer');
    $config->{csa} = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'coordsystem');

    # cache schema version
    $config->{mca}->get_schema_version if defined $config->{mca};

    # check we got slice adaptor - can't continue without a core DB
    die("ERROR: Could not connect to core database\n") unless defined $config->{sa};
}

# gets regulatory adaptors
sub get_reg_adaptors {
    my $config = shift;

    foreach my $type(@REG_FEAT_TYPES) {
        next if defined($config->{$type.'_adaptor'});

        my $adaptor = $config->{reg}->get_adaptor($config->{species}, 'funcgen', $type);
        if(defined($adaptor)) {
            $config->{$type.'_adaptor'} = $adaptor;
        }
        else {
            delete $config->{regulatory};
            last;
        }
    }
}

# gets file handle for input
sub get_in_file_handle {
    my $config = shift;

    # define the filehandle to read input from
    my $in_file_handle = new FileHandle;

    if(defined($config->{input_file})) {

        # check defined input file exists
        die("ERROR: Could not find input file ", $config->{input_file}, "\n") unless -e $config->{input_file};

        if($config->{input_file} =~ /\.gz$/){
            $in_file_handle->open($config->{compress}." ". $config->{input_file} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
        }
        else {
            $in_file_handle->open( $config->{input_file} ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
        }
    }

    # no file specified - try to read data off command line
    else {
        $in_file_handle = 'STDIN';
        debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...") unless defined $config->{quiet};
    }

    return $in_file_handle;
}

# gets file handle for output and adds header
sub get_out_file_handle {
    my $config = shift;

    # define filehandle to write to
    my $out_file_handle = new FileHandle;

    # check if file exists
    if(-e $config->{output_file} && !defined($config->{force_overwrite})) {
        die("ERROR: Output file ", $config->{output_file}, " already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite\n");
    }

    if($config->{output_file} =~ /stdout/i) {
        $out_file_handle = *STDOUT;
    }
    else {
        $out_file_handle->open(">".$config->{output_file}) or die("ERROR: Could not write to output file ", $config->{output_file}, "\n");
    }

    # file conversion, don't want to add normal headers
    if(defined($config->{convert})) {
        # header for VCF
        if($config->{convert} =~ /vcf/i) {
            print $out_file_handle join "\t", (
                '#CHROM',
                'POS',
                'ID',
                'REF',
                'ALT',
                'QUAL',
                'FILTER',
                'INFO'
            );
            print $out_file_handle "\n";
        }

        return $out_file_handle;
    }

    # make header
    my $time = &get_time;
    my $db_string = $config->{mca}->dbc->dbname." on ".$config->{mca}->dbc->host if defined $config->{mca};
    $db_string .= "\n## Using cache in ".$config->{dir} if defined($config->{cache});
    my $version_string =
        "Using API version ".$config->{reg}->software_version.
        ", DB version ".(defined $config->{mca} && $config->{mca}->get_schema_version ? $config->{mca}->get_schema_version : '?');

    my $header =<<HEAD;
## ENSEMBL VARIANT EFFECT PREDICTOR v$VERSION
## Output produced at $time
## Connected to $db_string
## $version_string
## Extra column keys:
## CANONICAL    : Indicates if transcript is canonical for this gene
## HGNC         : HGNC gene identifier
## ENSP         : Ensembl protein identifer
## HGVSc        : HGVS coding sequence name
## HGVSp        : HGVS protein sequence name
## SIFT         : SIFT prediction
## PolyPhen     : PolyPhen prediction
## Condel       : Condel SIFT/PolyPhen consensus prediction
## MATRIX       : The source and identifier of a transcription factor binding profile aligned at this position
## HIGH_INF_POS : A flag indicating if the variant falls in a high information position of a transcription factor binding profile
HEAD

    # add headers
    print $out_file_handle $header;

    # add column headers
    print $out_file_handle '#', (join "\t", @OUTPUT_COLS);
    print $out_file_handle "\n";

    return $out_file_handle;
}

# convert a variation feature to a line of output
sub convert_vf {
    my $config = shift;
    my $vf = shift;

    my $convert_method = 'convert_to_'.lc($config->{convert});
    my $method_ref   = \&$convert_method;

    my $line = &$method_ref($config, $vf);
    my $handle = $config->{out_file_handle};

    if(scalar @$line) {
        print $handle join "\t", @$line;
        print $handle "\n";
    }
}

# converts to Ensembl format
sub convert_to_ensembl {
    my $config = shift;
    my $vf = shift;

    return [
        $vf->{chr} || $vf->seq_region_name,
        $vf->start,
        $vf->end,
        $vf->allele_string,
        $vf->strand,
        $vf->variation_name
    ];
}

# converts to VCF format
sub convert_to_vcf {
    my $config = shift;
    my $vf = shift;

    # look for imbalance in the allele string
    my %allele_lengths;
    my @alleles = split /\//, $vf->allele_string;

    foreach my $allele(@alleles) {
        $allele =~ s/\-//g;
        $allele_lengths{length($allele)} = 1;
    }

    # in/del/unbalanced
    if(scalar keys %allele_lengths > 1) {

        # we need the ref base before the variation
        # default to N in case we can't get it
        my $prev_base = 'N';

        unless(defined($config->{cache})) {
            my $slice = $vf->slice->sub_Slice($vf->start - 1, $vf->start -1);
            $prev_base = $slice->seq if defined($slice);
        }

        for my $i(0..$#alleles) {
            $alleles[$i] =~ s/\-//g;
            $alleles[$i] = $prev_base.$alleles[$i];
        }

        return [
            $vf->{chr} || $vf->seq_region_name,
            $vf->start - 1,
            $vf->variation_name,
            shift @alleles,
            (join ",", @alleles),
            '.', '.', '.'
        ];

    }

    # balanced sub
    else {
        return [
            $vf->{chr} || $vf->seq_region_name,
            $vf->start,
            $vf->variation_name,
            shift @alleles,
            (join ",", @alleles),
            '.', '.', '.'
        ];
    }
}


# converts to pileup format
sub convert_to_pileup {
    my $config = shift;
    my $vf = shift;

    # look for imbalance in the allele string
    my %allele_lengths;
    my @alleles = split /\//, $vf->allele_string;

    foreach my $allele(@alleles) {
        $allele =~ s/\-//g;
        $allele_lengths{length($allele)} = 1;
    }

    # in/del
    if(scalar keys %allele_lengths > 1) {

        if($vf->allele_string =~ /\-/) {

            # insertion?
            if($alleles[0] eq '-') {
                shift @alleles;

                for my $i(0..$#alleles) {
                    $alleles[$i] =~ s/\-//g;
                    $alleles[$i] = '+'.$alleles[$i];
                }
            }

            else {
                @alleles = grep {$_ ne '-'} @alleles;

                for my $i(0..$#alleles) {
                    $alleles[$i] =~ s/\-//g;
                    $alleles[$i] = '-'.$alleles[$i];
                }
            }

            @alleles = grep {$_ ne '-' && $_ ne '+'} @alleles;

            return [
                $vf->{chr} || $vf->seq_region_name,
                $vf->start - 1,
                '*',
                (join "/", @alleles),
            ];
        }

        else {
            warn "WARNING: Unable to convert variant to pileup format on line number ", $config->{line_number} unless defined($config->{quiet});
            return [];
        }

    }

    # balanced sub
    else {
        return [
            $vf->{chr} || $vf->seq_region_name,
            $vf->start,
            shift @alleles,
            (join "/", @alleles),
        ];
    }
}

# prints a line of output from the hash
sub print_line {
    my $config = shift;
    my $line = shift;
    return unless defined($line);

    $line->{Extra} = join ';', map { $_.'='.$line->{Extra}->{$_} } keys %{ $line->{Extra} || {} };

    my $output = join "\t", map { $line->{$_} || '-' } @OUTPUT_COLS;

    my $fh = $config->{out_file_handle};

    print $fh "$output\n";
}

# outputs usage message
sub usage {
    my $usage =<<END;
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

version $VERSION

By Will McLaren (wm2\@ebi.ac.uk)

http://www.ensembl.org/info/docs/variation/vep/vep_script.html

Usage:
perl variant_effect_predictor.pl [arguments]

Options
=======

--help                 Display this message and quit
--verbose              Display verbose output as the script runs [default: off]
--quiet                Suppress status and warning messages [default: off]
--no_progress          Suppress progress bars [default: off]

--config               Load configuration from file. Any command line options
                       specified overwrite those in the file [default: off]

-i | --input_file      Input file - if not specified, reads from STDIN. Files
                       may be gzip compressed.
--format               Specify input file format - one of "ensembl", "pileup",
                       "vcf", "hgvs", "id" or "guess" to try and work out format.
-o | --output_file     Output file. Write to STDOUT by specifying -o STDOUT - this
                       will force --quiet [default: "variant_effect_output.txt"]
--force_overwrite      Force overwriting of output file [default: quit if file
                       exists]

--species [species]    Species to use [default: "human"]

-t | --terms           Type of consequence terms to output - one of "ensembl", "SO",
                       "NCBI" [default: ensembl]

--sift=[p|s|b]         Add SIFT [p]rediction, [s]core or [b]oth [default: off]
--polyphen=[p|s|b]     Add PolyPhen [p]rediction, [s]core or [b]oth [default: off]
--condel=[p|s|b]       Add Condel SIFT/PolyPhen consensus [p]rediction, [s]core or
                       [b]oth. Add 1 (e.g. b1) to option to use old Condel algorithm
                       [default: off]

NB: SIFT, PolyPhen and Condel predictions are currently available for human only

--regulatory           Look for overlaps with regulatory regions. The script can
                       also call if a variant falls in a high information position
                       within a transcription factor binding site. Output lines have
                       a Feature type of RegulatoryFeature or MotifFeature
                       [default: off]

NB: Regulatory consequences are currently available for human and mouse only

--hgnc                 If specified, HGNC gene identifiers are output alongside the
                       Ensembl Gene identifier [default: off]
--hgvs                 Output HGVS identifiers (coding and protein). Requires database
                       connection [default: off]
--protein              Output Ensembl protein identifer [default: off]
--gene                 Force output of Ensembl gene identifer - disabled by default
                       unless using --cache or --no_whole_genome [default: off]
--canonical            Indicate if the transcript for this consequence is the canonical
                       transcript for this gene [default: off]

--coding_only          Only return consequences that fall in the coding region of
                       transcripts [default: off]
--most_severe          Ouptut only the most severe consequence per variation.
                       Transcript-specific columns will be left blank. [default: off]
--summary              Output only a comma-separated list of all consequences per
                       variation. Transcript-specific columns will be left blank.
                       [default: off]
--per_gene             Output only the most severe consequence per gene. Where more
                       than one transcript has the same consequence, the transcript
                       chosen is arbitrary. [default: off]


--check_ref            If specified, checks supplied reference allele against stored
                       entry in Ensembl Core database [default: off]
--check_existing       If specified, checks for existing co-located variations in the
                       Ensembl Variation database [default: off]
--check_alleles        If specified, the alleles of existing co-located variations
                       are compared to the input; an existing variation will only
                       be reported if no novel allele is in the input (strand is
                       accounted for) [default: off]

--no_intergenic        Excludes intergenic consequences from the output [default: off]

--check_frequency      Turns on frequency filtering. Use this to include or exclude
                       variants based on the frequency of co-located existing
                       variants in the Ensembl Variation database. You must also
                       specify all of the following --freq flags [default: off]
--freq_pop [pop]       Name of the population to use e.g. hapmap_ceu for CEU HapMap,
                       1kg_yri for YRI 1000 genomes. See documentation for more
                       details
--freq_freq [freq]     Frequency to use in filter. Must be a number between 0 and 0.5
--freq_gt_lt [gt|lt]   Specify whether the frequency should be greater than (gt) or
                       less than (lt) --freq_freq
--freq_filter          Specify whether variants that pass the above should be included
  [exclude|include]    or excluded from analysis

--chr [list]           Select a subset of chromosomes to analyse from your file. Any
                       data not on this chromosome in the input will be skipped. The
                       list can be comma separated, with "-" characters representing
                       a range e.g. 1-5,8,15,X [default: off]
--gp                   If specified, tries to read GRCh37 position from GP field in the
                       INFO column of a VCF file. Only applies when VCF is the input
                       format and human is the species [default: off]

--convert              Convert the input file to the output format specified.
  [ensembl|vcf|pileup] Converted output is written to the file specified in
                       --output_file. No consequence calculation is carried out when
                       doing file conversion. [default: off]

--refseq               Use the otherfeatures database to retrieve transcripts - this
                       database contains RefSeq transcripts (as well as CCDS and
                       Ensembl EST alignments) [default: off]
--host                 Manually define database host [default: "ensembldb.ensembl.org"]
-u | --user            Database username [default: "anonymous"]
--port                 Database port [default: 5306]
--password             Database password [default: no password]
--genomes              Sets DB connection params for Ensembl Genomes [default: off]
--registry             Registry file to use defines DB connections [default: off]
                       Defining a registry file overrides above connection settings.
--db_version=[number]  Force script to load DBs from a specific Ensembl version. Not
                       advised due to likely incompatibilities between API and DB

--no_whole_genome      Run in old-style, non-whole genome mode [default: off]
--buffer_size          Sets the number of variants sent in each batch [default: 5000]
                       Increasing buffer size can retrieve results more quickly
                       but requires more memory. Only applies to whole genome mode.

--cache                Enables read-only use of cache [default: off]
--dir [directory]      Specify the base cache directory to use [default: "\$HOME/.vep/"]
--write_cache          Enable writing to cache [default: off]
--build [all|list]     Build a complete cache for the selected species. Build for all
                       chromosomes with --build all, or a list of chromosomes (see
                       --chr). DO NOT USE WHEN CONNECTED TO PUBLIC DB SERVERS AS THIS
                       VIOLATES OUR FAIR USAGE POLICY [default: off]

--compress             Specify utility to decompress cache files - may be "gzcat" or
                       "gzip -dc" Only use if default does not work [default: zcat]

--skip_db_check        ADVANCED! Force the script to use a cache built from a different
                       database than specified with --host. Only use this if you are
                       sure the hosts are compatible (e.g. ensembldb.ensembl.org and
                       useastdb.ensembl.org) [default: off]
--cache_region_size    ADVANCED! The size in base-pairs of the region covered by one
                       file in the cache. [default: 1MB]
END

    print $usage;
}
