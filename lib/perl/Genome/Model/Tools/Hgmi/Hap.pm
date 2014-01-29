package Genome::Model::Tools::Hgmi::Hap;

# bdericks: This is a command that in turn calls the sub-commands needed
# to execute the BAP pipeline. This process is being converted to a build
# process, so this script will soon be obsolete
#
# TODO Why not have the _execute_build method of the bacterial gene annotation
# processing profile do all this directly? Better yet, why not make this a
# workflow? Doing so would require fixing all the directory changing first. Also,
# probably won't be able to use an xml file for the workflow, since skipping
# a few steps is necessary at times

use strict;
use warnings;

# FIXME Remove this dependency
use lib "/gsc/scripts/opt/bacterial-bioperl";

use Genome;
use Command;

use Carp;
use English;
use Cwd;
use File::Path qw(mkpath);
use File::Spec;
use File::Temp qw/ tempfile tempdir /;
use YAML qw( LoadFile DumpFile );
use Data::Dumper;

class Genome::Model::Tools::Hgmi::Hap (
    is => 'Command',
    has => [
        config => {
            is  => 'String',
            doc => 'YAML file for reading',
        },
    ],
    has_optional => [
        gen_example => {
            is => 'Boolean',
            doc => 'Generate an example yaml config file',
        },
        internalhash => {
            is => 'HashRef',
            is_transient => 1,
            doc => 'internal',
        },
        dev => {
            is => 'Boolean',
            doc => 'development flag for testing',
        },
        skip_core_check => {
            is => 'Boolean',
            doc => 'skips core genes check',
            default => 0,
        },
        skip_ber => {
            is => 'Boolean',
            doc => 'skips the JCVI product naming tools',
            default => 0,
        },
        skip_protein_annotation => {
            is => 'Boolean',
            doc => 'skips running bap_finish, protein annotation and ber',
            default => 0,
        },
		keep_pep => {
            is => 'Boolean',
            doc => 'Keep temporary fasta file of gene proteins',
            default => 1
        },
        pep_file => {
            is => 'String',
            doc => 'Fasta file of gene proteins',
        },
        keggscan_version => {
            is => 'Number',
            doc => 'Version of KEGGScan to use',
        },
        interpro_version => {
            is => 'Text',
            doc => 'Version of iprscan and data to use',
        },
    ]
);

sub help_brief { return 'Runs the HGMI tools pipeline' }
sub help_synopsis { help_brief() }
sub help_detail {
    return 'Runs the entire HGMI pipeline, includes prediction, merging, finishing, PAP, and BER';
}

sub execute {
    my $self = shift;

    local $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
    $self->debug_message("Beginning execution of HAP command");

    # FIXME I'd like to get rid of the configuration file altogether
    confess "No configuration file found at " . $self->config unless -f $self->config;
    my $config = LoadFile($self->config);

    if($self->interpro_version and not exists $config->{iprpath}) {
        $config->{iprpath} = '/gsc/scripts/pkg/bio/iprscan/iprscan-' . $self->interpro_version .'/bin/iprscan';
    }

    # Core gene check only relevant to bacteria
    if($config->{cell_type} eq 'VIRAL') {
        $self->skip_core_check(1);
    }

    # FIXME This directory structure can really be simplified
    $self->debug_message("Creating directory structure.");

    # Directory structure is created
    my $dir_builder = Genome::Model::Tools::Hgmi::DirBuilder->create(
        path                  => $config->{path},
        org_dirname           => $config->{org_dirname},
        assembly_version_name => $config->{assembly_name},
        assembly_version      => $config->{assembly_version},
        pipe_version          => $config->{pipe_version},
        cell_type             => $config->{cell_type}
    );
    confess "Could not create directory builder object!" unless $dir_builder;
    my $dir_build_rv = $dir_builder->execute;
    confess "Trouble executing directory builder!" unless defined $dir_build_rv and $dir_build_rv == 1;

    $self->debug_message("Directories created, now collecting sequence!");

    # FIXME Should not be changing directories like this
    my $next_dir = $config->{path} . "/"
        . $config->{org_dirname} . "/"
        . $config->{assembly_name} . "/"
        . $config->{assembly_version} . "/"
        . "Sequence/Unmasked";
    confess "Directory does not exist: $next_dir" unless -d $next_dir;
    chdir($next_dir);
    $self->debug_message("Changed directory to $next_dir");

    # Collect sequence into directory
    # FIXME Add output directory param here so changing directories isn't necessary
    my $collect_sequence = Genome::Model::Tools::Hgmi::CollectSequence->create(
        seq_file_name  => $config->{seq_file_name},
        seq_file_dir   => $config->{seq_file_dir},
        minimum_length => $config->{minimum_length},
    );
    confess "Could not create collect sequence object!" unless $collect_sequence;
    my $collect_seq_rv = $collect_sequence->execute;
    confess "Trouble executing collect sequence!" unless defined $collect_seq_rv and $collect_seq_rv == 1;

    $self->debug_message("Sequence collection done, now naming sequence.");

    # Run sequence name tool (what does this do?) 
    my $sequence_name = Genome::Model::Tools::Hgmi::SequenceName->create(
        locus_tag        => $config->{locus_tag},
        fasta            => $collect_sequence->new_ctgs_out,
        analysis_version => $config->{pipe_version},
        #acedb_version    => $config->{acedb_version},
        acedb_version    => $self->acedb_version_lookup( $config->{acedb_version} ),
        project_type     => $config->{project_type},
        path             => $config->{path},
    );
    confess "Could not create sequence name object!" unless $sequence_name;
    my $seq_name_rv = $sequence_name->execute;
    confess "Trouble executing sequence name!" unless defined $seq_name_rv and $seq_name_rv == 1;

    # Once again, changing directory because a tool doesn't accept an output directory option
    # FIXME Remove this!
    $next_dir = $config->{path} . "/"
        . $config->{org_dirname} . "/"
        . $config->{assembly_name} . "/"
        . $config->{assembly_version} . "/" . "BAP" . "/"
        . $config->{pipe_version}
        . "/Sequence";
    confess "Directory does not exist: $next_dir" unless -d $next_dir;
    chdir($next_dir);
    $self->debug_message("Changed directory to $next_dir");

    $self->debug_message("Sequence naming complete, now making prediction models.");

    # Making prediction models, which are thrown into current working directory
    # FIXME Add an output directory parameter
    my $make_models = Genome::Model::Tools::Hgmi::MkPredictionModels->create(
        locus_tag  => $config->{locus_tag},
        fasta_file => $sequence_name->new_output,
    );
    confess "Could not create make models object!" unless $make_models;
    my $make_models_rv = $make_models->execute;
    confess "Trouble executing make models!" unless defined $make_models_rv and $make_models_rv == 1;

    # FIXME Remove directory hopping
    $next_dir = $config->{path} . "/"
        . $config->{org_dirname} . "/"
        . $config->{assembly_name} . "/"
        . $config->{assembly_version} . "/" . "BAP" . "/"
        . $config->{pipe_version};
    confess "Directory does not exist: $next_dir" unless -d $next_dir;
    chdir($next_dir);
    $self->debug_message("Changed directory to $next_dir");

    $self->debug_message("Prediction models created, now running gene prediction!");

    # Create prediction object and execute
    # FIXME Add output directory param, do not rely on current working directory
    my $predict = Genome::Model::Tools::Hgmi::Predict->create(
        organism_name    => $config->{organism_name},
        locus_tag        => $config->{locus_tag},
        project_type     => $config->{project_type},
        ncbi_taxonomy_id => $config->{ncbi_taxonomy_id},
        gram_stain       => $config->{gram_stain},
        locus_id         => $config->{locus_id},
        dev              => $self->dev,
        runner_count     => $config->{runner_count},
    );
    confess "Could not make prediction object!" unless $predict;

    # Skip execution if previous gene prediction execution was successful
    # TODO Can be removed when this is a workflow
    if ($predict->is_valid()) {
        $self->debug_message("Prediction has already been run successfully, continuing!");
    }
    else {
        my $predict_rv = $predict->execute;
        confess "Trouble executing gene prediction!" unless defined $predict_rv and $predict_rv == 1;
    }

    $self->debug_message("Gene prediction run, now starting gene merging!");

    # Create merge object and execute
    my $merge = Genome::Model::Tools::Hgmi::Merge->create(
        organism_name => $config->{organism_name},
        locus_tag     => $config->{locus_tag},
        project_type  => $config->{project_type},
        runner_count  => $config->{runner_count},
        dev           => $self->dev,
        nr_db         => $config->{nr_db},
        iprpath       => $config->{iprpath},
        ipr_version   => $self->interpro_version,
    );
    confess "Could not create gene merging object!" unless $merge;

    # Skip execution if previous gene merging run was successful
    # TODO Can be removed when this is a workflow
    if ($merge->is_valid) {
        $self->debug_message("Skipping gene merging step, previous execution was successful!");
    }
    else {
        my $merge_rv = $merge->execute;
        confess "Trouble executing gene merging!" unless defined $merge_rv and $merge_rv == 1;
    }

    $self->debug_message("Gene merging done, now tagging overlaps.");

    # Overlapping genes aren't properly tagged in merging, so an extra step is run to tag them
    my $tag_overlaps = Genome::Model::Tools::Bacterial::TagOverlaps->create(
        sequence_set_id => $merge->sequence_set_id,
        dev             => $self->dev,
    );
    confess "Could not create tag overlaps object!" unless $tag_overlaps;
    my $overlap_rv = $tag_overlaps->execute;
    confess "Trouble executing tag overlaps!" unless defined $overlap_rv and $overlap_rv == 1;

    # FIXME Do not change directories like this, have the tool accept an output directory param!
    $next_dir = $config->{path} . "/"
        . $config->{org_dirname} . "/"
        . $config->{assembly_name} . "/"
        . $config->{assembly_version} . "/"
        . "Genbank_submission/"
        . $config->{pipe_version}
        . "/Annotated_submission";
    confess "Directory does not exist: $next_dir" unless -d $next_dir;
    chdir($next_dir);
    $self->debug_message("Changed directory to $next_dir");

    $self->debug_message("Overlaps are totally tagged, running rrna screen.");

    # Running rrna screen
    my $rrna_screen = Genome::Model::Tools::Hgmi::RrnaScreen->create(
        sequence_set_id => $merge->sequence_set_id,
        dev => $self->dev,
        rrna_database => '/gscmnt/278/analysis/HGMI/rRNA_testing/16s_23srnadb',
    );
    confess "Could not create rrna screen object!" unless $rrna_screen;
    my $rrna_screen_rv = $rrna_screen->execute;
    confess "Trouble executing rrna screen!" unless defined $rrna_screen_rv and $rrna_screen_rv == 1;

    $self->debug_message("Rrna screen done, running finishing step.");

    # Run finish step
    my $finish = Genome::Model::Tools::Hgmi::Finish->create(
        sequence_set_id  => $merge->sequence_set_id,
        locus_tag        => $config->{locus_tag},
        organism_name    => $config->{organism_name},
        project_type     => $config->{project_type},
        acedb_version    => $config->{acedb_version},
        assembly_name    => $config->{assembly_name},
        org_dirname      => $config->{org_dirname},
        assembly_version => $config->{assembly_version},
        pipe_version     => $config->{pipe_version},
        path             => $config->{path},
        nr_db            => $config->{nr_db},
        dev              => $self->dev,
        skip_acedb_parse => $config->{skip_acedb_parse},
    );
    confess "Could not create finish object!" unless $finish;
    my $finish_rv = $finish->execute;
    confess "Trouble executing finish step!" unless defined $finish_rv and $finish_rv == 1;

    $self->debug_message("Done with finishing!");

    chdir($next_dir);
    $self->debug_message("Changed directory to $next_dir");

    # Run the core gene check (unless the user says we shouldn't!)
    unless ($self->skip_core_check)
    {
        $self->debug_message("Preparing to run core gene check!");

        my $core_gene = Genome::Model::Tools::Hgmi::CoreGenes->create(
            cell_type       => $config->{cell_type},
            sequence_set_id => $merge->sequence_set_id,
            dev             => $self->dev,
        );
        confess "Could not create core gene check object!" unless $core_gene;
        my $core_rv = $core_gene->execute;
        confess "Trouble executing core gene check!" unless defined $core_rv and $core_rv == 1;
    }
    else {
        $self->debug_message("Skipping core gene check.");
    }

    # If the user wants to skip protein annotation, we're done!
    if($self->skip_protein_annotation) {
        $self->debug_message("Skipping protein annotation!");
        return 1;
    }

    # Also skip PAP if there's no workflow file spepcifed in the config file! 
    unless (defined $config->{workflowxml}) {
        $self->debug_message("No workflow xml path provided in config file, so PAP/BER will not be run. Exiting...");
        return 1;
    }

    $self->debug_message("Moving data from mgap to biosql: Hap.pm");
    $self->mgap_to_biosql($config->{locus_tag}, $merge->sequence_set_id);

    $self->debug_message("Creating peptide file: Hap.pm");
    $self->get_gene_peps($config->{locus_tag});
	
    $self->debug_message("Getting ready to run PAP!");

    my $gram_stain = $config->{gram_stain};
    confess 'Gram stain not specified in config file, cannot start PAP workflow!' unless defined $gram_stain;
    unless ($gram_stain eq 'positive' or $gram_stain eq 'negative') {
        confess "Gram stain invalid: must be either positive or negative, not $gram_stain";
    }

    my $base_archive_dir = File::Spec->catdir(
        $config->{path},          $config->{org_dirname},
        $config->{assembly_name}, $config->{assembly_version},
    );

    my $blastp_archive_dir = File::Spec->catdir($base_archive_dir, 'Blastp', $config->{pipe_version}, 'Hybrid');
    my $interpro_archive_dir = File::Spec->catdir($base_archive_dir, 'Interpro', $config->{pipe_version}, 'Hybrid');
    my $keggscan_archive_dir = File::Spec->catfile($base_archive_dir, 'Kegg', $config->{pipe_version}, 
        'Hybrid', join('.', 'KS-OUTPUT', $config->{locus_tag}, 'CDS', 'pep'));
    my $psortb_archive_dir = File::Spec->catfile($base_archive_dir, 'psortB', $config->{pipe_version}, 'Hybrid');

    my $send = Genome::Model::Tools::Hgmi::SendToPap->create(
        skip_db_upload       => ((exists $config->{skip_db_upload}) ? $config->{skip_db_upload} : 0),
        locus_tag            => $config->{locus_tag},
        sequence_set_id      => $merge->sequence_set_id,
        sequence_name        => $config->{assembly_name},
        organism_name        => $config->{organism_name},
        workflow_xml         => $config->{workflowxml},
        gram_stain           => $config->{gram_stain},
        blastp_archive_dir   => $blastp_archive_dir,
        interpro_archive_dir => $interpro_archive_dir,
        keggscan_archive_dir => $keggscan_archive_dir,
        psortb_archive_dir   => $psortb_archive_dir,
        dev                  => $self->dev,
		pep_file			 => $self->pep_file,
        keggscan_version     => $self->keggscan_version,
        interpro_version     => $self->interpro_version,
    );
    confess 'Could not create command object for send to pap!' unless $send;

    $self->debug_message("Parameters for protein annotation:\n" . Data::Dumper::Dumper($send));
    confess 'Could not run pap!' unless $send->execute;

    # jcvi product naming goes here.
    # need tochdir
    # $path/Acedb/$acedb_version/ace_files/$locus_tag/$pipe_version
    unless($self->skip_ber) 
    {

        my $acedb_version
            = $self->acedb_version_lookup( $config->{acedb_version} );

        $next_dir = $config->{path}
            . "/Acedb/"
            . $acedb_version
            . "/ace_files/"
            . $config->{locus_tag} . "/"
            . $config->{pipe_version};

        warn qq{\n\nACEDB_Dir: $acedb_version\n\n};

        unless ( -d $next_dir )
        {
            croak
                qq{\n\nThe directory : '$next_dir', does not exit, from Hap.pm: $OS_ERROR\n\n};
        }
        my $cwd = getcwd();
        unless ( $cwd eq $next_dir )
        {
            chdir($next_dir)
                or croak
                "Failed to change to '$next_dir', from Hap.pm: $OS_ERROR\n\n";
            $self->debug_message("Changed directory to $next_dir");
        }

        #run /gsc/scripts/gsc/annotation/biosql2ace <locus_tag>
        # make sure files are not blank. croak if they are

        my @biosql2ace = (
            'biosql2ace',
            $config->{locus_tag},
        );
        if ( $self->dev )
        {
            push( @biosql2ace, '--dev' );
        }
        my ( $b2a_out, $b2a_err );
        IPC::Run::run( \@biosql2ace, '>', \$b2a_out, '2>', \$b2a_err )
            or croak "cant dump biosql to ace\n\n";

    # check that output files are not empty
        my @outputfiles = ( );
        if($config->{workflowxml} =~ /noblastp/)
        {
            @outputfiles = qw(merged.raw.sorted.ace merged.raw.sorted.ipr.ace REPORT-top_new.ks.ace)
        }
        else
        {
            @outputfiles = qw(briefID.fof.ace merged.raw.sorted.ace merged.raw.sorted.ipr.ace REPORT-top_new.ks.ace)
        }

        foreach my $outputfile ( @outputfiles )
        {
            my $size = -s $outputfile;
            if ( $size == 0 )
            {
                croak
                    "file from biosql2ace dump, $outputfile , is empty...from Hap.pm\n\n";
            }
        }

        my $ber_config = undef;

        my $jcvi = Genome::Model::Tools::Ber::AmgapBerProtName->create(
            sequence_set_id => $merge->sequence_set_id,
            config          => $self->config(),
        );

        if ( $self->dev )
        {
            $jcvi->dev(1);
        }

        if ($jcvi)
        {
            $jcvi->execute()
                or croak
                "can't run protein product naming step...from Hap.pm\n\n";
        }
        else
        {
            croak "can't set up product naming step...from Hap.pm\n\n";
        }
    } # ber skipping

    return 1;
}

sub read_config
{
    my $self = shift;

    my $conf = $self->config;
    unless ( -f $conf )
    {
        carp "no config file $conf ...Hap.pm\n\n";
        return undef;
    }

    my $confhash = LoadFile($conf);

    return 1;
}

sub build_empty_config
{
    my $self     = shift;
    my $dumpfile = $self->config;
    my $config   = {

        #dir builder stuff
        'path'          => "",
        'org_dirname'   => "<organism abbreviated name>",
        'assembly_name' =>
            "<full org name, locus_tag, finish assembly, pipeline>",
        'assembly_version' => "",
        'pipe_version'     => "",
        'cell_type'        => "<BACTERIA or ARCHEA or VIRAL>",

        #collect sequence stuff
        'seq_file_name'  => "",
        'seq_file_dir'   => "",
        'minimum_length' => "",

        # sequence name
        'assembly_version' => "",
        'locus_id'         => "<locus_tag w/o DFT/FNL>",
        'acedb_version'    => "",

        #mk prediction mods
        'locus_tag' => "",

        #predict
        'runner_count'            => "",
        'organism_name'           => "",
        'project_type'            => "",
        'gram_stain'              => "",
        'ncbi_taxonomy_id'        => "",
        'predict_script_location' => "<optional>",

        #merge
        # uses some of the same items from predict
        'merge_script_location' => "<optional>",

        #finish
        # uses some of the same items from predict
        'acedb_version'          => "",
        'skip_acedb_parse'       => "<optional>",
        'finish_script_location' => "<optional>",
        'workflowxml'            => "",
    };
    DumpFile( $dumpfile, $config );    # check return?
    return 1;
}

sub acedb_version_lookup
{
    my $self         = shift;
    my $v            = shift;
    my $acedb_lookup = undef;

    my %acedb_version_lookup = (
        'DEV'  => 'Development',
        'DEV1'  => 'Development_1',
        'DEV2'  => 'Development_2',
        'V1'  => 'Version_1.0',
        'V2'  => 'Version_2.0',
        'V3'  => 'Version_3.0',
        'V4'  => 'Version_4.0',
        'V5'  => 'Version_5.0',
        'V6'  => 'Version_6.0',
        'V7'  => 'Version_7.0',
        'V8'  => 'Version_8.0',
        'V9'  => 'Version_9.0',
        'V10' => 'Version_10.0',
    );

    if ( exists( $acedb_version_lookup{$v} ) )
    {
        $acedb_lookup = $acedb_version_lookup{$v};
    }
    
   else ##veena added 05/09
    {
        $v =~ s/V(\d+)/Version_$1.0/;
        $acedb_lookup = $v;
    }

    return $acedb_lookup;
}

sub get_gene_peps {
    my $self = shift;
    my ($locus_tag) = @_;

    my $dbadp;

    if($self->dev()) {
        $dbadp = Bio::DB::BioDB->new(
            -database => 'biosql',
            -user     => 'sg_user',
            -pass     => 'sgus3r',
            -dbname   => 'DWDEV',
            -driver   => 'Oracle'
        );
    }
    else {
        $dbadp = Bio::DB::BioDB->new(
            -database => 'biosql',
            -user     => 'sg_user',
            -pass     => 'sg_us3r',
            -dbname   => 'DWRAC',
            -driver   => 'Oracle'
        );
    }

    my $cleanup = $self->keep_pep ? 0 : 1;
    my $tempdir = tempdir(
        CLEANUP => $cleanup,
        DIR     => '/gscmnt/temp212/info/annotation/PAP_tmp', # FIXME Shouldn't have hard-coded paths for temp stuff
    );
    my ($fh, $file) = tempfile(
        "pap-XXXXXX",
        DIR    => $tempdir,
        SUFFIX => '.fa'
    );

    unless(defined($self->pep_file)) {
        $self->pep_file($file);
    }

    $file = $self->pep_file();
    my $seqout = new Bio::SeqIO(
        -file   => ">$file",
        -format => "fasta"
    );

    my $adp = $dbadp->get_object_adaptor("Bio::SeqI");
    $adp->dbh->{'LongTruncOk'} = 0;
    $adp->dbh->{'LongReadLen'} = 1000000000;

    my $query = Bio::DB::Query::BioQuery->new();
    $query->datacollections( [ "Bio::PrimarySeqI s", ] );

    $query->where( ["s.display_id like '$locus_tag%'"] );
    my $res = $adp->find_by_query($query);
    GENE: while (my $seq = $res->next_object()) {
        my $gene_name = $seq->display_name();
        my @feat = $seq->get_SeqFeatures();
        FEATURE: foreach my $f (@feat) {
            my $display_name = $f->display_name();
            next FEATURE if $f->primary_tag ne 'gene';
            next if $f->has_tag('Dead');
            my $ss;
            $ss = $seq->subseq( $f->start, $f->end );

            unless(defined($ss)) {
                die "failed to fetch sequence for '$display_name' from BioSQL";
            }

            my $pep = $self->make_into_pep( $f->strand,
                $ss,
                $display_name );
            $seqout->write_seq($pep);
        }
    }

    unless (-f $file) {
        print STDERR "the fasta file $file doesn't exist! SendToPap.pm\n";
        return 0;
    }
    unless(-s $file > 0) {
        print STDERR "the fasta file $file still empty! SendToPap.pm\n";
    }

    return 1;
}

sub mgap_to_biosql
{
    my $self      = shift;
	my ($locus_tag, $ssid, $testnorun) = @_;

    # FIXME Deployed code... NOOOOOOOOOOOOOO!
    my @command = (
        'bap_load_biosql', '--sequence-set-id', $ssid,
        #'--tax-id',        $taxid,
    );

    if($self->dev)
    {
        push(@command,'--dev');
        push(@command,'--biosql-dev');
    }
    my ($cmd_out,$cmd_err);

	if($testnorun)
	{
# just for testing.
		print join(" " , @command),"\n";;
		return 1;
	}

	IPC::Run::run(
			\@command,
			\undef,
			'>',
			\$cmd_out,
			'2>',
			\$cmd_err,
			) or croak "can't load biosql from mgap!!!\n$cmd_err SendToPap.pm";

	print STDERR $cmd_err,"\n";
	return 1;

}

sub make_into_pep
{
    my ( $self, $strand, $subseq, $name ) = @_;

    my $seq = new Bio::Seq(
        -display_id => $name,
        -seq        => $subseq
    );
    my $newseq;
    if ( $strand < 0 )
    {
        $newseq = $seq->revcom->translate->seq;
    }
    else
    {
        $newseq = $seq->translate->seq;
    }
    my $seqobj = new Bio::Seq(
        -display_id => $name,
        -seq        => $newseq
    );
    return $seqobj;
}

1;

# $Id$
