package Genome::Model::GenePrediction::Command::Bacterial::Merge;

use strict;
use warnings;
use Genome;

use lib "/gsc/scripts/opt/bacterial-bioperl";

use BAP::Config;
use BAP::DB::Organism;
use BAP::DB::CodingGene;
use BAP::DB::SequenceSet;
use BAP::JobSource::InterGenicBlastX;
use BAP::JobSource::Phase2BlastP;
use BAP::GeneMerge;

use PP;
use PP::LSF;
use PP::RPC;

use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;

use Carp;
use Data::Dumper;
use DateTime;
use DateTime::Duration;
use DBI;
use DBD::SQLite;
use English;
use File::stat;
use File::Slurp qw(slurp);
use File::Temp;

use IPC::Run;
use Genome::Utility::Email;
use XML::DOM::XPath;

class Genome::Model::GenePrediction::Command::Bacterial::Merge {
    is  => 'Command',
    has => [
        sequence_set_id => {
            is  => 'Integer',
            doc => "sequence set id",
        },
        nr_db => {
            is  => 'String',
            doc => "path to non-redundant protein database",
        },
    ],
    has_optional => [
        iprpath => {
            is => 'Scalar',
            doc => "alternative path to iprscan",
            default => "/gsc/scripts/bin/iprscan",
        },
		 ipr_version => {
			is => 'String',
            valid_values => ['4.5', '4.7', '4.8', '4.8-40'],
            default => '4.8',
            },
        runner_count => {
            is  => 'Integer',
            doc => "number of runner jobs to spawn off to LSF; default 10",
            default => 10,
        },
        job_stdout => {
            is  => 'String',
            doc => "file to write stdout",
        },
        job_stderr => {
            is  => 'String',
            doc => "file to write stderr",
        },
        only_phase => {
            is  => 'Integer',
            doc => "only run specified phase (1,2,3,4,or 5)",
            valid_values => ['1', '2', '3', '4', '5'],
        },
        debug_file => {
            is  => 'String',
            doc => "path to debug file",
        },
        dev => {
            is      => 'Boolean',
            doc     => "use development database",
            default => 0,
        },
        skip_blastx => {
            is      => 'Boolean',
            doc     => "skip blastx",
            default => 0,
        },
        no_mail => {
            is      => 'Boolean',
            doc     => "don't send email upon completion",
            default => 0,
        },
        tmp_usage => {
            is      => 'Integer',
            doc     => "tmp usage in megabytes; default 100Mbytes",
            default => 100,
        },
        overlap_percent => {
            is  => 'Integer',
            doc => "gene over lap percentage",
            default => 30,
        },
        overlap_bp => {
            is  => 'Integer',
            doc => "number of bp genes may over lap",
            default => 300,
        },
        rpc_queue => {
            is      => 'String',
            doc     => "LSF queue to use",
            default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        },
        rpc_core_number => {
            is      => 'String',
            doc     => "number of cores to use for LSF jobs; default 2",
            default => 2,
        },
        min_gene_length => {
            is      => 'Integer',
            doc     => "minimum gene length; default 60 bp",
            default => 60,
        },
        rfam_seed => {
            is => 'Scalar',
            doc => "path to Rfam.seed file",
            default => $ENV{GENOME_SW} . "/rfam/rfam-8.1/Rfam.seed",
        },
        _debug_fh => {
            is  => 'Scalar',
            doc => "",

        },
        _rpc_args => {
            is  => 'Hashref',
            doc => "",
        },
        _selected_genes => {
            is => 'Hashref',
            doc => "",
        },
        _network_temp_args => {
            is => 'Hashref',
            doc => "",
        },
        _rfam_special => {
            is => 'Hashref',
            doc => "",
        },
        _fetched_sequences => {
            is => 'Hashref',
            doc => "",
        },

    ],
};

#our %rpc_args;
#our $debug_fh;

$File::Temp::KEEP_ALL = 1;


sub sub_command_sort_position {30}

sub help_brief
{
"perform gene merging operations on predicted bacterial genes";
}

sub help_synopsis
{
return <<"EOS"
genome model gene-prediction merge [options]
EOS
}

sub help_detail
{
return <<"EOS"
runs gene merging on predicted genes for the specified sequence set id.
EOS
}

sub execute
{
    my $self = shift;
    my $user = Genome::Sys->username;

    $self->debug_message("iprpath : ". "/gsc/scripts/pkg/bio/iprscan/iprscan-". $self->ipr_version ."/bin/iprscan");

    my $runner_count = $self->runner_count;

    my %temp_args = (
                       'TEMPLATE' => 'mgap_XXXXXXXX',
                       'DIR'      => '/gscmnt/temp212/info/annotation/BAP_tmp',
                       'SUFFIX'   => '.temp',
                    );

    $self->_network_temp_args( \%temp_args );
    my $dt_started = mark_time();

    if ( $self->dev ) { $BAP::DB::DBI::db_env = 'dev'; }

    unless ( -e $self->nr_db )
    {
        croak "nr-db does not exist; " . $self->nr_db;
    }

    unless ( -r $self->nr_db )
    {
        croak "nr-db is not readable;" . $self->nr_db;
    }

    # Note the mtime of the blast database
    my $data_version
        = DateTime->from_epoch( epoch => stat( $self->nr_db )->mtime() )
        ->mdy();

    # Note BAP/MGAP version and Subversion revision
    my $svn_version = '$Revision$';
    $svn_version =~ s/\D+//g;
    my $bap_version = BAP::Config->new()->version();

    # Find all the riboswitches / leaders in Rfam and note the accessions
    my %rfam_special;
    my $rfam_io = Bio::AlignIO->new(
        -file => $self->rfam_seed 
    );

    # This loop emits thousands and thousands of bioperl warnings, which are the result
    # of Bio::Seq::Meta objects being created by the next_aln method below. Setting
    # verbose to -1 suppresses these warnings and keeps the error log from being
    # filled with hundreds of thousands of lines of warning messages.
    $rfam_io->verbose(-1);
    while (my $aln = $rfam_io->next_aln) {
        my $accession   = $aln->accession();
        my $description = $aln->description();

        if (($description =~ /riboswitch/i) || ($description =~ /leader/i)) {
            $rfam_special{$accession} = 1;
        }
    }
    $self->_rfam_special(\%rfam_special);

    # Blow up the files we're going to write stdout/stderr to if they exist.
    foreach my $out_file ( $self->job_stdout, $self->job_stderr )
    {

        if ( defined($out_file) && ( -e $out_file ) )
        {
            unlink($out_file) || carp "Failed to unlink '$out_file'!";
        }

    }

    my $debug_fh   = IO::File->new();
    my $debug_file = $self->debug_file;
    $debug_fh->open(">$debug_file")
        or die "Can't open '$debug_file': $OS_ERROR";
    $self->_debug_fh($debug_fh);

## Try not to hork the database file if we catch a signal.
    #local $SIG{INT} = \&{ $self->handle_sigint };
    my $rpc_queue    = $self->rpc_queue;
    my $rpc_core_num = $self->rpc_core_number;
    my $tmp_usage    = $self->tmp_usage;

    my %rpc_args = (
        runner_count => $runner_count,
        app_init     => 0,
        port         => 7654 + $PID,
        pp_type      => 'lsf',
        q            => "'$rpc_queue'",
        n            => "$rpc_core_num",
        R            => "'span[hosts=1] rusage[mem=4096, tmp=$tmp_usage]'",
        maxmessage   => 81920000,
        lib_paths    => [ UR::Util::used_libs ],
    );

    if ( defined( $self->job_stdout ) ) {
        $rpc_args{'o'} = $self->job_stdout;
    }
    if ( defined( $self->job_stderr ) ) {
        $rpc_args{'e'} = $self->job_stderr;
    }

    $self->_rpc_args( \%rpc_args );

## Run whatever chunks 'o code we're supposed to.
## WARNING! Stupid Tricks with Code Refs ahead.

    my %fetched_sequences = ();

    {
        my $sequence_set
            = BAP::DB::SequenceSet->retrieve( $self->sequence_set_id );

        my @sequences = $sequence_set->sequences();

        %fetched_sequences = map { $_->sequence_name() => 1 } @sequences;

    }
    $self->_fetched_sequences(\%fetched_sequences); 

    my %selected_genes = ();
    # is it necessary to have a local block?
    {

        my @run_phases = ();

        if ( defined( $self->only_phase ) )
        {
            @run_phases = ( $self->only_phase );
        }
        else
        {

            @run_phases = ( 1, 2, 3, 4, 5 );
            if($self->skip_blastx) {
                @run_phases = ( 1,2,4,5);
            }
        }
        foreach my $phase (@run_phases)
        {
            $self->debug_message("running phase ". $phase);
#            unless ( defined( $self->skip_blastx )
#                && ( $phase == 3 ) )
#            {

                $self->debug_message("before phase ". $phase);

                if($phase == 1) {
                    $self->phase1();
                }
                elsif( $phase == 2 ) {
                    $self->phase2();
                }
                elsif( $phase == 3 ) {
                    $self->phase3();
                }
                elsif( $phase == 4 ) {
                    $self->phase4();
                }
                elsif( $phase == 5 ) {
                    $self->phase5();
                }
                $self->debug_message("after phase ".$phase);
                BAP::DB::DBI->dbi_commit();
                BAP::DB::DBI->db_Main->disconnect();
#            }
            $self->debug_message("completed phase ". $phase);

        }

    }
    my $dt_finished = mark_time();

    my $sequence_set
        = BAP::DB::SequenceSet->retrieve( $self->sequence_set_id );

    $sequence_set->software_version($bap_version);
    $sequence_set->data_version($data_version);
    $sequence_set->update();
    BAP::DB::DBI->dbi_commit();

    my $sequence_set_name = $sequence_set->sequence_set_name();

    unless ( $self->dev )
    {
        $self->log_run( $dt_started, $dt_finished, $sequence_set );
    }

    unless ( $self->no_mail )
    {

        $self->send_mail( $self->sequence_set_id, $sequence_set_name, $dt_started,
            $dt_finished, $bap_version, $svn_version, $user, );

    }

    #exit 0;

    return 1;
}

sub handle_sigint
{
    carp "shutting down (caught SIGINT in the face)";

    BAP::DB::DBI->dbi_rollback();
    exit(1);
}

sub mark_time
{

    return DateTime->now( time_zone => 'America/Chicago' );

}

sub send_mail
{
    my $self = shift;
    my ( $ss_id, $ss_name, $started, $finished, $bap_version, $svn_version,
        $user )
        = @_;

    my $date_started = join( ' on ', $started->hms(':'), $started->ymd('/') );
    my $date_finished
        = join( ' on ', $finished->hms(':'), $finished->ymd('/') );
    my $duration      = DateTime::Duration->new( $finished - $started );
    my $hours_running = $duration->in_units('hours');
    my $from   = Genome::Utility::Email::construct_address($user);

    my $subject = "BMG script mail for MGAP SSID: $ss_id ($ss_name)";

    my $body = <<BODY;
The bap_merge_genes script ($bap_version/$svn_version) has finished running MGAP SSID: $ss_id ($ss_name).
The job began at $date_started and ended at $date_finished.
The total run time of the script is:  $hours_running  hours.
BODY

    Genome::Utility::Email::send(
        from    => $from,
        to      => [ $from, 'kpepin@watson.wustl.edu' ], #FIXME remove hardcoded address
        subject => $subject,
        body    => $body,
    );
}

sub check_failed_jobs
{

    my ($job_source) = @_;

    my $failed_jobs = $job_source->failed_jobs();

    if ( $#$failed_jobs > 0 )
    {

        foreach my $job ( @{$failed_jobs} )
        {

            my $job_id  = $job->job_id();
            my $class   = ref($job);
            my $seq_obj = $job->seq();
            my $seq_id  = $seq_obj->display_id();

            warn "job $job_id ($class) ($seq_id) failed";

        }

        warn "There were failed jobs - rolling back";

        BAP::DB::DBI->dbi_rollback();
        exit(1);

    }

}

sub phase1
{
    my $self = shift;
    $self->debug_message("in phase_1");
    $self->best_per_locus( 'phase_1', 'phase_0' );

}

sub phase2
{
    my $self = shift;
    $self->debug_message("in phase_2");
    $self->check_overlapping( 'phase_2', 'phase_1' );

}

sub phase3
{
    my $self       = shift;
    my $fasta_fh   = File::Temp->new();
    my $fasta_file = $fasta_fh->filename();
    my %selected_genes = %{$self->_selected_genes};
    my $debug_fh = $self->_debug_fh;
    #my %rpc_args = $self->_rpc_args; # not use this.
    my $rpc_queue    = $self->rpc_queue;
    my $rpc_core_num = $self->rpc_core_number;
    my $tmp_usage    = $self->tmp_usage;
    my %rpc_args = (
        runner_count => $self->runner_count,
        app_init     => 0,
        port         => 7654 + $PID,
        pp_type      => 'lsf',
        q            => "'$rpc_queue'",
        n            => "$rpc_core_num",
        R            => "'span[hosts=1] rusage[mem=4096, tmp=$tmp_usage]'",
        maxmessage   => 81920000,
        lib_paths    => [ UR::Util::used_libs ],
    );

    my $feature_fh   = File::Temp->new();
    my $feature_file = $feature_fh->filename();

    {

        my $seqstream
            = Bio::SeqIO->new( -fh => $fasta_fh, -format => 'Fasta' );
        my $sequence_set
            = BAP::DB::SequenceSet->retrieve( $self->sequence_set_id );

        my @sequences = $sequence_set->sequences();
        my @sequence_names = map { $_->sequence_name() } @sequences;

        $self->sanity_check_sequences( 'phase_3', \@sequence_names );

        my @fetched_gene_names = ();

        foreach my $sequence (@sequences)
        {

            my $seq = Bio::Seq->new(
                -seq => $sequence->sequence_string(),
                -id  => $sequence->sequence_name(),
            );

            my @coding_genes = $sequence->coding_genes( phase_2 => 1 );

            push @fetched_gene_names, map { $_->gene_name() } @coding_genes;

            my @rna_genes     = $sequence->rna_genes();
            my @gene_features = ();

            {
                my @blastx_genes
                    = $sequence->coding_genes( source => 'blastx' );

                foreach my $blastx_gene (@blastx_genes)
                {

                    my @proteins = $blastx_gene->protein();

                    foreach my $protein (@proteins)
                    {
                        $protein->delete();
                    }

                    $blastx_gene->delete();

                }
            }
        P3G: foreach my $gene ( @coding_genes, @rna_genes )
            {

                my $gene_name = $gene->gene_name();

                print $debug_fh
                    join( "\t", 'phase_3', $gene_name, "fetched from db" ),
                    "\n";

                if ( $gene->isa('BAP::DB::RNAGene') )
                {

                    if ( $gene->redundant() )
                    {
                        print $debug_fh join( "\t",
                            'phase_3', $gene_name,
                            "skipping redundant rRNA" ),
                            "\n";
                        next P3G;
                    }

                }

                if ( $gene->isa('BAP::DB::CodingGene') )
                {

                    $gene->phase_3(1);
                    $gene->update();

                    print $debug_fh
                        join( "\t", 'phase_3', $gene_name, "selected" ), "\n";
                    $selected_genes{'phase_3'}{$gene_name} = 1;

                }

                my $gene_feature
                    = $self->gene_to_feature( $sequence->sequence_name(), $gene );

                $gene_feature->attach_seq($seq);

                push @gene_features, $gene_feature;

            }

            foreach my $gene_feature (@gene_features)
            {
                print $feature_fh $gene_feature->gff_string(), "\n";
            }

            $seqstream->write_seq($seq);

        }

        $self->sanity_check_genes( 'phase_3', 'phase_2', \@fetched_gene_names );

        $seqstream->close();
        close($fasta_fh);
        close($feature_fh);

    }

    BAP::DB::DBI->dbi_commit();
    BAP::DB::DBI->db_Main->disconnect();

    my $seqstream
        = Bio::SeqIO->new( -file => $fasta_file, -format => 'Fasta' );
    my $featstream = Bio::Tools::GFF->new( -file => $feature_file );

    my $job_source
        = BAP::JobSource::InterGenicBlastX->new( $seqstream, $featstream,
        $self->nr_db, $self->rpc_core_number);

    local $rpc_args{job_source} = $job_source;

    my $server = PP::RPC->new(%rpc_args);

    $server->start;

    check_failed_jobs($job_source);

    my $collection   = $job_source->feature_collection();
    my @blastx_genes = $collection->get_all_features();

    my %blastx_seq_ids = map { $_->seq_id() => 1 } @blastx_genes;

    foreach my $blastx_seq_id ( keys %blastx_seq_ids )
    {

        my $sequence = BAP::DB::Sequence->retrieve(
            'sequence_name'   => $blastx_seq_id,
            'sequence_set_id' => $self->sequence_set_id,
        );
        my $sequence_obj = Bio::Seq->new(
            -id  => $sequence->sequence_name(),
            -seq => $sequence->sequence_string(),
        );

        my @blastx_genes
            = grep { $_->seq_id() eq $blastx_seq_id; } @blastx_genes;
        @blastx_genes = sort { $a->start() <=> $b->start() } @blastx_genes;

        foreach my $i ( 0 .. $#blastx_genes )
        {

            my $blastx_gene = $blastx_genes[$i];
            my $seq_id      = $blastx_gene->seq_id();
            my $gene_name   = join( '.', $seq_id, 'blastx', $i + 1 );

            $blastx_gene->display_name($gene_name);

        }

        foreach my $gene (@blastx_genes)
        {

            my $gene_name = $gene->display_name();
            my $gene_seq_obj
                = $sequence_obj->trunc( $gene->start, $gene->end );

            print $debug_fh
                join( "\t", 'phase_3', $gene_name, "gene created" ), "\n";
            print $debug_fh
                join( "\t", 'phase_3', $gene_name, "selected" ), "\n";

            $selected_genes{'phase_3'}{$gene_name} = 1;

            unless ( $gene->strand() > 0 )
            {
                $gene_seq_obj = $gene_seq_obj->revcom();
            }

            my $protein_seq_obj = $gene_seq_obj->translate();

            my $coding_gene = BAP::DB::CodingGene->insert(
                {   gene_name       => $gene_name,
                    sequence_id     => $sequence->sequence_id(),
                    start           => $gene->start,
                    end             => $gene->end,
                    strand          => $gene->strand(),
                    source          => 'blastx',
                    sequence_string => $gene_seq_obj->seq(),
                    phase_0         => 0,
                    phase_1         => 0,
                    phase_2         => 0,
                    phase_3         => 1,
                    phase_4         => 0,

                    phase_5 => 0,
                }
            );

            BAP::DB::Protein->insert(
                {   protein_name    => $gene_name,
                    gene_id         => $coding_gene->gene_id(),
                    sequence_string => $protein_seq_obj->seq(),
                    internal_stops  => 0,
                }
            );

        }

    }
    $self->_selected_genes(\%selected_genes);
    $self->_rpc_args(\%rpc_args);
    unlink $fasta_file, $feature_file;
    return 1;
}

sub phase4 {
    my $self = shift;
    if ($self->skip_blastx) {
        $self->best_per_locus('phase_4', 'phase_2');
    }
    else {
        $self->best_per_locus('phase_4', 'phase_3');
    }
}

sub phase5
{
    my $self = shift;
    $self->check_overlapping( 'phase_5', 'phase_4' );

}

sub tag_redundant
{
    my $self = shift;
    my ($feature_ref) = @_;
    #$self->debug_message("inside tag_redundant");
    my $is_redundant = sub {

        my ( $f, $g ) = @_;

        if ( $f->strand() == $g->strand() )
        {

            if (   ( $f->start() >= $g->start() )
                && ( $f->end() <= $g->end() )
                && ( $f->length() < $g->length() ) )
            {
                $f->add_tag_value( 'redundant', 1 );
            }

            elsif (( $g->start() >= $f->start() )
                && ( $g->end() <= $f->end() )
                && ( $g->length() < $f->length() ) )
            {
                $g->add_tag_value( 'redundant', 1 );
            }

        }

    };
    #$self->debug_message("sending code ref to find_overlaps");
    $self->find_overlaps( $is_redundant, $feature_ref );
    #$self->debug_message("find_overlaps has been run");

}

sub tag_preferred
{

    my ($feature_ref) = @_;

    my %locus = ();

    foreach my $feature ( @{$feature_ref} )
    {

        if ( $feature->strand() > 0 )
        {
            push @{ $locus{ $feature->strand() }{ $feature->end() } },
                $feature;
        }
        else
        {
            push @{ $locus{ $feature->strand() }{ $feature->start() } },
                $feature;
        }

    }

    foreach my $strand ( keys %locus )
    {

        foreach my $end ( keys %{ $locus{$strand} } )
        {

            my @features = @{ $locus{$strand}{$end} };

            my $preferred_feature;

            if ( @features > 1 )
            {
                $preferred_feature = pick_feature( \@features );
            }
            else
            {
                ($preferred_feature) = @features;
            }

            $preferred_feature->add_tag_value( 'preferred', 1 );

        }

    }

}

sub pick_feature
{

    my ($feature_ref) = @_;

    my %score = (
        ( map { $_->source_tag() => 50 } @{$feature_ref} ),
        'genemark' => 10,
        'glimmer3' => 20,
        'glimmer2' => 30,
        'blastx'   => 40,
    );

    @{$feature_ref} = sort {
               $score{ $a->source_tag() } <=> $score{ $b->source_tag() }
            || $b->length() <=> $a->length()
    } @{$feature_ref};

    return $feature_ref->[0];

}

sub tag_overlapping
{
    my $self = shift;

    my ($feature_ref) = @_;

    my %seen_overlap = ();

    my $overlap = sub {

        my ( $f, $g ) = @_;

        my $f_gene_name = $f->display_name();
        my $g_gene_name = $g->display_name();

        if (   exists( $seen_overlap{$f_gene_name}{$g_gene_name} )
            || exists( $seen_overlap{$g_gene_name}{$f_gene_name} ) )
        {
            return;
        }

        my $overlap_count;

        if ( $g->end() < $f->end() )
        {
            $overlap_count = $g->length();
        }
        else
        {
            $overlap_count = abs( $f->end() - $g->start() ) + 1;
        }

        my $overlap_percent_f
            = sprintf( "%.1f", ( $overlap_count / $f->length ) * 100 );
        my $overlap_percent_g
            = sprintf( "%.1f", ( $overlap_count / $g->length ) * 100 );

        if (   ( $overlap_count > $self->overlap_bp )
            || ( $overlap_percent_f > $self->overlap_percent )
            || ( $overlap_percent_g > $self->overlap_percent ) )
        {

            $f->add_tag_value( 'overlap' => $g_gene_name );
            $g->add_tag_value( 'overlap' => $f_gene_name );
            $seen_overlap{$f_gene_name}{$g_gene_name} = 1;
            $seen_overlap{$g_gene_name}{$f_gene_name} = 1;

        }

    };

    $self->find_overlaps( $overlap, $feature_ref );

}

sub find_overlaps
{
    my $self = shift;
    my ( $callback, $feature_ref ) = @_;

    unless ( ref($callback) eq 'CODE' )
    {
        croak 'callback arg is not a CODE ref';
    }

    foreach my $feature ( @{$feature_ref} )
    {

        unless ( $feature->isa('Bio::SeqFeatureI') )
        {
            croak 'feature does not implement Bio::SeqFeatureI';
        }

    }

    unless ( @{$feature_ref} > 1 )
    {
        return;
    }
    @{$feature_ref}
        = sort { $a->start() <=> $b->start() || $a->length() <=> $b->length() }
        @{$feature_ref};

OUTER: foreach my $o ( 0 .. $#$feature_ref )
    {

        my $f = $feature_ref->[$o];

        if ( $o < $#$feature_ref )
        {

        INNER1: foreach my $i ( $o + 1 .. $#$feature_ref )
            {

                my $g = $feature_ref->[$i];

                if ( $g->start() > $f->end() )
                {
                    last INNER1;
                }
                else
                {
                    $callback->( $f, $g );
                }

            }

        }

        if ( $o > 0 )
        {

        INNER2: foreach my $i ( reverse( 0 .. $o - 1 ) )
            {

                my $g = $feature_ref->[$i];

                if ( $g->end() < $f->start() )
                {
                    last INNER2;
                }
                else
                {
                    $callback->( $g, $f );
                }

            }

        }

    }

}

sub gene_to_feature
{
    my $self = shift;
    my ( $sequence_name, $coding_gene ) = @_;

    my $gene_feature = Bio::SeqFeature::Generic->new(
        -seq_id => $sequence_name,
        -start  => $coding_gene->start(),
        -end    => $coding_gene->end(),
        -strand => $coding_gene->strand(),
        -source => $coding_gene->source(),
    );

    $gene_feature->add_tag_value( 'gene_id', $coding_gene->gene_id() );
    $gene_feature->display_name( $coding_gene->gene_name() );
    $gene_feature->add_tag_value(
        'Sequence' => "Sequence " . $coding_gene->gene_name() );

    $gene_feature->primary_tag('Sequence');

    return $gene_feature;

}

sub blastp
{
    my $self = shift;
    my %rpc_args = %{$self->_rpc_args};
    my ( $fasta_file, $blast_db ) = @_;

    my $job_source = BAP::JobSource::Phase2BlastP->new( $blast_db, $fasta_file,
        $self->rpc_core_number);
    local $rpc_args{'job_source'} = $job_source;
    my $server = PP::RPC->new(%rpc_args);

    $server->start();

    check_failed_jobs($job_source);
    $self->_rpc_args(\%rpc_args);
    return $job_source->evidence();

}

sub iprscan
{
    my $self = shift;
    my ($fasta_file) = @_;
    #my %network_temp_args = %{$self->_network_temp_args};
    my %network_temp_args = %{$self->_network_temp_args};
    my %evidence = ();
    %network_temp_args = (
                       'TEMPLATE' => 'mgap_XXXXXXXX',
                       'DIR'      => '/gscmnt/temp212/info/annotation/BAP_tmp',
                       'SUFFIX'   => '.temp',
                    );

#???    $self->_
    #my $temp_fh = File::Temp->new( %network_temp_args, );
    my $temp_fh = File::Temp->new( 
                       'TEMPLATE' => 'mgap_XXXXXXXX',
                       'DIR'      => '/gscmnt/temp212/info/annotation/BAP_tmp',
                       'SUFFIX'   => '.temp',);
    my $temp_fn = $temp_fh->filename();
    $temp_fh->close();
    $self->debug_message("temp filename: $temp_fn");
    #my $err_fh = File::Temp->new( %network_temp_args, );
    my $err_fh = File::Temp->new( 
                       'TEMPLATE' => 'mgap_XXXXXXXX',
                       'DIR'      => '/gscmnt/temp212/info/annotation/BAP_tmp',
                       'SUFFIX'   => '.temp',);
    my $err_fn = $err_fh->filename();
    $err_fh->close();
    $self->debug_message("err fn: $err_fn");
    # originally hardcoded to these at various points in the past.
     #'/gscmnt/974/analysis/iprscan16.1/iprscan/bin/iprscan',
     #          '/gscmnt/974/analysis/iprscan16.1/iprscan/bin/iprscan.hacked',
     #   '/gsc/scripts/bin/iprscan',
     #   '/gscmnt/temp212/info/annotation/InterProScan/iprscan16.1/iprscan/bin/iprscan.hacked',    
    my @cmd = (
        #$self->iprpath,
		#"/gsc/scripts/pkg/bio/iprscan/iprscan-4.8/bin/iprscan",
		"/gsc/scripts/pkg/bio/iprscan/iprscan-". $self->ipr_version ."/bin/iprscan",
        '-cli',
        '-appl hmmpfam',
        '-goterms',
        '-verbose',
        '-iprlookup',
        '-seqtype p',
        '-format ebixml',
        "-i $fasta_file",
        "-o $temp_fn",
    );

    $self->debug_message("about to run interproscan");
    my $cmd = join( ' ', @cmd );
    $self->debug_message("ipr cmd: $cmd");
    my $pp = PP->run(
        pp_type => 'lsf',
        command => $cmd,
        q       => $self->rpc_queue,
        R       => "'rusage[mem=1024]'",
        e       => $err_fn,
#        o       => "iprscan-debug.out", # we should probably try to store this somewhere.
    );

    $pp->wait_on();
    $pp->update_stats();

    my $exit_status = $pp->stats->[2];
    $self->debug_message("iprscan done $exit_status");
    if ( $exit_status eq 'EXIT' )
    {
        carp "exit status returned as 'EXIT'. Sleeping";
        sleep(60);
    }
    $pp->update_stats();
    $exit_status = $pp->stats->[2];

    # should watch out for this?  it could be 'EXIT' as well
    unless ( $exit_status eq 'DONE' )
    {
        croak "iprscan failed with status '$exit_status', check $err_fn";
    }

    #unless (-z $err_fn) {
    #    my $err_text = slurp($err_fn);
    #    croak "iprscan is not happy: $err_text";
    #}

    my $seqio = Bio::SeqIO->new( -file => $temp_fn, -format => 'interpro' );

    while ( my $seq = $seqio->next_seq() )
    {

        my $accession = $seq->accession();

        unless ( defined($accession) ) { next; }

        my @features = $seq->get_SeqFeatures();

        if (@features) { $evidence{$accession} = 1; }

    }
    unlink $temp_fn;

    return \%evidence;

}

sub best_per_locus
{
    my $self = shift;
    #$self->debug_message("starting best per locus");
    my %selected_genes ;
    if(defined($self->_selected_genes)) {
        %selected_genes = %{$self->_selected_genes};
    }
    my $debug_fh = $self->_debug_fh;
    my ( $current_phase, $previous_phase ) = @_;

    my $sequence_set
        = BAP::DB::SequenceSet->retrieve( $self->sequence_set_id );

    my @sequences = $sequence_set->sequences();
    my @sequence_names = map { $_->sequence_name() } @sequences;
    #$self->debug_message("sequences retrieved, sanity checking...");
    $self->sanity_check_sequences( $current_phase, \@sequence_names );

    my @fetched_gene_names = ();

    #$self->debug_message("sanity checked, pulling genes out...");
    foreach my $sequence (@sequences)
    {

        my @coding_genes = ();

        if ( defined($previous_phase) )
        {
            @coding_genes = $sequence->coding_genes( $previous_phase => 1 );
            push @fetched_gene_names, map { $_->gene_name() } @coding_genes;
        }
        else
        {
            @coding_genes = $sequence->coding_genes();
        }

        my @gene_features = ();

        foreach my $coding_gene (@coding_genes)
        {

            my $gene_name = $coding_gene->gene_name();

            print $debug_fh
                join( "\t", $current_phase, $gene_name, "fetched from db" ),
                "\n";

            $coding_gene->$current_phase(0);
            $coding_gene->update();

            ## Ignore genes with bad translations.
            ## Ignore (nonfragment) genes with missing starts / stops
            if ( $coding_gene->internal_stops() )
            {
                print $debug_fh join( "\t",
                    $current_phase, $gene_name, "has internal stops" ),
                    "\n";
                next;
            }

            unless ( $coding_gene->fragment() )
            {
                if (   $coding_gene->missing_start()
                    || $coding_gene->missing_stop() )
                {
                    print $debug_fh join( "\t",
                        $current_phase, $gene_name, "missing start/stop" ),
                        "\n";
                    next;
                }
            }

            ## Ignore glimmer2 predictions, they're staying in the DB for possible
            ## data mining down the road, but at this point, we don't want them
            ## included in our gene sets
            if ( $coding_gene->source() eq 'glimmer2' )
            {
                print $debug_fh join( "\t",
                    $current_phase, $gene_name, "ignoring glimmer2" ),
                    "\n";
                next;
            }

            my $gene_feature
                = $self->gene_to_feature( $sequence->sequence_name(), $coding_gene );
            unless ( $gene_feature->length() >= $self->min_gene_length )
            {
                print $debug_fh
                    join( "\t", $current_phase, $gene_name, "length < 90bp" ),
                    "\n";
                next;
            }

            push @gene_features, $gene_feature;

        }

        my %predicted_by = ();

        foreach my $gene_feature (@gene_features)
        {

            my $source_tag = $gene_feature->source_tag();

            push @{ $predicted_by{$source_tag} }, $gene_feature;

        }

        foreach my $source ( keys %predicted_by )
        {
            #$self->debug_message("tagging redundant");
            $self->tag_redundant( \@{ $predicted_by{$source} } );
            #$self->debug_message("tagged redundant");

        }

        my @redundant_features
            = grep { $_->has_tag('redundant'); } @gene_features;

        foreach my $redundant_feature (@redundant_features)
        {
            my $gene_name = $redundant_feature->display_name();
            print $debug_fh
                join( "\t", $current_phase, $gene_name, "redundant" ), "\n";
        }

        @gene_features
            = grep { !( $_->has_tag('redundant') ); } @gene_features;
        #$self->debug_message("tagging preferred gene features");
        tag_preferred( \@gene_features );

        @gene_features = grep { $_->has_tag('preferred'); } @gene_features;

        my %this_phase = ();

        foreach my $gene_feature (@gene_features)
        {
            my $gene_name = $gene_feature->display_name();
            $this_phase{$gene_name} = 1;
        }

        foreach my $coding_gene (@coding_genes)
        {

            my $gene_name = $coding_gene->gene_name();

            if ( exists( $this_phase{ $coding_gene->gene_name() } ) )
            {
                $coding_gene->$current_phase(1);
                $coding_gene->update();
                print $debug_fh
                    join( "\t", $current_phase, $gene_name, "selected" ),
                    "\n";
                $selected_genes{$current_phase}{$gene_name} = 1;
            }

        }

    }
    $self->_selected_genes(\%selected_genes);
    #$self->_rfam_special(\%rfam_special);
    $self->sanity_check_genes( $current_phase, $previous_phase,
        \@fetched_gene_names );

    #$self->debug_message("finished best per locus");

    return 1;
}

sub check_overlapping
{
    my $self = shift;
    my $debug_fh = $self->_debug_fh;
    my %selected_genes = %{$self->_selected_genes};
    my %network_temp_args = %{$self->_network_temp_args};
    my %rfam_special = %{$self->_rfam_special};
    my ( $current_phase, $previous_phase ) = @_;
    %network_temp_args = (
                       'TEMPLATE' => 'mgap_XXXXXXXX',
                       'DIR'      => '/gscmnt/temp212/info/annotation/BAP_tmp',
                       'SUFFIX'   => '.temp',
                    );

    my $sequence_set
        = BAP::DB::SequenceSet->retrieve( $self->sequence_set_id );

    my @sequences         = $sequence_set->sequences();
    my @sequence_names    = map { $_->sequence_name() } @sequences;
    my @overlapping_genes = ();
    my @short_genes       = ();
    my @trna_overlaps     = ();

    my %delete         = ();
    my %rrna_delete    = ();
    my %trna_delete    = ();
    my %overlap_graph  = ();
    my %redundant_rrna = ();

    my @fetched_gene_names = ();

    $self->sanity_check_sequences( $current_phase, \@sequence_names );

    foreach my $sequence (@sequences)
    {

        my $seq = Bio::Seq->new(
            -seq => $sequence->sequence_string(),
            -id  => $sequence->sequence_name(),
        );

        my @coding_genes = $sequence->coding_genes( $previous_phase => 1 );

        push @fetched_gene_names, map { $_->gene_name() } @coding_genes;

        my @rna_genes  = $sequence->rna_genes();
        my @trna_genes = $sequence->trna_genes();

        my @gene_features      = ();
        my @rna_gene_features  = ();
        my @trna_gene_features = ();

        foreach my $coding_gene (@coding_genes)
        {
            my $gene_name = $coding_gene->gene_name();

            print $debug_fh
                join( "\t", $current_phase, $gene_name, "fetched from db" ),
                "\n";

            $coding_gene->$current_phase(0);
            $coding_gene->update();

            my $gene_feature
                = $self->gene_to_feature( $sequence->sequence_name(), $coding_gene );

            $gene_feature->attach_seq($seq);

            $gene_feature->add_tag_value( 'type' => 'coding' );

            push @gene_features, $gene_feature;

        }

        foreach my $rna_gene (@rna_genes)
        {

            my $acc = $rna_gene->acc();

            my $gene_name = $rna_gene->gene_name();

            if (   ( $acc eq 'RF00029' )
                || ( $acc eq 'RF00234' )
                || ( $acc eq 'RF00382' ) )
            {
                print $debug_fh join( "\t",
                    $current_phase, $gene_name, "skipping rRNA ($acc)" ),
                    "\n";
                next;
            }

            if ( $rna_gene->redundant() )
            {
                print $debug_fh join( "\t",
                    $current_phase, $gene_name, "skipping redundant rRNA" ),
                    "\n";
                next;
            }
            print $debug_fh join( "\t",
                $current_phase, $gene_name, "fetched from db (rRNA)" ),
                "\n";

            my $gene_feature
                = $self->gene_to_feature( $sequence->sequence_name(), $rna_gene );

            $gene_feature->attach_seq($seq);

            $gene_feature->add_tag_value( 'type' => 'rRNA' );

            if ( exists( $rfam_special{$acc} ) )
            {
                $gene_feature->add_tag_value( 'overlap_50' => 1 );
                print $debug_fh join( "\t",
                    $current_phase, $gene_name, "is a riboswitch/leader" ),
                    "\n";
            }

            push @rna_gene_features, $gene_feature;

        }
        BAP::GeneMerge::tag_redundant_rfam( \@rna_gene_features );

        my @nonredundant_rna_gene_features = ();

        foreach my $feature (@rna_gene_features)
        {
            if ( $feature->has_tag('redundant') )
            {
                $redundant_rrna{ $feature->display_name() } = 1;
            }
            else
            {
                push @nonredundant_rna_gene_features, $feature;
            }
        }

        @rna_gene_features = @nonredundant_rna_gene_features;

        foreach my $trna_gene (@trna_genes)
        {

            my $gene_name = $trna_gene->gene_name();

            print $debug_fh join( "\t",
                $current_phase, $gene_name, "fetched from db (tRNA)" ),
                "\n";

            my $gene_feature
                = $self->gene_to_feature( $sequence->sequence_name(), $trna_gene );

            $gene_feature->attach_seq($seq);

            $gene_feature->add_tag_value( 'type' => 'tRNA' );

            push @trna_gene_features, $gene_feature;

        }

        BAP::DB::DBI->dbi_commit();
        BAP::DB::DBI->db_Main->disconnect();

        BAP::GeneMerge::tag_rna_overlap(
            [ @gene_features, @rna_gene_features, @trna_gene_features, ] );
        foreach my $rrna_overlap ( grep { $_->has_tag('delete_rrna_overlap') }
            @gene_features )
        {
            my $gene_name = $rrna_overlap->display_name();
            print $debug_fh
                join( "\t", $current_phase, $gene_name, "overlaps rRNA" ),
                "\n";
            $rrna_delete{$gene_name} = 1;
        }

        foreach my $trna_overlap ( grep { $_->has_tag('delete_trna_overlap') }
            @gene_features )
        {
            my $gene_name = $trna_overlap->display_name();
            print $debug_fh
                join( "\t", $current_phase, $gene_name, "overlaps tRNA" ),
                "\n";
            $trna_delete{$gene_name} = 1;
        }

        @gene_features
            = grep { !$_->has_tag('delete_rrna_overlap'); } @gene_features;
        @gene_features
            = grep { !$_->has_tag('delete_trna_overlap'); } @gene_features;
        my @trna_overlap_features
            = grep { $_->has_tag('check_trna_overlap'); } @gene_features;
        my @short_gene_features = grep { $_->length() < 120; } @gene_features;

        if (@trna_overlap_features)
        {
            push @trna_overlaps, @trna_overlap_features;
        }

        if (@short_gene_features)
        {
            push @short_genes, @short_gene_features;
        }

        my $graph = BAP::GeneMerge::graph_overlaps( $self->overlap_bp,
            $self->overlap_percent, \@gene_features );

        my @overlapping_gene_features
            = grep { $graph->has_vertex( $_->display_name() ) }
            (@gene_features);

        $overlap_graph{ $sequence->sequence_name() } = $graph;

        if (@overlapping_gene_features)
        {
            push @overlapping_genes, @overlapping_gene_features;
        }

    }

    my $blastp_evidence  = {};
    my $iprscan_evidence = {};
    my %evidence         = ();

    if (   ( @overlapping_genes > 0 )
        || ( @short_genes > 0 )
        || (@trna_overlaps) )
    {

        #my $fasta_fh = File::Temp->new( %network_temp_args, );
        my $fasta_fh = File::Temp->new( 
                       'TEMPLATE' => 'mgap_XXXXXXXX',
                       'DIR'      => '/gscmnt/temp212/info/annotation/BAP_tmp',
                       'SUFFIX'   => '.temp',);
        my $fasta_fn = $fasta_fh->filename();
        my $seqio = Bio::SeqIO->new( -fh => $fasta_fh, -format => 'Fasta' );

        my %wrote_sequence = ();

        foreach my $gene ( @overlapping_genes, @short_genes, @trna_overlaps )
        {

            my $gene_name   = $gene->display_name();
            my $protein_seq = $gene->seq->translate();

            unless ( exists( $wrote_sequence{$gene_name} ) )
            {
                $protein_seq->display_id($gene_name);
                $seqio->write_seq($protein_seq);
                $wrote_sequence{$gene_name} = 1;
            }

            if ( $gene->has_tag('overlap') )
            {
                print $debug_fh join( "\t",
                    $current_phase, $gene_name, "overlaps another gene" ),
                    "\n";
            }

            if ( $gene->length() < 120 )
            {
                print $debug_fh join( "\t",
                    $current_phase, $gene_name, "shorter than 120bp" ),
                    "\n";
            }

            if ( $gene->has_tag('check_trna_overlap') )
            {
                print $debug_fh join( "\t",
                    $current_phase, $gene_name,
                    "overlaps tRNA on other strand" ),
                    "\n";
            }

        }

        $fasta_fh->close();
        $blastp_evidence = $self->blastp( $fasta_fn, $self->nr_db );
        $iprscan_evidence = $self->iprscan($fasta_fn);

        unlink $fasta_fn;

        %evidence = map { $_ => 1 }
            ( keys( %{$blastp_evidence} ), keys( %{$iprscan_evidence} ) );

    TRNAO: foreach my $gene (@trna_overlaps)
        {

            my $gene_name = $gene->display_name();

            if ( exists( $evidence{$gene_name} ) )
            {
                print $debug_fh join( "\t",
                    $current_phase, $gene_name,
                    "has evidence (trna_overlap)" ),
                    "\n";
                next TRNAO;
            }
            $trna_delete{$gene_name} = 1;

            my $graph = $overlap_graph{ $gene->seq_id() };

            if ( $graph->has_vertex($gene_name) )
            {
                $graph->delete_vertex($gene_name);
            }

        }

    SGENE: foreach my $gene (@short_genes)
        {

            my $gene_name   = $gene->display_name();
            my $gene_length = $gene->length();

            if (exists( $iprscan_evidence->{$gene_name} )
                || ( exists( $evidence{$gene_name} )
                    && ( $gene->length() >= 90 ) )
                )
            {
                print $debug_fh join( "\t",
                    $current_phase, $gene_name,
                    "has evidence ($gene_length bp)" ),
                    "\n";
                next SGENE;

            }

            $delete{$gene_name} = 1;

            my $graph = $overlap_graph{ $gene->seq_id() };

            if ( $graph->has_vertex($gene_name) )
            {
                $graph->delete_vertex($gene_name);
            }

        }

        my %overlaps_on = ();

        foreach my $gene (@overlapping_genes)
        {

            my $gene_name = $gene->display_name();
            my $graph     = $overlap_graph{ $gene->seq_id() };

            if ( $graph->has_vertex($gene_name) )
            {

                push @{ $overlaps_on{ $gene->seq_id } }, $gene;

                if ( exists( $iprscan_evidence->{$gene_name} ) )
                {
                    $gene->add_tag_value( 'pfam_evidence' => 1 );
                }

                if ( exists( $blastp_evidence->{$gene_name} ) )
                {
                    $gene->add_tag_value( 'blastp_evidence' => 1 );
                }

            }

        }
        @overlapping_genes = ();

        foreach my $seq_id ( keys %overlaps_on )
        {

            BAP::GeneMerge::fancy_tag_overlapping( $overlap_graph{$seq_id},
                $overlaps_on{$seq_id} );

            push @overlapping_genes, @{ $overlaps_on{$seq_id} };

        }

    OGENE: foreach my $gene (@overlapping_genes)
        {

            my $gene_name = $gene->display_name();

            if ( exists( $iprscan_evidence->{$gene_name} ) )
            {
                print $debug_fh join( "\t",
                    $current_phase, $gene_name, "has pfam evidence" ),
                    "\n";
            }
            elsif ( exists( $blastp_evidence->{$gene_name} ) )
            {
                print $debug_fh join( "\t",
                    $current_phase, $gene_name, "has blastp evidence" ),
                    "\n";
            }
            else
            {
                print $debug_fh join( "\t",
                    $current_phase, $gene_name, "has no evidence" ),
                    "\n";
            }

            if ( $gene->has_tag('delete_overlap') )
            {

                $delete{$gene_name} = 1;

                my ($delete_for) = $gene->each_tag_value('delete_for');

                my $msg = "discarded in favor of $delete_for";

                print $debug_fh
                    join( "\t", $current_phase, $gene_name, $msg ), "\n";

            }

        }

    }

    foreach my $sequence (@sequences)
    {

        my @coding_genes = $sequence->coding_genes( $previous_phase => 1 );

        foreach my $coding_gene (@coding_genes)
        {

            my $gene_name = $coding_gene->gene_name();

            if ( exists( $blastp_evidence->{$gene_name} ) )
            {
                $coding_gene->blastp_evidence(1);
                $coding_gene->update();
            }

            if ( exists( $iprscan_evidence->{$gene_name} ) )
            {
                $coding_gene->pfam_evidence(1);
                $coding_gene->update();
            }
            unless ( exists( $delete{$gene_name} )
                || exists( $rrna_delete{$gene_name} )
                || exists( $trna_delete{$gene_name} ) )
            {

                $coding_gene->$current_phase(1);
                $coding_gene->update();

                print $debug_fh
                    join( "\t", $current_phase, $gene_name, "selected" ),
                    "\n";

                $selected_genes{$current_phase}{$gene_name} = 1;

            }

        }

        my @rna_genes = $sequence->rna_genes();
        foreach my $rna_gene (@rna_genes)
        {

            my $gene_name = $rna_gene->gene_name();

            if ( exists( $redundant_rrna{$gene_name} ) )
            {

                $rna_gene->redundant(1);
                $rna_gene->update();

                print $debug_fh
                    join( "\t", $current_phase, $gene_name, "redundant" ),
                    "\n";
            }

        }

    }
    $self->sanity_check_genes( $current_phase, $previous_phase,
        \@fetched_gene_names );
    return 1
}

sub log_run
{
    my $self = shift;
    my ( $started, $finished, $sequence_set ) = @_;

    # should we just grab sequence_set via sequence_set_id?

    my $elapsed_seconds = ( $finished->epoch() - $started->epoch() );
    my $db_file         = BAP::Config->new()->activity_db_file();

    my $sequence_id   = $sequence_set->sequence_set_id();
    my $sequence_name = $sequence_set->sequence_set_name();

    my $organism = BAP::DB::Organism->retrieve(
        'organism_id' => $sequence_set->organism_id() )->organism_name();

    my $host = 'unknown';
    my $user = 'unknown';

    if ( exists( $ENV{LSB_HOSTS} ) )
    {
        $host = $ENV{LSB_HOSTS};
    }
    elsif ( exists( $ENV{HOST} ) )
    {
        $host = $ENV{HOST};
    }

    if (Genome::Sys->username)
    {
        $user = Genome::Sys->username;
    }
    elsif ( exists( $ENV{LOGIN} ) )
    {
        $user = $ENV{LOGIN};
    }
    elsif ( exists( $ENV{USERNAME} ) )
    {
        $user = $ENV{USERNAME};
    }

    my $dbh = DBI->connect( "dbi:SQLite:dbname=$db_file", '', '',
        { RaiseError => 1, AutoCommit => 1 } );
    my $sql = <<SQL;
INSERT INTO activity_log (activity, 
                          sequence_id, 
                          sequence_name, 
                          organism_name,
                          host,
                          user, 
                          started,
                          finished) 
                  VALUES (?,?,?,?,?,?,
                          strftime('%s', 'now') - $elapsed_seconds,
                          strftime('%s', 'now')
                         );
SQL

    $dbh->do( $sql, {}, 'merge', $sequence_id, $sequence_name, $organism,
        $host, $user );

}

sub sanity_check_genes
{
    my $self = shift;
    my %selected_genes = %{$self->_selected_genes};
    my ( $current_phase, $previous_phase, $gene_array_ref ) = @_;

    my %fetched = map { $_ => 1 } @{$gene_array_ref};

    foreach my $gene ( keys %{ $selected_genes{$previous_phase} } )
    {
        unless ( exists( $fetched{$gene} ) )
        {

            my $msg = <<MSG;
*******************************************************************************
* SANITY CHECK FAILED
*
* selected $gene in $previous_phase, but failed to fetch it in $current_phase
*******************************************************************************
MSG

            croak $msg;

        }

    }
    $self->_selected_genes(\%selected_genes);
    return 1;
}

sub sanity_check_sequences
{
    my $self = shift;
    my ( $current_phase, $sequence_ref ) = @_;
    my %fetched_sequences = %{$self->_fetched_sequences};

    my %fetched = map { $_ => 1 } @{$sequence_ref};

    foreach my $sequence ( keys %fetched_sequences )
    {

        unless ( exists( $fetched{$sequence} ) )
        {

            my $msg = <<MSG;
*******************************************************************************
* SANITY CHECK FAILED
* fetched $sequence at startup, but failed to fetch it in $current_phase
*******************************************************************************
MSG

            croak $msg;

        }

    }

    foreach my $sequence ( keys %fetched )
    {

        unless ( exists( $fetched_sequences{$sequence} ) )
        {
            my $msg = <<MSG;
*******************************************************************************
* SANITY CHECK FAILED
* fetched $sequence in $current_phase, but did not fetch it at startup
*******************************************************************************
MSG

            croak $msg;

        }

    }
    return 1;
}

1;
