package Genome::Model::GenePrediction::Command::Bacterial::Predict;

#use lib '/gsc/scripts/opt/bacterial-bioperl';

use strict;
use warnings;

use Genome;

use BAP::Config;
use BAP::DB::CodingGene;
use BAP::DB::Organism;
use BAP::DB::Sequence;
use BAP::DB::SequenceSet;
use BAP::DB::RNAGene;
use BAP::DB::tRNAGene;
use GAP::JobSource::Composite;
use BAP::JobSource::Genemark;
use BAP::JobSource::Glimmer3;
use GAP::JobSource::tRNAscan;
use GAP::JobSource::RfamScan;
use GAP::JobSource::RNAmmer;
use GAP::JobSource::MetaRna;

use Bio::Seq;
use Bio::SeqIO;

use PP::RPC;

use Carp;
use DateTime;
use DateTime::Duration;
use DBI;
use DBD::SQLite;
use English;
use File::Temp;
use Getopt::Long;
use Genome::Utility::Email;
use Pod::Usage;

class Genome::Model::GenePrediction::Command::Bacterial::Predict {
    is  => 'Command',
    doc => "",
    has => [
        sequence_set_id => {
            is  => 'Integer',
            doc => "sequence set id for the genome",
        },
        glimmer3_model => {
            is  => 'String',
            doc => "path to glimmer3 model file",
        },
        glimmer3_pwm => {
            is  => 'String',
            doc => "path to glimmer3 pwm",
        },
        genemark_model => {
            is  => 'String',
            doc => "path to genemark model file",
        },
        domain => {
            is  => 'String',
            doc => "domain of genome (either bacteria, archaea, virial)",
        },

    ],
    has_optional => [
        no_mail => {
            is      => 'Boolean',
            doc     => "turn off mailing completion report",
            default => 0,
        },
        dev => {
            is      => 'Boolean',
            doc     => "use development databases",
            default => 0,
        },
        circular_dna => {
            is      => 'Boolean',
            doc     => "genome's dna is circular",
            default => 0,
        },
        tmp_usage => {
            is      => 'Integer',
            doc     => "megabytes to reserved for tmp space",
            default => 100,
        },
        runner_count => {
            is      => 'Integer',
            doc     => "number of runner jobs to spawn to lsf",
            default => 10,
        },
        job_stdout => {
            is  => 'String',
            doc => "file to write out stdout to",
        },
        job_stderr => {
            is  => 'String',
            doc => "file to write out stderr to",
        },
    ],

};

sub help_brief {
    "run bacterial gene prediction";
}

sub help_synopsis {
return <<EOS
genome model gene-prediction predict [options]
EOS
}

sub help_detail {
    return <<EOS
Used to run Rfam, RNAmmer, GeneMark and Glimmer3 in BAP/MGAP

This script is intented to be used with the HGMI/bacterial annotation pipeline (bap/mgap).  This program runs  Rfam, RNAmmer, Glimmer3 and GeneMark gene predictors and stores the output into our oracle database.

EOS
}

sub execute
{
    my $self           = shift;
    my $glimmer3_model = $self->glimmer3_model;
    my $glimmer3_pwm   = $self->glimmer3_pwm;
    my $genemark_model = $self->genemark_model;
    my $sequence_set_id = $self->sequence_set_id;
    my $circular_dna_flag = $self->circular_dna;
    my $tmp_usage = $self->tmp_usage;
    my $runner_count = $self->runner_count;
    my $skip_mail_flag = $self->no_mail;
    my $user = Genome::Sys->username;
    my $dt_started = mark_time();

    if ($self->dev) { $BAP::DB::DBI::db_env = 'dev'; }

    unless ( -e $glimmer3_model )
    {
        croak "glimmer3-model '$glimmer3_model' is empty";
    }
    unless ( -e $glimmer3_pwm )
    {
        croak "glimmer3-model '$glimmer3_pwm' is empty";
    }
    unless ( -e $genemark_model )
    {
        croak "genemark-model '$genemark_model' is empty";
    }

    unless ( ( $self->domain eq 'archaea' )
        || ( $self->domain eq 'bacteria' )
        || ( $self->domain eq 'virial' ) )
    {
        croak "domain '".$self->domain."' is invalid";
    }

    foreach my $out_file ( $self->job_stdout, $self->job_stderr )
    {

        if ( defined($out_file) && ( -e $out_file ) )
        {
            unlink($out_file) || carp "Failed to unlink '$out_file'!";
        }

    }

    my %gene_count = ();

    my %source_fixup = (
        'glimmer3' => 'Glimmer3',
        'genemark' => 'GeneMark',
        'trnascan' => 'tRNAscan',
        'rfam'     => 'rfam',
        'meta_rna' => 'meta_rna',
        'rnammer'  => 'rnammer',
    );

    my $svn_version = '$Revision$';
    $svn_version =~ s/\D+//g;
    my $bap_version = BAP::Config->new()->version();

    my $allfasta_fh    = File::Temp->new();
    my $gt40fasta_fh   = File::Temp->new();
    my $allfasta_file  = $allfasta_fh->filename();
    my $gt40fasta_file = $gt40fasta_fh->filename();

    my $allseqstream
        = Bio::SeqIO->new( -fh => $allfasta_fh, -format => 'Fasta' );
    my $gt40seqstream
        = Bio::SeqIO->new( -fh => $gt40fasta_fh, -format => 'Fasta' );

    my $sequence_set = BAP::DB::SequenceSet->retrieve($sequence_set_id);

    $sequence_set->software_version($bap_version);
    $sequence_set->update();

    my $sequence_set_name = $sequence_set->sequence_set_name();
    my @sequences         = $sequence_set->sequences();

    foreach my $sequence (@sequences)
    {

        my $seq = Bio::Seq->new(
            -seq => $sequence->sequence_string(),
            -id  => $sequence->sequence_name(),
        );
        if ( $seq->length() >= 41 )
        {

            $gt40seqstream->write_seq($seq);

        }

        $allseqstream->write_seq($seq);

        $gene_count{ $sequence->sequence_name() } = {
            'glimmer3' => 0,
            'genemark' => 0,
            'trnascan' => 0,
            'rfam'     => 0,
            'meta_rna' => 0,
            'rnammer'  => 0,
        };

        # Blow up any predictions from previous runs
        # (proteins should be automagically deleted
        #  due to the has_many relationship from
        #  coding_gene to protein)
        my @coding_genes = $sequence->coding_genes();
        my @trna_genes   = $sequence->trna_genes();
        my @rna_genes    = $sequence->rna_genes();

        foreach my $gene ( @coding_genes, @trna_genes, @rna_genes )
        {
            $gene->delete();
        }

    }

    BAP::DB::DBI->dbi_commit();

    $allseqstream->close();
    $gt40seqstream->close();

    close($allfasta_fh);
    close($gt40fasta_fh);

    my $genemark_job_source
        = BAP::JobSource::Genemark->new(
        Bio::SeqIO->new( -file => $gt40fasta_file, -format => 'Fasta' ),
        $genemark_model, $circular_dna_flag, );
    my $glimmer3_job_source
        = BAP::JobSource::Glimmer3->new(
        Bio::SeqIO->new( -file => $gt40fasta_file, -format => 'Fasta' ),
        $glimmer3_model, $glimmer3_pwm, $circular_dna_flag, );
    my $trnascan_job_source
        = GAP::JobSource::tRNAscan->new(
        Bio::SeqIO->new( -file => $allfasta_file, -format => 'Fasta' ),
        $self->domain, );
    my $metarna_job_source
        = GAP::JobSource::MetaRna->new(
        Bio::SeqIO->new( -file => $allfasta_file, -format => 'Fasta' ),
        $self->domain, );
    my $rfamscan_job_source
        = GAP::JobSource::RfamScan->new(
        Bio::SeqIO->new( -file => $allfasta_file, -format => 'Fasta' ),
        );
    my $rnammer_job_source
        = GAP::JobSource::RNAmmer->new(
        Bio::SeqIO->new( -file => $allfasta_file, -format => 'Fasta' ),
        $self->domain, );

    my $job_source = GAP::JobSource::Composite->new(
        $genemark_job_source, $glimmer3_job_source, $trnascan_job_source,
        $metarna_job_source,  $rfamscan_job_source, $rnammer_job_source,
    );

    my %rpc_args = (
        job_source   => $job_source,
        runner_count => $runner_count,
        app_init     => 0,
        port         => 7654 + $PID,
        pp_type      => 'lsf',
        R            => "rusage[tmp=$tmp_usage]",
        q            => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        maxmessage   => 16384000,
        lib_paths    => [
            UR::Util::used_libs,
            '/gsc/scripts/opt/bacterial-bioperl',
        ],
    );

    if ( defined($self->job_stdout) ) { $rpc_args{'o'} = $self->job_stdout; }
    if ( defined($self->job_stderr) ) { $rpc_args{'e'} = $self->job_stderr; }

    my $server = PP::RPC->new(%rpc_args);

    $self->debug_message("Starting rpc server");
    $server->start();

    my $failed_jobs = $job_source->failed_jobs();

    if ($#$failed_jobs > 0) {
        $self->warning_message("Found failed jobs! Here are some details:");
        my $msg;
        foreach my $job ( @{$failed_jobs}) { 
            my $job_id  = $job->job_id();
            my $class   = ref($job);
            my $seq_obj = $job->seq();
            my $seq_id  = $seq_obj->display_id();
            $msg .= "job $job_id ($class) ($seq_id) failed";
        }
        $self->warning_message($msg . "\nRolling back database changes!");
        BAP::DB::DBI->dbi_rollback();
        confess;
    }

    $self->debug_message("RPC server has run all jobs, none appear to have failed!");

    my %features = ();

    foreach my $feature (@{$job_source->feature_ref()}){
        push @{ $features{ $feature->seq_id() } }, $feature;
    }

    $self->debug_message("Iterating through sequences produced by predictors...");

    SEQUENCE: foreach my $seq_id ( keys %features ) {

        $self->debug_message("Retrieving sequence with set id $sequence_set_id and name $seq_id");
        my $sequence = BAP::DB::Sequence->retrieve(
            sequence_set_id => $sequence_set_id,
            sequence_name   => $seq_id,
        );

        my $seq_obj = Bio::Seq->new(
            -seq => $sequence->sequence_string(),
            -id  => $sequence->sequence_name(),
        );

        unless ($sequence) {
            $self->warning_message("Failed to fetch sequence object from BAP DB with id $sequence_set_id and name $seq_id!");
            next SEQUENCE;
        }

        my @features = @{ $features{$seq_id} };

        FEATURE: foreach my $feature (@features)
        {
            unless ($feature->isa('Bio::SeqFeatureI')) {
                $self->warning_message("This feature does not implement Bio::SeqFeatureI!");
                next FEATURE;
            }

            my $source = $feature->source_tag();

            if ( $source eq 'Glimmer_3.X' ) { $source = 'glimmer3'; }
            elsif ( $source eq 'Genemark.hmm.pro' )
            {
                $source = 'genemark';
            }
            elsif ( $source eq 'tRNAscan-SE' ) { $source = 'trnascan'; }
            elsif ( $source eq 'Infernal' )    { $source = 'rfam'; }
            elsif ( $source =~ /rnammer/i ) { $source = 'rnammer'; }
            else
            {
                $self->warning_message("unknown source_tag '$source'");
                next FEATURE;
            }

            # GeneMark emits Bio::Tools::Prediction::Gene features, not
            # Bio::SeqFeature::Generic features like Glimmer3,
            # so, grab the Bio::Tools::Prediction::Exon and discard the
            # rest (the exon will have the fuzzy location, if there is one)
            if ( $source eq 'genemark' )
            {

                ($feature) = $feature->exons();

                unless ( defined($feature) )
                {
                    warn "genemark prediction had no exons";
                    next FEATURE;
                }

            }

            my $strand   = $feature->strand();
            my $start    = $feature->start();
            my $end      = $feature->end();
            my $location = $feature->location();

            unless ( ( $source eq 'glimmer3' )
                || ( $source eq 'genemark' ) )
            {

                # Coords for GeneMark/Glimmer will get fixed later
                if ( $start > $end ) {
                    ( $start, $end ) = ( $end, $start );
                }

            }

            unless ( defined($seq_id) ) { $seq_id = 'undef'; }

            my $feature_id
            = join( ':', $seq_id, $source, $strand, $start, $end );

            if ( $seq_id eq 'undef' )
            {
                $self->warning_message('($feature_id) - undefined seq_id!');
                next FEATURE;
            }

            $gene_count{$seq_id}{$source}++;

            my $gene_name;
            if ( $source eq 'trnascan' )
            {
                $gene_name = join( '.',
                    $seq_id,
                    ( join( '', 't', $gene_count{$seq_id}{$source} ) ) );
            }
            else
            {
                $gene_name = join '.', $seq_id, $source_fixup{$source},
                $gene_count{$seq_id}{$source};
            }

            if (   ( $source eq 'glimmer3' )
                || ( $source eq 'genemark' ) )
            {

                my $internal_stops = 0;
                my $fragment       = 0;
                my $missing_start  = 0;
                my $missing_stop   = 0;
                my $wraparound     = 0;
                if ( $feature->has_tag('wraparound') ) {
                    $wraparound = 1;
                }
                if ( $location->isa('Bio::Location::Fuzzy') )
                {
                    $fragment = 1;
                }

# Fractional codons mean the reading frame is horked, and presently the
# upstream parsers aren't passing along frame info, so we have to do our best to
# reverse engineer the correct one
                unless ( ( ( abs( $end - $start ) + 1 ) % 3 ) == 0 )
                {

# Without a Bio::Location::Fuzzy, we cannot determine the correct reading frame
                    unless ( $location->isa('Bio::Location::Fuzzy') )
                    {
                        warn
                        "($feature_id) - length is not a multiple of 3 and location is not fuzzy";
                        next FEATURE;
                    }
                    my $fuzzy_start = 0;
                    my $fuzzy_end   = 0;
                    my $extra_bases
                    = ( ( abs( $end - $start ) + 1 ) % 3 );
                    my $start_pos_type = $location->start_pos_type();
                    my $end_pos_type   = $location->end_pos_type();

                    if (   $start_pos_type eq 'BEFORE'
                        || $start_pos_type eq 'AFTER' )
                    {
                        $fuzzy_start = 1;
                    }

                    if (   $end_pos_type eq 'AFTER'
                        || $end_pos_type eq 'BEFORE' )
                    {
                        $fuzzy_end = 1;
                    }

                    # We have to have one inexact end to shave bases off
                    unless ( $fuzzy_start || $fuzzy_end )
                    {
                        warn
                        "($feature_id) - location is fuzzy, but start and stop are both exact";
                        next FEATURE;
                    }

# With two inexact ends, there is no hope of restoring the correct reading frame
# intended by the predictor (we don't have an exact anchor point)
                    if ( $fuzzy_start && $fuzzy_end )
                    {
                        warn
                        "($feature_id) - both start and end are fuzzy";
                        next FEATURE;
                    }

                    if ($fuzzy_start)
                    {
                        if ( $start_pos_type eq 'BEFORE' )
                        {
                            $start += $extra_bases;
                        }
                        if ( $start_pos_type eq 'AFTER' )
                        {
                            $start -= $extra_bases;
                        }
                    }

                    if ($fuzzy_end)
                    {
                        if ( $end_pos_type eq 'BEFORE' )
                        {
                            $end += $extra_bases;
                        }
                        if ( $end_pos_type eq 'AFTER' )
                        {
                            $end -= $extra_bases;
                        }
                    }

                }

# The GeneMark parser emits SeqFeatures for minus strand genes with start < end,
# the Glimmer parser emits SeqFeatures for minus strand genes with start > end,
# now that we're done screwing with them, flip them to be consistent
                if ( $start > $end ) {
                    ( $start, $end ) = ( $end, $start );
                }

                my $gene_seq_obj = $seq_obj->trunc( $start, $end );

                unless ( $strand > 0 )
                {
                    $gene_seq_obj = $gene_seq_obj->revcom();
                }

                my $protein_seq_obj = $gene_seq_obj->translate();

                if ( $protein_seq_obj->seq() =~ /\*.+/ )
                {
                    $internal_stops = 1;
                }

                unless ( $protein_seq_obj->seq() =~ /\*$/ )
                {
                    $missing_stop = 1;
                }

                my $first_codon = substr( $gene_seq_obj->seq(), 0, 3 );
                unless ( $first_codon =~ /tg$/i ) { $missing_start = 1; }

                my $coding_gene_obj = BAP::DB::CodingGene->insert(
                    {   gene_name       => $gene_name,
                        sequence_id     => $sequence->sequence_id(),
                        start           => $start,
                        end             => $end,
                        strand          => $strand,
                        source          => $source,
                        sequence_string => $gene_seq_obj->seq(),
                        internal_stops  => $internal_stops,
                        missing_start   => $missing_start,
                        missing_stop    => $missing_stop,
                        fragment        => $fragment,
                        wraparound      => $wraparound,
                        phase_0         => 1,
                        phase_1         => 0,
                        phase_2         => 0,
                        phase_3         => 0,
                        phase_4         => 0,
                        phase_5         => 0,
                    }
                );

                my $protein_name = $gene_name;
                my $gene_id      = $coding_gene_obj->gene_id();

                BAP::DB::Protein->insert(
                    {   protein_name    => $protein_name,
                        gene_id         => $gene_id,
                        sequence_string => $protein_seq_obj->seq(),
                        internal_stops  => $internal_stops,
                    }
                );

            }
            elsif ( $source eq 'trnascan' )
            {

                my $score = $feature->score();

                my ($codon) = $feature->each_tag_value('Codon');
                my ($aa)    = $feature->each_tag_value('AminoAcid');

                BAP::DB::tRNAGene->insert(
                    {   gene_name   => $gene_name,
                        sequence_id => $sequence->sequence_id(),
                        start       => $start,
                        end         => $end,
                        strand      => $strand,
                        source      => $source,
                        score       => $score,
                        codon       => $codon,
                        aa          => $aa,
                    }
                );

            }
            elsif ( $source eq 'rfam' )
            {

                my $score              = $feature->score();
                my ($rfam_accession)   = $feature->each_tag_value('acc');
                my ($rfam_description) = $feature->each_tag_value('id');
                my ($rfam_product)
                = $feature->each_tag_value('rfam_prod');

                if ( $rfam_description =~ /tRNA/i ) { next FEATURE; }
                if ( $score <= 50 ) { next FEATURE; }

                BAP::DB::RNAGene->insert(
                    {   gene_name   => $gene_name,
                        sequence_id => $sequence->sequence_id(),
                        start       => $start,
                        end         => $end,
                        acc         => $rfam_accession,
                        desc        => $rfam_description,
                        strand      => $strand,
                        source      => $source,
                        score       => $score,
                        rfam_prod   => $rfam_product,
                    }
                );

            }
            elsif ( $source eq 'rnammer' )
            {

                my $score = $feature->score();
                my ($description) = $feature->each_tag_value('group');

                BAP::DB::RNAGene->insert(
                    {   gene_name   => $gene_name,
                        sequence_id => $sequence->sequence_id(),
                        start       => $start,
                        end         => $end,
                        acc         => 'RNAmmer',
                        desc        => $description,
                        strand      => $strand,
                        source      => $source,
                        score       => $score,
                    }
                );
            }

        }

    }

    $self->debug_message("Committing to BAP database!");
    BAP::DB::DBI->dbi_commit();
    my $dt_finished = $self->mark_time();

    unless ($self->dev) { $self->log_run( $dt_started, $dt_finished, $sequence_set ); }

    unless ($skip_mail_flag)
    {
        $self->send_mail( $sequence_set_id, $sequence_set_name, $dt_started,
            $dt_finished, $bap_version, $svn_version, $user, );
    }

    return 1;
}


sub handle_sigint
{
    carp "aborting - caught SIGINT in the face";

    BAP::DB::DBI->dbi_rollback();
    exit(1);
}

sub mark_time
{
    my $self = shift;
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
    my $duration        = DateTime::Duration->new( $finished - $started );
    my $minutes_running = $duration->in_units('minutes');
    my $from     = Genome::Utility::Email::construct_address($user);

    my $subject = "BPG script mail for MGAP SSID: $ss_id ($ss_name)";

    my $body = <<BODY;
The bap_predict_genes script ($bap_version/$svn_version) has finished running for MGAP SSID: $ss_id ($ss_name).
The job began at $date_started and ended at $date_finished.
The total run time of the script is:  $minutes_running  minutes.
BODY

    Genome::Utility::Email::send(
        from    => $from,
        to      => [ $from, 'kpepin@watson.wustl.edu'], #FIXME remove hardcoded address
        subject => $subject,
        body    => $body,
    );

    return 1;
}

sub log_run
{
    my $self = shift;
    $self->debug_message("Inserting logging information into database.");
    my ( $started, $finished , $sequence_set) = @_;
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

    $dbh->do( $sql, {}, 'predict', $sequence_id, $sequence_name, $organism,
        $host, $user );

    $self->debug_message("Done logging!");
    return 1;
}

1;
