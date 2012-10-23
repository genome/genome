package Genome::Model::Tools::Hgmi::RrnaScreen;

use strict;
use warnings;

use Genome;
use Command;
use Carp;
use IPC::Run;
use File::Temp qw( tempdir tempfile );
use File::Slurp;

use BAP::DB::DBI;
use BAP::DB::CodingGene;
use BAP::DB::GeneTag;
use BAP::DB::Tag;

use BPlite;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is         => 'Command',
    has        => [

        'sequence_set_id' => {
            is  => 'String',
            doc => "sequence set id of organism",
        },
        'dev' => {
            is          => 'Boolean',
            doc         => "development flag",
            is_optional => 1,
            default     => 0,
        },
        'rrna_database' => {
            is          => 'String',
            doc         => "rrna database",
            is_optional => 1,
            default => '/gscmnt/278/analysis/HGMI/rRNA_testing/16s_23srnadb',
        },
        'script_location' => {
            is          => 'String',
            doc         => "location of bap_rrna_screen",
            is_optional => 1,
            default     => "/gsc/scripts/bin/bap_rrna_screen",
        }
    ],
);

sub help_brief
{
    "Runs the rrna screen after everything has been finished.";
}

sub help_synopsis
{
    my $self = shift;
    return <<"EOS"
need to put help synopsis here
EOS
}

sub help_detail
{
    my $self = shift;
    return <<"EOS"
RNAMMER usually has trouble detecting rrna gene fragments in assemblies,
thus this tool for screening gene sets for genes that are really rrna 
fragments.  These genes will get tagged in mgap as 'Dead'.
EOS
}


sub execute
{
    my $self = shift;
    my $rrnascreen_script = $self->script_location();
    $self->run_screen();
    return 1;
}

sub run_screen
{
    my $self = shift;

    my $sequence_set_id = $self->sequence_set_id;
    my $dev             = $self->dev;

    # changed to default.
    #my $rrnadb = '/gscmnt/278/analysis/HGMI/rRNA_testing/16s_23srnadb';
    my $rrnadb = $self->rrna_database;

    if ($dev)
    {
        $BAP::DB::DBI::db_env = 'dev';
    }

    # dump peps out,

    my @pep_export = (
        'gmt', 'bacterial', 'export-transcripts', '--sequence-set-id',
        $sequence_set_id,
    );

    if ($dev)
    {
        push( @pep_export, '--dev' );
    }

    my ( $fasta_data, $stderr );
    IPC::Run::run( \@pep_export, '>', \$fasta_data, '2>', \$stderr ) or croak;

    my ( $tmpfh, $tmpfasta )
        = tempfile( "rrna_screen_fastaXXXXXX", TEMPDIR => './' );
    write_file( $tmpfasta, $fasta_data );

    my ( $tmpfh2, $tmpblast )
        = tempfile( "rrna_screen_blastXXXXXX", TEMPDIR => './' );

    # run rrna screen
    # this needs to goto the blades, -K option?
    #my @screen = ('blastn',
    #              $rrnadb,
    #              $tmpfasta,
    #              'B=1',
    #              'V=1',
    #              '-o',
    #              $tmpblast);
    my ( $bsubh, $bsub )
        = tempfile( "rrna_screen_bsubXXXXXX", TEMPDIR => './' );
    my $bsubout = $bsub . ".out";
    my $bsuberr = $bsub . ".err";
    my @screen  = (
        'bsub',   '-K',    '-o',      $bsubout, '-e',  $bsuberr,
        'blastn', $rrnadb, $tmpfasta, 'B=1',    'V=1', '-o',
        $tmpblast
    );

    my $stdout;
    $stderr = undef;

    IPC::Run::run( \@screen, '>', \$stdout, '2>', \$stderr, )
        or croak "problem running screen step\n$stderr";

    # parse output
    #    my @parse = ('bap_parse_rrna_screen',
    #                 '-input',
    #                 $tmpblast,
    #                 '-output',
    #                 $tmpblast."parsed",
    #                 '-num_hits',
    #                 1,
    #                 '-percent',
    #                 70,
    #                 '-fol',
    #                 0.7,
    #                 '-query',
    #                 $tmpfasta
    #                 );
    #    $stderr = undef;
    #
##    IPC::Run::run(\@parse,
    #                  '>',
    #                  \$stdout,
    #                  '2>',
    #                  \$stderr,
    #                  ) or croak "can't parse blast output";
    #    #                     input, output, num_hits, percent, fol, query
    $self->parse_screen_results( $tmpblast, $tmpblast . "parsed",
        1, 70, 0.7, $tmpfasta );

    my $rnahits;
    $stderr = undef;
    my @parsed = read_file( $tmpblast . "parsed" );
    my @hits = grep {/====/} @parsed;
    my @totag;
    foreach my $hit (@hits)
    {
        my @f = split( " ", $hit );
        push( @totag, $f[1] );

    }

    print STDERR "found ", scalar @totag, " genes to tag dead...\n";

    # tag hits.

    my ($tag) = BAP::DB::Tag->search(
        {   tag_name  => 'Dead',
            tag_value => 'rrna hit'
        }
    );
    foreach my $gene2tag (@totag)
    {

        # get gene id from coding_gene
        # create row in gene_tag
        my ($coding_gene)
            = BAP::DB::CodingGene->search( { gene_name => $gene2tag } );
        if ($coding_gene)
        {
            my $genetag = BAP::DB::GeneTag->find_or_create(
                {   gene_id => $coding_gene,
                    tag_id  => $tag
                }
            );
        }

        my $phase5_gene = $gene2tag;
        $phase5_gene =~ s/\.(\d+)$/\.p5_hybrid\.$1/;
        write_file(
            'delete.rrnahits.ace',
            { append => 1 },
            "Sequence $phase5_gene\nDead \"hits rRNA gene\"\n\n"
        );
    }
    unlink($tmpblast);
    unlink($tmpblast."parsed"); # keep for now?
    unlink($tmpfasta);
    unlink($bsubout);
    unlink($bsuberr);
    unlink($bsub);

    BAP::DB::DBI->dbi_commit();
    return 1;
}

sub parse_screen_results
{
    my $self        = shift;
    my $input       = shift;
    my $output_file = shift;
    my $num_hits    = shift;
    my $percent     = shift;
    my $fol         = shift;
    my $query_file  = shift;

    #my @parse = ('bap_parse_rrna_screen',
    #             '-input',
    #             $tmpblast,
    #             '-output',
    #             $tmpblast."parsed",
    #             '-num_hits',
    #             1,
    #             '-percent',
    #             70,
    #             '-fol',
    #             0.7,
    #             '-query',
    #             $tmpfasta
    #             );
    # below gutted from bap_parse_rrna_screen
    #Parse query file and get lengths of query
    my $query = new FileHandle "$query_file";
    my %query_sizes;
    my $current_name;
    while (<$query>)
    {
        chomp;
        my $line = $_;
        next if ( $line =~ /^\s*$/ );
        if ( $line =~ /^>/ )
        {
            ####$line =~ s/^>gi\|\d+\|//; #Get rid of the gi # part of name, since blast (without -gi switch) will not have it either
            $line =~ s/^\s+//;
            $line =~ s/\s+$//;
            $current_name = $line;
            $current_name =~ s/^>//;
            my @current_name = split( /\s+/, $current_name );
            $current_name = $current_name[0];
        }
        else
        {
            $query_sizes{$current_name} += length($line);
        }
    }

    #Parse blast file
    my $data_hash;
    my $blast           = new FileHandle($input);
    my $multiple_report = new BPlite::Multi( \*$blast );
    while ( my $report = $multiple_report->nextReport )
    {
        my $name = $report->query;
        while ( my $subject = $report->nextSbjct )
        {
            my $subject_name = $subject->name;
            my $best_length  = 0;
            while ( my $hsp = $subject->nextHSP )
            {

#           Collect data in %data_hash ----> NOTE: WE'RE ONLY GRABBING THE LONGEST HSP HERE IF MORE THAN ONE
#           But for Makedonka's purposes it shouldn't matter
                my $current_length = $hsp->length;
                if ( $current_length > $best_length )
                {
                    $best_length = $current_length;
                    $data_hash->{$name}->{$subject_name}->{pvalue} = $hsp->P;
                    $data_hash->{$name}->{$subject_name}->{match}
                        = $hsp->percent;
                    $data_hash->{$name}->{$subject_name}->{sstart}
                        = $hsp->sbjctBegin;
                    $data_hash->{$name}->{$subject_name}->{send}
                        = $hsp->sbjctEnd;
                    $data_hash->{$name}->{$subject_name}->{qstart}
                        = $hsp->queryBegin;
                    $data_hash->{$name}->{$subject_name}->{qend}
                        = $hsp->queryEnd;
                    $data_hash->{$name}->{$subject_name}->{bitscore}
                        = $hsp->bits;
                }
            }
        }
    }
    $blast->close;

    #Prepare output file
    my $output = new FileHandle(">>$output_file");

#Sort data_hash and output the number of top hits requested by user (with min e-value constraint)
    my $i;
    print $output
        "<p_value>  <bitscore>  <%match>  <query start>  <query end>  <subject start>  <subject end>  <subject name>\n";
    print $output
        "(NOTE:displaying only values for longest HSP in each subject)\n";
    foreach my $hit ( keys( %{$data_hash} ) )
    {
        $i = 0;
        foreach my $sbjct (
            sort
            {
                $data_hash->{$hit}->{$a}->{pvalue} <=> $data_hash->{$hit}
                    ->{$b}->{pvalue}
            } keys( %{ $data_hash->{$hit} } )
            )
        {
            if ( $i < $num_hits )
            {
                $i++;
                my $pvalue = $data_hash->{$hit}->{$sbjct}->{pvalue};
                my $match  = $data_hash->{$hit}->{$sbjct}->{match};
                my $qstart = $data_hash->{$hit}->{$sbjct}->{qstart};
                my $qend   = $data_hash->{$hit}->{$sbjct}->{qend};
                my $sstart = $data_hash->{$hit}->{$sbjct}->{sstart};
                my $send   = $data_hash->{$hit}->{$sbjct}->{send};
                my $bits   = $data_hash->{$hit}->{$sbjct}->{bitscore};
                my $query_aligned_length = abs( $qend - $qstart );
                my $clean_hit            = $hit;

                if ( $clean_hit =~ /WARNING/ )
                {
                    my $stupidvar = 1;
                }

                my @clean_hit = split( /\(/, $clean_hit );
                $clean_hit = shift(@clean_hit);

                #$clean_hit =~ s/\(\d+\s+letters\;\s+record\s+\d+\)//;
                #$clean_hit =~ s/\(\d+\s+letters\)//;
                $clean_hit =~ s/^\s+//;
                $clean_hit =~ s/\s+$//;
                my @c_h = split( /\s+/, $clean_hit );
                $clean_hit = $c_h[0];

                unless (( defined $query_sizes{$clean_hit} )
                    and ( defined $fol ) )
                {
                    my $stupidvar = 1;
                }

                my $min_needed_alignment_length
                    = int( ( $query_sizes{$clean_hit} * $fol ) + .5 );
                $sbjct =~ s/^>//;
                if ( $query_aligned_length >= $min_needed_alignment_length )
                {
                    if ( $percent <= $match )
                    {
                        if ( $i == 1 )
                        {
                            print $output "========== $hit ==========";
                        }
                        print $output
                            "$pvalue\t$bits\t$match\t$qstart\t$qend\t$sstart\t$send\t$sbjct\n";
                    }
                }
            }
        }
    }
    $output->close;

    return 1;
}

1;
