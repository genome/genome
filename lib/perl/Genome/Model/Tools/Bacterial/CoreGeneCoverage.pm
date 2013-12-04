package Genome::Model::Tools::Bacterial::CoreGeneCoverage;

# bdericks: This entire module needs a major overhaul.

use strict;
use warnings;

use Genome;
use Carp;
use IPC::Run;
use File::Slurp;
use PP::LSF;

class Genome::Model::Tools::Bacterial::CoreGeneCoverage (
    is => 'Command',
    has => [
        fasta_file => { 
            is => 'Filepath',
            doc => 'Fasta containing protein sequence on which core gene coverage is checked',
        },
        percent_id => { 
            is => 'Number',
            doc => 'Acceptable percent identity',
        },
        fraction_of_length => { 
            is => 'Number',
            doc => 'Minimum fraction of the total length of a gene that needs to be present to be counted as covered',
        },
        cell_type => { 
            is => 'String',
            doc => 'Determines which core gene set/group file to use',
            valid_values => ['ARCHAEA', 'BACTERIA'],
        },
        output_file => {
            is => 'FilePath',
            doc => 'File to write coverage results to',
        },
    ],
    has_optional => [
        _passed => {
            is => 'Boolean',
            is_transient => 1,
            doc => 'Set by command module during execution, if true the core gene check passed',
        },
        minimum_percent_coverage => {
            is => 'Number',
            default => 90,
            doc => 'Minimum percentage of core genes that need to be covered for the core gene check to pass',
        },
        survey_type => {
            is => 'Text',
            default => 'geneset',
            doc => 'To create blast db, run blastp for geneset and tblastn for assembly',
            valid_values => ['geneset', 'assembly'],
        },
        # TODO Need to replace these paths with reference sequence builds of some sort
        bacterial_query_file => {
             is => 'FilePath',
             doc => 'File containing bacterial core gene sequences',
             default => '/gscmnt/ams1102/info/core_genes/bacteria/CoreGenes.faa',
        },
        bacterial_group_file => {
             is => 'FilePath',
             doc => 'File containing bacterial core group sets',
             example_values=> ['/gscmnt/ams1102/info/core_genes/bacteria/CoreGroups_66.cgf'],
        },
        archaea_query_file => {
             is => 'FilePath',
             doc => 'File containing archaeal core gene sequences',
             example_values => ['/gscmnt/ams1102/info/core_genes/archaea/Archaea_coreset_104.gi.faa']
        },
        archaea_group_file => {
            is => 'FilePath',
            doc => 'File containing archaeal core group sets',
            default => '/gscmnt/ams1102/info/core_genes/archaea/Archaea_coreset_104.gi.cgf',
        },
        run_locally => {
            is => 'Boolean',
            doc => 'Make blast run locally instead of submitting to LSF',
        },
    ],
);

sub help_brief {
    return 'Detects the presence of core genes/groups in the given set of predicted genes';
}
sub help_synopsis { return help_brief() }

sub help_detail {
    return 'Given a file with a bunch of protein sequences, determines the total coverage of core genes';
}

my %XD_FORMAT_PARAMS_FOR_SURVEY_TYPE = (
    'geneset'  => ['-p', '-I'],
    'assembly' => ['-n', '-I'],
);

my %BLASTER_FOR_SURVEY_TYPE = (
    'geneset'  => 'blastp',
    'assembly' => 'tblastn',
);

sub execute {
    my $self = shift;
    $self->status_message("Running core gene coverage command");

    my $params = $XD_FORMAT_PARAMS_FOR_SURVEY_TYPE{$self->survey_type};
    # TODO Make sure this is a tool that others can install, otherwise turn this into a genome tool
    my @xdformat = ('xdformat', @$params, $self->fasta_file);
    my ($xdf_out,$xdf_err);
    my $xdf_rv = IPC::Run::run(\@xdformat,
        '>',
        \$xdf_out,
        '2>',
        \$xdf_err, );
    unless($xdf_rv) {
        $self->error_message("failed to format fasta file ".$self->fasta_file."\n".$xdf_err."\n");
        return 0; # or should we exit(1)?
    }

    # TODO Is the bsub necessary here?
    # bsub a blastp query on that,
    my ($bsubout,$bsuberr) = ($self->fasta_file.".blastout", $self->fasta_file.".blasterr");
    my $blastresults = $self->fasta_file.".blastp_results";
    my $blast_query_file;
    if($self->cell_type eq 'BACTERIA' ) {
        $blast_query_file = $self->bacterial_query_file;
    }
    elsif($self->cell_type eq 'ARCHAEA' ) 
    {
        $blast_query_file = $self->archaea_query_file;
    }

    my $blaster = $BLASTER_FOR_SURVEY_TYPE{$self->survey_type};
    my $blastp_cmd = join(' ', $blaster, $self->fasta_file, $blast_query_file, '-o', $blastresults);

    if ( $self->run_locally ) {
        my $rv = eval { Genome::Sys->shellcmd( cmd => $blastp_cmd ) };
        $self->error_message("Failed to run blast with command: $blastp_cmd") and return
            if not $rv;
    }
    else {
    $self->status_message("bsubbing blastp command: $blastp_cmd");
    my $blastp_job = PP::LSF->create(
        command => $blastp_cmd,
        q => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        o => $bsubout,
        e => $bsuberr,
    );
    my $start_rv = $blastp_job->start;
    confess "Trouble starting LSF job for blastp ($blastp_cmd)" unless defined $start_rv and $start_rv;

    my $wait_rv = $blastp_job->wait_on;
    confess "Trouble while waiting for LSF job for blastp ($blastp_cmd) to complete!" unless defined $wait_rv and $wait_rv;
    }


    $self->status_message("Blastp done, parsing");

    # TODO Change into a module call
    # run parse_blast_results_percid_fraction_oflength on the output

    my $run_dir = File::Basename::dirname($self->output_file);
    my $cov_pid_out = $run_dir.'/Cov_30_PID_30';

    my @parse = (
        'gmt','bacterial','parse-blast-results',
        '--input', $blastresults,
        '--output', $cov_pid_out,
        '--num-hits', 1,
        '--percent', $self->percent_id,
        '--fol', $self->fraction_of_length,
        '--blast-query', $blast_query_file,
    );

    my ($parse_stdout, $parse_stderr);
    my $parse_rv = IPC::Run::run(\@parse,
        '>',
        \$parse_stdout,
        '2>',
        \$parse_stderr, );
    unless($parse_rv) {
        $self->error_message("failed to parse output ".$blastresults."\n".$parse_stderr);
        return 0;
    }

    $self->status_message("Done parsing, calculating core gene coverage percentage");

    my @core_gene_lines = read_file( $cov_pid_out );
    #@core_gene_lines = grep /====/
    my $core_groups_ref_arry = $self->get_core_groups_coverage(
        @core_gene_lines,
    );

    # the easy way, but we will have to change it eventually
    my $cmd1 = "grep \"====\" $cov_pid_out | awk '{print \$2}' | sort | uniq | wc -l";
    my $core_gene_count = `$cmd1`;
    chomp $core_gene_count;
    my $cmd2 = "grep \">\" ".$blast_query_file." | wc -l"; # counting seqs
    my $query_count = `$cmd2`;
    chomp $query_count;
    my $core_groups = scalar(@$core_groups_ref_arry);

    # this doesn't make alot of sense yet.
    # TODO Seriously? Instead of illuminating what's going on here, you're just gonna say you're confused too?
    my $core_pct = 100 * $core_gene_count / $query_count;
    my $coregene_pct = sprintf("%.02f",100 * $core_gene_count / $query_count);

    my $output_fh;
    if ($self->output_file eq 'STDOUT') {
        $output_fh = $self->output_file;
    }
    else {
        if (-e $self->output_file) {
            unlink $self->output_file;
            $self->status_message('Removed existing output file at ' . $self->output_file);
        }
        $output_fh = IO::File->new($self->output_file, 'w');
    }
    confess 'Could not get file handle for output file ' . $self->output_file unless $output_fh;

    $output_fh->print("Perc of Coregenes present in this assembly: $coregene_pct \%\n");
    $output_fh->print("Number of Core Groups present in this assembly: $core_groups\n");    
    #$output_fh->print("Core gene count: $core_gene_count\n");
    #$output_fh->print("Query count: $query_count\n");
    
    $self->status_message("Perc of Coregenes present in this assembly: $coregene_pct \%");
    $self->status_message("Number of Core Groups present in this assembly: $core_groups");
    #$self->status_message("Core gene count: $core_gene_count");
    #$self->status_message("Query count: $query_count");

    if($core_pct <= $self->minimum_percent_coverage) {
        $output_fh->print("Core gene test FAILED!\n");
        $self->status_message("Core gene test FAILED!");
        $self->_passed(0);
    }
    else {
        $output_fh->print("Core gene test PASSED\n");
        $self->status_message("Core gene test PASSED");
        $self->_passed(1);
    }

    $output_fh->close;

    # the below replicates 'cat Cov_30_PID_30 CoregeneTest_result >Cov_30_PID_30.out
    # TODO This can be replaced with Genome::Sys->cat I think
    my $covdata = read_file( $cov_pid_out );
    my $cgtest_result = read_file($self->output_file);
    write_file( $cov_pid_out.'.out',$covdata,"\n",$cgtest_result,"\n");

    # unlink temp files. these really should be absolute path
    # and should go to a writable directory....
    unlink($blastresults);
    unlink($bsubout);
    unlink($bsuberr);
    unlink($cov_pid_out);
    unlink($self->fasta_file) if $self->survey_type eq 'geneset';
    unlink($self->fasta_file.".xpd");
    unlink($self->fasta_file.".xpi");
    unlink($self->fasta_file.".xps");
    unlink($self->fasta_file.".xpt");
    # tblastn excess files
    unlink($self->fasta_file.'.xni');
    unlink($self->fasta_file.'.xns');
    unlink($self->fasta_file.'.xnd');
    unlink($self->fasta_file.'.xnt');
	#unlink($self->output_file);

	## Gzip Cov_30_PID_30.out
	## We will remove the existing Cov_30_PID_30.out.gz from earlier run before gzip'ing the Cov_30_PID_30.out file
        if ( -e $cov_pid_out.'.out.gz') {
                unlink $cov_pid_out.'.out.gz';
		$self->status_message('Removed existing Cov_30_PID_30.out.gz file');
	}
	my $gzip_rv = IPC::Run::run (
			'gzip',
                        $cov_pid_out.'.out',
			); 
	unless($gzip_rv) {
		$self->error_message("failed to gzip Cov_30_PID_30.out\n");
		return 0;
	}

    return 1;
}

sub get_core_groups_coverage {
    my ( $self, @lines ) = @_;

    my $core_groups ;
    if($self->cell_type eq 'BACTERIA') {
        $core_groups = $self->bacterial_group_file;
    }
    elsif($self->cell_type eq 'ARCHAEA') 
    {
        $core_groups = $self->archaea_group_file;
    }

    use IO::String;
    # stolen from  get_core_groups_coverage script
    # altered to 
    #read list to make an array of genes covered
    #my $gene_list = `grep "===="  Cov_30_PID_30  | awk '{print \$2}' | sort | uniq`;
    my $gene_list;
    for my $line ( @lines ) {
        if ( $line =~ /^==/ ) {
            my ( $id ) = $line =~ /==\s+(\d+)\s+/;
            $gene_list .= $id."\n";
        } 
    }

    my $gl = IO::String->new($gene_list);
    my %gene_hash;
    while (<$gl>){
        chomp;
        my $line=$_;
        $gene_hash{$line}=1;
    }

    #Output file
    my @results = ( );

    #Read ortholog data
    my $cgf=new FileHandle("$core_groups");
    while(<$cgf>){
        chomp;
        my $line=$_;
        $line =~ s/\s+$//;
        my @array=split(/\s/,$line);
        my $flag=0;
        foreach my $a (@array){
            next if ($a eq "-");
            my $gi=(split(/\(/,$a))[0];
            if ($gene_hash{$gi}){
                #print $o "$line\n";
                push(@results, $line);
                last;
            }
        }
    }

#    $o->close;
    return \@results;

}


1;
