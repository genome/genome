package Genome::Model::Tools::Annotate::ImportInterpro::ExecuteIprscan;

use strict;
use warnings;
use Genome;
use IO::File;
use Benchmark qw(:all) ;
use File::Temp;

class Genome::Model::Tools::Annotate::ImportInterpro::ExecuteIprscan{
    is => 'Genome::Model::Tools::Annotate',
    has => [
        interpro_version => { 
            is => 'Number',
            is_input => 1,
            is_optional => 1,
            default => 4.5,
            doc => 'Version of Interpro used.  This option is currently nonfunctional  The default is 4.5',
        },
        scratch_dir => { 
            is => 'Path',
            is_input => 1,
            doc => 'files for fasta generation, iprscan output, etc. are written to this directory'
        },
        benchmark => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            is_input => 1,
            doc => 'if set, run times are displayed as status messages after certain steps are completed (x, y, z, etc)',
        },
    ]
};

#TODO: Write me
sub help_synopsis {
    return <<EOS
TODO
EOS
}

#TODO: Write me
sub help_detail{ 
    return <<EOS
TODO
EOS
}

sub execute{
    my $self = shift;

    my $scratch_dir = $self->scratch_dir;
    die "Could not get tmp directory $scratch_dir" unless $scratch_dir; #TODO: Sanity check this
    my $iprscan_dir = '/gsc/scripts/pkg/bio/iprscan/iprscan-'.$self->interpro_version; #defaults to 4.5; 
    die "Could not find interpro version ".$self->interpro_version unless -d $iprscan_dir; 

    # db disconnect to avoid Oracle failures killing our long running stuff
    if (Genome::DataSource::GMSchema->has_default_handle) {
        $self->debug_message("Disconnecting GMSchema default handle.");
        Genome::DataSource::GMSchema->disconnect_default_dbh();
    }

    #converter.pl requires this environment variable to be set to the iprscan directory.  It refuses to run if this isn't set
    my $old_iprscan_home = $ENV{'IPRSCAN_HOME'}; #save this value so it can be reset at the end of the script
    $ENV{'IPRSCAN_HOME'} = $iprscan_dir; 
    #Run each .fasta through iprscan and throw the results into a tab delimited temp file
    my $pre_iprscan = Benchmark->new;
    my %iprscan;
    my %fastas = $self->get_fastas();
    my $output_text = File::Temp->new(UNLINK => 0,
                                      DIR => $scratch_dir,
                                      TEMPLATE => 'import-interpro_iprscan-output-text_XXXXX');
    for my $fasta_file (keys %fastas){
        $self->debug_message("Dealing with fasta $fasta_file");
        my $iprscan_temp = File::Temp->new(UNLINK => 0,
                                           DIR => $scratch_dir,
                                           TEMPLATE => 'import-interpro_iprscan-result_XXXXXX');
        my $iprscan_output = $iprscan_temp->filename;
        #This will run the iprscan, appending STDOUT and STDERR to the $ouput_file
        #We will eventually parse this $output_file for failure messages and run the restart jobs commands (which are normally printed to STDERR)
        Genome::Sys->shellcmd(cmd => $iprscan_dir.'/bin/iprscan -cli -i ' . $fasta_file . ' -o ' . $iprscan_output . ' -seqtype p -appl hmmpfam -appl superfamily -appl hmmsmart -appl patternscan -iprlookup -goterms -verbose -format raw >> ' . $output_text . ' 2>&1' ,) or die "iprscan failed: $!"; 
        $iprscan{$iprscan_output} = $iprscan_temp;
        print "Fetched 25000 transcripts" . "\n";
    }
    my $post_iprscan = Benchmark->new;
    my $iprscan_time = timediff($post_iprscan, $pre_iprscan);
    $self->status_message('iprscan: ' . timestr($iprscan_time, 'noc')) if $self->benchmark;

    #TODO: Restart failed jobs correctly
    my @restart_commands = $self->_find_restart_commands($output_text->filename, $iprscan_dir); 
    for my $cmd (@restart_commands){
        $self->debug_message("Restarting job with command: $cmd");
        Genome::Sys->shellcmd(cmd => $cmd) or die "iprscan restart failed for:\n $cmd\n $!"; 
    }

    #reset environment variables
    $ENV{'IPRSCAN_HOME'} = $old_iprscan_home;
}

sub get_fastas{
    my $self = shift;
    my $scratch_dir = $self->scratch_dir;
    my %fastas;

    while (my $file = glob ("$scratch_dir/import-interpro_fasta_*")) {
        $fastas{$file}++;
    }
    return %fastas;
}

#This parses the output for notifications of failed jobs and 
#collects the commands into a list.
sub _find_restart_commands{
#TODO: this clearly doesn't work.  The report file in question required for
#the restart lives someplace in $iprscan_dir/tmp/.  This also generally 
#requires some sort of manual intervention to fix the error before the restart will work...
    my ($self, $output_file, $iprscan_dir) = @_;

    my $iprscan_bin = $iprscan_dir . '/bin/';
    my $output = IO::File->new($output_file, "r");
    $output or die "Could not open output_file $output_file: $!";

    my @commands;
    for my $line (<$output>){
        if ($line =~ /^\.\//){  
            $line =~ s/^\.\//$iprscan_bin/; 
            push @commands, $line;
        }
    }
    $output->close or die "Could not close output file $output_file: $!";
    return @commands;
}

1;
