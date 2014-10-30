package Genome::Model::Tools::Annotate::ImportInterpro::Run;

use strict;
use warnings;
use Genome;
use Benchmark qw(:all) ;

my $low  = 1000;
my $high = 20000;
UR::Context->object_cache_size_lowwater($low);
UR::Context->object_cache_size_highwater($high);

class Genome::Model::Tools::Annotate::ImportInterpro::Run{
    is => 'Genome::Model::Tools::Annotate',
    has => [
        reference_transcripts => {
            is => 'String',
            is_input => 1, 
            is_optional => 0,
            doc => 'provide name/version number of the reference transcripts set you would like to use ("NCBI-human.combined-annotation/0").',
        },
    ],
    has_optional => [
        interpro_version => { 
            is => 'Number',
            is_input => 1,
            is_optional => 1,
            default => 4.5,
            doc => 'Version of Interpro used.  This option is currently nonfunctional  The default is 4.5',
        },
        chunk_size => {
            is => 'Number',
            is_input => 1,
            is_optional => 1,
            default => 25000,
            doc => 'Number of sequences submitted to interpro at a time.  Defaults to 25000',
        },
        commit_size => {
            is => 'Number',
            is_input => 1,
            is_optional => 1,
            default => 100,
            doc => 'Number of Interpro results saved at a time.  Defaults to 100',
        },
        benchmark => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            is_input => 1,
            doc => 'if set, run times are displayed as status messages after certain steps are completed (x, y, z, etc)',
        },
        log_file => { 
            is => 'Path',
            is_optional => 1,
            default => '/dev/null',
            is_input => 1,
            doc => 'if set, STDOUT and STDERR are directed to this file path for logging purposes.  Defaults to /dev/null'
        },
        scratch_dir => { 
            is => 'Path',
            is_input => 1,
            doc => 'files for fasta generation, iprscan output, etc. are written to this directory'
        },
    ],
};

sub help_synopsis {
    return <<EOS
gmt annotate import-interpro --reference-transcripts NCBI-human.combined-annotation/54_36p
EOS
}

sub help_detail{
    return <<EOS
This runs the interpro import.  It takes the name and version number of the set of reference transcripts.

The current version uses IPRscan version 4.5.  Work in in progress to add support for other versions of IPRscan.

This tool runs all of the transcripts in the set through IPRscan, parses the results, and creates a Genome::InterproResults object for each result
EOS
}

sub execute {
    my $self = shift;

    print "Starting ImportInterpro Run" ."\n";
   
    my $total_start = Benchmark->new;
   
    my $log_file = $self->log_file; #TODO: sanity check this
    open (OLDOUT, ">&STDOUT");
    open (OLDERR, ">&STDERR");

    close(STDOUT);
    close(STDERR);

    open(STDOUT, "> $log_file") or die "Can't redirect STDOUT: $!";
    open(STDERR, "> $log_file") or die "Can't redirect STDERR: $!";
   
    my ($model_name, $build_version) = split("/", $self->reference_transcripts);
    my $model = Genome::Model->get(name => $model_name);
    die "Could not get model $model_name" unless $model;
    my $build = $model->build_by_version($build_version);
    die "Could not get imported annotation build version $build_version" unless $build;
    my $chunk_size = $self->chunk_size;
    die "Could not get chunk-size $chunk_size" unless $chunk_size; 
    die "chunk-size of $chunk_size is invalid.  Must be between 1 and 50000" if($chunk_size > 50000 or $chunk_size < 1);
    my $commit_size = $self->commit_size;
    die "Could not get commit-size $commit_size" unless $commit_size; 
    die "commit-size of $commit_size is invalid.  Must be greater than 1" if($commit_size < 1);
    my $iprscan_dir = '/gsc/scripts/pkg/bio/iprscan/iprscan-'.$self->interpro_version; #defaults to 4.5; 
    die "Could not find interpro version ".$self->interpro_version unless -d $iprscan_dir; 

    my $scratch_dir = $self->scratch_dir;
    die "Could not get tmp directory $scratch_dir" unless $scratch_dir; #TODO: Sanity check this
    print "Starting GenerateTranscript Fastas" . "\n";
    my $fasta_success = Genome::Model::Tools::Annotate::ImportInterpro::GenerateTranscriptFastas->execute(
        build => $build,
        chunk_size => $chunk_size,
        benchmark => $self->benchmark,
        scratch_dir => $scratch_dir,
    );
    die "Could not generate .fasta files: $!" unless $fasta_success->result;
    print "Finished GenerateTranscriptFastas" . "\n";
    my $interpro_success = Genome::Model::Tools::Annotate::ImportInterpro::ExecuteIprscan->execute(
        benchmark => $self->benchmark,
        scratch_dir => $scratch_dir,
    );
    die "Could not complete Interpro scan: $!" unless $interpro_success->result;
    print "Finished Iprscan" . "\n";

    #Here, we are manually unloading UR objects that are sticking around in
    #the cache after they should be deleted.  This prevents fatal, out of
    #memory errors. Its a bit of a hack, but we were pressed for time.
    for my $class (qw/Genome::Transcript Genome::Protein Genome::DataSource::Proteins Genome::DataSource::Transcripts UR::DataSource::File/) {
        my @o = $class->is_loaded;
        for my $o (@o) {
            $o->unload;
        }
    }   

    print "Finished UR Object Unload" . "\n";
    my $results_success = Genome::Model::Tools::Annotate::ImportInterpro::GenerateInterproResults->execute(
        build => $build,
        benchmark => $self->benchmark,
        scratch_dir => $scratch_dir,
        commit_size => $commit_size,
        reference_transcripts => $self->reference_transcripts,
    );
    die "Could not generate Interpro results: $!" unless $results_success->result;
    print "Finished GenerateInterproResuts" . "\n";
    
    close(STDOUT);
    close(STDERR);

    open (STDOUT, ">&OLDOUT");
    open (STDERR, ">&OLDERR");

    close (OLDOUT);
    close (OLDERR);
    
    my $total_finish = Benchmark->new;
    my $total_time = timediff($total_finish, $total_start);
    $self->status_message('Total: ' . timestr($total_time, 'noc')) if $self->benchmark;

    return 1;
}
1;

=pod

=head1 Name

Genome::Model::Tools::Annotate::ImportInterpro::Run

=head1 Synopsis

Gets every transcript for a given build, runs them through Interpro, and creates Genome::InterproResult objects from the results

=head1 Usage

 in the shell:

     gmt annotate import-interpro run --reference-transcripts NCBI-human.combined-annotation/54_36p --scratch-dir my_dir

 in Perl:

     Genome::Model::Tools::Annotate::ImportInterpro::Run->execute(
         reference_transcripts => 'NCBI-human.combined-annotation/54_36p',
         interpro_version => '4.1', #default 4.5
         chunk_size => 40000, #default 25000
         log_file => mylog.txt, #default /dev/null
         scratch_dir =>  'myDir',
     );

=head1 Methods

=over

=item variant_file

A string containing the name and version number of the reference transcripts to use as input.  The format is:
name/version_number

=item 

=back

=head1 See Also

B<Genome::InterproResult>, 

=head1 Disclaimer

Copyright (C) 2010 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

 B<Jim Weible> I<jweible@genome.wustl.edu>

=cut


#$Id: 
