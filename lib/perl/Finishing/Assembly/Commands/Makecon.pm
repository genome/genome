package Finishing::Assembly::Commands::Makecon;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;

my ($whole, $project, $dest, $out_file, $test, $proj_dir, $organism);

my %name :name(name:r) :isa(string) :desc('Project name');
my %dest :name(destination:o) :isa(dir_rw) :desc('Destination directory for consensus');
my %whole :name(whole:o) :isa(string) :desc('Whole');

sub START
{
    my $self = shift;

    $project = $self->name;
    $dest = $self->destination || Cwd::cwd();
    $out_file = "$dest/$project.con";
    $whole = $self->whole;

    return 1;
}

sub execute
{
    my $self = shift;

    use Cwd;
    use GetSeq;
    use GSCApp;
    use Findid::Utility;
    use ProjectWorkBench::Model::Ace::Dir;

    #prepare();

    App->init;

    my $proj = GSC::Project->get(name => $project);
    print "$project doesn't have a GSC::Project object\n" unless $proj;

    $organism = Findid::Utility->convert_GSC_to_ana($proj->get_species_name) if $proj;
    $proj_dir = Findid::Utility->get_proj_dir($project);
    my $edit_dir = $proj_dir.'/edit_dir' if $proj_dir;
    my $pstatus  = $proj->project_status if $proj;

    my $ace = ProjectWorkBench::Model::Ace::Dir->new(dir => $edit_dir) 
    if $edit_dir and -d $edit_dir;
    my $recent_ace = $ace->recent_ace if $ace;

    print "\n    Project   :  $project\n";
    print "    Status    :  $pstatus\n"      if $pstatus;
    print "    Organism  :  $organism\n"     if $organism;
    print "    Newest Ace:  $recent_ace\n\n" if $recent_ace;

    #get_test();
    $self->get_test($pstatus);

    my %params = (
        type   => 'project', 
        name   => $project,
        prefix => 1
    );

    if ($test) {
        print "Get seq from the most recent ace file of $project\n\n";   
    }
    else { 
        my $org_dir = Findid::Utility->get_org_finished_seq_dir($organism);

        unless ($org_dir) {
            print "No analysis finished seq dir available for $organism\n";
            exit 1;
        }

        print "Get finished seq of $project from analysis storage\n\n";

        $params{finished} = 1;
        $params{whole} = 1 if $whole;
    }

    GetSeq->new(%params)->get_file($out_file);  

    exit 1 unless -s $out_file;

    print "Done.\n\n";

    return 1;
}

sub prepare {
    my $usage = 'makecon [-whole] project_name destination';
    my $error = "Wrong commandline\nUsage is:\n$usage\n";

    die $error unless @ARGV > 0 and @ARGV < 4;

    if (@ARGV == 3) {
        die $error unless $ARGV[0] eq '-whole';
        ($whole, $project, $dest) = @ARGV;
    }
    elsif (@ARGV == 2) {
        if ($ARGV[0] eq '-whole') {
            ($whole, $project) = @ARGV;
        }
        else {
            ($project, $dest) = @ARGV;
        }
    }
    else {
        ($project) = @ARGV;
    }

    $dest = cwd unless $dest;

    $dest    =~ s/^\s*(\S+?)\s*$/$1/;
    $project =~ s/^\s*(\S+?)\s*$/$1/;

    $project  = uc($project);
    $out_file = "$dest/$project.con";

    @ARGV = ();

    return;
}

sub get_test {

    my ($self, $pstatus) = @_;

    if ($pstatus) {

        if ($pstatus eq 'redundant') {
            if (-d $proj_dir) {
                my $statusfile = "$project.submit.status";

                if (-e "$proj_dir/$statusfile") {
                    my $check = `cat $proj_dir/$statusfile 2>/dev/null`;
                    $test++ unless $check =~ /\sSUBMIT\s/;
                }
                else {
                    print "No $statusfile\nAssume $project never presubmitted\n";
                    $test++;
                }
            }
            else {
                print "Can't locate $project in file system\n";
                print "Assume it is finished and offline. Try get seq from analysis\n";
            }
        }
        elsif ($pstatus ne 'submitted' && $pstatus ne 'resubmitted' && $pstatus ne 'submitted_redundant'){
            $test++;
        }
    }
    else {
        print "$project doesn't have project status.\n";
        print "Still try to get the seq from seqmgr in case it exists\n\n";

        my $match = "$proj_dir/edit_dir/*ace*";

        if (glob($match)) {
            $test++;
        }
        else {
            print "$project doesn't have valid acefiles.\n";
            print "Check spelling, make seqmgr link or create project in DB\n";
            exit 1;
        }
    }

    $test = 1 if $organism and $organism eq 'maize';

    return;
}

1;

=pod

=head1 Name

Finishing::Assembly::Commands::Makecon

=head1 Synopsis

=head1 Usage

=head1 Methods

=head2 

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head1 See Also

=head1 Disclaimer

Copyright (C) 2005 - 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$

