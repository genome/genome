package Genome::Model::Tools::BacterialContaminationScreen;

use strict;
use warnings;

use Genome;

use IO::File;
class Genome::Model::Tools::BacterialContaminationScreen
{
   is => 'Command',
            doc => "THIS IS A TEST",
   has => [
            type => {
                        is => 'String',
                        valid_values => ['read','contig','16S']
                    },
            e => {
                    calculate_from => [ 'type' ],
                         calculate => q( 
                                                if (not defined $type) 
                                                {
                                                    # this should not happen
                                                    return;
                                                }
                                                elsif ($type eq 'read') 
                                                {
                                                    return '1e-20';
                                                }
                                                elsif ($type eq 'contig') 
                                                {
                                                    return '1e-100';
                                                }
                                                elsif ($type eq '16S') 
                                                {
                                                    return '1e-100';
                                                }
                                                else 
                                                {
                                                    die "unkown type $type!";
                                                }
                                            ),
                    doc => '-e param for running megablast',
                 },

            minimum_size => {
                    calculate_from => [ 'type' ],
                         calculate => q( 
                                                if (not defined $type) 
                                                {
                                                    # this should not happen
                                                    return;
                                                }
                                                elsif ($type eq 'read') 
                                                {
                                                    return 120;
                                                }
                                                elsif ($type eq 'contig') 
                                                {
                                                    return 150; 
                                                }
                                                elsif ($type eq '16S') 
                                                {
                                                    return 150;
                                                }
                                                else 
                                                {
                                                    die "unkown type $type!";
                                                }
                                            ),
                    doc => '-e param for running megablast',
                 },
            input_file => {
                is => 'String',
                doc => 'list of inputs for queries (ex: /gscmnt/temp100/research/jxu/HGMI/bin/in.txt)', 
            },
            output_file => {
                is => 'String',
                doc => 'file for commands to be output to (ex: /gscmnt/temp100/research/jxu/HGMI/bin/commandLines.txt)', 
            },
            use_cluster => {
                                is => 'String',
                                default => 'N',
                                doc => 'in legacy code, use_cluster was determined by running the cmd `/gscmnt/temp100/research/jxu/bin/get_fasta_stats orfs.faa |wc|awk \'\{print \$1\}\'` and checking if the value returned is > 1000.  However, this line was commented out, so defaulted to N'
                            }
          ],
};    
sub create 
{ 
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;
    
    unless ( defined $self->type) {
        $self->error_message("A type is expected");
        return;
    }
    unless (defined $self->input_file)
    {
        $self->error_message("An input file is expected");
    }
    unless (defined $self->output_file)
    {
        $self->error_message("An output file is expected");
    }
    return $self;
}

sub execute()
{
   my $self = shift;
   my $fh = new IO::File;
   my $output_file = $self->output_file;
   if ($fh->open("> $output_file"))
   {
        my $cmd;
        my$query = 'query.fna';
        $cmd = '/gsc/scripts/bin/auto-GC-check.pl &';
        $fh->print("$cmd\n");
        $cmd = '/gsc/scripts/bin/auto-GCSeqsize-check.pl &';
        $fh->print("$cmd\n");
        $cmd = 'getorf -sequence ' . $query . ' -out orfs.faa.orig -minsize '.$self->minimum_size;
        $fh->print("$cmd\n");
        $cmd = '/gsc/scripts/bin/removeOverlappingORFsFromFasta.pl orfs.faa.orig > orfs.faa'; 
        $fh->print("$cmd\n");
        $fh->close; 
        $self->_process_input;
   } 
   else
   {
        die("couldn't open output file " . $self->output_file);
   }
}

sub _process_input
{ 
    my ($self) = @_; 
    my $fh = new IO::File;
    my $in_file = $self->input_file;
    if ($fh->open("< $in_file"))
    {
        while(<$fh>)
        {
            my (@t);
            my ($tool, $query, $db, $path); 
            chomp($_);
            next if $_ =~ /^\#/;
            @t = split(/[\s\t]+/, $_);
            $tool = $t[0];
            $query = $t[1];
            $db = $t[2];
            $path = $t[3];
            $self->_run_query($tool, $query, $db, $path);
        }
        $fh->close;
    }
    else
    {
       die("No input file $in_file"); 
    }
}
sub _run_query
{
  my ($self, $tool, $query, $db, $path, $cmd) = @_;
  my ($outFile, $outOutFile, $annotFile, $headerFile, $currPath, $qscript, $origQuery, $firstKnown);
  my $output_file = $self->output_file;
  my $fh = new IO::File;

  if ($fh->open(">> $output_file"))
  {    
      if ($tool =~ /mega/) 
      {
        $outFile = $query.'_'.$tool.'_'.$db.'_'.$self->e.'.out';
        $cmd ='megablast -d '.$path.'/'.$db.' -i '.$query.' -p 0.9 -e '.$self->e.' -D 3 >'.$outFile; 
        $fh->print("$cmd\n");
        $outOutFile = $outFile.'.bestHit';
        $cmd = 'cat '.$outFile.' |sort -k 1 | /gsc/scripts/bin/sortHitsInBlock_1.pl | /gsc/scripts/bin/getBestScoreOneHit.pl > '.$outOutFile;
        $fh->print("$cmd\n");
        $annotFile = $outOutFile.'.annot';
        if ($db eq 'nt') 
        { 
            $headerFile = '/gsc/scripts/share/bacterial_contamination_screening/nt.header.index';
        } 
        elsif ($db eq 'SSU') 
        { 
            $headerFile = '/gsc/scripts/share/bacterial_contamination_screening/SSU.header.index';
        } 
        else 
        {
          $fh->print("Error in $headerFile\n:");
       }
       $cmd = '/gsc/scripts/bin/crossTables.pl '.$headerFile.' '.$outOutFile.' | sort -k 12 -n -r > '.$annotFile;
       $fh->print ("$cmd\n");
     } 
     elsif (($tool =~ /blast/) && ($query ne 'orfs.faa') ) 
     {
       $outFile = $query.'_'.$tool.'_'.$db.'_'.$self->e.'.out';
       $cmd = "blastall -p ".$tool." -d ".$path.'/'.$db." -i ".$query." -e ".$self->e." -a 2 > ".$outFile;
       $fh->print("$cmd\n");
       $firstKnown = $outFile.'.firstKnownHit';
        $cmd = '/gsc/scripts/bin/blast2firstKnownHit.pl '.$outFile.' |grep -v hypothe > '.$firstKnown; 
        $fh->print("$cmd\n");
        $cmd = '/gsc/scripts/bin/sortColumn.pl '.$firstKnown.' 5 > '.$firstKnown.'.sorted';
        $fh->print("$cmd\n");

      } 
      elsif ($query eq 'orfs.faa')  
      {
        if ($self->use_cluster eq 'N') 
        {
          $outFile = $query.'_'.$tool.'_'.$db.'_'.$self->e.'.out';
          $cmd = "blastall -p ".$tool." -d ".$path.'/'.$db." -i ".$query." -e ".$self->e." -a 2 > ".$outFile;  
          $fh->print("$cmd\n");
        } 
        else 
        {
          open (FOF, "< fasta.fof") or die "couldn't open fasta.fof: $!\n";
          while (<FOF>) 
          {
            chomp($_);
            $query = $_;
            $qscript = 'qscript.'.$query; 
            $currPath = `pwd`;
            chomp($currPath);
            $origQuery = $query;
            $query = $currPath.'/'.$query;
            $outFile = $currPath.'/'.$origQuery.'_'.$tool.'_'.$db.'_'.$self->e.'.out';
            open (OUT, "> $qscript");
            $cmd = "cd ".$currPath;
            print (OUT "$cmd\n");
            $cmd = "blastall -p ".$tool." -d ".$path.'/'.$db." -i ".$query." -e ".$self->e." -a 2 > ".$outFile;
            print (OUT "$cmd\n"); 
          }
          close(FOF);
          $cmd = '';
          $fh->print("echo Please run cluster version of orfs.faa.*\n");
        }
      } 
      else 
      {
        print "$tool unknown\n"; exit; 
      }
      $fh->close;
  }
  else
  {
    die ("unable to open output file $output_file");
  }
}#run_tool

sub help_detail {                           # This is what the user will see with --help <---
    return <<EOS
gmt breakpoint-pal --breakpoint-id chr11:36287905-36288124


EOS
}
1;
