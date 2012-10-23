package Genome::Model::Tools::Hgmi::SequenceName;

use strict;
use warnings;

use IO::File;
use File::Slurp;
use File::Path qw( make_path );
use Cwd;
require "pwd.pl";
use Bio::SeqIO;
use Bio::Seq;
use Carp;


UR::Object::Type->define(
class_name => __PACKAGE__,
is => 'Command',
has => [
        'fasta'            => { is => 'String',
				doc => "fasta file" },
        'analysis_version' => { is => 'String',
                                doc => "Analysis version" },
        'locus_tag'        => { is => 'String',
				doc => "Locus ID string" },
        'acedb_version'    => { is => 'String',
				doc => "Ace DB version" },
        'new_output'       => { is => 'String',
				doc => "new fasta output file",
				is_optional => 1,
			      },
        'path' => { is => 'String',
                    doc => "base path for output files",
                   }, 
        'project_type' => { is => 'String',
                            doc => "project type for pipeline [HGMI|CORE]...",
                          },
       ]


			);

sub help_brief
{
    "tool for renaming sequences (HGMI Projects only!)";
}

sub help_synopsis
{
    my $self = shift;
    return <<"EOS"
This command is intended for use with HGMI projects only!
EOS
}

sub help_detail
{
    my $self = shift;
    return <<"EOS"
This command is intended for use with HGMI projects only!
EOS
}

sub execute
{
    my $self = shift;

    my @line = split(/\./,$self->fasta);

    my $revised_fasta = $line[0] = 'contigs';
    my $fasta = $line[1] = 'fasta';
    my $short_ver_num = $self->version_lookup($self->analysis_version);
    my $new = "newname";
    my $new_output_file = join(".",$self->locus_tag,
                                   $short_ver_num,
                                   $line[0],
                                   $new,
                                   $line[1]);

    $self->new_output($new_output_file);

    my $instream = new Bio::SeqIO(-file => $self->fasta, -format => 'fasta');
    my $outstream = new Bio::SeqIO(-file => ">$new_output_file",
                                   -format => 'fasta');

    my @seq_gcfile;

    while( my $seqobj = $instream->next_seq)
    {
        my $seq_id = $seqobj->primary_id();
        my $new_seq_id = join("_",$self->locus_tag,$seq_id);
        push(@seq_gcfile, $new_seq_id);
        my $seq = $seqobj->seq;
        if ($seq =~ m/x/ix)
        {carp "\n\nWARNING this contig contains X's : $new_seq_id\n\n";}
        $seq = uc($seq);
        my $newseq = new Bio::Seq(-seq => $seq,
                                  -id => $new_seq_id
                                  );
        $outstream->write_seq($newseq);

    }

    # this dependence on being in the current directory 
    # needs to be fixed.
    my $cwd = getcwd();
    my @cwd = split(/\//x,$cwd);
    if($#cwd < 9)
    {
        croak "the current working directory seems short,\nare you in the right place?\n$cwd\n$#cwd";
    }
    # this is getting hard coded here.  this is not good.
    my $hgmi_acedb_patha;

    my $acedb_version = $self->acedb_version;
    #$acedb_version =~ s/V(\d)/Version_$1\.0/;
    #$hgmi_acedb_patha = $self->path. "/".$self->project_type ."/Acedb/". $acedb_version .
    $hgmi_acedb_patha = $self->path. "/Acedb/". $acedb_version .
                        "/ace_files/". $self->locus_tag;
    $self->status_message("acedb path will be: ".$hgmi_acedb_patha);
#    $hgmi_acedb_patha = "/gscmnt/278/analysis/HGMI/Acedb/". $acedb_version . 
#                        "/ace_files/". $self->locus_tag;

    unless (-e $hgmi_acedb_patha)
    {
        make_path($hgmi_acedb_patha, { verbose => 1, mode => 0775 })
            or croak "Can't make path $hgmi_acedb_patha: $!";
#        mkdir qq{$hgmi_acedb_patha} 
#	or croak "Can't make $hgmi_acedb_patha: $!\n";
    }

    my $newHGMIpath = $hgmi_acedb_patha."/".$self->analysis_version;

    unless (-e $newHGMIpath)
    {
        make_path($newHGMIpath, { verbose => 1, mode => 0775 })
            or croak "Can't make path $newHGMIpath: $!";
#        mkdir qq{$newHGMIpath}
#	or croak "Can't make $newHGMIpath: $!\n";
    }

    my $hgmi_acedb_path = $newHGMIpath;

    my ($Intergenic,$BAPseq,$Ensemblseq,$Rfamseq);
    # FIXME: the @cwd[..] is bad, need to get rid of it.
    # need to change directory construction here
    # @cwd[0..8] for CORE projects
    if($self->project_type eq 'HGMI')
    {
        # $self->path / org_dirname / assembly name /assembly version / Sequence / Unmasked
        my $dirbase = $cwd;
        $dirbase =~ s/Sequence\/Unmasked//;
        my $testval = join("\/",@cwd[0..7]);
        $self->status_message("is this $testval better\nthan $dirbase ?");
        #$Intergenic = join("\/",@cwd[0..7],'Gene_merging',$self->analysis_version,'Hybrid','intergenic');
        $Intergenic = join("\/",$dirbase,'Gene_merging',$self->analysis_version,'Hybrid','intergenic');
        #$BAPseq = join("\/", @cwd[0..7],'BAP',$self->analysis_version,'Sequence');
        $BAPseq = join("\/", $dirbase ,'BAP',$self->analysis_version,'Sequence');
        #$Ensemblseq = join("\/", @cwd[0..7],'Ensembl_pipeline',$self->analysis_version,'Sequence');
        $Ensemblseq = join("\/", $dirbase,'Ensembl_pipeline',$self->analysis_version,'Sequence');
        #$Rfamseq = join("\/", @cwd[0..7],'Rfam',$self->analysis_version);
        $Rfamseq = join("\/", $dirbase,'Rfam',$self->analysis_version);
    }
    elsif($self->project_type eq 'CORE')
    {
        $Intergenic = join("\/",@cwd[0..8],'Gene_merging',$self->analysis_version,'Hybrid','intergenic');
        $BAPseq = join("\/", @cwd[0..8],'BAP',$self->analysis_version,'Sequence');
        $Ensemblseq = join("\/", @cwd[0..8],'Ensembl_pipeline',$self->analysis_version,'Sequence');
        $Rfamseq = join("\/", @cwd[0..8],'Rfam',$self->analysis_version);

    }
    elsif($self->project_type eq 'VIRAL')
    {
        $Intergenic = join("\/",@cwd[0..8],'Gene_merging',$self->analysis_version,'Hybrid','intergenic');
        $BAPseq = join("\/", @cwd[0..8],'BAP',$self->analysis_version,'Sequence');
        $Ensemblseq = join("\/", @cwd[0..8],'Ensembl_pipeline',$self->analysis_version,'Sequence');
        $Rfamseq = join("\/", @cwd[0..8],'Rfam',$self->analysis_version);

    }
    else
    {
        $self->error_message("unknown project type: ".$self->project_type);
        croak;
    }

    # the presence of the symlink target should be
    # tested before creating the symlink
    # symlink OLDFILE,NEWFILE
    my $hgmiOUTPUTfilelink = qq{$hgmi_acedb_path/$new_output_file};

    unless ( -e $hgmiOUTPUTfilelink) {
      symlink "$cwd/$new_output_file","$hgmi_acedb_path/$new_output_file" ||
	die "Can't make symlink path: $!\nsource: $cwd/$new_output_file\ndest: $hgmi_acedb_path/$new_output_file\n ";
    }

    my $intergenicNEWfilelink = qq{$Intergenic/$new_output_file};

    unless ( -e $intergenicNEWfilelink) {
      symlink "$cwd/$new_output_file","$Intergenic/$new_output_file" or 
	die "Can't make symlink path Intergenic: $!\nsource: $cwd/$new_output_file\ndest: $Intergenic/$new_output_file";
    }

    my $bapNEWfilelink = qq{$BAPseq/$new_output_file};
    
    unless ( -e $bapNEWfilelink ) {
      symlink "$cwd/$new_output_file","$BAPseq/$new_output_file" or 
	die "Can't make symlink path BAP: $!\nsource: $cwd/$new_output_file\ndest: $BAPseq/$new_output_file";
    }
   
    my $ensemblNEWfilelink = qq{$Ensemblseq/$new_output_file};

    unless (-e $ensemblNEWfilelink ) {
      symlink "$cwd/$new_output_file","$Ensemblseq/$new_output_file" or 
	die "Can't make symlink path Ensembl: $!\nsource: $cwd/$new_output_file\ndest: $Ensemblseq/$new_output_file";
    }
    my $rfamNEWfilelink = qq{$Rfamseq/$new_output_file};
    
    unless ( -e $rfamNEWfilelink) {
      symlink "$cwd/$new_output_file","$Rfamseq/$new_output_file" or 
	die "Can't make symlink path Rfam: $!\nsource: $cwd/$new_output_file\ndest: $Rfamseq/$new_output_file\n" ;
    }

    chdir($Rfamseq);
    my $cwd2 = getcwd();
    if(exists($ENV{HGMI_DEBUG}))
    {
        print "\n$cwd2\n\n";
    }

    chdir($hgmi_acedb_path);
    my $gc_file = "genomic_canonical.ace";
    my $new_gc_outfile = join(".",$self->locus_tag,$gc_file);
    my @gc_outlines = ( );
    foreach my $gc_name (@seq_gcfile)
    {
        push(@gc_outlines,"Sequence $gc_name\n");
        push(@gc_outlines,"Genomic_canonical\n\n");
    }
    write_file($new_gc_outfile,@gc_outlines) or 
        croak "can't write to $new_gc_outfile: $!";;


    return 1;
}


sub version_lookup
{
    my $self = shift;
    my $v = shift;
    my $lookup = undef;
    # FIXME : have this be a regex substitution, otherwise we will
    # have to constantly update this.
    my %version_lookup = ( 
		       'Version_1.0' => 'v1', 'Version_2.0' => 'v2',
		       'Version_3.0' => 'v3', 'Version_4.0' => 'v4',
		       'Version_5.0' => 'v5', 'Version_6.0' => 'v6',
		       'Version_7.0' => 'v7', 'Version_8.0' => 'v8',
		       );

    if(exists($version_lookup{$v}))
    {
        $lookup = $version_lookup{$v};
    }

    return $lookup;
}

# for the future
#sub directory_create
#{
#    my $self = shift;
#
#    return 1;
#}

1;
