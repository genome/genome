package Genome::Model::Tools::Predictor::Ber;

use strict;
use warnings;
use Genome;
use Command;
use Getopt::Long;

use Carp;
use File::Path qw(mkpath);
use File::Spec;
use Cwd;
use YAML qw( LoadFile );
use English;

class Genome::Model::Tools::Predictor::Ber {
    is => 'Genome::Model::Tools::Predictor::Base',
    doc => 'execute BER gene predictor',
    
};

sub help_brief
{
    "Runs the entire BER product naming pipeline";
}

sub help_synopsis
{
    return <<"EOS"
Runs the entire BER product naming pipeline. This tool pulls data from biosql, runs blastp, hmmpfam, btab, hmmtab and anno-sqlite.bash.
EOS
}

sub help_detail
{
    return <<"EOS"
Runs the entire BER product naming pipeline. This tool pulls data from biosql, runs blastp, hmmpfam, btab, hmmtab and anno-sqlite.bash.  In addition, this script will also create the BER_product_naming.ace file, writes the tace script, parses all ace files into acedb, gathers QC stats from the phase5 file and queries acedb to make sure gene counts match, then writes the rt ticket blurb and mails the user when finished.
EOS
}

sub requires_chunking {
    return 0;
}

sub run_predictor {
    my $self = shift;

    $self->status_message("Starting BER");
    $self->status_message("Creating kegg output directory at " . $self->output_directory);

    if (-d $self->output_directory) {
        $self->warning_message("Existing output directory found at " . $self->output_directory . ", removing!");
        my $rv = system("rm -rf " . $self->output_directory);
        confess "Could not remove " . $self->output_directory unless defined $rv and $rv == 0;
    }

    make_path($self->output_directory);
    chmod(0775, $self->output_directory);



#Requirements for re-run: Not in any modules
#-If we rerun BER we usually delete the sqlite db that exists
#-If we run BER more than one time on same genome and have diff data sets then can’t use same config files. Can’t change the locus tag becuase gene name is the same for those that still exist.
#-Can overwrite config files with new name list.


    $self->_create_config_files();
    $self->_setup_dirs();
    $self->_prep_input_files();
    $self->_run_anno_sqlite();
    $self->_finish();
    


    return 1;
}

sub _get_locus_from_fasta_header {
    my $self = shift;
    
    $self->input_fasta_file;
    
    #get headers
	my $output_fh = Genome::Sys->open_file_for_reading($self->raw_output_path) or die "Could not get file handle for " . $self->raw_output_path;
    
    my $fasta_header = $output_fh->getline;
   	chomp $fasta_header;
    
    #check header format 
    my ($locus_id) = $fasta_header =~ /^\>\w*?\-?(\w+)_Contig/;

    #get locus id and ensure only one exists
    return $locus_id;
}


sub _create_config_files {

    my $self = shift;

    my $locus_id = $self->_get_locus_from_fasta_header;
    
    my $locus=-1;
    
    my $OUT1 = Genome::Sys->open_file_for_writing($locus_id.'_asm_feature') or die "failed to open _asm_feature";
	my $OUT2 = Genome::Sys->open_file_for_writing($locus_id.'_asmbl_data') or die "failed to open _asmbl_data";
	my $OUT3 = Genome::Sys->open_file_for_writing($locus_id.'_ident2') or die "failed to open _ident2";
	my $OUT4 = Genome::Sys->open_file_for_writing($locus_id.'_stan') or die "failed to open _stan";

	
	$OUT1->print("asmbl_id\tend3\tend5\tfeat_name\tfeat_type\n");
	$OUT2->print("id\tname\ttype\n");
	$OUT3->print("complete\tfeat_name\tlocus\n");
	$OUT4->print("asmbl_data_id\tasmbl_id\tiscurrent\n");
	
	
	my $count=0;
	my $asmblid=0;

	my $file = 'FULLPATH Final_BER_Naming.05-15-13.fof';
	my $fh = Genome::Sys->open_file_for_reading($file) or die "failed to open $file";
	
	while (my $line = $fh->getline)  {
	  my @fea_part;
	  my @fea_name;
	  chomp $line;
	  my @arr=split('\t',$line);
	  if($arr[1]!~/^$locus_id/)
	  {
	    @fea_part=split('-',$arr[1]);
	    if($fea_part[1]=~/$locus_id/)
	    {
	      @fea_name=split("_",$fea_part[1]);
	    }
	    elsif($fea_part[2]=~/$locus_id/)
	    {
	      @fea_name=split("_",$fea_part[2]);
	    }
	    $fea_name[1]=~s/\D+//i;
	    if($fea_name[1] != $asmblid)
	    {
	      $asmblid=$fea_name[1];
	      $count++;
	      
	      $OUT2->print($count,"\tContig\tcontig\n");
	      $OUT4->print($count,"\t",$count,"\t","1\n");
	    }
	  }
	  else
	  {
	    @fea_part=split('\.',$arr[1]);
	    @fea_name=split("_",$fea_part[0]);
	    $fea_name[1]=~s/\D+//i;
	    if($fea_name[1] !~ $asmblid)
	    {
	      $asmblid=$fea_name[1];
	      $count++;
	      
	      $OUT2->print($count,"\tContig\tcontig\n");
	      $OUT4->print($count,"\t",$count,"\t","1\n");
	    }
	  }  
	
	#  if($#fea_part==3){$locus=$fea_part[3];}
	#  else{$locus=$fea_part[2];}
	  $locus++;
	  
	  $OUT1->print($count,"\t",$arr[3],"\t",$arr[2],"\t",$arr[1],"\tORF\n");
	  $OUT3->print(" \t",$arr[1],"\t",$fea_name[0]);
	  
	  $OUT3->printf("%05d\n",$locus);
	}	
	
    return 1;

#    replace with real code to produce config files in a config dir in our output directory
#    {
#        Genome/Model/Tools/Ber/AmgapPrepareBER.pm
#            
#            Manually running details:
#            > perl ~/git/xu_script/BER_config_from_txtfile.pl -locus TELCIRDFT Final_BER_Naming.05-15-13.fof
#                
#                This command creates 4 files
#                 TELCIRDFT_asm_feature
#                 TELCIRDFT_asmbl_data
#                 TELCIRDFT_ident2
#                 TELCIRDFT_stan
#                                    
#                                    run the command todos on the 4 files
#                                        > todos TELCIRDFT_*
#    }
}

sub _setup_dirs {
	
	my $self = shift;

    my $locus_id = $self->_get_locus_from_fasta_header;
    
    
    Genome::Sys->create_directory($locus_id);
     
    map{Genome::Sys->create_directory($locus_id."/$_")}qw(fasta hmm ber);
     
    Genome::Sys->create_symlink("/gscmnt/gc6125/info/annotation/worm_analysis/T_circumcincta/T_circumcincta.14.0.ec.cg.pg/Version_1.0/PAP/Version_1.0/BER/*pep", $locus_id);
    
    chdir($locus_id."/fasta");
    
	Genome::Sys->shellcmd(cmd => "dbshatter ../*pep");
	
	chdir("../../../db/CSV");
	
	#Revisit - Multiple file copying?
	Genome::Sys->copy_file("/gscmnt/gc6125/info/annotation/worm_analysis/T_circumcincta/T_circumcincta.14.0.ec.cg.pg/Version_1.0/PAP/Version_1.0/BER/Version_1.0/($locus_id)_*", .);

	chdir("/gscmnt/gc9002/info/annotation/BER/autoannotate_v2.5/data/genomes/($locus_id)/fasta");
	
	
	#Genome/Model/Tools/Ber/AmgapBerProtName.pm
	#
	#Manually running details:
	#location is (this is taken from the MGAP pipeline config file that’s made manually )
	#/gscmnt/gc9002/info/annotation/BER/autoannotate_v2.5/data/genomes
	#
	#mkdir TELCIRDFT
	#cd TELCIRDFT
	#mkdir fasta hmm ber
	#ln -s /gscmnt/gc6125/info/annotation/worm_analysis/T_circumcincta/T_circumcincta.14.0.ec.cg.pg/Version_1.0/PAP/Version_1.0/BER/*pep .
	#
	#cd fasta
	#dbshatter ../*pep
	#
	#cd ../../../db/CSV
	#
	#cp /gscmnt/gc6125/info/annotation/worm_analysis/T_circumcincta/T_circumcincta.14.0.ec.cg.pg/Version_1.0/PAP/Version_1.0/BER/Version_1.0/TELCIRDFT_* .
	#
	#cd /gscmnt/gc9002/info/annotation/BER/autoannotate_v2.5/data/genomes/TELCIRDFT/fasta

}

sub _prep_input_files {

#3. Prep input files for BER
#
#	my $cmd = Genome::Model::Tools::Ber::N|Berxxxxx->create();
#  	my $rv = $cmd->execute;
#
#   i) Running blastp & hmmpfam
#       Genome/Model/Tools/Ber/BerRunBlastphmmpfam.pm
#
#       Manully running details:
#       > for file in `cat dbshatter.fof`; do ln -s $file $file.fasta; done
#       > ~kpepin/git/staging/run_ber_prep.csh TELCIRDFT
#
#         This step takes some time, overnight, it blasts each gene vs an nr db and runs hmmpfam vs an hmm db.
#
#         Check output and that the hmm and ber fof files created by the script are linked to src dir
#         ls -lt  /gscmnt/gc9002/info/annotation/BER/autoannotate_v2.5/src/*TELCIR*fof
#
#         For the next step you can cd src_dir
#        /gscmnt/gc9002/info/annotation/BER/autoannotate_v2.5/src, but script does it as well.
#
#    ii) Converting blastp & hmmpfam output to btab & htab files respectively
#        
#       Genome/Model/Tools/Ber/BerRunBtabhmmtab.pm
#
#        Details:
#        > ~kpepin/git/staging/run_ber_prep.step2.csh TELCIRDFT
#
#       This will take some time, overnight, as well to parse each of the blastp/hmmpfam files to btab and htab files.

}

sub _run_anno_sqlite {
	
	
#   Run anno-sqlite.bash
#   
#   Genome/Model/Tools/Ber/BerRunAnnoSqlite.pm
#   
#   Manually running details:
#  /gscmnt/gc9002/info/annotation/BER/autoannotate_v2.5/src
#  > bsub -o TELCIRDFT.out -e TELCIRDFT.err -R 'select[type=LINUX64]' ./anno-sqlite.bash  
#  TELCIRDFT 130521 gram-
#
#  This will take X hours to run, depends on how many genes to process. For TELCIRDFT  
#  Processed 25,567 genes
#  Started 08:01am
#  Ended
#
#  To check that the process is running ok
#  1. bjobs | grep gram
#  2. wc /gscmnt/gc9002/info/annotation/BER/autoannotate_v2.5/out/sqlite-locusID-date.dat
#
#Number of sequences in the pep file should be equal to the number of lines in the sqlite dat file.Should do this check in PAP to make sure all is complete.

}

sub _finish {
Genome/Model/Tools/Ber/BerRunFinish.pm
    
   This above module has a lot of stuff that’s very specific to MGAP pipeline. To me, this may not be of much use except the ace file generation part.

   Additional info:
   In house script to do the above:
   ~kpepin/scripts/parse_dat2ace.pl sqlite-locusID-date.dat > sqlite-locusID-date.dat.ace

   Parse into acedb and redump new .tbl file to get naming issues that need to be resolved

}

# Should contain all code necessary to parse the raw output of the predictor.
sub parse_output {
    die "Override in subclasses of " . __PACKAGE__;
}

# Any filtering logic should go here.
sub filter_results {
    die "Override in subclasses of " . __PACKAGE__;
}

# Should return a path to the file that should be executed using the current value of version.
sub tool_path_for_version {
	#possibly move version control into here
    die "Override in subclasses of " . __PACKAGE__;
}
    
# Should create an ace file from the raw output of the predictor
sub create_ace_file {
	#~kpepin/scripts/parse_dat2ace.pl [sample.dat]>[sample.dat.ace]
	#find out which variable represents "sample"
    die "Override in subclasses of " . __PACKAGE__;
}

sub read_config
{
    my $self = shift;
    
    my $conf = $self->config;
    unless(-f $conf)
    {
        carp "no config file $conf ... AmgapBerProtName.pm \n\n";
        return undef;
    }

    my $confhash = LoadFile($conf);

    return 1;
}

1;

