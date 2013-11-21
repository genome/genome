package Genome::Model::Tools::ManualReview::CreateDirectory;
use strict;
use warnings;
use Command;
use Data::Dumper;
use IO::File;
use PP::LSF;
use File::Temp;
use File::Basename;
class Genome::Model::Tools::ManualReview::CreateDirectory
{
    is => 'Command',                       
    has => 
    [ 
        map_list => 
        {
            type => 'String',
            is_optional => 0,
            doc => "File of input maps",
        },
        snp_file =>
        {
            type => 'String',
            is_optional => 0,
            doc => "File of variants",    
        },
        output_dir =>
        {
            type => 'String',
            is_optional => 0,
            doc => "Directory to generate Manual Review reports",    
        }
    ], 
};
############################################################
sub help_brief {   
    return;
}
sub help_synopsis { 
    return;
}
sub help_detail {
    return <<EOS 
    Creates a manual review directory from a given map list, snp file, and output dir.  
EOS
}
############################################################
sub execute { 
    my $self = shift;
    my $out_dir = $self->output_dir;
    
    my $snps = $self->snp_file;
    my $maplist = $self->map_list;
    #unless(`uname -m` =~ /x86_64/)
    #{
    #    $self->error_message( "manual-review must be run on a x64 system.\n");
    #    return 0;
    #}
    
    if(!-e $out_dir) {`mkdir -p $out_dir`;}
    my $gmt = `which gmt`;
    chomp($gmt);
    my $fh = IO::File->new($maplist);
    my @lines = <$fh>;
    chomp @lines;
    my @jobs;

    if(!-e "$out_dir/all.map")
    {
        `mkdir -p $out_dir/temp_map`;
        foreach my $line (@lines)
        {
            my $l = basename($line);
            my $f = $line;
            my $o = $f;
            
            $o = $out_dir."/temp_map/".$l;
            print $o,"\n";
            die "File $o alreadys exists!\n" if(-e $o);
            #system("bsub -q aml -oo $l.log gmt maq get-intersect --input=$f --snpfile=$snps --output=$o");
            my %job_params = (
                pp_type => 'lsf',
                q => $ENV{GENOME_LSF_QUEUE_SHORT},
                command => "$gmt maq get-intersect --input=$f --snpfile=$snps --output=$o",
                o => "$l.log",
            );
            my $job = PP::LSF->create(%job_params);
            $self->error_message("Can't create job: $!")
                and return unless $job;
            push @jobs, $job;            
        }
        
        foreach(@jobs)
        {
            $_->start;
        }
        while(1)
        {
            foreach(@jobs)
            {
                if(defined $_ && $_->has_ended){
                    if($_->is_successful) {$_ = undef;}
                    else { print Dumper $_; print "\n"; die "Job failed.\n"}                    
                }
            }
            foreach(@jobs)
            {
                if(defined $_) { goto SLEEP;}
            }
            last; #if we're here then we're done
SLEEP:      sleep 30;
        }
    }
    if(!-e "$out_dir/all.map")
    {
        @jobs = ();
        foreach (@lines) 
        {
            $_ = $out_dir."/temp_map/".basename($_);
        }
        my $command_string;
        my @temp_lines = @lines;
        my $final_merge = $ENV{GENOME_SW} . "/maq/maq-0.6.8_x86_64-linux/maq mapmerge $out_dir/all.map";
        for(my $i=0;$i<@lines/1000;$i++)
        {
            my @temp = splice(@temp_lines,0,1000);
            my $maps = join ' ',@temp;
            $command_string .=$ENV{GENOME_SW} . "/maq/maq-0.6.8_x86_64-linux/maq mapmerge $out_dir/temp_map/all.map.$i $maps\n";
            $final_merge .= " $out_dir/temp_map/all.map.$i";
        }
        $command_string .= "$final_merge\n";
        my $fh = IO::File->new(">$out_dir/temp_map/command");
        print $fh $command_string;
        $fh->close;
        `chmod 755 $out_dir/temp_map/command`;
        my $maps = join ' ',@lines;       
        #system("bsub -q aml -R 'select[type=LINUX64]'-oo mapmerge.log maq mapmerge $out_dir/all.map $maps");
        my %job_params = (
                pp_type => 'lsf',
                q => $ENV{GENOME_LSF_QUEUE_SHORT},
                R => 'select[type=LINUX64]',
                command => "$out_dir/temp_map/command",#"maq mapmerge $out_dir/all.map $maps",
                oo => "mapmerge.log",
            );
        my $job = PP::LSF->create(%job_params);
        $self->error_message("Can't create job: $!")
            and return unless $job;
        $job->start;
        while(1)
        {
            if($job->has_ended){
                if($job->is_successful) {last;}
                else {die "Job failed.\n";}                    
            }
            sleep 30;
        }            
    }
    if(1)
    {
        @jobs = ();
        $fh = IO::File->new($snps);    
        foreach my $line (<$fh>)
        {
            chomp $line;
            my ($seq, $pos) = split /\s+/,$line; 
            my $seqpos = $seq.'_'.$pos;
            if(!-e "$out_dir/$seqpos"){`mkdir -p $out_dir/$seqpos`;}
            my $temp_fh = IO::File->new(">$out_dir/$seqpos/annotation.tsv");
            print $temp_fh $line,"\n"; 
        }
        $fh->seek(0,0);
        foreach my $line (<$fh>)
        {   
            chomp $line;
            my ($seq, $pos) = split /\s+/,$line; 
            my $seqpos = $seq.'_'.$pos;
            #system("bsub -q aml -W 50 -oo $seqpos.log gmt maq get-intersect --input=$out_dir/all.map --snpfile=$out_dir/$seqpos/annotation.tsv --output=$out_dir/$seqpos/$seqpos --justname=2");
            my %job_params = (
                pp_type => 'lsf',
                q => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
                #W => 5,
                command => "$gmt maq get-intersect --input=$out_dir/all.map --snpfile=$out_dir/$seqpos/annotation.tsv --output=$out_dir/$seqpos/$seqpos --justname=2",
                oo => "$seqpos.log",
            );
            my $job = PP::LSF->create(%job_params);
            $self->error_message("Can't create job: $!")
                and return unless $job;
            push @jobs, $job; 
        }
        foreach(@jobs)
        {
            $_->start;
        }
        
        while(1)
        {
            foreach(@jobs)
            {
                if(defined $_ && $_->has_ended){
                    if($_->is_successful) {$_ = undef;}
                    else {$self->error_message( "gmt maq get-intersects failed.\n"); return;}                    
                }
            }
            foreach(@jobs)
            {
                if(defined $_) { goto SLEEEP;}
            }
            last; #if we're here then we're done
SLEEEP:      sleep 30;
        }
        #`rm -rf $out_dir/temp_map`;        
    }
    @jobs = ();
    
    my $proj_fof = File::Temp->new(UNLINK => 1);
    my @projects = `\\ls -d -1 $out_dir/*/`;
    @projects = grep { !($_ =~ /temp_map/) } @projects;
    print $proj_fof @projects;
    #print $out_dir,"\n",$proj_fof,"\n";
    
    $proj_fof->close;  
    my $gm = `which genome`;  
    chomp $gm;$gm .= " model";
    #return Genome::Model::Tools::PrepareNextgenAce->execute(fof => $proj_fof->filename, basedir => $out_dir);
    foreach my $line (@projects)
    {   
        chomp $line;
        my $seqposline = `cat $line/annotation.tsv`;
        chomp $seqposline;
        my ($seq, $pos) = $seqposline =~ /(\S+)\t(\S+)/;
        my $start = $pos-300;
        my $end = $pos +300;
        mkdir ("$line/edit_dir");
        mkdir("$line/phd_dir");

        my $ref_seq_dir = Genome::Config::reference_sequence_directory();

        my $command = "$gm write ace maq  --ace=$line/edit_dir/$seq"."_$pos.ace --chromosome=$seq --start=$start --end=$end --refseq=$ref_seq_dir/NCBI-human-build36/$seq.fasta --maq=$line/$seq"."_$pos.map";
        print $command,"\n";
        my %job_params = (
            pp_type => 'lsf',
            q => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
            R => 'select[type=LINUX64]',
            #W => 15,
            command =>$command,#"$gmt manual-review prepare-nextgen-ace --project-dir=$line --basedir=$out_dir",
            oo => "$line.log",
        );
        my $job = PP::LSF->create(%job_params);
        $self->error_message("Can't create job: $!")
            and return unless $job;
        push @jobs, $job; 
    }
    foreach(@jobs)
    {
        $_->start;
    }

    while(1)
    {
        foreach(@jobs)
        {
            if(defined $_ && $_->has_ended){
                if($_->is_successful) {$_ = undef;}
                else {$self->error_message( "$_->{command} failed\n"); return;}                    
            }
        }
        foreach(@jobs)
        {
            if(defined $_) { goto SLEEEEP;}
        }
        last; #if we're here then we're done
SLEEEEP:      sleep 30;
    }
    return 1;
}
############################################################
1;
