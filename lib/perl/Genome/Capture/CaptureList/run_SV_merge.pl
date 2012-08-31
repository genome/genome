
#use strict;
use Cwd;
use Config::IniFiles;
use File::Path qw(make_path);



my $ini_file = $ARGV[0];
if(!-e $ini_file || -z $ini_file) {
    print STDERR "Cannot find $ini_file.  Exiting\n";
    exit 5;
}

my $number;
$number = $ARGV[1] if(defined $ARGV[1]);

my $no_filter = 0;
$no_filter = $ARGV[2] if(defined $ARGV[2]);

my $cfg = Config::IniFiles->new( -file => $ini_file);

my $current_dir = getcwd();
my $project = rem_white_space($cfg->val('REQUIRED','project'));
$project =~ s/#NUMBER/$number/g;
my $out_dir = rem_white_space($cfg->val('REQUIRED','output_dir'));
$out_dir =~ s/#NUMBER/$number/g;
my ($proj,$num) = ($project =~ /([^\d]+)(\d+)/); 

#my $output_dir = "${out_dir}/${project}/no_filters";
my $output_dir = $out_dir;
#if(!-d "$output_dir") {
#    make_path("$output_dir");
#}

#default configuration file
my $config = "merge_config_$number.txt";
my $config_ctx = "merge_config_ctx_$number.txt";


make_config_file($ini_file,$project,$output_dir);
make_comparison();



sub make_comparison{


    #run intrachromosomal comparison
    print STDERR "run intrachromosomal comparison\n";
    my $cmd = "perl bigComparison_final_simple_more.pl -h -n -r 0.75 ${output_dir}/${config} > ${output_dir}/$project";
    &run_system($cmd);

    #run interchromosomal comparison
    print STDERR "run interchromosomal comparison\n";
    my $sv_file = $cfg->val('BreakDancer','file'); #grab a random SV file
    $cmd = "perl bigComparison_final_simple_more.pl -h -n -x 1 ${output_dir}/${config_ctx} ${sv_file} > ${output_dir}/${project}_ctx";
    &run_system($cmd);
if($no_filter == 1){
    return;
}
    #
    $cmd = "grep BD ${output_dir}/${project} | grep -v AS | perl -ane '\$a=0; foreach ( \@F[1..\$#F] ) { ( \$ncn,\$tcn ) = ( \$_=~/Ncn(\\S+)\:Tcn(\\S+)/ ) ; ( \$tp ) = (\$_=~/tp(\\S+):/); \$a = 1 if(\$ncn - \$tcn > 0.5 && \$tp =~ /DEL/ || \$tcn-\$ncn>0.5 && \$tp =~ /ITX/);} print \"\$_\" if(\$a);\' > ${output_dir}/${project}.CR2.tmp";
    &run_system($cmd);

    #filter step
    print STDERR "Run Filter\n";
    $cmd = "perl filters.pl ${output_dir}/${project}.CR2.tmp ${output_dir}/${project}.CR2 ${output_dir}/${project}.CR2.filteredout";
    &run_system($cmd);

    print STDERR "Run Count\n";
    
    # create a hash to do the filter work
    my $f = "${output_dir}/${project}";
    
    #$cmd = "./count_CRs.sh ${out_dir}/${project} $proj $num > $output_dir/${proj}.count1";
    
    my $output_capture = "${output_dir}/${project}.capture";

    if(-e $output_capture){
        `rm -f ${output_capture}`;
    }

    open FILE_out, ">", $output_capture || die $!;
    
    open FILE, "<", $f || die $!;
    my %filter;
    while(<FILE>){
        my ($ID,) = split("\t", $_);
        if($_ =~ /AS/ || $_ =~ /SD/ || $_ =~ /PD/ && $_ =~ /CN/ || $_ =~ /PD/ && $_ !~ /BD/ && $_ !~ /AS/ && $_ !~ /PDN/ || $_ =~ /HD/ || $_ =~ /CBS/){
            if($_ =~ /AS/){
                $filter->{AS} ++;
            }  elsif($_ =~ /SD/){
                $filter->{SD} ++;
            } elsif($_ =~ /PD/ && $_ =~ /CN/){
                $filter->{PDCN} ++;
            } elsif($_ =~ /PD/ && $_ !~ /BD/ && $_ !~ /AS/ && $_ !~ /PDN/){
                $filter->{PDT} ++;
            } elsif($_ =~ /HD/){
                $filter->{HD} ++;
            } elsif($_ =~ /CBS/){
                $filter->{CBS} ++;
            }
            $filter->{$ID} = 1;
            print FILE_out $_;
        }
    }

    close FILE;

    open FILE, "<", "$f.CR2" or die $!;

    while(<FILE>){
        my ($ID,) = split("\t", $_);
        if(!$filter->{$ID}){
            $filter->{BD} ++;
            $filter->{$ID} = 1;
            print FILE_out $_;
        }
    }

    close FILE;

    open FILE, "<", "${f}_ctx" or die $!;

    while(<FILE>){
        $filter->{CTX} ++;
        print FILE_out $_;
    }

    close FILE_out;

    my $output_count = "${output_dir}/${project}.count1";
    open FILE_OUT, ">", $output_count || die $!;
    print FILE_OUT "AS: " . $filter->{AS} . "\n";
    print FILE_OUT "SD: " . $filter->{SD} . "\n";
    print FILE_OUT "PD & CN: " . $filter->{PDCN} . "\n";
    print FILE_OUT "PDT: " . $filter->{PDT} . "\n";
    print FILE_OUT "HD: " . $filter->{HD} . "\n";
    print FILE_OUT "BD: " . $filter->{BD} . "\n";
    print FILE_OUT "CTX: " . $filter->{CTX} . "\n";
    print FILE_OUT "CBS: " . $filter->{CBS} . "\n";
    close FILE_out;


#    print STDERR "Running nimblegen from gmt...\n";          
    my $output_nimblegen = "${output_dir}/${project}.nimblegen";
    my $output_count2  = "${output_dir}/${project}.count2";
    my $filtered_out = "${output_dir}/${project}.nimblegen.filteredout";
    $cmd = "gmt nimblegen design-from-sv --sv-file ${output_capture}  --output-file ${output_nimblegen} --filtered-out-file ${output_filteredout} --count-file ${output_count2}"
#    &run_system($cmd);  

#    $cmd = "grep AS ${f} > ${output_capture}";
#    #&run_system($cmd);
#    system($cmd);
#    print STDERR "    Run Count on SD (CR1)\n";
#    $cmd = "grep SD ${f} | grep -v AS >> ${output_capture}";
#    #&run_system($cmd);
#    system($cmd);
#    #print "${f}.CR2\n";
#    print STDERR "    Run Count on filtered unassembled breakdancer (CR2)\n";
#    $cmd = "cat ${f}.CR2 >> ${output_capture}";
#    #$run_system($cmd);
#    system($cmd);

#    print STDERR "    Run Count on Pindel and Copy Number (CR3)\n";
#    $cmd = "grep PD ${f} | grep CN >> ${output_capture}";
#    #&run_system($cmd);
#    system($cmd);
#    print STDERR "    Run Count on Pindel Tumor only (CR4)\n";
#    $cmd = "grep PD ${f} | grep -v CN | grep -v BD | grep -v AS | grep -v PDN >> ${output_capture}";
#    #&run_system($cmd);
#    system($cmd);
#    print STDERR "    Run Count on Copy Number only (CR5)\n";
#    $cmd = "grep CN ${f} | grep -v BD | grep -v PD | grep -v AS | perl -ane '(\${s})=(\$F[1]=~/sc(\\S+)\\:/); print \"\$_\" if(\$s>=100)' >> ${output_capture}";
#    #&run_system($cmd);
#    system($cmd);
#    print STDERR "    Run Count on Hydra only\n";
#    $cmd = "grep HD ${f} | grep -v BD | grep -v PD | grep -v AS | grep -v CN >> ${output_capture}";
#    #&run_system($cmd);
#    system($cmd);
#    print STDERR "    Run Count on CBS only\n";
#    $cmd = "grep CBS ${f} | grep -v CN >> ${output_capture}";
#    system($cmd);
#    print STDERR "Run Count on CTX only\n";
#    $cmd = "cat ${f}_ctx >> ${output_capture}";
#    #&run_system($cmd);
#    system($cmd);

}


sub run_system {
    my $cmd = shift;
    my $exit_code = 0;
    $exit_code = system($cmd);
    while($exit_code) {
	print STDERR "rerunnning $cmd\n";
	$exit_code = system($cmd);
    }
  


}


sub inter_chromosomal_comparison {

    
    my $sv_file = $cfg->val('BreakDancer_FILE','chr1'); #grab a random SV file

    my $cmd = "perl bigComparison_final_simple_more.pl -h -n -x 1 ${output_dir}/${config_ctx} ${sv_file} > ${output_dir}/${project}_ctx";
    system($cmd);


}


sub make_config_file {

    my $ini_file = shift;
    my $project = shift;
    my $out_dir = shift;

    my $cfg = Config::IniFiles->new( -file => $ini_file );

    open(OUT, "> $output_dir/$config" ) or die "Unable to write to config.txt\n";
    
    my ($header,$loc,$line_skip,$rule);
      
    #assembly_file
    if($cfg->val('Assembly','all_file')){
    $header = "${project}.AS";
    $loc = rem_white_space($cfg->val('Assembly','all_file'));  
    $loc =~ s/#NUMBER/$number/g;
    $line_skip = 0;
    $rule = rem_white_space($cfg->val('Rules','AS'));
    print OUT "$header\t$loc\t$line_skip\t$rule\n";
}

    ($header,$loc,$line_skip,$rule)=();

    if($cfg->val('SquareDancer','file')){
    #squaredancer_file
    $header = "${project}.SD";
    $loc = rem_white_space($cfg->val('SquareDancer','file'));  
    $loc =~ s/#NUMBER/$number/g;
    $line_skip = 0;
    $rule = $cfg->val('Rules','SD');
    print OUT "$header\t$loc\t$line_skip\t$rule\n";
}
    
    ($header,$loc,$line_skip,$rule)=();
   if($cfg->val('Pindel','tumor_file')){
      $header="${project}.PDT";
     $loc = rem_white_space($cfg->val('Pindel','tumor_file'));
    $loc =~ s/#NUMBER/$number/g;
   $line_skip = 0;
   $rule = rem_white_space($cfg->val('Rules','PDT'));
   print OUT "$header\t$loc\t$line_skip\t$rule\n";
  } 

    ($header,$loc,$line_skip,$rule)=();
   if($cfg->val('Pindel','normal_file')){
      $header="${project}.PDN";
     $loc = rem_white_space($cfg->val('Pindel','normal_file'));
    $loc =~ s/#NUMBER/$number/g;
   $line_skip = 0;
   $rule = rem_white_space($cfg->val('Rules','PDN'));
   print OUT "$header\t$loc\t$line_skip\t$rule\n";
  } 

    ($header,$loc,$line_skip,$rule)=();
   #breakdancer_files
   if($cfg->val('BreakDancer','dir')){
    my @chrom = (1 .. 22);
    my @chrom = (@chrom, 'X', 'Y');

    my $header_base = "${project}.BD";
    foreach my $x (@chrom) {
	my $suffix = "chr${x}";
	$header = "$header_base.${suffix}"; 
	my $loc_dir = rem_white_space($cfg->val('BreakDancer','dir'));
    $loc_dir =~ s/#NUMBER/$number/g;
	$loc = "$loc_dir/${project}.${suffix}.sv";
	$line_skip = 0;
	$rule = $cfg->val('Rules','BD');
	print OUT "$header\t$loc\t$line_skip\t$rule\n";
    }
}
    ($header,$loc,$line_skip,$rule)=();
    #Copy-number_file
    if($cfg->val('CNA','file')){
    $header = "${project}.CNA";
    $loc = rem_white_space($cfg->val('CNA','file'));  
    $loc =~ s/#NUMBER/$number/g;
    $line_skip = 0;
    $rule = rem_white_space($cfg->val('Rules','CN'));
    print OUT "$header\t$loc\t$line_skip\t$rule\n";
}

    
    ($header,$loc,$line_skip,$rule)=();
   if($cfg->val('HYDRA','file')){
      $header="${project}.HD";
     $loc = rem_white_space($cfg->val('HYDRA','file'));
    $loc =~ s/#NUMBER/$number/g;
   $line_skip = 0;
   $rule = rem_white_space($cfg->val('Rules','HD'));
   print OUT "$header\t$loc\t$line_skip\t$rule\n";
  } 

    ($header,$loc,$line_skip,$rule)=();
    #cbs_file
    if($cfg->val('CBS','file')){
    $header = "${project}.CBS";
    $loc = rem_white_space($cfg->val('CBS','file'));  
    $loc =~ s/#NUMBER/$number/g;
    $line_skip = 0;
    $rule = rem_white_space($cfg->val('Rules','CBS'));
    print OUT "$header\t$loc\t$line_skip\t$rule\n";
}

    close OUT;

    open(OUT, "> $output_dir/$config_ctx" ) or die "Unable to write to config_ctx.txt\n";
    ($header,$loc,$line_skip,$rule)=(); 
    #assembly_ctx_file
    if($cfg->val('Assembly','CTX_file')){
    $header = "${project}.AS";
    $loc = rem_white_space($cfg->val('Assembly','CTX_file'));  
    $loc =~ s/#NUMBER/$number/g;
    $line_skip = 0;
    $rule = rem_white_space($cfg->val('Rules','CTX'));
    print OUT "$header\t$loc\t$line_skip\t$rule\n";
}

    if($cfg->val('SquareDancer','file')){
    $header = "${project}.SD";
    $loc = rem_white_space($cfg->val('SquareDancer','file'));  
    $loc =~ s/#NUMBER/$number/g;
    $line_skip = 0;
    $rule = rem_white_space($cfg->val('Rules','SD'));
    print OUT "$header\t$loc\t$line_skip\t$rule\n";
}

    close OUT;

}


sub rem_white_space {

    my $string = shift;

    $string =~ s/^\s+//g;
    $string =~ s/\s+$//g;

    return $string;

}
