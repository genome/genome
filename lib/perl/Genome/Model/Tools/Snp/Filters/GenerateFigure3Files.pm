package Genome::Model::Tools::Snp::Filters::GenerateFigure3Files;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Statistics::R;
use Workflow;

class Genome::Model::Tools::Snp::Filters::GenerateFigure3Files
{
    is => 'Command',
    has => [
        basedir                     => { is => 'String', 
                                         is_input => 1,
                                         doc => 'Put some docmumentation here' },
        snp_report_file             => { is => 'String',
                                         is_optional=> 1, 
                                         doc => 'Put some docmumentation here' },
        ref_seq_id                  => { is => 'String', 
                                         is_input => 1,
                                         doc => 'Put some docmumentation here' },
        binomial_output_file        => { is => 'String', 
                                         is_input => 1,
                                         doc => 'Put some docmumentation here' },
        parent_event                => { is => 'Genome::Model::Event',
                                         is_optional=>1,
                                         is_input => 1,
                                         doc => 'The parent event' },
        validated_somatic_variants  => { is => 'String', 
                                         is_optional=>1,
                                         is_output => 1,
                                         doc => 'Put some docmumentation here' },
        log_metric       =>             { is => 'Boolean',
                                         is_optional=>1,
                                         default=>1,
                                         doc=>'should this run write to the db about its statistics or not',
                                       },
    ], 
};

#----------------------------------
sub execute {
    my $self = shift;
    $DB::single = $DB::stopper;
    if(!defined $self->log_metric) {
        unless($self->snp_report_file) {
            $self->error_message("If you are not supplying a parent event and do not wish to log metrics, you must supply a relevant annotation file.");
            return;
        }
    }
    elsif($self->log_metric) {
        unless($self->parent_event) {
            $self->error_message("If you wish to log metrics, you must supply a parent event to tie them to and get annotation from.");
            return;
        }
    }
    
    unless($self->generate_figure_3_files()){
        return;
    }
    return 1;
}
#----------------------------------
 
sub generate_figure_3_files {
    my $self = shift;  
    
    my $somatic_file = $self->binomial_output_file;
    my $dir= $self->basedir;
    my $snp_file=$self->snp_report_file;
    unless($snp_file) {
        my $prior = Genome::Model::Event->get(id => $self->parent_event->prior_event);
        $snp_file= $prior->snp_report_file;
    }
    #Hash for ranking somatic statuses
    my %rank = ('S' => 0,
        'G' => 1,
        'WT' => 2,
        'A' => 3,
        'V' => 3,
        'O' => 3,
        'LQ' => 4,
        'X' => 4,
        'NC' => 5,
    );
    #END HASH
    
    #my $snp_file = $self->_report_file('snp');
      my $dbsnp_fh = IO::File->new(">$dir" . "/somatic_variants_in_d_v_w_" . $self->ref_seq_id .
            ".csv");
       my $dbsnp_count=0;      
       my $non_coding_fh = IO::File->new(">$dir" .
           "/non_coding_tumor_only_variants_" . $self->ref_seq_id .
            ".csv");
       my $non_coding_count=0;      
       my $novel_tumor_fh = IO::File->new(">$dir" .
           "/novel_tumor_only_variants_" . $self->ref_seq_id .
            ".csv");
       my $novel_tumor_count=0;      
       my $silent_fh = IO::File->new(">$dir" .
           "/silent_tumor_only_variants_" . $self->ref_seq_id .
            ".csv");
       my $silent_count=0;     
       my $nonsynonymous_fh = IO::File->new(">$dir" .
           "/non_synonymous_splice_site_variants_" . $self->ref_seq_id .
            ".csv");
       my $nonsynonymous_count=0;      
       my $var_never_manreview_fh = IO::File->new(">$dir" .
           "/var_never_manreview_" . $self->ref_seq_id . ".csv");
       my $never_manreview_count=0;
       my $var_pass_manreview_fh = IO::File->new(">$dir" .
           "/var_pass_manreview_" . $self->ref_seq_id .
            ".csv");
       my $var_pass_manreview_count=0;     
       my $var_fail_manreview_fh = IO::File->new(">$dir" .
           "/var_fail_manreview_" . $self->ref_seq_id .
            ".csv");
       my $var_fail_manreview_count=0;     
       my $var_fail_valid_assay_fh = IO::File->new(">$dir" .
           "/var_fail_valid_assay_" . $self->ref_seq_id .
            ".csv");
       my $var_fail_valid_assay_count=0;      
       my $var_complete_validation_fh = IO::File->new(">$dir" .
           "/var_complete_validation_" . $self->ref_seq_id .
            ".csv");
       my $var_complete_validation_count=0;      
       my $validated_snps_fh = IO::File->new(">$dir" .
           "/validated_snps_" . $self->ref_seq_id .
            ".csv");
       my $validated_snps_count=0;     
       my $false_positives_fh = IO::File->new(">$dir" .
           "/false_positives_" . $self->ref_seq_id .
            ".csv");
       my $false_positives_count=0;     
       my $validated_somatic_var_fh = IO::File->new(">$dir" .
           "/validated_somatic_variants_" . $self->ref_seq_id .
            ".csv");
       my $validated_somatic_var_count=0;     
        #added to track things that were passed through manual review but
        #don't have a validation status. Could be pending or could be missing
        #from db
       my $passed_but_no_status_count=0;     
       my $passed_but_no_status_fh = IO::File->new(">$dir" .
           "/var_pass_manreview_but_no_val_status_" . $self->ref_seq_id .
            ".csv");
         my $annotation_fh = IO::File->new($snp_file);
        if(!defined($annotation_fh)) {
            $self->error_message("Could not open report file: $snp_file");
            return;
        }
        my $somatic_fh = IO::File->new($somatic_file);
        if(!defined($somatic_fh)) {
            $self->error_message("Could not open file of somatic mutations: $somatic_file.");
            return;
        }
        my @cur_somatic_snp;
        my @cur_anno_snp;
        my $anno_line;
        my $somatic_line;
        #throw away header
        $somatic_fh->getline;
        #throw away the header, but preload anno_line so our loop gets off the ground. i rate this hack:medium special
        $anno_line = $annotation_fh->getline;
        #end throw away header section
      while(($somatic_line=$somatic_fh->getline) && defined $anno_line) {
          chomp $somatic_line;
          if (!defined $somatic_line || !defined $anno_line) {
              #the filter file must be over if we're here
              last;
              
          }
          @cur_somatic_snp = split(/,\s*/, $somatic_line);
          
          while(!@cur_anno_snp || ($cur_anno_snp[1] < $cur_somatic_snp[1]) ) {
              #we hit this block because a) this is our first time through
              #or b) the last annotation position is smaller than the current somatic snp
              # snp value   
              $anno_line = $annotation_fh->getline;
              chomp $anno_line;
              if(!defined $anno_line) {
                  $self->error_message("Annotation file has ended before somatic file. This may be ok.");
                  $self->error_message("Last somatic snp was\n " . join (" ", @cur_somatic_snp) . "last anno snp was\n " . join(" ", @cur_anno_snp) );
                  last;
              }
              @cur_anno_snp= split (/,\s*/, $anno_line);
          }
         
          
          while($cur_anno_snp[1] == $cur_somatic_snp[1] ) {
          #if we get here then the idea is we have a somatic line with the same position as a snp line.
     
             #call in Brian's somatic file and Eddie's report

             #For Eddie's output we need to know the type of variant and also the
             #dbSNP and Watson/Venter status
             my @report_indexes = (0,1,2,3,5,8,13,18,19,20);
             #larson erroneously claimed this was correct
             #my @report_indexes = (0,1,2,3,5,8,13,16,17,18); 

             #this is taken care of implicitly by the loop actually...damn pair programming
             if(defined($cur_somatic_snp[0]) && defined($cur_anno_snp[0])) {
                 #it's genic and in Eddie's report and passed Brian's filters
                 my ($chromosome, $begin, $end,
                     $variant_allele, $reference_allele, $gene, $variant_type,$dbsnp,
                     $watson, $venter) = @cur_anno_snp[@report_indexes];    

                 #Test if seen in dbSNP or Watson/Venter
                 if((defined($dbsnp) && $dbsnp ne '0') || (defined($watson) && $watson ne '0' ) || (defined($venter) && $venter ne '0')) {
                     #previously identified
                     $self->_write_array_to_file(\@cur_anno_snp, $dbsnp_fh);
                     $dbsnp_count++;
                 }
                 else {
                     $self->_write_array_to_file(\@cur_anno_snp,
                         $novel_tumor_fh);    
                     $novel_tumor_count++;    
                     #nonsynonymous
                     if( $variant_type =~ /missense|nonsense|nonstop|splice_site/i) {
                         #output those that are coding
                         $self->_write_array_to_file(\@cur_anno_snp,
                             $nonsynonymous_fh);
                             $nonsynonymous_count++;

                         my $variant_detail = Genome::VariantReviewDetail->get(
                             chromosome => $chromosome, 
                             begin_position => $begin, 
                             #end_position => $end,
                             variant_type => 'S',
                             #insert_sequence_allele1 => $variant_allele,
                             #delete_sequence => $reference_allele,
                         );
                         #this if is the same thing with the alleles commented out...but if you uncomment them then they're different. 
                         #leave in in case we go back to this.
                         if(!defined($variant_detail)) {
                             #try and see if it was a biallelic variant site
                             $variant_detail = Genome::VariantReviewDetail->get(
                                 chromosome => $chromosome, 
                                 begin_position => $begin, 
                                 variant_type => 'S',
                                 #end_position => $end,
                                 #insert_sequence_allele2 => $variant_allele,
                                 #delete_sequence => $reference_allele,
                             );
                         }
                         if(defined($variant_detail)) {
                             #it's been sent to manual review
                             my $decision = $variant_detail->pass_manual_review; 
                             #if there is a somatic status then it probably
                             #passed manual review but wasn't documented.
                             #Or it was directly passed along.
                             my $db_status = $variant_detail->somatic_status;
                             my $status;
                             #take care of the possibility of discrepancy
                             #NOTE: this only checks for two levels of discrepancy
                             #This should really be recursive
                             my ($level1_discrep, $level2_discrep);
                                 ($level1_discrep) = $db_status =~ /^DISCREPANCY\((.*)/ if $db_status;
                                 if(defined($level1_discrep)) {
                                     ($level2_discrep) = $level1_discrep =~ /^DISCREPANCY\((.*)/;
                                 }
                                 if(defined($level2_discrep)) {
                                     my @status = split /:/, $level2_discrep;
                                     map {s/^\s*(\S*)\s*$/$1/} @status;
                                     my $new_status = 'NC';
                                     foreach my $discr_status (@status) {
                                         if($rank{$discr_status} < $rank{$new_status} ) {
                                             $new_status = $discr_status;
                                         }
                                     }
                                     $status = $new_status;
                                 }
                             elsif(defined($level1_discrep)) {
                                 my @status = split /:/, $level1_discrep;
                                 map {s/^\s*(\S*)\s*$/$1/} @status;
                                 my $new_status = 'NC';
                                 foreach my $discr_status (@status) {
                                     unless(defined($rank{$discr_status})) {
                                         $self->error_message("No entry for |$discr_status|");
                                     }
                                     if($rank{$discr_status} < $rank{$new_status} ) {
                                         $new_status = $discr_status;
                                     }
                                 }
                                 $status = $new_status;
                             }
                             else {
                                 $status = $db_status;
                             }


                             if(defined($decision) && lc($decision) eq 'yes'
                                 || defined($status)) {
                                 $self->_write_array_to_file(\@cur_anno_snp,
                                     $var_pass_manreview_fh);
                                 $var_pass_manreview_count++;

                                 if(defined($status)) {
                                     $status = uc($status);
                                     if($status eq 'S') {
                                         $self->_write_array_to_file(\@cur_anno_snp,
                                             $validated_somatic_var_fh);
                                         $validated_somatic_var_count++;
                                     }
                                     elsif($status eq 'WT') {
                                         $self->_write_array_to_file(\@cur_anno_snp,
                                             $false_positives_fh);
                                         $false_positives_count++;
                                     }
                                     elsif($status eq 'G') {
                                         $self->_write_array_to_file(\@cur_anno_snp,
                                             $validated_snps_fh);
                                         $validated_snps_count++;
                                     }
                                     else {
                                         $self->_write_array_to_file(\@cur_anno_snp, $var_fail_valid_assay_fh);
                                         $var_fail_valid_assay_count++;
                                     }

                                 }
                                 else {
                                     #else no validation status
                                     $passed_but_no_status_count++;     
                                     $self->_write_array_to_file(\@cur_anno_snp, $passed_but_no_status_fh);  
                                 }


                             }
                             else {
                                 #TODO This may not be valid if for instance
                                 #maybe's were passed along to validation
                                 $self->_write_array_to_file(\@cur_anno_snp,
                                     $var_fail_manreview_fh);
                                 $var_fail_manreview_count++;
                             }

                         }
                         else {
                             #we're actually not tracking this case in the
                             #figure
                             #Not in database, but should be!
                             #Either novel or missing from database
                             $never_manreview_count++;
                             $self->_write_array_to_file(\@cur_anno_snp,
                                 $var_never_manreview_fh);
                         }

                     }
                     elsif( $variant_type eq 'silent') {
                         $self->_write_array_to_file(\@cur_anno_snp,
                             $silent_fh);
                             $silent_count++;
                     }
                     else {
                         $self->_write_array_to_file(\@cur_anno_snp,
                             $non_coding_fh);
                             $non_coding_count++;
                     }
                 }
             }

       
     
             $anno_line = $annotation_fh->getline;
             chomp $anno_line;
              if(!defined $anno_line) {
                  #annotation file has ended before somatic one has...this is bad
                  $self->error_message("Annotation file has ended before somatic file. This is probably ok.\n");
                  $self->error_message("Last somatic snp was\n " . @cur_somatic_snp . "last anno snp was\n " . @cur_anno_snp );
                  last;
              }
              @cur_anno_snp= split (/,\s*/, $anno_line);
             
          }
      }
      #close all filehandles.

      if($self->log_metric) {
          my $metric;
          my $event = $self->parent_event;
          $metric = $event->add_metric(
              name=> "somatic_variants_in_d_v_w",
              value=> $dbsnp_count,
          );

          $metric = $event->add_metric(
              name=> "non_coding_tumor_only_variants",
              value=> $non_coding_count,
          );


          $metric = $event->add_metric(
              name=> "novel_tumor_only_variants",
              value=> $novel_tumor_count,
          );


          $metric = $event->add_metric(
              name=> "silent_tumor_only_variants",
              value=> $silent_count,
          );
          $metric = $event->add_metric(
              name=> "non_synonymous_splice_site_variants",
              value=> $nonsynonymous_count,
          );


          $metric = $event->add_metric(
              name=> "var_pass_manreview",
              value=> $var_pass_manreview_count,
          );


          $metric = $event->add_metric(
              name=> "var_fail_manreview",
              value=> $var_fail_manreview_count,
          );
          $metric = $event->add_metric(
              name=> "var_fail_valid_assay",
              value=> $var_fail_valid_assay_count,
          );
          $metric = $event->add_metric(
              name=> "var_complete_validation",
              value=> $var_complete_validation_count,
          );

          $metric = $event->add_metric(
              name=> "validated_snps",
              value=> $validated_snps_count,
          );
          $metric = $event->add_metric(
              name=> "false_positives",
              value=> $false_positives_count,
          );
          $metric = $event->add_metric(
              name=> "validated_somatic_variants",
              value=> $validated_somatic_var_count,
          );
          $metric = $event->add_metric(
              name=> "var_never_sent_to_manual_review",
              value=> $never_manreview_count,
          );
          $metric = $event->add_metric(
              name=> "var_pass_manreview_but_no_val_status",
              value=> $passed_but_no_status_count,
          );
      }




      $dbsnp_fh->close;
      $non_coding_fh->close;
      $novel_tumor_fh->close;
      $silent_fh->close;
      $nonsynonymous_fh->close;
      $var_pass_manreview_fh->close;
      $var_fail_manreview_fh->close;
      $var_fail_valid_assay_fh->close;
      $var_complete_validation_fh->close;
      $validated_snps_fh->close;
      $false_positives_fh->close;
      $validated_somatic_var_fh->close;
      $annotation_fh->close;
      $somatic_fh->close;

  }

sub _write_array_to_file {
    my $self=shift;
    my $array_ref_to_write=shift;
    my $file_handle=shift;

    my $line_to_write = join (" ", @{$array_ref_to_write});
    $file_handle->print($line_to_write . "\n");

    return 1;
}
