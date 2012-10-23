# -*-Perl-*-

######################################
# Copyright (C) 2002 Shiaw-Pyng YAng #
# Washington University, St. Louis   #
# All Rights Reserved.               #
######################################

package ChimpaceObjects;

use strict;

##########################
# Methods to be exported #
##########################

require Exporter;

@aceObjects::ISA    = qw (Exporter);
@aceObjects::EXPORT = qw (get_ace_contigDNA;
			  get_ace_contigQual;
			  get_all_contigNames;
			  Total_Contigs;
			  Total_Reads;
			  Total_Bases;
			  Contig_Orientation;
			  Contig_Length;
			  Read_Orientation;
			  Read_StartPos;
			  Read_EndPos;
			  ReadsInContig;
			  Number_ReadsInContig;
			  getAlignUnpadedClipStart;
			  getAlignUnpadedClipEnd);
			  
##########################
# Methods to be imported #
##########################

########################################
# PACKAGE DEFINITION FOR "ace" OBJECTS #
########################################

sub new {    # file handle type inputs
    my ( $class, %INPUT ) = @_;

    my ( $FH, $parse_qual, $parse_seq ) =
      ( $INPUT{'-acefilehandle'}, $INPUT{'-parse_qual'}, $INPUT{'-parse_seq'} );
    
    # Output: object(instance) of aceObjects class
    my $ace;

    # Create hash pointer
    $ace = {};

    # Bless it into current aceObjects class
    bless $ace, $class;

    $ace->{'ContigName'} = [];

    #parrsing the file handel and store what is needed
    #set input separator to be paragraph
    local $/ = "";

    ###########################################
    #modify
    #############################################
    my ($pleft,$pright,$upleft,$upright,$strand);

    while (<$FH>) {
        if ( length($_) >= 3 ) {
            if ( substr( $_, 0, 3 ) eq "AS " ) {
                my ( $ctg_num, $read_num ) = $_ =~ /^AS (\d+)\s+(\d+)/;
                $ace->{'CONTIGS_NUMBER'} = $ctg_num;
                $ace->{'READS_NUMBER'}   = $read_num;
            }

            #CO <contig name> <# of bases> <# of reads in contig> <# of base segments in contig> <U or C>
            if ( substr( $_, 0, 3 ) eq "CO " ) {

                # found a CO line
                my ( $ContigName, 
		     $contig_bases_num, 
		     $reads_in_contig_num,
		     $ori )  = $_ =~ /^CO (\S+)\s+(\d+)\s+(\d+)\s+\d+\s+([CU])/;
                $::Contig_Name = $ContigName;

                #delete the CO line and all new lines and whit space, remove any pads
                s/CO.*\n//;
                s/\n//g;
                s/\s+//g;
                

##################################################################
#Count pads before depad
##################################################################
		my $pcnt=0;
		my @bases=split //,$_;
		my @pvec=(0) x $contig_bases_num;
		for(my $n=0;$n<=$#bases;$n++) {
		    ++$pcnt if($bases[$n] eq '*');
		    $pvec[$n]=$pcnt;
		    #$ace->{'pvec'}->{$ContigName}->[$n]=$pcnt;
		}
		$ace->{'pvec'}->{$ContigName}=\@pvec;
##################################################################
##################################################################
		s/\*//g;
                tr/xnacgt/XNACGT/;

                #store the data
                my $contig_length = length($_);
                push ( @{ $ace->{'ContigName'} }, $::Contig_Name );
                $ace->{'DNA'}->{$ContigName}              = $_;
                $ace->{'ContigPadLength'}->{$ContigName}  = $contig_bases_num;
                $ace->{'ContigOrient'}->{$ContigName}     = $ori;
                $ace->{'NumReadsInContig'}->{$ContigName} = $reads_in_contig_num;
                $ace->{'ContigLength'}->{$ContigName}     = $contig_length;
            }

            # BQ list of base qualities for the unpadded consensus bases
            if ( substr( $_, 0, 3 ) eq "BQ\n" ) {

                # found a BQ line
                #delete the BO line and all new lines and white space, remove any pads

            }

            #AF BS block
            if ( substr( $_, 0, 3 ) eq "AF " ) {

                #find AF BS block
                my @AF_lines = grep ( /^AF /, split ( /\n/, $_ ) );

                foreach my $af_line (@AF_lines) {
                    chomp $af_line;
                    my ( $af, $read, $read_ori, $start_pos ) =  split ( /\s+/, $af_line );
                    my ( $template, $ext ) = split ( /\./, $read );
                    my ($rtype) = $ext =~ /^([sfrxyzbigca])\d+/;

                    push ( @{ $ace->{'ReadsBelongTo'}->{$::Contig_Name} },$read );
                    $ace->{'ReadOri'}->{$read}      = $read_ori;
		    my $changed_start_pos ;
		    if ($start_pos < 0) {
			 $changed_start_pos = 1;
			 $ace->{'MatchSense'}->{$read} = 1;
			 $ace->{'OrigReadStartPos'}->{$read} = $start_pos ;
		     } else {
			 $changed_start_pos = $start_pos ;
			 $ace->{'MatchSense'}->{$read} = 0;
			 $ace->{'OrigReadStartPos'}->{$read} = $start_pos ;
		     }
                    $ace->{'ReadStartPos'}->{$read} = $changed_start_pos;
		}
            }

            #--- RD Block ---#
            if ( substr( $_, 0, 3 ) eq "RD " ) {

                #RD <read name> <# of padded bases> <# of whole read info items> <# of readtags>
                my ( $read, $numofpadbases, $numreadinfo, $numreadtags ) =
                  $_ =~ /^RD (\S+)\s+(\d+)\s+(\d+)\s+(\d+)/;
                $::Read = $read;

                #delete the RD line and all new lines and whit space, remove any pads
                s/^RD.*\n//;
                s/\n//g;
                s/\s+//g;
                s/\*//g;
                $ace->{'ReadDNA'}->{$read}             = $_;
                $ace->{'ReadPaddedLength'}->{$read}    = $numofpadbases;
                $ace->{'ReadNonPaddedLength'}->{$read} = length($_);
                $ace->{'NumWholeReadInfo'}->{$read}    = $numreadinfo;
                $ace->{'NumReadTags'}->{$read}         = $numreadtags;
		$ace->{'ReadEndPos'}->{$read}  =  $ace->{'OrigReadStartPos'}->{$read}+ $numofpadbases - 1; 
#################################################################################
#MODIFY 
#################################################################################
		
		
		$pleft=$ace->{'OrigReadStartPos'}->{$read};
		$pright=$ace->{'ReadEndPos'}->{$read};
		$strand=$ace->{'ReadOri'}->{$read};
		
		
		if($pleft >= 1) {
		    $upleft=$pleft - $ace->{'pvec'}->{$::Contig_Name}->[$pleft-1];
		    
		}
		elsif($strand eq 'U') {
		    $upleft=1 - $ace->{'pvec'}->{$::Contig_Name}->[0];
		}
		else { $upleft= $pleft; }
		
		$ace->{'UnPaddedReadStartPos'}->{$read} = $upleft;
		if($pright <= $ace->{'ContigPadLength'}->{$::Contig_Name}) {
		    
		    $upright= $pright - $ace->{'pvec'}->{$::Contig_Name}->[$pright-1];
		}
		elsif($strand eq 'C') {
		    $upright=$ace->{'ContigLength'}->{$::Contig_Name};
		}
		else {
		    
		    $upright=$pright - $ace->{'pvec'}->{$::Contig_Name}->[$#{$ace->{'pvec'}->{$::Contig_Name}}];
		}
		$ace->{'UnPaddedReadEndPos'}->{$read} = $upright;
		#warn "pleft = $pleft,$pright,$strand,$read, $upleft,$upright \n";
		 
#####################################################################
#####################################################################



      
	    }
            if ( substr( $_, 0, 3 ) eq "QA " ) {

                #QA <qual clipping start> <qual clipping end> <align clipping start> <alignclipping end>
                #If the entire read is low quality, then <qual clipping start> and <qual clipping end> will both be -1
                my ($qualclip_start,  $qualclip_end,$alignclip_start, $alignclip_end
                  ) = $_ =~ /^QA ([-]?\d+)\s+([-]?\d+)\s+([-]?\d+)\s+([-]?\d+)/;
                #$alignclip_start = 0 if (! defined $alignclip_start);
		#$alignclip_end = 1
		next if ($alignclip_start == -1 && $alignclip_end == -1);
                $ace->{'QualClipStart'}->{$::Read} = $qualclip_start;
                $ace->{'QualClipEnd'}->{$::Read}   = $qualclip_end;
                $ace->{'AlignClipStart'}->{$::Read} = $alignclip_start;
                $ace->{'AlignClipend'}->{$::Read}  = $alignclip_end;


		###################################################
		#######################################################
		#firstly, defined padded $alignclipleft,$alignclipright
		#secondly, use the same formula for unpad position
		my ($alignclipleft,$alignclipright);
		#if ($strand eq 'U'){
		    $alignclipleft = $pleft+ $alignclip_start - 1;
		    $alignclipright = $pleft + $alignclip_end - 1; 
		#}
		#else
		#{
		    #$alignclipright = $pright - $alignclip_start + 1;
		    #$alignclipleft = $pright - $alignclip_end + 1;
		#}
		
		
		if($alignclipleft >= 1) {
		    $ace->{'AlignUnpadedClipStart'}->{$::Read} = $alignclipleft - $ace->{'pvec'}->{$::Contig_Name}->[$alignclipleft-1];
		    
		}
		elsif($strand eq 'U') {
		    $ace->{'AlignUnpadedClipStart'}->{$::Read} = 1 - $ace->{'pvec'}->{$::Contig_Name}->[0];
		}
		else { $ace->{'AlignUnpadedClipStart'}->{$::Read} = $alignclipleft; }
		
		
		if($alignclipright <= $ace->{'ContigPadLength'}->{$::Contig_Name} && $alignclipright > 0) {
		    
		    $ace->{'AlignUnpadedClipEnd'}->{$::Read} = $alignclipright - ($ace->{'pvec'}->{$::Contig_Name}->[$alignclipright-1]);
		}
		elsif($strand eq 'C') {
		    $ace->{'AlignUnpadedClipEnd'}->{$::Read} = $ace->{'ContigLength'}->{$::Contig_Name};
		}
		else {
		    
		    $ace->{'AlignUnpadedClipEnd'}->{$::Read} = $alignclipright - $ace->{'pvec'}->{$::Contig_Name}->[$#{$ace->{'pvec'}->{$::Contig_Name}}];
		}
		
		###################################
		######################################
                s/^QA.*\n//;
            }
            if ( substr( $_, 0, 16 ) eq "DS CHROMAT_FILE:" ) {

                #DS CHROMAT_FILE: <name of chromat file> <name of PHD file> <chemistry> <dye> <date/time of the phd file> <template> 
                #my ( $junk, $chromat_file, $phd_file, $chemistry, $dye,
                #    $date_time_phd, $template, $direction )
                #  = split ( /[A-Z_]+\:/, $_ );
                #$ace->{'CHROMAT_FILE'} = chomp($chromat_file) if (defined $chromat_file);
                #$ace->{'PHD_FILE'}     = chomp($phd_file) if (defined $phd_file);
                #$ace->{'CHEM'}         = chomp($chemistry) if (defined $chemistry);
                #$ace->{'DYE'}          = chomp($dye) if (defined $chemistry);
                #$ace->{'TIME'}         = chomp($date_time_phd)  if (defined $date_time_phd);
                #$ace->{'TEMPLATE'}     = chomp($template) if (defined $template);
                #$ace->{'DIRECTION'}    = chomp($direction) if (defined $direction);

            }
            if ( substr( $_, 0, 3 ) eq "WR{" ) {
                #s/^WR.*\n//;
                #my ( $read, $tag_type, $program_created, $yymmdd ) =
                #  split ( /\s+/, $_ );
                #$ace->{'TAG_TYPE'}->{$read}        = $tag_type;
                #$ace->{'PROGRAM_CREATED'}->{$read} = $program_created;
                #$ace->{'YYMMDD'}->{$read}          = $yymmdd;
                #s/^\{.*\n//;
            }
        }
    }

    # Return newly created object
    return $ace;
}

#given contig name, this method return the contig DNA string
#usage $dna = $ace -> get_ace_contigDNA(-contig => contigname);
sub get_ace_contigDNA {

    # Input
    my ( $self, %INPUT ) = @_;
    my ($contigname) = $INPUT{'-contig'};
    return $self->{'DNA'}->{$contigname};

}

#given contig name, this method return the contig quality array
#usage @qual = $ace -> get_ace_contigQual(-contig => contigname);
sub get_ace_contigQual {

    # Input
    my ( $self, %INPUT ) = @_;
    my ($contigname) = $INPUT{'-contig'};
    return @{ $self->{'QUAL'}->{$contigname} };

}

#this method return all contig names
#usage @names = $ace -> get_all_contigNames();
sub get_all_contigNames {

    # Input
    my ($self) = @_;
    return @{ $self->{'ContigName'} };
}

#--- Find the total number of contigs ---#

sub Total_Contigs {
    my ($self) = @_;
    return $self->{'CONTIGS_NUMBER'};
}

#--- Find the total number of Reads ---#
sub Total_Reads {
    my ($self) = @_;
    return $self->{'READS_NUMBER'};
}

#--- Find the total bases in a contig ---#
sub Total_Bases {
    my ( $self, %INPUT ) = @_;
    my ($contigname) = $INPUT{'-contig'};
    return $self->{'ContigPadLength'}{$contigname};
}

#--- Find the Orientation of the contig ---#
sub Contig_Orientation {
    my ( $self, %INPUT ) = @_;
    my ($contigname) = $INPUT{'-contig'};
    return $self->{'ContigOrient'}->{$contigname};
}

#--- Find the length of a contig ---#
sub Contig_Length {
    my ( $self, %INPUT ) = @_;
    my ($contigname) = $INPUT{'-contig'};
    return $self->{'ContigLength'}->{$contigname};
}

sub Read_Orientation {
    my ( $self, %INPUT ) = @_;
    my ($readname) = $INPUT{'-read'};
    return $self->{'ReadOri'}->{$readname};
}

sub Read_StartPos {
    my ( $self, %INPUT ) = @_;
    my ($readname) = $INPUT{'-read'};
    return $self->{'ReadStartPos'}->{$readname};
}

sub Read_EndPos {
    my ( $self, %INPUT ) = @_;
    my ($readname) = $INPUT{'-read'};
    return $self->{'ReadEndPos'}->{$readname};
}

#--- Find all the Reads in a contig ---#
sub ReadsInContig {
    my ( $self, %INPUT ) = @_;
    my ($contigname) = $INPUT{'-contig'};
    return @{ $self->{'ReadsBelongTo'}->{$contigname} };
}


sub Number_ReadsInContig {
    my ( $self, %INPUT ) = @_;
    my ($contigname) = $INPUT{'-contig'};
    return $self->{'NumReadsInContig'}->{$contigname};
}


sub FindReadType {
    my ($self,%INPUT) = @_;
    my ($readname) = $INPUT{'-read'};
    return $self->{'ReadType'}->{$readname} ;
}

sub getUnpadedStart{
my ($self,%INPUT) = @_;
my ($readname) = $INPUT{'-read'};
return $self->{'UnPaddedReadStartPos'}->{$readname};
}

sub getUnpadedEnd{
my ($self,%INPUT) = @_;
my ($readname) = $INPUT{'-read'};
return $self->{'UnPaddedReadEndPos'}->{$readname};
}

sub getAlignUnpadedClipStart {
    my ($self,%INPUT) = @_;
    my ($readname) = $INPUT{'-read'};
    return $self->{'AlignUnpadedClipStart'}->{$readname} ;
}

sub getAlignUnpadedClipEnd {
    my ($self,%INPUT) = @_;
    my ($readname) = $INPUT{'-read'};
    return $self->{'AlignUnpadedClipEnd'}->{$readname} ;
}

sub Find_UnComplimentryReads {
    my ( $self, %INPUT ) = @_;
    my @UnComplimentryReads ;
    my ($contigname) = $INPUT{'-contig'};
    foreach my $read (@{$self->{'ReadsBelongTo'}->{$contigname}}) {
	if ($self->{'ReadOri'}->{$read} eq "U") {
	    push (@UnComplimentryReads,$read);
	}
    }
    return @UnComplimentryReads ;
}

sub find_complimentryreads {
    my ( $self, %INPUT ) = @_;
    my @ComplimentryReads ;
    my ($contigname) = $INPUT{'-contig'};
    foreach my $read (@{$self->{'ReadsBelongTo'}->{$contigname}}) {
	if ($self->{'ReadOri'}->{$read} eq "C") {
	    push (@ComplimentryReads,$read);
	}
    }
    return @ComplimentryReads ;
}





#sub parse_ace_reads {
#
#}
#sub parse_ace_readsqual {
#
#}

