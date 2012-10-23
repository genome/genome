package PostData;
#----------------------------------------------------------#
# Copyright (C) 1997 Washington University, St. Louis, MO. #
# All Rights Reserved.                                     #
#                                                          #
# Author: Ian Korf                                         #
# Send all comments to ikorf@sapiens.wustl.edu             #
#                                                          #
# DISCLAIMER: THIS SOFTWARE IS PROVIDED "AS IS"            #
#             WITHOUT WARRANTY OF ANY KIND.                #
#----------------------------------------------------------#
require Exporter;
@PostData::ISA = qw(Exporter);
@EXPORT = qw(PostData);

#----------------------------------------------------------------------
=head1 NAME PostData

=cut

#----------------------------------------------
# PostData - output simple data structures
#
# Prints out a data structure in regular format
# Only follows ARRAY, HASH, and SCALAR
# Does not terminate recursive data structures
#----------------------------------------------
sub PostData {
	my ($data,$max,$level) = @_;
	$level = 1 unless defined $level;
	$max = 100 unless defined $max;
	my ($i,$key,$value);
	my $tab = "  " x $level;

	my $ref = ref($data);
	if (not $ref) {
		if (defined $data) {print "$data\n"}
		else {print "UNDEFINED VALUE\n"}
	}
	elsif ($ref eq 'SCALAR') {
		print "$ref\n";
		return if $max == $level;
		print $tab,"$$data\n";
	}
	elsif ($ref eq 'ARRAY') {
		print "$ref\n";
		$level++;
		return if $max == $level;
		for($i=0;$i<@$data;$i++) {
			print $tab,"[$i] = ";
			PostData($data->[$i],$max,$level);
		}
	}
	elsif ($ref eq 'REF') {
		print "$ref\n";
		$level++;
		return if $max == $level;
		PostData($$data,$max,$level);
	}
	elsif ($ref eq 'HASH') {
		print "$ref\n";
		$level++;
		return if $max == $level;
		foreach $key (sort mySort keys %$data) {
			print $tab,"{$key} = ";
			PostData($data->{$key},$max,$level)
		}
	}
	elsif ($ref ne 'CODE') {
		print "$ref\n";
		$level++;
		return if $max == $level;
		foreach $key (sort mySort keys %$data) {
			print $tab,"{$key} = ";
			PostData($data->{$key},$max,$level)
		}
	}
	else {
		print "$ref not handled by PostData\n";
	}
}

sub mySort {
	if ($a =~ /^\d+$/ and $b =~ /^\d+$/) {$a <=> $b}
	else                             {$a cmp $b}
}

1;
