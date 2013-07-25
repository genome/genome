# DBD::Oracle Oraperl emulation. This file is not relevant to the
# emulation but is included for completeness only.
# I have updated %ora_types in case it's used. Tim Bunce.

# oraperl.ph
#
# Various constants which may be useful in oraperl programs
#
# Author:	Kevin Stock
# Date:		28th October 1991
# Last Change:	8th April 1992


# Oraperl error codes, set in $ora_errno

$ORAP_NOMEM	= 100001;	# out of memory
$ORAP_INVCSR	= 100002;	# invalid cursor supplied
$ORAP_INVLDA	= 100003;	# invalid lda supplied
$ORAP_NOSID	= 100004;	# couldn't set ORACLE_SID
$ORAP_BADVAR	= 100005;	# bad colon variable sequence
$ORAP_NUMVARS	= 100006;	# wrong number of colon variables
$ORAP_NODATA	= 100007;	# statement does not return data


# Oraperl debugging codes for $ora_debug
# From version 2, you shouldn't really use these.

$ODBG_EXEC	=   8;		# program execution
$ODBG_STRNUM	=  32;		# string/numeric conversions
$ODBG_MALLOC	= 128;		# memory allocation/release

# Oracle datatypes
# I don't know whether these are valid for all versions.

%ora_types =
(
	 1,	'character array',
	 2,	'number',
	 3,	'signed integer',
	 4,	'float',
	 7,	'packed decimal',
	 8,	'long string',
	 9,	'varchar',
	11,	'rowid',
	12,	'date',
	15,	'varraw',
	23,	'raw',
	24,	'long raw',
	96,	'char',
	106,'mlslabel',
);

1;
