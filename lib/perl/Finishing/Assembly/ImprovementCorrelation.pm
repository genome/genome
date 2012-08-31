package Finishing::Assembly::ImprovementCorrelation;

use strict;
use warnings;

use base 'Finishing::Assembly::Item';


sub contig_count{
    return shift->contigs->count;
}
1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/ImprovementCorrelation.pm $
#$Id: ImprovementCorrelation.pm 30959 2007-12-13 20:43:51Z ebelter $
