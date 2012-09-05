#############################################
### Functions for testing significance of ###
### per-gene categorized mutation rates   ###
#############################################

# Fetch command line arguments
args = commandArgs();
input_file = as.character(args[4]);
output_file = as.character(args[5]);
run_type = as.character(args[6]);
processors = as.numeric(args[7]);
skip_low_mr_genes = as.numeric(args[8]);

# See if we have the necessary packages installed to run in parallel
is.installed <- function( mypkg ) is.element( mypkg, installed.packages()[,1] );
parallel = FALSE;
if( processors > 1 & is.installed( 'doMC' ) & is.installed( 'foreach' )) {
    parallel = TRUE;
}

gethist <- function( xmax, n, p, ptype = "positive_log" ) {
    dbinom( 0:xmax, n, p ) -> ps;
    ps = ps[ps > 0];
    lastp = 1 - sum( ps );
    if( lastp > 0 ) ps = c( ps, lastp );
    if( ptype == "positive_log" ) ps = -log( ps );
    return( ps );
}

binit <- function( x, hmax, bin, dropbin = T ) {
    bs = as.integer( x / bin );
    bs[bs > hmax/bin] = hmax / bin;
    bs[is.na( bs )] = hmax / bin;
    tapply( exp(-x), as.factor( bs ), sum ) -> bs;
    bs = bs[bs>0];
    bs = -log( bs );
    if( dropbin ) bs = as.numeric( bs );
    return( bs );
}

convolute_b <- function( a, b ) {
    tt = NULL;
    for( j in b ) { tt = c( tt, ( a + j )); }
    return( tt );
}

mut_class_test <- function( x, xmax = 100, hmax = 25, bin = 0.001 ) {
    x = as.data.frame( x );
    colnames( x ) = c( "Class", "n", "x", "e" );
    x$p = NA; x$lh0 = NA; x$lh1 = NA;
    tot_muts = x[( x$Class == "Overall" ),]$x;
    tot_bps = x[( x$Class == "Overall" ),]$n;
    overall_bmr = x[( x$Class == "Overall" ),]$e;

    # Remove the row containing overall MR and BMR because we don't want it to be a tested category
    x = x[( x$Class != "Overall" ),];

    # If user wants to skip testing genes with low MRs, measure the relevant MRs of this gene
    gene_mr = 0; indel_mr = 0; indel_bmr = 0; trunc_mr = 0; trunc_bmr = 0;
    if( skip_low_mr_genes == 1 ) {
        if( tot_bps > 0 ) { gene_mr = tot_muts / tot_bps; }
        if( x[( grep( "Indels", x$Class )),]$n > 0 ) { indel_mr = x[( grep( "Indels", x$Class )),]$x / x[( grep( "Indels", x$Class )),]$n; }
        indel_bmr = x[( grep( "Indels", x$Class )),]$e;
        if( nrow( x[( grep( "Truncations", x$Class )),] ) > 0 ) {
            if( x[( grep( "Truncations", x$Class )),]$n > 0 ) { trunc_mr = x[( grep( "Truncations", x$Class )),]$x / x[( grep( "Truncations", x$Class )),]$n; }
            trunc_bmr = x[( grep( "Truncations", x$Class )),]$e;
        }
    }

    # Set pvals of 1 for genes with zero mutations, zero covered bps, or zero overall BMR
    if( tot_muts <= 0 | tot_bps <= 0 | overall_bmr <= 0 ) {
        p.fisher = 1; p.lr = 1; p.convol = 1; qc = 1;
    }
    # If user wants to skip testing genes with low MRs, give them pvals of 1
    else if( skip_low_mr_genes == 1 & gene_mr < overall_bmr & indel_mr <= indel_bmr & trunc_mr <= trunc_bmr ) {
        p.fisher = 1; p.lr = 1; p.convol = 1; qc = 1;
    }
    else {
        # Skip testing mutation categories that have zero BMR, or if this gene has #muts >= #covd bps
        x = x[( x$n > 0 & x$n > x$x & x$e > 0 ),];
        rounded_mut_cnts = round(x$x);
        for( i in 1:nrow(x) ) {
            x$p[i] = binom.test( rounded_mut_cnts[i], x$n[i], x$e[i], alternative = "greater" )$p.value;
            x$lh0[i] = dbinom( rounded_mut_cnts[i], x$n[i], x$e[i], log = T );
            x$lh1[i] = dbinom( rounded_mut_cnts[i], x$n[i], x$x[i] / x$n[i], log = T );
            ni = x$n[i]; ei = x$e[i];
            gethist( xmax, ni, ei, ptype = "positive_log" ) -> bi;
            binit( bi, hmax, bin ) -> bi;
            if( i == 1 ) { hist0 = bi; }
            if( i > 1 & i < nrow(x) ) { hist0 = convolute_b( hist0, bi ); binit( hist0, hmax, bin ) -> hist0; }
            if( i == nrow(x)) { hist0 = convolute_b( hist0, bi ); }
        }

        # Fisher combined p-value
        q = ( -2 ) * sum( log( x$p ));
        df = 2 * length( x$p );
        p.fisher = 1 - pchisq( q, df );

        # Likelihood ratio test
        q = 2 * ( sum( x$lh1 ) - sum( x$lh0 ));
        df = sum( x$lh1 != 0 );
        if( df > 0 ) p.lr = 1 - pchisq( q, df );
        if( df == 0 ) p.lr = 1;

        # Convolution test
        bx = -sum( x[,"lh0"] );
        p.convol = sum( exp( -hist0[hist0>=bx] ));
        qc = sum( exp( -hist0 ));
    }

    # Return results
    rst = list( x = cbind( x, tot_muts, p.fisher, p.lr, p.convol, qc ));
    return( rst );
}

dotest <- function( idx, mut, zgenes ) {
    step = round( length( zgenes ) / processors );
    start = step * ( idx - 1 ) + 1;
    stop = step * idx;
    if( idx == processors ) { stop = length( zgenes ); }
    tt = NULL;
    for( Gene in zgenes[start:stop] ) {
        mutgi = mut[mut$Gene==Gene,];
        mut_class_test( mutgi[,2:5], hmax = 25, bin = 0.001 ) -> z;
        tt = rbind( tt, cbind( Gene, unique( z$x[,(9:11)] )));
    }
    return( tt );
}

combineresults <- function( a, b ) {
    return( rbind( a, b ));
}

smg_test <- function( gene_mr_file, pval_file ) {
    read.delim( gene_mr_file ) -> mut;
    colnames( mut ) = c( "Gene", "Class", "Bases", "Mutations", "BMR" );
    mut$BMR = as.numeric( as.character( mut$BMR ));
    tt = NULL;

    # Run in parallel if we have the needed packages, or fall back to the old way
    if( parallel ) {
        library( 'doMC' );
        library( 'foreach' );
        registerDoMC();
        cat( "Parallel backend installed - splitting across", processors, "cores\n" );

        options( cores = processors );
        mcoptions <- list( preschedule = TRUE );
        zgenes = unique( as.character( mut$Gene ));
        tt = foreach( idx = 1:processors, .combine="combineresults", .options.multicore = mcoptions ) %dopar% {
            dotest( idx, mut, zgenes );
        }
        write.table( tt, file = pval_file, quote = FALSE, row.names = F, sep = "\t" );
    }
    else {
        for( Gene in unique( as.character( mut$Gene ))) {
            mutgi = mut[mut$Gene==Gene,];
            mut_class_test( mutgi[,2:5], hmax = 25, bin = 0.001 ) -> z;
            tt = rbind( tt, cbind( Gene, unique( z$x[,(9:11)] )));
        }
        write.table( tt, file = pval_file, quote = FALSE, row.names = F, sep = "\t" );
    }
}

smg_fdr <- function( pval_file, fdr_file ) {
    read.table( pval_file, header = T, sep = "\t" ) -> x;

    #Calculate FDR measure and write FDR output
    p.adjust( x[,2], method="BH" ) -> fdr.fisher;
    p.adjust( x[,3], method="BH" ) -> fdr.lr;
    p.adjust( x[,4], method="BH" ) -> fdr.convol;
    x = cbind( x, fdr.fisher, fdr.lr, fdr.convol );
    #Rank SMGs starting with lowest convolution test FDR, and then by Likelihood Ratio FDR
    x = x[order( fdr.convol, fdr.lr ),];
    write.table( x, file = fdr_file, quote = FALSE, row.names = F, sep = "\t" );
}

# Figure out which function needs to be invoked and call it
if( run_type == "smg_test" ) { smg_test( input_file, output_file ); }
if( run_type == "calc_fdr" ) { smg_fdr( input_file, output_file ); }
