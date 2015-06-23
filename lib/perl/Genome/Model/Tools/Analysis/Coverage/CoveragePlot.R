
# loads add-on packages
#library(graphics);
library(grid);


cs_scale <- function(data, breaks=5, range=NULL, fraction=0.0, order=FALSE)
{
    # if there is no data points
    if (length(data) == 0) stop("No data for scaling");

    # creates a scale vector
    if(! is.null(range))
    {
        # extends a numerical range by a small percentage, i.e., fraction, on both sides.
        data <- extendrange(data, r=range(range, na.rm=T), f=fraction);
    }
    else
    {
        # extends a numerical range by a small percentage, i.e., fraction, on both sides.
        data <- extendrange(data, f=fraction);
    }
    
    # gets pretty breakpoints
    p <- pretty(data, n=breaks);
    
    if (order)
    {
        # gets the data scope
        s <- cs_range(data);
    
        if (s[3]) return(rev(p))
        else return(p);
    }
    else
    {
        return(p);
    }
}



cs_scale_scope <- function(scale, factor=0.05)
{
    # returns a vector containing the minimum and maximum of all the given arguments in the given order
    s <- cs_range(scale);
    
    # extends a numerical range by a small percentage, i.e., fraction, on both sides
    r <- extendrange(scale, f=factor);
    
    if (s[3]) return(c(r[2], r[1]))
    else return(c(r[1], r[2]));
}



cs_range <- function(x)
{
    # gets the minimum and maximum
    r <- range(x);
    
    if (x[1] <= x[length(x)]) reverse <- FALSE
    else reverse <- TRUE;
    
    if (reverse) return(c(r[2], r[1], reverse))
    else return(c(r[1], r[2], reverse));
}



cs_range_limit_xy <- function(x, y, xlim=NULL, ylim=NULL)
{
    if (is.null(xlim)) stop("xlim not defined");
    
    # checks the data region
    minx = min(x);
    maxx = max(x);
    miny = min(y);
    maxy = max(y);
    
    # if all x and y data points are within the given range
    if (xlim[1] <= minx & maxx <= xlim[2])
    {
        if (is.null(ylim))
        {
            return(list(x=x, y=y));
        }
        else if (ylim[1] <= miny & maxy <= ylim[2])
        {
            return(list(x=x, y=y));
        }
    }
    
    # selects the data points that are within the given range
    data <- list(x=vector(mode="numeric"), y=vector(mode="numeric"));
    
    # adds into a new dataset
    idx <- 1;
    for(i in 1:length(x))
    {
        if (is.null(ylim))
        {
            if (xlim[1] <= x[i] & x[i] <= xlim[2])
            {
                # a new dataset
                data$x[idx] <- x[i];
                data$y[idx] <- y[i];
                
                # increases the index number
                idx <- idx + 1;
            }
        }
        else
        {
            if (xlim[1] <= x[i] & x[i] <= xlim[2] & ylim[1] <= y[i] & y[i] <= ylim[2])
            {
                # a new dataset
                data$x[idx] <- x[i];
                data$y[idx] <- y[i];
                
                # increases the index number
                idx <- idx + 1;
            }
        }
    }
    
    return(data);
}


cs_vmargin <- function(x, visible=T, rot=0, label=FALSE, default=unit(1.5, "lines"))
{
    margin <- default;
    
    # computes the histogram of x to infer the best height of the marks on x-axis
    for (i in 1:length(x))
    {
        if (i == 1)
        {
            max = toString(x[i]);
        }
        else
        {
            # converts a scale value into a string
            xstr <- toString(x[i]);
            
            if (nchar(xstr) > nchar(max))
            {
                # to determine the maximum character length among scale elements
                max = xstr;
            }
        }
    }
    
    # determines the margin between the borders of plot viewport and parent viewport
    if (visible)
    {
        if (rot == 0)
        {
            #margin <- margin + unit(1, "strheight", max);
            margin <- margin + unit(1, "lines");
        }
        else if (rot == 90)
        {
            margin <- margin + unit(1, "strwidth", max);
        }
        else
        {
            stop("rot is either 0 or 90");
        }
    }
    
    # determines the margin where the X label is placed
    if (label) margin <- margin + unit(3, "strheight", "O");
    
    return (margin);
}



cs_hmargin <- function(y, visible=T, rot=0, label=FALSE, default=unit(1.5, "lines"))
{
    margin <- default;
    
    # computes the histogram of y to infer the best height of the marks on y-axis
    for (i in 1:length(y))
    {
        if (i == 1)
        {
            max = toString(y[i]);
        }
        else
        {
            # converts a scale value into a string
            ystr <- toString(y[i]);
            
            if (nchar(ystr) > nchar(max))
            {
                # to determine the maximum character length among scale elements
                max = ystr;
            }
        }
    }
    
    # determines the margin between the borders of plot viewport and parent viewport
    if (visible)
    {
        if (rot == 0)
        {
            margin <- margin + unit(1, "strwidth", max);
        }
        else if (rot == 90)
        {
            #margin <- margin + unit(1, "strheight", max);
            margin <- margin + unit(1, "lines");
        }
        else
        {
            stop("rot is either 0 or 90");
        }
    }
    
    # determines the margin where the Y label is placed
    if (label) margin <- margin + unit(3, "strheight", "O");
    
    return (margin);
}



cs_viewport2 <- function(left, bottom, xscale, xfactor=0.05, yscale, yfactor=0.05, lines=0, default=unit(1.5, "lines"), name=NULL)
{
    # sets up the width and height considering a space for a title
    width = unit(1.0, "npc") - left - default;
    height = unit(1.0, "npc") - bottom - default - unit(lines, "lines");
    
    # creates a viewport for a plot region
    vp <- viewport(
        x = left,
        y = bottom,
        width = width,
        height = height,
        xscale = cs_scale_scope(xscale, factor=xfactor),
        yscale = cs_scale_scope(yscale, factor=yfactor),
        just = c("left", "bottom"),
        #clip = TRUE,

        name = name
    );
    
    return(vp);
}



cs_base <- function(title=NULL, lines=2, xlab=NULL, xmargin=2, ylab=NULL, ymargin=2, xscale=NULL, xfactor=0.05, xvisible=T, xrot=0, xtick=-1, xspace=-1, yscale=NULL, yfactor=0.05, yvisible=T, yrot=0, ytick=-1, yspace=-1, default = unit(1.5, "lines"), gp=gpar(), name=NULL, parent=NULL, vp=NULL, xmar = NULL, ymar = NULL)
{
    # the chracter identifiers for the grid components
    if (is.null(name)) nv <- "cs_base" else nv = name;        # the grid viewport
    nt = "title";  # main title
    nx = "xaxis";   # x axis
    ny = "yaxis";   # y axis
    nxl = "xlabel"; # x label
    nyl = "ylabel"; # x label
    nb = "border";   # border
    
    # checks the x and y scale for a plot
    if(is.null(xscale)) stop("xscale not defined");
    if(! is.vector(xscale, mode="numeric")) stop("wrong data type for xscale");
    if(is.null(yscale)) stop("yscale not defined");
    if(! is.vector(yscale, mode="numeric")) stop("wrong data type for yscale");
    
    # creates and sets up a viewport for a plot
    if (is.null(vp))
    {
        # sets up the margin of a plot region as a horizontal space from the left
        xmar <- cs_hmargin(yscale, visible=yvisible, rot=yrot, label=(! is.null(ylab)), default=default);
        
        # sets up the margin of a plot region as a vertical space from the bottom
        ymar <- cs_vmargin(xscale, visible=xvisible, rot=xrot, label=(! is.null(xlab)), default=default);
        
        # the lines for a title
        if (is.null(title)) titlelines <- 0
        else                titlelines <- lines;
        
        # creates a viewport
        # old line: cv <- cs_viewport2(left=xmar, bottom=ymar, xscale=cs_range(xscale), xfactor=xfactor, yscale=cs_range(yscale), yfactor=yfactor, lines=titlelines, default=default, name=nv);
        cv <- cs_viewport2(left=xmar, bottom=ymar, xscale=xscale, xfactor=xfactor, yscale=yscale, yfactor=yfactor, lines=titlelines, default=default, name=nv);
    }
    else
    {
        # uses the given viewport
        cv <- vp;
    }


    # draws the title
    if (is.null(title)) gt <- NULL
    else                gt <- textGrob(title, y=unit(1, "npc") + unit(lines, "lines"), name=nt, vp=cv);
    
    # draws the x axis on the plot
    gx <- xaxisGrob(at=xscale, name=nx, label=xvisible, vp=cv);
    gx$children$ticks$y1 <- unit(xtick / 2, "lines");
    if (xrot == 90)
    {
        gx <- editGrob(gx, gPath = "labels", rot=xrot, y=unit(abs(xtick)*xspace, "lines"), just=c("right", "center"));
    }
    else
    {
        gx <- editGrob(gx, gPath = "labels", y=unit(abs(xtick)*xspace, "lines"), just=c("center", "top"));
    }
    
    # draws the x axis on the plot
    gy <- yaxisGrob(at=yscale, name=ny, label=yvisible, vp=cv);
    gy$children$ticks$x1 <- unit(ytick / 2, "lines");
    if (yrot == 90)
    {
        gy <- editGrob(gy, gPath = "labels", rot=90, x=unit(abs(ytick)*yspace, "lines"), just=c("center", "bottom"));
    }
    else
    {
        gy <- editGrob(gy, gPath = "labels", x=unit(abs(ytick)*yspace, "lines"));
    }
    
    # draws the x label
    if (is.null(xlab)) gxl <- NULL
    else        gxl <- textGrob(xlab, y=unit(xmargin, "strheight", xlab) - ymar, name=nxl, vp=cv);

    # draws the y label
    if (is.null(ylab)) gyl <- NULL
    else        gyl <- textGrob(ylab, x=unit(ymargin, "strheight", ylab) - xmar, rot=90, name=nyl, vp=cv);
    
    # draws the border line
    gb <- rectGrob(name=nb, vp=cv);
    

    # creates a grid graphical objects
    g <- gTree(
         name=paste("g", nv, sep=""), 
         children = gList(gt, gxl, gyl, gx, gy, gb),
         childrenvp = cv,   # the current viewport into the gTree
         gp = gp,
         vp = parent
        );    # cs_base

    return(g);
}



grid.cs_lines <- function(x, y, p=c(0,0.2,0.4,0.6,0.8,1.0), pch=1, size=unit(1,"char"), pcol=0, default.units="npc", arrow=NULL, gp=gpar(), draw=TRUE, name=NULL, vp=NULL)
{
    # TODO: declaration of class identifier and grid components
    # sets with a string giving the class attribute for the list.struct
    cl <- "cs_lines";
    
    # lists up the chracter identifiers for the grid components
    nl = "lines";    # the lines object
    np = "points";     # the points object
    
    
    # checks the input data
    if (missing(x)) stop("Missing x");
    if (missing(y)) stop("Missing y");
    if (missing(p)) stop("Missing p");
    if (length(x) != length(y)) stop("x and y numeric vectors with different size");
    if (length(x) < length(p)) stop("p numeric vectors with larger size");
    if (length(p) == 0) stop("p required");
    
    # TODO: data custimization (optional)
    # calculates the index number
    i <- as.integer(p * (length(x) - 1)) + 1;
    
    
    # TODO: creation of a child gTree object 
    # creates an area gTree
    g <- gTree(
        name=name, 
        children = gList(),
        childrenvp = vp,
        gp = gp,
        vp = NULL,
        cl = cl);
    
    # the default area option
    slot.names <- names(gp);
    if ("col" %in% slot.names) gp$col <- NULL;
    
    # TODO: adds grid objects as components
    # adds the lines
    g <- addGrob(g, linesGrob(x, y, default.units=default.units, arrow=arrow, name=nl, gp=gp, vp=vp));
    
    # adds the points
    g <- addGrob(g, pointsGrob(x[i], y[i], pch=pch, size=size, default.units=default.units, name=np, gp=gpar(col=pcol), vp=vp));

    if (draw) 
        grid.draw(g);
    
    invisible(g);
}


cs_linesGrob <- function(x, y, p=c(0,0.2,0.4,0.6,0.8,1.0), pch=1, size=unit(1,"char"), pcol=0, default.units="npc", arrow=NULL, gp=gpar(), name=NULL, vp=NULL)
{
    grid.cs_lines(x=x, y=y, p=p, pch=pch, size=size, pcol=pcol, default.units=default.units, arrow=arrow, gp=gp, draw=FALSE, name=name, vp=vp);
}



grid.cs_vruler <- function(scale=NULL, breaks=NULL, type=c("top", "bottom"), inner=TRUE, size=0.02, lty=NULL, lwd=1, gp=gpar(), draw=TRUE, name=NULL, vp=NULL)
{
    # TODO: declaration of class identifier and grid components
    # sets with a string giving the class attribute for the list.struct
    cl <- "cs_vruler";
    
    # lists up the chracter identifiers for the grid components
    nm = "major";    # the lines object
    ns = "minor";    # the lines object
    nl = "lines";    # the lines object
    
    # checks the input data
    if (missing(scale)) stop("scale argument required");
    #if (missing(breaks)) stop("breaks argument required");
    if (length(scale) < 2) stop("Invalid scale");
    
    if (is.null(breaks))
    {
        breaks <- length(pretty(c(scale[1], scale[2]))) - 1;
    }
    
    # the size factor for major and minor marks
    msize = size;      # msize = 0.03;
    ssize = size / 2;  # ssize = 0.015;
    
    # the type of the vertical ruler
    top = FALSE;
    bottom = FALSE;
    if (length(type) == 2)
    {
        if (type[1] == "bottom" || type[2] == "bottom") bottom = TRUE;
        if (type[1] == "top" || type[2] == "top") top = TRUE;
    }
    else
    {
        if (type == "top") top = TRUE;
        if (type == "bottom") bottom = TRUE;
    }
    
    
    # TODO: data custimization (optional)
    
    # the default area option
    #slot.names <- names(gp);
    #if (! ("col" %in% slot.names)) gp$col <- "grey";
    
    # TODO: adds grid objects as components
    # creates a child gTree for major marks
    gcm <- gTree(
         name=nm, 
         children = gList(),
         childrenvp = vp,
         gp = gp,
         vp = NULL,
         cl = NULL);

    # creates a child gTree for minor marks
    gcs <- gTree(
         name=ns, 
         children = gList(),
         childrenvp = vp,
         gp = gp,
         vp = NULL,
         cl = NULL);

    # creates a child gTree for gridlines
    gcl <- gTree(
         name=nl, 
         children = gList(),
         childrenvp = vp,
         gp = gp,
         vp = NULL,
         cl = NULL);
    
    # adds major and minor marks
    cur <- min(scale);
    if (breaks > 0)
        width <- abs(scale[1] - scale[2]) / breaks
    else
        width <- abs(scale[1] - scale[2]);
    
    i = 0;
    while(cur <= max(scale))
    {
        # create a major and minor mark
        if (i == 0)
        {
            # as a major mark
            if (top)
            {
                if(inner) m <- linesGrob(x=unit(c(cur, cur), "native"), y=unit(c(1 - msize, 1), "npc"), name=paste("t", cur, sep=""), gp=gp, vp=vp)
                else m <- linesGrob(x=unit(c(cur, cur), "native"), y=unit(c(1, 1 + msize), "npc"), name=paste("t", cur, sep=""), gp=gp, vp=vp)
                
                gcm <- addGrob(gcm, m);
            }
            
            if (bottom)
            {
                if(inner) m <- linesGrob(x=unit(c(cur, cur), "native"), y=unit(c(0, msize), "npc"), name=paste("b", cur, sep=""), gp=gp, vp=vp)
                else m <- linesGrob(x=unit(c(cur, cur), "native"), y=unit(c(-msize, 0), "npc"), name=paste("b", cur, sep=""), gp=gp, vp=vp);
                
                gcm <- addGrob(gcm, m);
            }
        }
        else
        {
            # as a minor mark
            if (top)
            {
                if(inner) m <- linesGrob(x=unit(c(cur, cur), "native"), y=unit(c(1 - ssize, 1), "npc"), name=paste("t", cur, sep=""), gp=gp, vp=vp)
                else m <- linesGrob(x=unit(c(cur, cur), "native"), y=unit(c(1, 1 + ssize), "npc"), name=paste("t", cur, sep=""), gp=gp, vp=vp)
                
                gcs <- addGrob(gcs, m);
            }
            
            if (bottom)
            {
                if(inner) m <- linesGrob(x=unit(c(cur, cur), "native"), y=unit(c(0, ssize), "npc"), name=paste("b", cur, sep=""), gp=gp, vp=vp)
                else m <- linesGrob(x=unit(c(cur, cur), "native"), y=unit(c(-ssize, 0), "npc"), name=paste("b", cur, sep=""), gp=gp, vp=vp);
                
                gcs <- addGrob(gcs, m);
            }
        }
        
        # creates a gridline
        if (! is.null(lty))
        {
            l <- linesGrob(x=unit(c(cur, cur), "native"), y=unit(c(0, 1), "npc"), name=paste("l", cur, sep=""), gp=gpar(lty=lty, lwd=lwd), vp=vp);
            gcl <- addGrob(gcl, l);
        }
        
        # moves the mark position
        cur <- cur + width;
        i <- i + 1;
        
        if (i == breaks) i <- 0;
    }

    
    # TODO: creation of a child gTree object 
    # creates an area gTree
    g <- gTree(
         name=name, 
         children = gList(gcl, gcs, gcm),
         childrenvp = vp,
         gp = NULL,
         vp = NULL,
         cl = cl);
    
    if (draw) 
        grid.draw(g);
    
    invisible(g);
}


cs_vrulerGrob <- function(scale=NULL, breaks=NULL, type=c("top", "bottom"), inner=TRUE, size=0.02, lty=NULL, lwd=1, gp=gpar(), name=NULL, vp=NULL)
{
    grid.cs_vruler(scale=scale, breaks=breaks, type=type, inner=inner, size=size, lty=lty, lwd=lwd, gp=gp, draw=FALSE, name=name, vp=vp);
}



grid.cs_hruler <- function(scale=NULL, breaks=NULL, type=c("left", "right"), inner=TRUE, size=0.02, lty=NULL, lwd=1, gp=gpar(), draw=TRUE, name=NULL, vp=NULL)
{
    # TODO: declaration of class identifier and grid components
    # sets with a string giving the class attribute for the list.struct
    cl <- "cs_hruler";
    
    # lists up the chracter identifiers for the grid components
    nm = "major";    # the lines object
    ns = "minor";    # the lines object
    nl = "lines";    # the lines object
    
    # checks the input data
    if (missing(scale)) stop("scale argument required");
    #if (missing(breaks)) stop("breaks argument required");
    if (length(scale) < 2) stop("Invalid scale");
    
    if (is.null(breaks))
    {
        breaks <- length(pretty(c(scale[1], scale[2]))) - 1;
    }
    
    # the size factor for major and minor marks
    msize = size;      # msize = 0.03;
    ssize = size / 2;  # ssize = 0.015;
    
    # the type of the vertical ruler
    right = FALSE;
    left = FALSE;
    if (length(type) == 2)
    {
        if (type[1] == "left" || type[2] == "left") left = TRUE;
        if (type[1] == "right" || type[2] == "right") right = TRUE;
    }
    else
    {
        if (type == "right") right = TRUE;
        if (type == "left") left = TRUE;
    }
    
    
    # TODO: data custimization (optional)
    
    # the default area option
    #slot.names <- names(gp);
    #if (! ("col" %in% slot.names)) gp$col <- "grey";
    
    # TODO: adds grid objects as components
    # creates a child gTree for major marks
    gcm <- gTree(
         name=nm, 
         children = gList(),
         childrenvp = vp,
         gp = gp,
         vp = NULL,
         cl = NULL);

    # creates a child gTree for minor marks
    gcs <- gTree(
         name=ns, 
         children = gList(),
         childrenvp = vp,
         gp = gp,
         vp = NULL,
         cl = NULL);

    # creates a child gTree for gridlines
    gcl <- gTree(
         name=nl, 
         children = gList(),
         childrenvp = vp,
         gp = gp,
         vp = NULL,
         cl = NULL);
    
    # adds major and minor marks
    cur <- min(scale);
    if (breaks > 0)
        width <- abs(scale[1] - scale[2]) / breaks
    else
        width <- abs(scale[1] - scale[2]);
        
    i = 0;
    while(cur <= max(scale))
    {
        # create a major and minor mark
        if (i == 0)
        {
            # as a major mark
            if (right)
            {
                if(inner) m <- linesGrob(x=unit(c(1 - msize, 1), "npc"), y=unit(c(cur, cur), "native"), name=paste("t", cur, sep=""), gp=gp, vp=vp)
                else m <- linesGrob(x=unit(c(1, 1 + msize), "npc"), y=unit(c(cur, cur), "native"), name=paste("t", cur, sep=""), gp=gp, vp=vp);
                
                gcm <- addGrob(gcm, m);
            }
            
            if (left)
            {
                if(inner) m <- linesGrob(x=unit(c(0, msize), "npc"), y=unit(c(cur, cur), "native"), name=paste("b", cur, sep=""), gp=gp, vp=vp)
                else m <- linesGrob(x=unit(c(-msize, 0), "npc"), y=unit(c(cur, cur), "native"), name=paste("b", cur, sep=""), gp=gp, vp=vp);
                
                gcm <- addGrob(gcm, m);
            }
        }
        else
        {
            # as a minor mark
            if (right)
            {
                if(inner) m <- linesGrob(x=unit(c(1 - ssize, 1), "npc"), y=unit(c(cur, cur), "native"), name=paste("t", cur, sep=""), gp=gp, vp=vp)
                else m <- linesGrob(x=unit(c(1, 1 + ssize), "npc"), y=unit(c(cur, cur), "native"), name=paste("t", cur, sep=""), gp=gp, vp=vp);
                
                gcs <- addGrob(gcs, m);
            }
            
            if (left)
            {
                if(inner) m <- linesGrob(x=unit(c(0, ssize), "npc"), y=unit(c(cur, cur), "native"), name=paste("b", cur, sep=""), gp=gp, vp=vp)
                else m <- linesGrob(x=unit(c(-ssize, 0), "npc"), y=unit(c(cur, cur), "native"), name=paste("b", cur, sep=""), gp=gp, vp=vp);
                
                gcs <- addGrob(gcs, m);
            }
        }
        
        # creates a gridline
        if (! is.null(lty))
        {
            l <- linesGrob(x=unit(c(0, 1), "npc"), y=unit(c(cur, cur), "native"), name=paste("l", cur, sep=""), gp=gpar(lty=lty, lwd=lwd), vp=vp);
            gcl <- addGrob(gcl, l);
        }
        
        # moves the mark position
        cur <- cur + width;
        i <- i + 1;
        
        if (i == breaks) i <- 0;
    }

    
    # TODO: creation of a child gTree object 
    # creates an area gTree
    g <- gTree(
         name=name, 
         children = gList(gcl, gcs, gcm),
         childrenvp = vp,
         gp = NULL,
         vp = NULL,
         cl = cl);
    
    if (draw) 
        grid.draw(g);
    
    invisible(g);
}


cs_hrulerGrob <- function(scale=NULL, breaks=NULL, type=c("left", "right"), inner=TRUE, size=0.02, lty=NULL, lwd=1, gp=gpar(), name=NULL, vp=NULL)
{
    grid.cs_hruler(scale=scale, breaks=breaks, type=type, inner=inner, size=size, lty=lty, lwd=lwd, gp=gp, draw=FALSE, name=name, vp=vp);
}



grid.cs_colorLegend <- function(x, y, width=3, height=unit(0.5, "char"), title=NULL, cols=NULL, labels=NULL, gp=gpar(), draw=TRUE, name=NULL, vp=NULL)
{
    # TODO: declaration of class identifier and grid components
    # sets with a string giving the class attribute for the list.struct
    cl <- "cs_colorLegend";
    
    # lists up the chracter identifiers for the grid components
    nt = "title";
    nc = "legend";
    nl = "labels";
    
    # checks the input data
    if (is.null(cols)) stop("cols is not specified for color");
    if (is.null(labels)) stop("labels is not specified for label");
    if (length(cols) != length(labels)) stop("cols and labels are vectors of different size");
        
    # TODO: data custimization (optional)
    # makes a numeric vector specifying x and y locations to draw a polygon
    
    # the default area option
    if (is.null(x)) x <- unit(0.5, "npc")    else x <- unit(x, "native");
    if (is.null(y)) y <- unit(0.5, "npc")    else y <- unit(y, "native");
    if (is.null(width)) w <- unit(3, "char") else w <- unit(width, "char");
    
    # TODO: creation of a child gTree object 
    # for the color legend and labels
    gc <- gTree(
         name=nc, 
         children = gList(),
         childrenvp = vp,
         gp = gp,
         vp = NULL,
         cl = NULL);

    gl <- gTree(
         name=nl, 
         children = gList(),
         childrenvp = vp,
         gp = gp,
         vp = NULL,
         cl = NULL);
    
    # creates each color box
    xc <- x;
    for (i in 1:length(cols))
    {
        r <- rectGrob(x=xc, y=y, width=w, height=height, just=c("left"), gp=gpar(col="white", fill=cols[i]), name=paste("color", i, sep=""), vp=vp);
        l <- textGrob(labels[i], x=xc+unit(width/2, "char"), y=y-height, just="top", name=paste("label", i, sep=""), vp=vp); # gp=gpar(cex=0.8)
        xc <- xc + w;
                
        # adds to the child gTree
        gc <- addGrob(gc, r);
        gl <- addGrob(gl, l);
    }
    
    # creates a border    
    r <- rectGrob(x=x, y=y, width=w*length(cols), height=height, just=c("left"), name="border", vp=vp)
    gc <- addGrob(gc, r);

    # TODO: adds grid objects as components
    gt <- textGrob(title, x=x+w*(length(cols)/2), y=y+height, just="bottom", name=nt, vp=vp);

    
    # TODO: creation of a child gTree object 
    # creates an area gTree
    g <- gTree(
         name=name, 
         children = gList(gc, gt, gl),
         childrenvp = vp,
         gp = NULL,
         vp = NULL,
         cl = cl);
    
    if (draw) 
        grid.draw(g);
    
    invisible(g);
}


cs_colorLegendGrob <- function(x, y, width=unit(0.5, "char"), height=unit(3, "char"), title=NULL, cols=NULL, labels=NULL, gp=gpar(), name=NULL, vp=NULL)
{
    grid.cs_colorLegend(x=x, y=y, width=width, height=height, title=title, cols=cols, labels=labels, gp=gp, draw=FALSE, name=name, vp=vp);
}



grid.cs_vlegend <- function(just="center", x=unit(0.5, "npc"), y=unit(0.5, "npc"), pch=NULL, lty=NULL, col, labels, fontsize=12, cex=0.5, lwd=1, border=TRUE, gp=gpar(), draw=TRUE, name=NULL, vp=NULL)
{
    # TODO: declaration of class identifier and grid components
    # sets with a string giving the class attribute for the list.struct
    cl <- "cs_vlegend";
    
    # lists up the chracter identifiers for the grid components
    nv = "vlegend";
    np = "point";    # the points object
    nl = "line";    # the lines object
    nt = "text";    # the text object
    nb = "border";  # the rec object

    
    # checks the input data
    if (missing(col)) stop("Missing col");
    if (! is.null(pch) && length(pch) != length(col)) stop("pch and col with different size");
    if (! is.null(lty) && length(lty) != length(col)) stop("lty and col with different size");
    

    # TODO: data custimization (optional)
    # makes a numeric vector specifying x and y locations to draw a polygon
    
    # the default area option
    labels <- as.character(labels);
    nkeys <- length(labels);
    
    # TODO: creates viewports
    # creates a Grid layout
    legend.layout <- grid.layout(nkeys, 3,
                widths = unit.c(unit(2, "lines"), max(unit(rep(1, nkeys), "strwidth", as.list(labels))), unit(0.5, "lines")),
                heights = unit.pmin(unit(1.05, "lines"), unit(0.7, "lines") + unit(rep(1, nkeys), "strheight", as.list(labels))),
                # upgraded from    heights = unit.pmin(unit(2, "lines"), unit(0.7, "lines") + unit(rep(1, nkeys), "strheight", as.list(labels))),
                just=just);         # just=just


    vplist <- vector("list", nkeys * 2 + 1);
    
    for (i in 1:nkeys)
    {
        for (j in 1:2)
        {
            vplist[[(i - 1) * 2 + j]] <- viewport(layout.pos.row=i, layout.pos.col=j, name=paste("C", i, j, sep=""));
        }
    }
    vplist[[nkeys * 2 + 1]] <- viewport(layout.pos.row=1:nkeys, layout.pos.col=1:3, name="CA");
                
    # creates a viewport for a legend region
    v <- viewport(
        x = x,
        y = y,
        default.units = "npc",          # default.units = "native",
        just = "center",                  # do not change this parameter
        width=unit(1.0, "npc"),
        height=unit(1.0, "npc"),
        layout = legend.layout,
        name = nv,
    );
    
    
    # creates a viewport tree
    vt <- vpTree(v, do.call("vpList", vplist));
    
    
    # TODO: adds grid objects as components
    glist <- vector("list");
    
    # creates the border line
    if (border) glist[[1]] <- rectGrob(x=unit(0.5, "npc"), y=unit(0.5, "npc"), width=unit(1.0, "npc"), height=unit(1.0, "npc"), just="centre", name=nb, gp=gpar(col="black"), vp=vpPath(nv, "CA"));
    
    # creates symbols and their labels
    for (i in 1:nkeys)
    {
        if (!is.null(pch))
            glist[[length(glist) + 1]] <- pointsGrob(0.5, 0.5, pch=pch[i], name=paste(np, i, sep=""), gp=gpar(col=col[i], cex=cex), vp=vpPath(nv, paste("C", i, 1, sep="")));
                     
        if (!is.null(lty))
            glist[[length(glist) + 1]] <- linesGrob(x=c(0.2, 0.8), y=0.5, name=paste(nl, i, sep=""), gp=gpar(col=col[i], lty=lty[i], lwd=lwd), vp=vpPath(nv, paste("C", i, 1, sep="")));
                     
        glist[[length(glist) + 1]] <- textGrob(labels[i], x=0, y=0.5, just= c("left", "centre"), name=paste(nt, i, sep=""), gp=gpar(fontsize=fontsize), vp=vpPath(nv, paste("C", i, 2, sep="")));
    }
    
    # TODO: creation of a child gTree object 
    # creates an area gTree
    g <- gTree(
         name=name, 
         children = do.call("gList", glist),
         childrenvp = vt,
         
         gp = NULL,
         vp = vp,
         cl = cl);
    
    
    if (draw) 
        grid.draw(g);
    
    invisible(g);
}


cs_vlegendGrob <- function(just="center", x=unit(0.5, "npc"), y=unit(0.5, "npc"), pch=NULL, lty=NULL, col, labels, fontsize=12, cex=0.5, lwd=1, border=TRUE, gp=gpar(), name=NULL, vp=NULL)
{
    grid.cs_vlegend(just=just, x=x, y=y, pch=pch, lty=lty, col=col, labels=labels, fontsize=fontsize, cex=cex, lwd=lwd, border=border, gp=gp, draw=FALSE, name=name, vp=vp);
}



grid.cs_boxplot <- function(d, at=0, width=NULL, outline=FALSE, horizon=FALSE, gp=gpar(), draw=TRUE, name=NULL, vp=NULL)
{
    # TODO: declaration of class identifier and grid components
    # sets with a string giving the class attribute for the list.struct
    cl <- "cs_boxplot";
    
    # lists up the chracter identifiers for the grid components
    nl = "line";    # the lines object
    nb = "hingebox";    # the rect object
    nx = "max";    # the lines object
    nm = "median";    # the lines object
    nn = "min";    # the lines object
    no = "outline";   # the child object

    
    # checks the input data
    if (missing(d)) stop("No boxplot.stats object found");
    #if (missing(at)) stop("No at found");
    if (is.null(width)) width = d$n;  # assigns the width as the number of data

    # TODO: data custimization (optional)
    # makes a numeric vector specifying x and y locations to draw a polygon
    
    # the default area option
    boxcol <- "white";
    
    # TODO: adds grid objects as components
    if (horizon)
    {
        lg <- polylineGrob(x=unit(c(d$stats[1], d$stats[5]), "native"), y=c(at, at), default.units="native", name=nl, gp=gpar(lty=5, fill=boxcol), vp=vp);
        bg <- rectGrob(x=d$stats[2], y=at, width=(d$stats[4]-d$stats[2]), height=width, just=c("left", "center"), default.units="native", name=nb, gp=gpar(fill=boxcol), vp=vp);
        ng <- linesGrob(x=c(d$stats[1], d$stats[1]), y=c(at-width/3, at+width/3), default.units="native", name=nn, gp=gp, vp=vp);
        mg <- linesGrob(x=c(d$stats[3], d$stats[3]), y=c(at-width/2, at+width/2), default.units="native", name=nm, gp=gpar(lwd=1.5), vp=vp);
        xg <- linesGrob(x=c(d$stats[5], d$stats[5]), y=c(at-width/3, at+width/3), default.units="native", name=nx, gp=gp, vp=vp);
    }
    else
    {
        lg <- polylineGrob(x=c(at, at), y=unit(c(d$stats[1], d$stats[5]), "native"), default.units="native", name=nl, gp=gpar(lty=5, fill=boxcol), vp=vp);
        bg <- rectGrob(x=at, y=d$stats[2], width=width, height=(d$stats[4]-d$stats[2]), just=c("center", "bottom"), default.units="native", name=nb, gp=gpar(fill=boxcol), vp=vp);
        ng <- linesGrob(x=c(at-width/3, at+width/3), y=c(d$stats[1], d$stats[1]), default.units="native", name=nn, gp=gp, vp=vp);
        mg <- linesGrob(x=c(at-width/2, at+width/2), y=c(d$stats[3], d$stats[3]), default.units="native", name=nm, gp=gpar(lwd=1.5), vp=vp);
        xg <- linesGrob(x=c(at-width/3, at+width/3), y=c(d$stats[5], d$stats[5]), default.units="native", name=nx, gp=gp, vp=vp);
    }
    
    # TODO: creation of a child gTree object 
    # creates an area gTree
    g <- gTree(
         name=name, 
         children = gList(lg, bg, mg, ng, xg),
         childrenvp = vp,
         gp = gp,
         vp = NULL,
         cl = cl);
    
    # draws the outliers
    if (outline)
    {
        if (horizon)
        {
            # adds outliers to the parent gTree
            og <- pointsGrob(x=d$out, y=rep(at, length(d$out)), pch=21, size = unit(0.3, "char"), default.units="native", name=no, gp=gp, vp=vp);
            g <- addGrob(g, og);
        }
        else
        {
            # adds outliers to the parent gTree
            og <- pointsGrob(x=rep(at, length(d$out)), y=d$out, pch=21, size = unit(0.3, "char"), default.units="native", name=no, gp=gp, vp=vp);
            g <- addGrob(g, og);
        }
    }
    
    if (draw) 
        grid.draw(g);
    
    invisible(g);
}


cs_boxplotGrob <- function(d, at=0, width=NULL, outline=FALSE, horizon=FALSE, gp=gpar(), name=NULL, vp=NULL)
{
    grid.cs_boxplot(d=d, at=at, width=width, outline=outline, horizon=horizon, gp=gp, draw=FALSE, name=name, vp=vp);
}



cs_plot <- function(title=NULL, xlab=NULL, ylab=NULL, xscale=NULL, xfactor=0.05, yscale=NULL, yfactor=0.05, gp=gpar(), name=NULL, parent=NULL)
{
    # lists up the chracter identifiers for the grid components
    if (is.null(name)) nv = "cs_plot" else nv <- name;
    nvr = "vruler";    # the vertical ruler
    nhr = "hruler";    # the horizontal ruler

    # the default gpar option
    names <- names(gp);
    if (! ("fontsize" %in% names)) gp$fontsize <- 8;

    # creates a xy scatter plot region and organizes components for publication purpose
    g <- cs_base(title=title, xlab=xlab, xmargin=3, ylab=ylab, ymargin=2, xscale=xscale, xfactor=xfactor, xvisible=T, xrot=0, xtick=0.5, xspace=-1, yscale=yscale, yfactor=yfactor, yvisible=T, yrot=0, ytick=0.5, yspace=-1, default=unit(1, "lines"), gp=gp, name=nv, parent=parent);
    
    # draws vertical ruler gridlines on a plot
    g <- addGrob(g, cs_vrulerGrob(scale=g$children$xaxis$at, breaks=NULL, type="bottom", inner=T, size=0.02, lty=NULL, name=nvr, vp=g$childrenvp));
    g <- addGrob(g, cs_hrulerGrob(scale=g$children$yaxis$at, breaks=NULL, type="left", inner=T, size=0.02, lty=NULL, name=nhr, vp=g$childrenvp));
    
    return(g);
}




# --------------------------------------
# gets the command line arguments when run by Rscript
args <- commandArgs(trailingOnly=TRUE);
if (length(args) != 2) stop("Two parameters required");

# sets the workding directory
setwd(args[2]);



# -------------------------------------
# starts the graphic device driver for producing PDF graphics
pdf("./coverage-plot.pdf", w=10, h=11);   # a4, a4r, letter, legal, USr, etc



# -------------------------------------
# reads a file in table format and prints out basic information
b <- read.table("./boxplot.tab", header=TRUE, sep="\t", comment.char="", stringsAsFactors=FALSE);
c <- read.table("./coverage.tab", header=TRUE, sep="\t", comment.char="", stringsAsFactors=FALSE);
d <- read.table("./density.tab", header=TRUE, sep="\t", comment.char="", stringsAsFactors=FALSE);

# gets the sample list
samples <- b$sample;

# gets the number of ROIs and total bases
nroi <- c$nroi[1];
nbase <- c$total_ref_bases[1];



# -------------------------------------
# plotting parameters

# defines titles
title1 = "Sequencing coverage depth";
title2 = "Sequencing coverage by depth";
title3 = "Distribution of the coverage depths"
title4 = "";        # blank

# a title for the x axis and the y axis
xlab1 = NULL;
ylab1 = "Average coverage depth";

xlab2 = "Minimum depth filtration";
ylab2 = "Coverage (%)";

xlab3 = "Average coverage depth of each ROI";
ylab3 = "Density";

xlab4 = NULL;
ylab4 = "Coverage (%)";###


# color and symbol set
colors = c('black', 'red', 'green3', 'blue', 'orange', 'cyan', 'magenta', 'gold');
pchars = c(0, 16,   1, 17,   2, 18,   5, 15,   4, 6);

# defines the symbol position
posx <- c(0,0.2,0.4,0.6,0.8,1.0);



# -------------------------------------
# Main code for R graphics
## TODO: divides the graphic regions
grid.newpage();
vt1 <- viewport(x=unit(0.17, "npc"), y=unit(1.0, "npc"), width=unit(0.35, "npc"), height=unit(0.25, "npc"), just="top", name="panel1");
vt2 <- viewport(x=unit(0.57, "npc"), y=unit(1.0, "npc"), width=unit(0.45, "npc"), height=unit(0.35, "npc"), just="top", name="panel2");
vt3 <- viewport(x=unit(0.25, "npc"), y=unit(0.55, "npc"), width=unit(0.5, "npc"), height=unit(0.3, "npc"), just="top", name="panel3");
vt4 <- viewport(x=unit(0.75, "npc"), y=unit(0.55, "npc"), width=unit(0.5, "npc"), height=unit(0.3, "npc"), just="top", name="panel4");



## TODO: adds plots

# for plot1

# for the title
title1 <- paste(title1, " of ", nroi, " ROIs (total ", nbase, ' bases)', sep="");

# defines the scale of a plot
max1 <- 0;
for (i in 1:length(b$stats))
    max1 <- max(max1, as.numeric(unlist(strsplit(b$stats[i], ","))));

# if with extreme distribution
if (max1 < 5) max1 <- 5;

yscale1 <- cs_scale(c(0, max1), breaks=5);

# adjusts plot parameters
xfactor1 = -0.025 * length(samples) + 0.35;
if (xfactor1 < 0.05) xfactor1 = 0.05;
width1 = 0.025 * length(samples) + 0.3;
if (width1 > 0.6) width1 = 0.6;


# draws plot 1
g1 <- cs_plot(title=title1, xlab=xlab1, ylab=ylab1, xscale=c(1:length(samples)), xfactor=xfactor1, yscale=yscale1, gp=gpar(), name=NULL, parent=vt1);

#
g1 <- addGrob(g1, cs_hrulerGrob(scale=yscale1, breaks=NULL, type=c("left", "right"), inner=TRUE, lty=3, lwd=0.5, gp=gpar(col="gray", alpha=0.7), name="hruler", vp=g1$childrenvp));
#g1 <- removeGrob(g1, gPath=gPath("hruler", "minor"));

#
g1 <- removeGrob(g1, gPath=gPath("vruler", "minor"));
g1 <- editGrob(g1, gPath="xaxis", label=samples, gp=gpar());
g1$children$xaxis <- editGrob(g1$children$xaxis, gPath="labels", rot=90, just="right", gp=gpar(col="black"));


for (i in 1:length(samples))
{
    #
    s <- samples[i];
    
    #
    bs <- b[b$sample == s, ];
    
    #
    v <- list(stats=vector(mode="numeric"), n=0, conf=vector(mode="numeric"), out=vector(mode="numeric"));
    
    #
    v$stats <- as.numeric(unlist(strsplit(bs$stats, ",")));
    v$n <- as.numeric(bs$n);
    v$conf <- as.numeric(unlist(strsplit(bs$conf, ",")));
    if (grepl(",", bs$out)) {
        v$out <- as.numeric(unlist(strsplit(bs$out, ",")));
    }
    
    # gp=gpar(col=cols3[i]), gp=gpar(col="navy")
    name <- paste("boxplot", i, sep="");
    g1 <- addGrob(g1, cs_boxplotGrob(v, at=i, width=width1, outline=FALSE, horizon=FALSE, gp=gpar(col="black"), name=name, vp=g1$childrenvp));
    
    #
    g1 <- editGrob(g1, gPath=gPath(name, "median"), gp=gpar(lwd=1.5));
}



# for plot 2

# defines the scale of a plot
xscale2 <- as.numeric(unlist(strsplit(args[1], ",")));    # NULL: automatically setup
yscale2 <- c(0, 20, 40, 60, 80, 100);


# draws plot 2
g2 <- cs_plot(title=title2, xlab=xlab2, ylab=ylab2, xscale=xscale2, yscale=yscale2, gp=gpar(), name=NULL, parent=vt2);

#
g2 <- addGrob(g2, cs_hrulerGrob(scale=yscale2, breaks=NULL, type=c("left", "right"), inner=TRUE, lty=3, lwd=0.5, gp=gpar(col="gray", alpha=0.7), name="hruler", vp=g2$childrenvp));
#g2 <- removeGrob(g2, gPath=gPath("hruler", "minor"));

#
g2 <- removeGrob(g2, gPath=gPath("vruler"));
#g2 <- removeGrob(g2, gPath=gPath("vruler", "minor"));


# coverage
lpchs <- c();
lcols <- c();
for (i in 1:length(samples))
{
    #
    s <- samples[i];
    
    #
    ds <- c[c$sample == s, ];
    
    # defines the color and symbol
    col <- colors[as.integer(i / length(pchars)) + 1];
    pch <- pchars[(i - 1) %% length(pchars) + 1];
    
    # min_depth_filter    percent_coverage
    g2 <- addGrob(g2, linesGrob(ds$min_depth_filter, ds$percent_coverage, default.units="native", name=paste("lines-", s, sep=""), gp=gpar(lty=1, lwd=1, col=col), vp=g2$childrenvp));
    g2 <- addGrob(g2, pointsGrob(ds$min_depth_filter, ds$percent_coverage, pch=pch, size=unit(0.7,"char"), default.units="native", name=paste("points", s, sep=""), gp=gpar(col=col), vp=g2$childrenvp));
    
    #
    lpchs = c(lpchs, pch);
    lcols = c(lcols, col);
}


# for the legend
g2 <- addGrob(g2, cs_vlegendGrob(just="right", x=unit(1.0, "npc")+unit(1.0, "char"), y=unit(0.5, "npc"), pch=lpchs, lty=rep(1, length(lpchs)), col=lcols, labels=samples, border=F, fontsize=8, cex=0.8, gp=gpar(), name="vlegend", vp=g2$childrenvp));
#g2 <- addGrob(g2, cs_vlegendGrob(just=c("right", "top"), x=unit(1.0, "npc")+unit(1.0, "char"), y=unit(0.0, "npc"), pch=lpchs, lty=rep(1, length(lpchs)), col=lcols, labels=samples, border=F, fontsize=8, cex=0.8, gp=gpar(), name="vlegend", vp=g2$childrenvp));



# for plot 3

# for the title
title3 <- paste(title3, " of ", nroi, " ROIs", sep="");

# defines the scale of a plot
if (max1 <= 5) {
    xscale3 <- cs_scale(d$x, breaks=5, range=c(0, max(d$x)));
} else {
    xscale3 <- cs_scale(d$x, breaks=5, range=c(0, max1));
    xscale3 <- c(xscale3, xscale3[length(xscale3)] + xscale3[length(xscale3)] - xscale3[length(xscale3) - 1], xscale3[length(xscale3)] + (xscale3[length(xscale3)] - xscale3[length(xscale3) - 1]) *2);
}

yscale3 <- cs_scale(d$y, breaks=5);


# draws plot 2
g3 <- cs_plot(title=title3, xlab=xlab3, ylab=ylab3, xscale=xscale3, yscale=yscale3, xfactor=0.05, gp=gpar(), name=NULL, parent=vt3);

# density
for (i in 1:length(samples))
{
    #
    s <- samples[i];
    
    #
    ds <- d[d$sample == s, ];
    
    # filters out data within the given range by the new x and y scale
    p <- cs_range_limit_xy(ds$x, ds$y, xlim=range(xscale3), ylim=range(yscale3));
    
    # defines the color and symbol
    col <- colors[as.integer(i / length(pchars)) + 1];
    pch <- pchars[(i - 1) %% length(pchars) + 1];
    
    # min_depth_filter    percent_coverage
    if (length(p$x) > 1)
        g3 <- addGrob(g3, cs_linesGrob(p$x, p$y, p=c(0.0,0.2,0.4,0.6,0.8,1.0), pch=pch, size=unit(0.7, "char"), pcol=col, default.units="native", name=NULL, gp=gpar(col=col), vp=g3$childrenvp));
}


# for the legend
g3 <- addGrob(g3, cs_vlegendGrob(just="right", x=unit(0.0, "npc"), y=unit(-0.6, "npc"), pch=lpchs, lty=rep(1, length(lpchs)), col=lcols, labels=samples, border=F, fontsize=8, cex=0.8, gp=gpar(), name="vlegend", vp=g3$childrenvp));



# for plot 4
# for the title
#title4 <- paste(title4, " of ", nroi, " ROIs (total ", nbase, ' bases)', sep="");

# defines the scale of a plot
xscale4 <- seq.int(1, length(samples), by=1);
yscale4 <- c(0, 20, 40, 60, 80, 100);


# draws plot 4
g4 <- cs_plot(title=title4, xlab=xlab4, ylab=ylab4, xscale=xscale4, xfactor=xfactor1, yscale=yscale4, gp=gpar(), name=NULL, parent=vt4);

#
g4 <- addGrob(g4, cs_hrulerGrob(scale=yscale4, breaks=0, type=c("left", "right"), inner=TRUE, lty=3, lwd=0.5, gp=gpar(col="gray", alpha=0.7), name="hruler", vp=g4$childrenvp));
#g2 <- removeGrob(g2, gPath=gPath("hruler", "minor"));

#
g4 <- removeGrob(g4, gPath=gPath("vruler"));
#g4 <- removeGrob(g4, gPath=gPath("vruler", "minor"));

#
g4 <- editGrob(g4, gPath="xaxis", label=samples, gp=gpar());
g4$children$xaxis <- editGrob(g4$children$xaxis, gPath="labels", rot=90, just="right", gp=gpar(col="black"));


# coverage
lcols <- c("red");
lfilter <- c(0);
for (i in 1:length(samples))
{
    #
    s <- samples[i];
    
    #
    ds <- c[c$sample == s, ];
    
    #
    g4 <- addGrob(g4, rectGrob(x=i, y=0, width=width1, height=100, just="bottom", default.units="native", name=NULL, gp=gpar(col="red", fill="red", lty="blank"), vp=g4$childrenvp));
    
    # makes arrays
    g4 <- addGrob(g4, rectGrob(x=i, y=0, width=width1, height=ds$percent_coverage, just="bottom", default.units="native", name=NULL, gp=gpar(col=ds$col, fill=ds$col, lty="blank"), vp=g4$childrenvp));
    
    # for the legend
    if (length(lcols) <= 1)
    {
        lcols <- c(lcols, ds$col);
        lfilter <- c(lfilter, ds$min_depth_filter);
    }
}


# for the legend
g4 <- addGrob(g4, cs_colorLegendGrob(1, 120, width=1.8, height=unit(1.5, "char"), title="Minimum depth", cols=lcols, labels=lfilter, gp=gpar(), name="clegend", vp=g4$childrenvp));



# draws plots
grid.draw(g1);
grid.draw(g2);
grid.draw(g3);
grid.draw(g4);



# -------------------------------------
# shuts down the current active device
dev.off();

# prints out warning messages if any and exits
warnings();
quit(save="no");



