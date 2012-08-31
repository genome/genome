/*
 * File:        FixedHeader.js
 * Version:     1.0.2
 * CVS:         $Id$
 * Description: "Fix" a header at the top of the table, so it scrolls with the table
 * Author:      Allan Jardine (www.sprymedia.co.uk)
 * Created:     Wed 16 Sep 2009 19:46:30 BST
 * Modified:    $Date$ by $Author$
 * Language:    Javascript
 * License:     LGPL
 * Project:     Just a little bit of fun :-)
 * Contact:     www.sprymedia.co.uk/contact
 * 
 * Copyright 2009 Allan Jardine, all rights reserved.
 *
 */


(function($) {

/*
 * Function: $.fn.dataTableExt.FixedHeader
 * Purpose:  FixedHeader "class"
 * Returns:  same as _fnInit
 * Inputs:   same as _fnInit
 */
$.fn.dataTableExt.FixedHeader = function ( oTable )
{
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Public functions
	 */
	
	/*
	 * Function: fnUpdate
	 * Purpose:  Update the floating header from the current state of the DataTable
	 * Returns:  -
	 * Inputs:   bool:bUpdateDom - Update offset information (provided for speed options) - default true
	 * Notes:    This would be used when the DataTables state is changed. For example using 
	 *   fnSetColumnVis() to change the number of visible columns.
	 */
	this.fnUpdate = function ( bUpdateDom )
	{
		if ( typeof bUpdateDom == 'undefined' )
		{
			bUpdateDom = true;
		}
		
		_fnCloneThead();
		
		if ( bUpdateDom )
		{
			_iStart = $(_oSettings.nTable).offset().top;
			_iStartLeft = $(_oSettings.nTable).offset().left;
		}
	}
	
	
	
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Private variables
	 */
	
	/* The DataTables object */
	var _oTable;
	
	/* The DataTables settings object - easy access */
	var _oSettings;
	
	/* The cloned table node */
	var _nCTable;
	
	/* The starting x-position of the table on the document */
	var _iStart;
	
	/* The starting x-position of the table relative to it's parent */
	var _iOffset;
	
	/* Current display information, cached so we don't have to query the DOM */
	var _oCache = {
		"sPosition": "",
		"sTop": "",
		"sLeft": ""
	};
	
	var _bIsIE6;
	
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Private functions
	 */
	
	/*
	 * Function: _fnCloneTable
	 * Purpose:  Clone the table node and do basic initialisation
	 * Returns:  -
	 * Inputs:   -
	 */
	function _fnCloneTable ()
	{
		var nOrigTable = _oSettings.nTable;
		
		/* We know that the table _MUST_ has a DIV wrapped around it, because this is simply how
		 * DataTables works. Therefore, we can set this to be relatively position (if it is not
		 * alreadu absolute, and use this as the base point for the cloned header
		 */
		if ( $(nOrigTable.parentNode).css('position') != "absolute" )
		{
			nOrigTable.parentNode.style.position = "relative";
		}
		
		/* Need to know the table's position relative to the other elements */
		_iOffset = nOrigTable.offsetTop;
		
		/* Just a shallow clone will do - we only want the table node */
		_nCTable = nOrigTable.cloneNode( false );
		_nCTable.style.position = "absolute";
		_nCTable.style.top = _iOffset+"px";
		_nCTable.style.left = nOrigTable.offsetLeft+"px";
		_nCTable.className += " FixedHeader_Cloned";
		_nCTable.id += "_Cloned";
		
		/* Insert the newly cloned table into the DOM, on top of the "real" header */
		nOrigTable.parentNode.insertBefore( _nCTable, nOrigTable );
		
		/* Dev note: for some mental reason we can't use the offset of '_nCTable' in IE. The original
		 * table will do us nicely though
		 */
		_iStart = $(_oSettings.nTable).offset().top;
		_iStartLeft = $(_oSettings.nTable).offset().left;
		
		/* Add the scroll event handler to move the table header */
		$(window).scroll( function () {
			var iWindow = $(window).scrollTop();
			
			if ( _bIsIE6 )
			{
				if ( _iStart < iWindow )
				{
					var iNew = iWindow-_iStart+_iOffset;
					var iTbodyHeight = _oSettings.nTable.getElementsByTagName('tbody')[0].offsetHeight;
					
					if ( iNew < _iOffset+iTbodyHeight )
					{
						/* In the middle of the table */
						_fnUpdateCache( 'sTop', iNew+"px", 'top', _nCTable.style );
					}
					else
					{
						/* At the bottom of the table */
						_fnUpdateCache( 'sTop', (_iOffset+iTbodyHeight)+"px", 'top', _nCTable.style );
					}
				}
				else
				{
					/* Above the table */
					_fnUpdateCache( 'sTop', _iOffset+"px", 'top', _nCTable.style );
				}
			}
			else
			{
				if ( _iStart < iWindow )
				{
					var iNew = iWindow-_iStart+_iOffset;
					var iTbodyHeight = _oSettings.nTable.getElementsByTagName('tbody')[0].offsetHeight;
					
					if ( iNew < _iOffset+iTbodyHeight )
					{
						/* In the middle of the table */
						_fnUpdateCache( 'sPosition', 'fixed',          'position', _nCTable.style );
						_fnUpdateCache( 'sTop',      "0px",            'top',      _nCTable.style );
						_fnUpdateCache( 'sLeft',     _iStartLeft+"px", 'left',     _nCTable.style );
					}
					else
					{
						/* At the bottom of the table */
						_fnUpdateCache( 'sPosition', 'absolute',                 'position', _nCTable.style );
						_fnUpdateCache( 'sTop',      _iOffset+iTbodyHeight+"px", 'top',      _nCTable.style );
						_fnUpdateCache( 'sLeft',     "0px",                      'left',     _nCTable.style );
					}
				}
				else
				{
					/* Above the table */
					_fnUpdateCache( 'sPosition', 'absolute',    'position', _nCTable.style );
					_fnUpdateCache( 'sTop',      _iOffset+"px", 'top',      _nCTable.style );
					_fnUpdateCache( 'sLeft',     "0px",         'left',     _nCTable.style );
				}
			}
		} );
	}
	
	/*
	 * Function: _fnUpdateCache
	 * Purpose:  Check the cache and update cache and value if needed
	 * Returns:  -
	 * Inputs:   string:sCache - cache property
	 *           string:sSet - value to set
	 *           string:sProperty - object property to set
	 *           object:oObj - object to update
	 */
	function _fnUpdateCache ( sCache, sSet, sProperty, oObj )
	{
		if ( _oCache[sCache] != sSet )
		{
			oObj[sProperty] = sSet;
			_oCache[sCache] = sSet;
		}
	}
	
	/*
	 * Function: _fnCloneThead
	 * Purpose:  Clone the THEAD element used in the DataTable and add required event listeners
	 * Returns:  -
	 * Inputs:   -
	 */
	function _fnCloneThead ()
	{
		/* Remove any children the cloned table has */
		while ( _nCTable.childNodes.length > 0 )
		{
			$('thead th', _nCTable).unbind( 'click' );
			_nCTable.removeChild( _nCTable.childNodes[0] );
		}
		
		/* Clone the DataTables header */
		var nThead = $('thead', _oSettings.nTable).clone(false)[0];
		_nCTable.appendChild( nThead );
		
		/* Copy the widths across - apparently a clone isn't good enough for this */
		$("thead:eq(0)>tr th", _oSettings.nTable).each( function (i) {
			$("thead:eq(0)>tr th:eq("+i+")", _nCTable).width( $(this).width() );
		} );
		
		$("thead:eq(0)>tr td", _oSettings.nTable).each( function (i) {
			$("thead:eq(0)>tr th:eq("+i+")", _nCTable)[0].style.width( $(this).width() );
		} );
		
		/* Add the event handlers for sorting */
		$('thead th', _nCTable).click( function (e) {
			/* Don't try and do the sort ourselves - let DataTables take care of the logic */
			var iTrigger = $('thead th', _nCTable).index(this);
			$('thead th:eq('+iTrigger+')', _oSettings.nTable).click();
			_fnCloneThead();
		} );
	}
	
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Initialisation
	 */
	_oTable = oTable;
	_oSettings = _oTable.fnSettings();
	
	//_bIsIE6 = ($.browser.msie && $.browser.version=="6.0");
	_bIsIE6 = ($.browser.msie && ($.browser.version=="6.0"||$.browser.version=="7.0"));
	
	_fnCloneTable();
	_fnCloneThead();
}

})(jQuery);
