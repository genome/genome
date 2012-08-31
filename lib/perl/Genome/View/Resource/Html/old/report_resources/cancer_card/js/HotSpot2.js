/*
	Created: Jim O'Donnell, National Maritime Museum
	Date: 02/2007 
	Description: Controller script for zoomable, annotated pictures.
	Usage: 	Create a new zoomable image - HotSpotController.init(picID, size, imagePath, buttonID, callback)
			Add a note (optional) - HotSpotController.setHotspotSize(spotSize); - only need to call this once.
									HotSpotController.addHotspot(x,y,heading,note); - repeat as necessary.
	Input:	picID		id of img element to which the overlay will be added.
			size		magnifier size (pixels).
			imagePath	image path to use as the source of the magnified view.
			buttonID	id of the element to use as an on/off button.
			callback	optional callback function to execute once the magnifier has finished loading.
			spotSize	hotspot size (pixels).
			x,y			x,y coordinates (on the larger image) of a note.
			heading		note heading.
			note		note text (one para max). May contain inline HTML tags (<em>, <a>, <img> etc.)
						but not block level HTML (<p>, <div>, <h1> etc.)
	Usage example:
	<p id="ZTthumbnail">
		<img src="/collections/images/700/PY/32/PY3270.jpg" id="fullscreen"  alt="London, from Greenwich. Picture in the possession of Walter Fawkes Esqr. of Farnley" /> 
	</p>
	<script language="JavaScript" type="text/javascript">
		HotSpotController.init("fullscreen",300, '/collections/images/2000/PY3270.jpg', 'magnifyButton');
		HotSpotController.setHotspotSize(150);
		HotSpotController.addHotspot(0, 750, 'St Alfege Church', 'St Alfege, the parish church of Greenwich, commemorates Alfege, Saxon Archbishop of Canterbury. Having been captured by Viking raiders in 1012, he was taken to their camp at Greenwich and murdered after a drunken feast. A church was built on the supposed site of his death in the 12th century, succeeded by a later medieval one whose roof collapsed in a storm in 1710. The present church was completed by 1714 except for the tower, which is that of the earlier church faced in 1730 with Portland stone and with a spire also added at that time.');
		HotSpotController.addHotspot(887, 450, 'St. Paul’s Cathedral', 'A cathedral dedicated to St Paul has overlooked the City of London since AD 604. The present cathedral, depicted here, is the fourth to occupy the same site. Designed by Sir Christopher Wren, it was built between 1675 and 1710 after its predecessor was destroyed in the Great Fire of London, 1666.');
	</script>
*/
var Moveable = {
	obj: null,
	init: function(obj) {
		obj.move = this.move;
		this.x = parseInt(obj.style.left);
		this.y = parseInt(obj.style.top);
	},
	move: function(x,y) {
		Moveable.obj = this;
		var o = Moveable.obj;
        o.style.top = y +"px";
        o.style.left = x + "px";		
	}
};

var HotSpotPic = {
	obj: null,
	hotspots: new Array(),
	hotspotSize: null,
	activeHotspot: null,
	init: function(obj) {
		obj.hotspots = this.hotspots;
		obj.hotspotSize = this.hotspotSize;
		obj.activeHotspot = this.activeHotspot;
		obj.addHotspot = this.addHotspot;
		obj.checkMessage = this.checkMessage;
		obj.getHotspot = this.getHotspot;
	},
	addHotspot: function(x,y,node) {
		HotSpotPic.obj = this;
		var o = HotSpotPic.obj;
		var hotspot = {x: x, y: y, node: node}
		o.hotspots.push(hotspot);
	},
	getHotspot: function(x,y) {
//return any stored note at x,y

		HotSpotPic.obj = this;
		var o = HotSpotPic.obj;
		var tmpHs = null;
		for (i=0;i<o.hotspots.length;i++) {
			var hotspot = o.hotspots[i];
			var node = hotspot.node;
			if (node != null) {
				if (x > (hotspot.x-o.hotspotSize) && x < (hotspot.x+o.hotspotSize) && y > (hotspot.y-o.hotspotSize) && y < (hotspot.y+o.hotspotSize)) {
					tmpHs = hotspot;
				}
			}
		}
		return tmpHs;
	},
	checkMessage: function(x,y) {
//is there a note at position (x,y)? If so, make it visible.
		HotSpotPic.obj = this;
		var o = HotSpotPic.obj;
		var className = '';
		
		var hotspot = this.getHotspot(x,y);
//Has the active hotspot changed? If so, hide the old one and make the new one active.
		if (hotspot != this.activeHotspot) {
			
			if (this.activeHotspot) {
				fadeOut(this.activeHotspot.node, 100);
				this.activeHotspot.node.style.zIndex = 1;
			}
			
			this.activeHotspot = hotspot;
			
			if (hotspot) {
				var tmpNode = hotspot.node;
				if (hotspot.x > HotSpotController.xMax/(2*HotSpotController.scaleFactor)) className = 'leftAlign';
		 		tmpNode.style.display = "block";
				tmpNode.style.zIndex = 100;
				fadeIn(tmpNode, 0);
				document.getElementById('captionContainer').className = className;
			}
		}
		
	}

};

var HotSpotController = {
	init: function(thumbnailID, viewPortSize, filePath, toggleID, callback) {
	
//reduced size image
		var thumbImage = document.getElementById(thumbnailID);
		
// control to toggle the overlay on/off;
		var toggleControl = document.getElementById(toggleID);
		
		this.hasRun = true;
		
//Container needs to wrap around image.
		thumbImage.parentNode.style.width = thumbImage.width+'px';

//create the overlay container.
		this.overlay = document.createElement('div');
		this.overlay.setAttribute('id', 'ZToverlay');
		this.overlay.style.visibility = 'hidden';
		//this.overlay.style.width = viewPortSize+'px';
		//this.overlay.style.height = viewPortSize+'px';
		//this.overlay.style.overflow = 'auto';
		thumbImage.parentNode.insertBefore(this.overlay, thumbImage);
// link so overlay can have keyboard focus
		var accessLink = document.createElement('a');
		accessLink.setAttribute('href', 'javascript:void(0);'); // slightly icky. oh well.
		addEvent(accessLink, 'keydown', keyPress);
		addEvent(accessLink, 'keypress', keyPress);
		this.overlay.appendChild(accessLink);
//draggable square
		this.magnifier = document.createElement('div');
		this.magnifier.setAttribute('id', 'magnifier');
		accessLink.appendChild(this.magnifier);
//container for the full size image - defines the viewable region.
		var viewPort = document.createElement('div');
		viewPort.setAttribute('id', 'viewPort');
//full size image.
		var image = document.createElement('img');
			image.setAttribute('src', filePath);
			image.setAttribute('id', 'ZTview');
			image.style.top = "0";
			image.style.left = "0";
		viewPort.appendChild(image);
// cancel right clicks and context menu events on the hi-res image.
		addEvent(image, 'mousedown', cancelRightClick);
		addEvent(image, 'contextmenu', cancelEvent);
		accessLink.appendChild(viewPort);
		this.zoomedView = image;
//initialise the view port height and width.
		this.setSize(viewPort, viewPortSize, viewPortSize);
		
		this.viewPortSize = viewPortSize;
		
		viewPort.style.top = (-1*viewPortSize)-10 + "px";
		viewPort.style.left = (-1*viewPortSize)/2.4 + "px";

//make various bits and pieces moveable etc. These are all decorator classes which add methods to existing DOM objects.
		Moveable.init(this.zoomedView);
		Moveable.init(this.overlay);
		HotSpotPic.init(this.zoomedView);
		
// add a handler to finish setup once the full size image has finished loading.
// stupid opera doesn't fire onload when a cached image loads - use img.complete instead.
		if (image.complete) {
			this.finishSetup(thumbImage, toggleControl, callback);
		} else {
			addEvent(image, 'load', function() {HotSpotController.finishSetup(thumbImage, toggleControl, callback)});
//if the image doesn't existing, remove the loading message and quit without finishing the set-up.
			addEvent(image, 'error', function() {
											  deleteLoadingMessage('loadingMessage');
											  });
			setLoadingMessage(thumbImage, 'Downloading magnified view', 'loadingMessage');
		}
	},
	finishSetup: function(thumbImage, toggleControl, callback) {
		deleteLoadingMessage('loadingMessage');
//scale factor for scaling between motion of the small image and motion of the large image.
// don't try this until zoomedView.width > 0 ie. until zoomedView has loaded.
		HotSpotController.scaleFactor = thumbImage.width/this.zoomedView.width;
//rescale magnifier to match our viewport size.
		HotSpotController.setSize(this.magnifier, this.viewPortSize*this.scaleFactor, this.viewPortSize*this.scaleFactor);
//set drag boundaries for the magnifier
		HotSpotController.xMax = thumbImage.width - (this.viewPortSize * this.scaleFactor);
		HotSpotController.yMax = thumbImage.height - (this.viewPortSize * this.scaleFactor);
		

// make the overlay draggable and add a handler to respond to drag events.
		Drag.init(this.overlay, null, 0, this.xMax, 0, this.yMax);
		HotSpotController.overlay.onDrag = function(x, y) {
	//When the magnifier is dragged, move the large image back by an equivalent amount to bring the selected region into view.
	//Check to see if we should display a message at the new position.
			curr_x = x/HotSpotController.scaleFactor;
			curr_y = y/HotSpotController.scaleFactor;
	//Large image moves in the opposite direction to the magnifier.
			HotSpotController.zoomedView.move(-1*curr_x, -1*curr_y);
			HotSpotController.zoomedView.checkMessage(curr_x,curr_y);
		}
		
//Suggestion from Pat Lauke - start in the center of the picture.
		HotSpotController.jump(thumbImage.width/2,thumbImage.height/2);
// add event handlers for jumping around the picture and toggling the overlay on and off.
		addEvent(thumbImage, 'click', function(e) {
											   toggleControl.focus();
											   HotSpotController.jumpToMouseClick(e);
											   cancelEvent(e);
											   });
// Why doesn't the following line work for the enter key in Safari?
		addEvent(toggleControl, 'click', function(e) {
												  toggleControl.focus();
												  toggleVisibility(HotSpotController.overlay);
												  HotSpotController.toggleHotspots();
												  cancelEvent(e);
												  return false;
												  });
// Add support for keyboard control
		addEvent(toggleControl, 'keydown', keyPress);
		addEvent(toggleControl, 'keypress', keyPress);
		toggleControl.className = 'active';
//controls for notes. Show notes if query param present in URL.
		var noteButton = document.getElementById('ZTnotes');
		if (noteButton) {
			addEvent(noteButton, 'click', HotSpotController.toggleHotspots);
			noteButton.className = 'active';
		}
		if(callback) callback();
	},
	setHotspotSize: function(size) {
		this.zoomedView.hotspotSize = size;
	},
	addHotspot: function(x,y,heading,caption) {
		var captionContainer = document.getElementById('captionContainer');
		var node = document.createElement('li');
		node.className = 'HScaption';
		node.appendChild(getHTMLElement('h2',heading));
		node.appendChild(getHTMLElement('p',caption));
		if (!captionContainer) {
			captionContainer = document.createElement('ul');
			captionContainer.setAttribute('id', 'captionContainer');
			this.overlay.appendChild(captionContainer);
		}
		captionContainer.appendChild(node);
		this.zoomedView.addHotspot(x,y,node);
	},
	toggleHotspots: function() {
		if (HotSpotController.showNotes) {
			HotSpotController.hideHotspots();
		} else {
			HotSpotController.showHotspots();
		}
	},
	showHotspots: function() {
		var len = this.zoomedView.hotspots.length;
		for ( var i=0; i<len; ++i ){
			if(typeof(console) != 'undefined') console.log(i);
			var marker = this.getHotspotMarker(this.zoomedView.hotspots[i]);
			marker.style.display = 'block';
		}
		this.showNotes = true;
	},
	hideHotspots: function() {
		var len = this.zoomedView.hotspots.length;
		for ( var i=0; i<len; ++i ){
			var marker = this.getHotspotMarker(this.zoomedView.hotspots[i]);
			marker.style.display = 'none';
		}
		this.showNotes = false;
	},
	getHotspotMarker: function(hotspot) {
/*
create a new DOM node to display a hotspot on the screen.
The new node is a sibling of the thumbnail image.
We will use absolute positioning to place the new node relative to the 
top left corner of the thumbnail.
*/
		var id = 'hsx' + hotspot.x + 'y' + hotspot.y;
		if (!document.getElementById(id)) {
			var newNode = document.createElement('a');
			var image = document.createElement('img');
			var size = this.viewPortSize * this.scaleFactor;
			var parent = this.overlay.parentNode;
			
			this.setSize(image, size, size);
			this.setSize(newNode, size, size);

//without this, the link won't be clickable in IE6
//			image.setAttribute('src', '/collections/shared/images/transparent.gif');
//			image.setAttribute('alt', '');
//			newNode.appendChild(image);
			
			newNode.setAttribute('href', 'javascript:void(0);');
			newNode.style.position = 'absolute';
			newNode.style.display = 'none';
			newNode.style.zIndex = 10;
			newNode.style.left = hotspot.x*this.scaleFactor+'px';
			newNode.style.top = hotspot.y*this.scaleFactor+'px';
			newNode.className = 'translucent';
			newNode.id = id;
			addEvent(newNode, 'focus', function() {
					if (HotSpotController.overlay.style.visibility == 'hidden') {
						HotSpotController.overlay.style.visibility = 'visible';
					}
					HotSpotController.jump(parseInt(hotspot.x*HotSpotController.scaleFactor)+size/2,parseInt(hotspot.y*HotSpotController.scaleFactor)+size/2)
				});
			addEvent(newNode, 'click', function() {
					if (HotSpotController.overlay.style.visibility == 'hidden') {
						HotSpotController.overlay.style.visibility = 'visible';
					}
					HotSpotController.jump(parseInt(hotspot.x*HotSpotController.scaleFactor)+size/2,parseInt(hotspot.y*HotSpotController.scaleFactor)+size/2)
				});
			parent.appendChild(newNode);
			return newNode;
		} else {
			return document.getElementById(id);
		}
	},
	setSize: function(obj, height, width) {
//set the size of obj to height and width in pixels.
    	obj.style.width = width + "px";
		obj.style.height = height + "px";			
	},
	jumpToMouseClick: function(e) {
		var x = 0;
		var y = 0;
		var coords = getXY(e);
		
		x = coords.x;
		y = coords.y;
		HotSpotController.jump(x,y);
	},
	jump: function(x,y) {
		
		var size = parseInt(HotSpotController.magnifier.style.width)*0.5;
		x = parseInt(x);
		y = parseInt(y);
//		if (typeof console != 'undefined') console.log(size);
// suggestion from MCG list - keep x,y inside drag boundaries.
		x = x-size;
		y = y-size;
		x = Math.max(0,x);
		y = Math.max(0,y);
		x = Math.min(x, HotSpotController.xMax);
		y = Math.min(y, HotSpotController.yMax);
//		if (typeof console != 'undefined') console.log('x:'+x+' y:'+y);
		HotSpotController.overlay.move(x,y);
		HotSpotController.overlay.onDrag(x,y);
	},
	move: function(keyCode, step) {
	// experimental support for keyboard navigation.
	// move the overlay in response to the cursor keys.
	// only does up, down, left, right at present.
		var curr_x = parseInt(HotSpotController.overlay.style.left);
		var curr_y = parseInt(HotSpotController.overlay.style.top);
		var new_x = curr_x;
		var new_y = curr_y;
		
		switch(keyCode) {
			case 37: //left
				new_x = Math.max(curr_x - step, 0);
				new_y = curr_y;
			break;
			case 38: //up
				new_x = curr_x;
				new_y = Math.max(curr_y - step, 0);
			break;
			case 39: //right
				new_x = Math.min(curr_x + step, this.xMax);
				new_y = curr_y;
			break;
			case 40: //down
				new_x = curr_x;
				new_y = Math.min(curr_y + step, this.yMax);
			break;
		}
		HotSpotController.overlay.move(new_x, new_y);
		HotSpotController.zoomedView.move(-1*new_x/this.scaleFactor, -1*new_y/HotSpotController.scaleFactor);
		HotSpotController.zoomedView.checkMessage(new_x/this.scaleFactor, new_y/HotSpotController.scaleFactor);
	}

};

function getTextElement(elementName, text) {
	var node = document.createElement(elementName);
	var textNode = document.createTextNode(text);
	node.appendChild(textNode);
	return node;
}

function getHTMLElement(elementName, html) {
	var node = document.createElement(elementName);
	var d = document.createElement('div');
	d.innerHTML = html;

	while (d.firstChild) {
		node.appendChild(d.firstChild);
	}
	return node;
}

function toggleVisibility(node) {
	if(node.style.display == "none" || node.style.visibility == "hidden") {
		node.style.display = "block";
		node.style.visibility = "visible";
		standardsFadeIn(node,0);
		
	} else {
		
		node.style.display = "none";
	}
}

function getXY(e) {
// Return x,y position of mouse clicks relative to target element.
// Works in IE, Firefox and Safari. Haven't tried Opera.
		if (e.layerX) {
			x = e.layerX;
			y = e.layerY;
		} else if (e.offsetX) {
			x = e.offsetX;
			y = e.offsetY;
		} else {
			x = e.x;
			y = e.y;
		}
		
		return {x:x,y:y};
}

function cancelRightClick(e) {
	var rightclick;
	if (e.which) rightclick = (e.which == 3);
	else if (e.button) rightclick = (e.button == 2);
	if (rightclick) {
		cancelEvent(e);
	}
}

function cancelEvent(e) {
	if (e && e.preventDefault) e.preventDefault(); // DOM event cancel
	return false; // IE event cancel
}

function setLoadingMessage(node, message, id) {
	var lm = getTextElement('span', message);
	lm.setAttribute('id', id);
	node.parentNode.insertBefore(lm,node);
}

function deleteLoadingMessage(id) {
	var node = document.getElementById(id);
	if (node) node.parentNode.removeChild(node);
}

	function keyPress(e) {
		var step = 2;
		var keyCode = e.keyCode;
		
		if (e.shiftKey) step = 20;
		
	//Map safari keycodes to codes generated by other browsers.
		switch (keyCode) {
		  case 63234:
			keyCode = 37;
			break;
		  case 63232:
			keyCode = 38;
			break;
		  case 63235:
			keyCode = 39;
			break;
		  case 63233:
			keyCode = 40;
			break;
		  default:;
		}
	//pass the code for the key that was pressed to the controller's keyboard function
		if (keyCode >= 37 && keyCode <= 40) {
			HotSpotController.move(keyCode, step);
			if (e && e.preventDefault) e.preventDefault(); // DOM event cancel
			return false; // IE event cancel
		} else {
			return true;
		}
	}

	function setOpacity(obj, opacity) {
	  opacity = (opacity == 100)?99.999:opacity;
	  
	  // IE/Win
	  obj.style.filter = "alpha(opacity:"+opacity+")";
	  
	  // Safari<1.2, Konqueror
	  obj.style.KHTMLOpacity = opacity/100;
	  
	  // Older Mozilla and Firefox
	  obj.style.MozOpacity = opacity/100;
	  
	  // Safari 1.2, newer Firefox and Mozilla, CSS3
	  obj.style.opacity = opacity/100;
	}
	
	function standardsFadeIn(obj,opacity) {
		if (opacity <= 100) {
		  obj.style.opacity = opacity/100;
		  opacity += 10;
		  window.setTimeout(function() {
							standardsFadeIn(obj, opacity);
									 }, 
									 50);
		}
	}

	function fadeIn(obj,opacity) {
		if (opacity <= 100) {
		  setOpacity(obj, opacity);
		  opacity += 10;
		  window.setTimeout(function() {
							fadeIn(obj, opacity);
									 }, 
									 50);
		}
	}
	function fadeOut(obj,opacity) {
		if (opacity >= 0) {
		  setOpacity(obj, opacity);
		  opacity -= 10;
		  window.setTimeout(function() {
							fadeOut(obj, opacity);
									 }, 
									 50);
		} else {
			obj.style.display = 'none';
		}
	}
	function getQueryVariable(variable) {
  		var query = window.location.search.substring(1);
  		var vars = query.split("&");
  		for (var i=0;i<vars.length;i++) {
    		var pair = vars[i].split("=");
    		if (pair[0] == variable) {
      			return pair[1];
    		}
 		 }
  		return false;
	} 