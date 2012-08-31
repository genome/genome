//
// tjpzoom 2005.11.29 - 2005.12.02.
// details/usage: http://valid.tjp.hu/zoom/
//

var zoomw=160;
var zoomh=120;
var defzoomamount=4;

var zoomamountstep=1.2;
var zoomsizemin=1000;
var zoomsizemax=100000;
var zoomsizestep=1.2;
var zoomamountmin=1;
var zoomamountmax=12;

function zoom_set(evt) {
 var evt = evt?evt:window.event?window.event:null; if(!evt){ return;}
 if(zoomid=='' || parseInt(document.getElementById(zoomid+'_container').style.width) == 0) {return true;}
 if(evt.keyCode==37 || evt.keyCode==100) {zoomw/=zoomsizestep; zoomh/=zoomsizestep; if(zoomw*zoomh<zoomsizemin) {zoomh=Math.sqrt(zoomsizemin/zoomratio); zoomw=zoomh*zoomratio;} zoom_init(); zoom_move(); return;} //left
 if(evt.keyCode==39 || evt.keyCode==102) {
  zoomw*=zoomsizestep; zoomh*=zoomsizestep;
  if(zoomw*zoomh>zoomsizemax) {zoomh=Math.sqrt(zoomsizemax/zoomratio); zoomw=zoomh*zoomratio;}
  if(zoomw>objw) {zoomw=objw; zoomh=objw/zoomratio;}
  else if(zoomh>objh) {zoomh=objh; zoomw=objh*zoomratio}
  zoom_init(); zoom_move(); return;
 } //right
 if(evt.keyCode==40 || evt.keyCode==98)  {zoomamount/=zoomamountstep; if(zoomamount<zoomamountmin) {zoomamount=zoomamountmin;} zoom_init(); zoom_move(); return;} //down
 if(evt.keyCode==38 || evt.keyCode==104) {zoomamount*=zoomamountstep; if(zoomamount>zoomamountmax) {zoomamount=zoomamountmax;} zoom_init(); zoom_move(); return;} //up
 return;
}

function zoom_init() {
 document.getElementById(zoomid+'_clip').style.width=objw*zoomamount+'px';
 document.getElementById(zoomid+'_clip').style.height=objh*zoomamount+'px';
 setTimeout("document.getElementById('"+zoomid+"_trigger').focus();",0);
}

function zoom_move(evt) {
 if(typeof(evt) == 'object') {
  var evt = evt?evt:window.event?window.event:null; if(!evt){ return;}
  if(evt.pageX) {
   xx=evt.pageX - ffox;
   yy=evt.pageY - ffoy;
  } else {
   if(typeof(document.getElementById(zoomid)+1) == 'number') {return true;} //mert az ie egy ocska kurva
//   xx=evt.x - ieox;
//   yy=evt.y - ieoy;
   xx=evt.clientX - ieox;
   yy=evt.clientY - ieoy;
  }
 } else {
  xx=lastxx; yy=lastyy;
 }
 lastxx=xx; lastyy=yy;
 document.getElementById(zoomid+'_clip').style.margin=((yy-zoomh/2 > 0)?(zoomh/2-yy*zoomamount):(yy*(1-zoomamount)))+'px 0px 0px '+((xx-zoomw/2 > 0)?(zoomw/2-xx*zoomamount):(xx*(1-zoomamount)))+'px';
 document.getElementById(zoomid+'_container').style.margin=((yy-zoomh/2 > 0)?(yy-zoomh/2):(0))+'px 0px 0px '+((xx-zoomw/2 > 0)?(xx-zoomw/2):(0))+'px';

 w2=((xx+zoomw/2<objw)?((zoomw/2<xx)?(zoomw):(zoomw/2+xx)):(zoomw/2+objw-xx)); if(w2<0) {w2=0;} document.getElementById(zoomid+'_container').style.width=w2+'px';
 h2=((yy+zoomh/2<objh)?((zoomh/2<yy)?(zoomh):(zoomh/2+yy)):(zoomh/2+objh-yy)); if(h2<0) {h2=0;} document.getElementById(zoomid+'_container').style.height=h2+'px';
 return false;
}

function zoom_off() {
 if(zoomid != '') {
  document.getElementById(zoomid+'_container').style.width='0px';
  document.getElementById(zoomid+'_container').style.height='0px';
 }
 zoomid='';
}

function countoffset() {
 ieox=0;ieoy=0;
 for(zmi=0;zmi<50;zmi++) {
  zme='document.getElementById(zoomid)';
  for(zmj=1;zmj<=zmi;zmj++) {
   zme+='.offsetParent';
  }
  if(eval(zme)+1 == 1) {zmi=100} else {ieox+=eval(zme+'.offsetLeft'); ieoy+=eval(zme+'.offsetTop');}
 }
 ffox=ieox;
 ffoy=ieoy;
 if(document.documentElement && document.documentElement.scrollTop) {
  ieox-=document.documentElement.scrollLeft;
  ieoy-=document.documentElement.scrollTop;
 } else {
  ieox-=document.body.scrollLeft;
  ieoy-=document.body.scrollTop;
 }
}

function zoom_on(evt,ow,oh,lowres,highres) {
 thisid=lowres+highres+ow+oh;
 thisid='zoom_'+thisid.replace(/[^a-z0-9]/ig,'')
 if(zoomid==thisid) {return;}
 if(typeof(highres) == typeof(nemistudom)) {highres=lowres;}
 var evt = evt?evt:window.event?window.event:null; if(!evt){ return;}
 zoomamount=defzoomamount;
 objw=ow;
 objh=oh;
 zoomid=thisid;
 if(evt.pageX) {
  evt.target.parentNode.id=thisid;
  countoffset();
  lastxx=evt.pageX - ffox; //mert a ff egy ocska kurva
  lastyy=evt.pageY - ffoy;
 } else {
  evt.srcElement.parentNode.id=thisid;
  countoffset();
//  lastxx=evt.x - ieox; //mert az ie egy ocska kurva
//  lastyy=evt.y - ieoy;
  lastxx=evt.clientX - ieox;
  lastyy=evt.clientY - ieoy; 
 }
 document.getElementById(thisid).style.width=ow+'px';
 document.getElementById(thisid).style.height=oh+'px';
 document.getElementById(thisid).style.overflow='hidden';
 document.getElementById(thisid).innerHTML='<div style="position:absolute; overflow:hidden; width:0; height:0;" id="'+thisid+'_container" onmousemove="zoom_move(event);" onmouseout="zoom_off()"><img src="'+highres+'" alt="" id="'+thisid+'_clip" style="padding:0;margin:0;border:0;" /></div><img src="'+lowres+'" id="'+thisid+'_pic" alt="" style="padding:0;margin:0;border:0;width:'+ow+'px; height: '+oh+'px"/><input type="text" id="'+thisid+'_trigger"  style="width:0;height:0;border:0;margin:0; padding:0;overflow:hidden;"/>';
 if(zoomw>objw) {zoomw=objw; zoomh=objw/zoomratio;}
 else if(zoomh>objh) {zoomh=objh; zoomw=objh*zoomratio}
 zoom_init();
 zoom_move('');
 if(document.all) {
  document.onkeydown=zoom_set;
 } else {
  window.captureEvents(Event.KEYDOWN);
  window.onkeydown = zoom_set;
 }
 return false;
}

var zoomamount=defzoomamount; var objw;var objh;var zoomid=''; var zoomratio=zoomw/zoomh; var ieox=0; var ieoy=0; var ffox=0; var ffoy=0;