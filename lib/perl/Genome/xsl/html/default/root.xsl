<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="html"/>
  <xsl:output doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"/>
  <xsl:output doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"/>

  <xsl:strip-space elements="*"/>

  <xsl:template match="/">
    <html>
      <head>
        <script src="http://ajax.googleapis.com/ajax/libs/jquery/1.4.2/jquery.min.js"></script>

        <script language="javascript" type="text/javascript">
          <![CDATA[
                   $(document).ready(function(){

//Hide (Collapse) the toggle containers on load
$(".toggle_container").hide();

//Slide up and down on click
$("span.trigger").click(function(){

$(this).parent().next(".toggle_container").slideToggle("slow");
});

});
          ]]>
        </script>
        <link rel="shortcut icon" href="/static/images/gc_favicon.png" type="image/png" />
        <link rel="stylesheet" href="/static/css/master.css" type="text/css" media="screen" />
        <style type="text/css" media="screen">
          <![CDATA[
                   body {
                   padding: 0;
                   margin: 0;
                   background-color:#efefef;
                   font-family: Helvetica, Arial, sans-serif;
                   color: #333;
                   }

div.container {
width: 800px;
margin: 0 auto;
font-size: 13px;
}

div.background {
float: left;
background-color: #FFF;
width: 800px;
border-bottom: 5px solid #AAA;
margin-bottom: 30px;
}

div.content {
width: 760px;
padding: 15px 20px;
}

div.header_object {
background: #CCC;
border-bottom: 5px solid #AAA;
}

div.header_object h1 {
margin: 0;
padding: 0;
line-height: 35px;
font-size: 22px;
padding-left: 15px;
}

span.id {
font-weight: normal;
}

span.display_name {

}

p.exception {
font-weight: bold;
color: #C33;
}

span.trigger {
color: #999;
font-weight: normal;
text-decoration: none;
}

span.trigger:hover {
cursor: pointer;
}

div.toggle_container {
border: 1px solid #EFEFEF;
padding: 5px 10px;
}

div.toggle_container p {
font-size: 10px;
font-family: "Courier New", Courier, Consolas, "Lucida Sans Typewriter", monospace;
max-width: 380px;
}

table.aspects {
width:100%;
border-collapse: collapse;
}

table.aspects td {
padding: 6px 5px 6px 5px;
border-bottom: 1px solid #EFEFEF;
}

table.aspects td.name {
text-align: right;
vertical-align: top;
}

table.aspects td.value p {
margin: 0 0 5px 0;
padding: 0;
}

table.hash td {
border: none;
}

table.hash td.name {
font-weight: bold;
}

          ]]>
        </style>
      </head>
      <body>
        <div class="container">
          <div class="background">
            <xsl:apply-templates/>
          </div>
        </div>
      </body>
    </html>
  </xsl:template>

</xsl:stylesheet>
