<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template name="genome" match="object[@id='Genome']" priority="30">
          <div class="background">
            <div class="search_form_container">
              <div class="search_form">
                <form method="get" action="/cgi-bin/search/index.cgi">
                  <table cellspacing="0" cellpadding="0" border="0" class="form">
                    <tr>
                      <td style="white-space: nowrap; font-weight: bold;">
                        GC Search:
                      </td>
                      <td>
                        <input type="text" size="30" name="query" style="background-color: #FFF; font-size: 120%;"/><br/>
                      </td>
                      <td>
                        <input type="submit" class="search_button" value="Search"/>
                      </td>
                    </tr>
                  </table>
                </form>
              </div>
              <div class="search_help">
                <p></p>
              </div>
            </div>
            <div class="page_padding">
              <br/>
              <br/>
              <h2 class="form_group" style="padding-top: 15px;">Search for Models and Builds</h2>
              <form action="status.cgi" method="GET">
                <input type="hidden" name="search_type" value="model_name" />
                <table cellpadding="0" cellspacing="0" border="0" class="form" width="100%">
                  <colgroup>
                    <col width="25%"/>
                    <col width="30%"/>
                    <col width="100%"/>
                  </colgroup>
                  <tbody>
                    <tr>
                      <td class="label">Model Name:</td>
                      <td class="input">
                        <input type="text" name="genome-model-name" value="" />
                      </td>
                      <td>
                        <input type="submit" name="Search" value="Search for Model Name" />
                      </td>
                    </tr>
                  </tbody>
                </table>
              </form>
              <form action="status.cgi" method="GET">
                <input type="hidden" name="search_type" value="model_id" />
                <table cellpadding="0" cellspacing="0" border="0" class="form" width="100%">
                  <colgroup>
                    <col width="25%"/>
                    <col width="30%"/>
                    <col width="100%"/>
                  </colgroup>
                  <tbody>
                    <tr>
                      <td class="label">Model ID:</td>
                      <td class="input">
                        <input type="text" name="genome-model-id" value="" />
                      </td>
                      <td>
                        <input type="submit" name="Search" value="Search for Model ID" />
                      </td>
                    </tr>
                  </tbody>
                </table>
              </form>
              <form action="status.cgi" method="GET">
                <input type="hidden" name="search_type" value="model_user" />
                <table cellpadding="0" cellspacing="0" border="0" class="form" width="100%">
                  <colgroup>
                    <col width="25%"/>
                    <col width="30%"/>
                    <col width="100%"/>
                  </colgroup>
                  <tbody>
                    <tr>
                      <td class="label">Models Owned by User:</td>
                      <td class="input">
                        <select name="user_name">
                          <option value="" selected="selected">Select a User</option>
                          <xsl:for-each select="//users/user">
                            <option><xsl:attribute name="value"><xsl:value-of select="."/></xsl:attribute><xsl:value-of select="."/></option>
                          </xsl:for-each>
                        </select>
                      </td>
                      <td>
                        <input type="submit" name="Search" value="Search for Models" />
                      </td>
                    </tr>
                  </tbody>
                </table>
              </form>
              <form action="/rest/Genome/Model/Build/status.html" method="GET" name="buildform">
                <table cellpadding="0" cellspacing="0" border="0" class="form" width="100%">
                  <colgroup>
                    <col width="25%"/>
                    <col width="30%"/>
                    <col width="100%"/>
                  </colgroup>
                  <tbody>
                    <tr>
                      <td class="label">Build ID:</td>
                      <td class="input">
                        <input type="text" name="id" value="" />
                      </td>
                      <td>
                        <input type="button" name="Search" value="Search for Build ID" onClick="document.buildform.submit()" />
                      </td>
                    </tr>
                  </tbody>
                </table>
              </form>
              <form action="status.cgi" method="GET">
                <input type="hidden" name="search_type" value="build_status" />
                <table cellpadding="0" cellspacing="0" border="0" class="form" width="100%">
                  <colgroup>
                    <col width="25%"/>
                    <col width="30%"/>
                    <col width="100%"/>
                  </colgroup>
                  <tbody>
                    <tr>
                      <td class="label">Builds with Status:</td>
                      <td class="input">
                        <select name="event_status">
                          <option value="" selected="selected">Select a Status</option>
                          <xsl:for-each select="//event-statuses/event-status">
                            <option><xsl:attribute name="value"><xsl:value-of select="."/></xsl:attribute><xsl:value-of select="."/></option>
                          </xsl:for-each>
                        </select>
                      </td>
                      <td>
                        <input type="submit" name="Search" value="Search for Builds" />
                      </td>
                    </tr>
                  </tbody>
                </table>
              </form>
            </div>
          </div>
  </xsl:template>

</xsl:stylesheet> 
