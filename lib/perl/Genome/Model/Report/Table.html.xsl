<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
<xsl:output method="html"/>
<xsl:output doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"/>
<xsl:output doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"/>
<xsl:template match="/">
<xsl:variable name="newline">
<xsl:text>
</xsl:text>
</xsl:variable>
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
	<head>
		<meta http-equiv="Content-Type" content="text/html; charset=utf-8"/>
		<title><xsl:value-of select="//report-meta/description"/></title>
	</head>
	<body style="margin:0;padding:0;">
		<div align="center">
      <table border="2px" bordercolor="#f2f2f1" cellpadding="0" cellspacing="0" width="600" summary="report title and framework">
				<tbody>
					<tr>
						<td>
              <h1 style="font-family:Helvetica,Arial,sans-serif;font-size:115%;line-height:40px;font-weight:bold;color:#FFF;margin:0;padding:0 15px;background-color: #5f604c;border-bottom:5px solid #8dc643">
                <xsl:value-of select="//report-meta/name"/>
							</h1>
						</td>
					</tr>
          <tr>
						<td style="padding:5px 5px;">
              <table border="0" cellspacing="0" cellpadding="0" width="100%">
								  <tbody>
									<tr>
										<td style="padding:10px;" colspan="2">
                      <table border="1" cellspacing="0" cellpadding="0" width="100%" style="margin:5px 0 5px 5px;padding:0;" summary="report key value data table">
                        <xsl:for-each select="//headers/*">
                        <th style="text-align:left;font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 5px;">
                          <xsl:value-of select="."/>
                        </th>
                        </xsl:for-each>
												  <colgroup>
													  <col/>
													  <col width="100%"/>
												  </colgroup>
                        <xsl:for-each select="//datasets/*/*">
                        <tr>
                          <xsl:for-each select="./*">
                          <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 5px;">
                            <xsl:value-of select="."/>
                          </td>
                          </xsl:for-each>
                        </tr>
                        <xsl:value-of select="$newline"/>
                        </xsl:for-each>
                      </table>
                    </td>
                  </tr>
                </tbody>
              </table>
            </td>
          </tr>
					<tr>
						<td align="right" style="background-color: #f2f2f1;">
							<img src="http://genome.wustl.edu/images/uploads/genome_center_logo.png" width="106" height="50" alt="The Genome Center at Washington University" style="margin: 10px;" />
						</td>
					</tr>
				</tbody>
			</table>
		</div>
	</body>
</html>
</xsl:template>
</xsl:stylesheet>
