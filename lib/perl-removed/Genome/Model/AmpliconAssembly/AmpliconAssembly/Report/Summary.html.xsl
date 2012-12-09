<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="html"/>
  <xsl:output doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"/>
  <xsl:output doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"/>
  <xsl:template match="/">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
	<head>
		<meta http-equiv="Content-Type" content="text/html; charset=utf-8"/>
		<title>
          Summary for Amplicon Assembly Build #<xsl:value-of select="//report-meta/generator-params/build-id"/>
		</title>
	</head>
	<body style="margin:0;padding:0;">
		<div align="center">
			<table border="0" cellpadding="0" cellspacing="0" width="600" summary="report title and framework">
				<tbody>
					<tr>
						<td>
							<h1 style="font-family:Helvetica,Arial,sans-serif;font-size:105%;line-height:45px;font-weight:bold;color:#FFF;margin:0;padding:0 15px;background-color: #5f604c;border-bottom:5px solid #8dc643">
                              Summary for Amplicon Assembly Build #<xsl:value-of select="//report-meta/generator-params/build-id"/>
							</h1>
						</td>
					</tr>
					<tr>
						<td style="padding:15px 25px;">
                          <table border="0" cellspacing="0" cellpadding="0" width="100%" style="margin:20px 0 0 0;padding:0;" summary="report subsection header">
                            <tr>
                              <td style="border-bottom: 2px solid #CCC;">
                                <h2 style="font-family:Helvetica,Arial,sans-serif;font-size:100%;color:#36372c;margin:0 0 0 0;padding:5px 0 3px 5px;">
                                  Build Information
                                </h2>
                              </td>
                              <td style="border-bottom: 2px solid #CCC;" align="right">
                                <a style="text-decoration:none;font-size:70%;font-weight:bold;font-family:Helvetica,Arial,sans-serif;line-height:1.2;"><xsl:attribute name="href">https://gscweb.gsc.wustl.edu<xsl:value-of select="//model-info/data-directory"/>/reports/<xsl:value-of select="//report-meta/name"/>/report.html</xsl:attribute>view full report &gt;</a>
                              </td>
                            </tr>
                          </table>

							<!-- report section -->
							<table border="0" cellspacing="0" cellpadding="0" width="100%" summary="report section header">
								<tbody>
									<tr>
										<td style="padding:10px;" colspan="2">
											<!-- report subsection-->

											<table border="0" cellspacing="0" cellpadding="0" width="100%" style="margin:5px 0 5px 5px;padding:0;" summary="report key value data table">
												<colgroup>
													<col/>
													<col width="100%"/>
												</colgroup>
												<tr>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      build id: 
													</td>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                    <xsl:value-of select="//report-meta/generator-params/build-id"/> 													</td>
												</tr>
												<tr>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      data directory:
													</td>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      <a><xsl:attribute name="href">https://<xsl:value-of select="//model-info/data-directory"/></xsl:attribute><xsl:text>View Directory</xsl:text></a>
													</td>
												</tr>
												<tr>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      report generated:
													</td>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      <xsl:value-of select="//report-meta/date"/>
													</td>
												</tr>
												<tr>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      generator:
													</td>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      <xsl:value-of select="//report-meta/generator"/>
													</td>
												</tr>

											</table>
                                            
                                            <!-- report subsection-->
											<table border="0" cellspacing="0" cellpadding="0" width="100%" style="margin:20px 0 0 0;padding:0;" summary="report subsection header">
												<tr>
													<td style="border-bottom: 2px solid #CCC;">
														<h2 style="font-family:Helvetica,Arial,sans-serif;font-size:100%;color:#36372c;margin:0 0 0 0;padding:5px 0 3px 5px;">
															 Model Information
														</h2>
													</td>
												</tr>
											</table>
											<table border="0" cellspacing="0" cellpadding="0" width="100%" style="margin:5px 0 5px 5px;padding:0;" summary="report subsection data table">
												<colgroup>
													<col/>
													<col width="100%"/>
												</colgroup>
												<tr>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      model name:
													</td>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      <xsl:value-of select="//model-info/name"/>
													</td>
												</tr>
												<tr>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
														subject name:
													</td>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      <xsl:value-of select="//model-info/subject-name"/>
													</td>
												</tr>
												<tr>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
														subject type:
													</td>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      <xsl:value-of select="//model-info/subject-type"/>
													</td>
												</tr>
												<tr>
                                                  <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                    purpose:
                                                  </td>
                                                  <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                    <xsl:value-of select="//model-info/purpose"/>
                                                  </td>
												</tr>
												<tr>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      assembler:
													</td>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      <xsl:value-of select="//model-info/assembler"/>
													</td>
												</tr>
												<tr>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
													sequencing platform:
													</td>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      <xsl:value-of select="//model-info/sequencing-platform"/>
													</td>
												</tr>
												<tr>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
														sequencing center:
													</td>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      <xsl:value-of select="//model-info/sequencing-center"/>
													</td>
												</tr>
												<tr>
                                                  <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                    region of interest:
                                                  </td>
                                                  <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                    <xsl:value-of select="//model-info/region-of-interest"/>
                                                  </td>
												</tr>
												<tr>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      assembly size:
													</td>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      <xsl:value-of select="//model-info/assembly-size"/>
													</td>
												</tr>
												<tr>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      reverse seq primer:
													</td>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      <xsl:value-of select="//model-info/primer-seq-reverse"/>
													</td>
												</tr>
												<tr>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
													reverse amp primer:
													</td>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      <xsl:value-of select="//model-info/primer-amp-reverse"/>
													</td>
												</tr>
												<tr>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      forward seq primer:
													</td>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      <xsl:value-of select="//model-info/primer-seq-forward"/>
													</td>
												</tr>
												<tr>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      forward amp primer:
													</td>
													<td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                      <xsl:value-of select="//model-info/primer-amp-forward"/>
													</td>
												</tr>

											</table>
										

                                            <!-- report subsection-->
											<table border="0" cellspacing="0" cellpadding="0" width="100%" style="margin:20px 0 0 0;padding:0;" summary="report subsection header">
                                              <tr>
                                                <td style="border-bottom: 2px solid #CCC;">
                                                  <h2 style="font-family:Helvetica,Arial,sans-serif;font-size:100%;color:#36372c;margin:0 0 0 0;padding:5px 0 3px 5px;">
                                                    Assembly Stats
                                                  </h2>
                                                </td>
                                              </tr>
											</table>
											<table border="0" cellspacing="0" cellpadding="0" width="100%" style="margin:5px 0 5px 5px;padding:0;" summary="report subsection data table">
                                              <colgroup>
                                                <col/>
                                                <col width="100%"/>
                                              </colgroup>
                                              <tr>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  attempted:
                                                </td>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  <xsl:value-of select="//stat/attempted"/>
                                                </td>
                                              </tr>
                                              <tr>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  assembled:
                                                </td>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  <xsl:value-of select="//stat/assembled"/>
                                                </td>
                                              </tr>
                                              <tr>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  assembly success:
                                                </td>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  <xsl:value-of select="//stat/assembly-success"/>
                                                </td>
                                              </tr>
                                              <tr>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  average length:
                                                </td>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  <xsl:value-of select="//stat/length-average"/>
                                                </td>
                                              </tr>
                                              <tr>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  max. length:
                                                </td>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  <xsl:value-of select="//stat/length-maximum"/>
                                                </td>
                                              </tr>
                                              <tr>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  median length:
                                                </td>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  <xsl:value-of select="//stat/length-median"/>
                                                </td>
                                              </tr>
                                              <tr>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  minimum length:
                                                </td>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  <xsl:value-of select="//stat/length-minimum"/>
                                                </td>
                                              </tr>
                                              <tr>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  average base quality:
                                                </td>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  <xsl:value-of select="//stat/quality-base-average"/>
                                                </td>
                                              </tr>
                                              <tr>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  less than 20 bases<br/> per assembly quality
                                                </td>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  <xsl:value-of select="//stat/quality-less-than-20-bases-per-assembly"/>
                                                </td>
                                              </tr>
                                              <tr>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  assembled reads:
                                                </td>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  <xsl:value-of select="//stat/reads-assembled"/>
                                                </td>
                                              </tr>
                                              <tr>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  average assembled reads:
                                                </td>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  <xsl:value-of select="//stat/reads-assembled-average"/>
                                                </td>
                                              </tr>
                                              <tr>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  maximum assembled reads:
                                                </td>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  <xsl:value-of select="//stat/reads-assembled-maximum"/>
                                                </td>
                                              </tr>
                                              <tr>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  median assembled reads:
                                                </td>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  <xsl:value-of select="//stat/reads-assembled-median"/>
                                                </td>
                                              </tr>
                                              <tr>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;font-weight:bold;text-align:left;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  successful assembled reads
                                                </td>
                                                <td style="font-family:Helvetica,Arial,sans-serif;font-size:90%;color:#36372c;white-space:nowrap;padding:5px 5px 5px 0;">
                                                  <xsl:value-of select="//stat/reads-assembled-success"/>
                                                </td>
                                              </tr>
											</table>
										</td>
									</tr>
								</tbody>
							</table>
						</td>
					</tr>
					<tr>
						<td align="right" style="background-color: #f2f2f1;">
							<img src="cid:img1" width="106" height="50" alt="The Genome Center at Washington University" style="margin: 10px;" />
						</td>
					</tr>
				</tbody>
			</table>
		</div>
	</body>
</html>
  </xsl:template>
</xsl:stylesheet>
