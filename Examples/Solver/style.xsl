<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
    
    <xsl:template match="/">
		<html>
			
			
			<body style="font-family:Arial;font-size:12pt;background-color:#FFFFFF">
				<h2>Debug:</h2>
				<xsl:apply-templates/>
			</body>
			<script type="text/javascript">
				function toggle(id) 
				{
					
					if (document.getElementById(id).style.display == 'none') 
					{
						document.getElementById(id).style.display = 'block';
					} 
					else 
					{
						document.getElementById(id).style.display = 'none';
					}
				}
			</script>
		</html>
	</xsl:template>
    
    
    <xsl:template match="MPI">
    <div style="color:#dc143c;margin-top:5px; margin-bottom:5px; spacing: 5px;">
		<span style="font-weight:bold"><xsl:value-of select="."/></span>
    </div>
    </xsl:template>
    
    
    <xsl:template match="VALUE">
		<div style="margin-top:2px; margin-bottom:5px;">
		<table width="100%" border="0" cellpadding="1" cellspacing="1">
			<tr>
				<td align="left">
					<span style="font-weight:bold; color:#141099"><xsl:value-of select="@name"/> = </span>
					<span style="font-weight:bold; color:#991199"><xsl:value-of select="CONTENT"/></span>
				</td>
				<td align="center" width="60%">
					<span style="color:#808000"><xsl:value-of select="CODE"/></span>
				</td>
			</tr>
		</table>
		</div>
    </xsl:template>
    
    <xsl:template match="TEXT">
    <div style="background-color:white;padding:4px;">
		<span style="font-weight:bold"><xsl:value-of select="."/></span>
    </div>
    </xsl:template>
    
    <xsl:template match="FUNCTION">
	<span onclick="javascript:toggle('{@id}')">
    <div style="background-color:white;padding:2px; border: solid 1px black;">
		<table width="100%"><tr>
		<td align="left"><span style="font-weight:bold"><xsl:value-of select="@name"/></span></td>
		<td align="right"><span style="font-weight:bold;"><xsl:value-of select="TIME"/></span></td>
		</tr></table>
    </div>
    </span>
    <div style="margin-left:6px; margin-right:0px; margin-bottom:2px; border: solid 1px black; padding: 5px" id="{@id}">
		<xsl:apply-templates select="MPI | VALUE | TEXT | FUNCTION"/> 
		<!-- <xsl:apply-templates select="FUNCTION"/> -->
	</div>
	
    </xsl:template>

</xsl:stylesheet> 

