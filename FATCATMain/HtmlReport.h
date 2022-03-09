#ifndef WRITEHTML_H
#define WRITEHTML_H

#include "basic.h"

void WriteHtmlHead(char *title, ofstream &html)
{
	html<<"<html>\n<title>"<<title<<"</title>\n";
	html<<"<body bgcolor=\"#ffffff\" bgcolor=\"#ffffff\" text=\"#407040\" link=\"#333399\" alink=\"#9966ff\" vlink=\"#9966ff\" style=\"margin:1\">"<<endl;
}

void CenterCaption(int font, char *str, ofstream &html)
{
	html<<"<h"<<font<<"><center>"<<str<<"</center></h"<<font<<">\n";
}

void WriteCaption(int font, char *str, ofstream &html)
{
	html<<"<h"<<font<<">"<<str<<"</h"<<font<<">\n";
}

void WriteHtmlEnd(ofstream &html)
{
	html<<"<hr>"<<endl;
	html<<"<center>"<<endl;
	html<<"<a href=\"http://bioinformatics.ljcrf.edu\">Godzik Laboratory</a>,<br>"<<endl;
	html<<"<a href=\"http://www.burnham.org\">The Burnham Institute</a>,<br>"<<endl;
	html<<"Contact: <A href=\"mailto:yye@burnham.org\">yye@burnham.org</a><br>"<<endl;
	html<<"</body>"<<endl;
	html<<"</html>"<<endl;
}

void BegTable(ostream &html)
{
	html<<"<table border=1 align=center>"<<endl;
}

void EndTable(ostream &html)
{
	html<<"</table>\n";
}

void BegRow(ostream &html)
{
	html<<"<tr>\n";
}

void EndRow(ostream &html)
{
	html<<"</tr>\n";
}

template <class Tdata>
void WriteACell(Tdata data, ofstream &html)
{
	html<<"  <td align=center>"<<data<<"</td>"<<endl;
}

void WriteACell(char *str, ofstream &html)
{
	html<<"  <td align=center>"<<str<<"</td>"<<endl;
}

void WriteACell(int span, char *str, ofstream &html)
{
	html<<"  <td colspan="<<span<<" align=center>"<<str<<"</td>"<<endl;
}

void WriteACell_Link(char *link, char *str, ofstream &html)
{
	html<<"  <td align=center><a href=\""<<link<<"\">"<<str<<"</a></td>"<<endl;
}

//especially useful for SCOP et al, where link = link0+link1
//eg "http://pfam.wustl.edu/cgi-bin/getdesc?name=RNA_pol_A"
//link = "http://pfam.wustl.edu/cgi-bin/getdesc?name=" and link1 = "RNA_pol_A"
void WriteACell_Link(char *link0, char *link1, char *str, ofstream &html)
{
	html<<"  <td align=center><a href=\""<<link0<<link1<<"\">"<<str<<"</a></td>"<<endl;
}

//link = link0+link1, especially useful for SCOP
void WriteACell_Link(char *link0, char *link1, char *str1, char *link2, char *str2, ofstream &html)
{
	html<<"  <td align=center><a href=\""<<link0<<link1<<"\">"<<str1<<"</a>"<<" -- <a href=\""<<link0<<link2<<"\">"<<str2<<"</a></td>"<<endl;
}

#endif
