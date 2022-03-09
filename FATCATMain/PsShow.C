#include <cstring>
#include <cstdlib>

#include "basic.h"
#include "PsShow.h"

PsShow::PsShow(char *filename)
{
	Initial(filename, "%!PS-Adobe-3.0","Results in postscript");
}

PsShow::PsShow(char *filename, char *title)
{
	Initial(filename, "%!PS-Adobe-3.0", title);
}

PsShow::PsShow(char *filename, char *head, char *title)
{
	Initial(filename, head, title); //the head has boundingbox information
}

void PsShow::Initial(char *filename, char *head, char *title)
{
        char    *str = new char[1000];
        out.open(filename);
        if(!out) {
                cout<<"open file error "<<filename<<endl;
                exit(0);
        }

	out<<head<<endl;
        time_t lt = time(NULL);
        out<<"%%Creator: Yuzhen Ye"<<endl;
        out<<"%%Title: "<<title<<endl;
        out<<"%%Creation date: "<<ctime(&lt);
        out<<"%%EndComments"<<endl;
        out<<"%%BeginProlog"<<endl;
        out<<"/M { moveto } def"<<endl;
        out<<"/m { moveto } def"<<endl;
        out<<"/L { lineto } def"<<endl;
        out<<"/l { lineto } def"<<endl;
	out<<"/R { setrgbcolor } bind def"<<endl;
	out<<"/G { setgray } bind def"<<endl;
	out<<"/H { sethsbcolor } bind def"<<endl;
        out<<"/LW { setlinewidth } bind def"<<endl;
	out<<"/D { 0 setdash 0 setlinecap } bind def"<<endl;
	out<<"/ND { [] 0 setdash 1 setlinecap } bind def"<<endl;
	out<<"/box {\n\ttranslate\n\tnewpath\n\t-0 -3 moveto\n\t-0 7 lineto\n\t6 7 lineto\n\t6 -3 lineto\n\tclosepath\n\tfill\n\t} def\n";
	out<<"/cross {\n\ttranslate\n\tnewpath\n\t-4 -4 moveto\n\t4 4 lineto\n\t-4 4 moveto\n\t4 -4 lineto\n\tstroke\n\t} def\n";
	out<<"/arrow {\n\tnewpath\n\t0 0 moveto\n\t-8 -8 lineto\n\t-4 -8 lineto\n\t-4 -16 lineto\n\t4 -16 lineto\n\t4 -8 lineto\n\t8 -8 lineto\n\t0 0 lineto\n\tclosepath\n\tfill\n\t} def\n";
	out<<"/uparrow {\n\ttranslate\n\tarrow\n\t} def\n";
	out<<"/downarrow {\n\ttranslate\n\t180 rotate\n\tarrow\n\t} def\n";
	out<<"/leftarrow {\n\ttranslate\n\t90 rotate\n\tarrow\n\t} def\n";
	out<<"/rightarrow {\n\ttranslate\n\t-90 rotate\n\tarrow\n\t} def\n";
	out<<"/upleftarrow {\n\ttranslate\n\t45 rotate\n\tarrow\n\t} def\n";
	out<<"/uprightarrow {\n\ttranslate\n\t-45 rotate\n\tarrow\n\t} def\n";
	out<<"/downleftarrow {\n\ttranslate\n\t135 rotate\n\tarrow\n\t} def\n";
	out<<"/downrightarrow {\n\ttranslate\n\t-135 rotate\n\tarrow\n\t} def\n";
        out<<"/color_char { gsave\n\tsetrgbcolor\n\tmoveto\n\tshow\n\tgrestore\n\t} def\n";
        out<<"/color_inv { gsave\n\tsetrgbcolor\n\ttranslate\n\tbox fill\n\tgrestore\n\tmoveto\n\tshow\n\t} def\n";
        out<<"/white_inv { gsave\n\tsetrgbcolor\n\ttranslate\n\tbox fill\n\tgrestore\n\tgsave\n\tsetrgbcolor\n\tmoveto\n\tshow\n\tgrestore\n\t} def\n";
        out<<"/FR { /Times-Roman findfont } bind def"<<endl;
        out<<"/TI { /Times-Italic findfont } bind def"<<endl;
        out<<"/TB { /Times-Bold findfont } bind def"<<endl;
        out<<"/TBI { /Times-BoldItalic findfont } bind def"<<endl;
        out<<"/FC { /Courier findfont } bind def"<<endl;
        out<<"/CI { /Courier-Italic findfont } bind def"<<endl;
        out<<"/CB { /Courier-Bold findfont } bind def"<<endl;
        out<<"/CBI { /Courier-BoldItalic findfont } bind def"<<endl;
        out<<"/FG { /Symbol findfont } bind def"<<endl;
        out<<"/PR { scalefont setfont show } bind def"<<endl;
        out<<"/SR { scalefont setfont } bind def"<<endl;
        out<<"/SC { 1.0 G } def"<<endl;
        out<<"/SCR { 1.0 0 0 R } def"<<endl;
        out<<"/circle { newpath 0 360 arc gsave SC grestore stroke } bind def"<<endl;
        out<<"/SS { newpath 0 360 arc gsave SC grestore stroke } bind def"<<endl;
        out<<"/SSF { newpath 0 360 arc gsave SC fill grestore stroke } bind def"<<endl;
        out<<"/SSR { newpath 0 360 arc gsave SCR fill grestore stroke } bind def"<<endl;
        out<<"/C {"<<endl;
        out<<"  1 index scalefont setfont"<<endl;
        out<<"  exch stringwidth pop -2 div exch -3 div rmoveto"<<endl;
        out<<" } bind def"<<endl;
        out<<"/CR90 {"<<endl;
        out<<"  1 index scalefont setfont"<<endl;
        out<<"  exch stringwidth pop -2 div exch 3 div exch rmoveto"<<endl;
        out<<" } bind def"<<endl;
        out<<"%%EndProlog"<<endl;
        out<<"%%BeginSetup"<<endl;
        out<<"%%IncludeResource: font Times-Roman"<<endl;
        out<<"%%IncludeResource: font Times-Italic"<<endl;
        out<<"%%IncludeResource: font Times-Bold"<<endl;
        out<<"%%IncludeResource: font Times-BoldItalic"<<endl;
        out<<"%%IncludeResource: font Courier"<<endl;
        out<<"%%IncludeResource: font Courier-Italic"<<endl;
        out<<"%%IncludeResource: font Courier-Bold"<<endl;
        out<<"%%IncludeResource: font Courier-BoldItalic"<<endl;
        out<<"%%IncludeResource: font Symbol"<<endl;
        out<<"1 setlinecap 1 setlinejoin 1 LW 0 G"<<endl;
        out<<"%%EndSetup"<<endl;

	out<<"%%Page 1"<<endl;
        out<<"gsave"<<endl;

        lineWidth = 0.5;
        lineGray = 0.0;

        int     i, j;
        itemDescript = new char* [MAXITEM];
        itemColor = new double* [MAXITEM];
        for(i = 0; i < MAXITEM; i ++)   {
                itemDescript[i] = new char[100];
                itemDescript[i][0] = 0;
                itemColor[i] = new double[3];
        }

        RGBValue = new double* [MAXCOLOR];
        for(i = 0; i < MAXCOLOR; i ++)    {
                RGBValue[i] = new double[3];
                for(j = 0; j < 3; j ++) RGBValue[i][j] = 0.0;
        }
        RGBValue[0][0] = 0.0;
        RGBValue[0][1] = 0.0;
        RGBValue[0][2] = 0.0;
        int     add = 1;
        for(i = 0; i < 3; i ++) {
                for(j = 0; j < 3; j ++) {
                        if(i == j) RGBValue[add + i][j] = 1.0;
                        else    RGBValue[add + i][j] = 0.0;
                }
        }
        add = 4;
        for(i = 0; i < 3; i ++) {
                for(j = 0; j < 3; j ++) {
                        if(i == j) RGBValue[add + i][j] = 0.0;
                        else    RGBValue[add + i][j] = 1.0;
                }
        }
        add = 7;
        for(i = 0; i < 3; i ++) {
                for(j = 0; j < 3; j ++) {
                        if(i == j) RGBValue[add + i][j] = 1.0;
                        else    RGBValue[add + i][j] = 0.2 * j;
                }
        }
        add = 10;
        totalColor = add;

        itemNum = 0;
        scaleX = 1.0;
        scaleY = 1.0;
        translateX = 0.0;
        translateY = 0.0;
        markFontSize = 8;
        textFontSize = 10;
        titleFontSize = 12;
        textAngle = 0;
        ifpointSymble = 0;
        pointSymble[0] = 0;
	symbleRad = 2.0;
	ifPrtLine = 1;
        nowColorIndex = 0;
        itemNum = 0;

        fonttypenum = 9;
        fontDescription = new char*[fonttypenum];
        for(i = 0; i < fonttypenum; i ++)      {
                fontDescription[i] = new char[10];
        }
        strcpy(fontDescription[0], "FR"); //Times
        strcpy(fontDescription[1], "TB"); //Times-Bold
        strcpy(fontDescription[2], "TI"); //Times-Italic
        strcpy(fontDescription[3], "TBI"); //Times-BoldItalic
        strcpy(fontDescription[4], "FC"); //Courier
        strcpy(fontDescription[5], "CB"); //Courier-Bold
        strcpy(fontDescription[6], "CI"); //Courier-Italic
        strcpy(fontDescription[7], "CBI"); //Courier-BoldItalic
        strcpy(fontDescription[8], "FG"); //Symblom

        axisXMin = 80;
        axisXMax = 500;
        axisYMin = 250;
        axisYMax = 700; //default page, small
	axisYMaxReal = axisYMax;
	axisYMinReal = axisYMin;
        pointxTitle = (axisXMax - axisXMin) / 2.0 + axisXMin;
        pointyTitle = axisYMin - 60;
	xAxisDescript = 1;
	yAxisDescript = 1;
	ifDrawAxisBox = 1;
	upsidedown = 0;
        delete[] str;
}


PsShow::~PsShow()
{
        out<<"grestore"<<endl;
        out<<"showpage"<<endl;
        out<<"%%Trailer"<<endl;
        out.close();
        int     i;
        for(i = 0; i < MAXITEM; i ++)   {
                delete[] itemDescript[i];
                delete[] itemColor[i];
        }
        delete[] itemDescript;
        delete[] itemColor;
        for(i = 0; i < MAXCOLOR; i ++)  {
                delete[] RGBValue[i];
        }
        for(i = 0; i < fonttypenum; i ++)       {
                delete[] fontDescription[i];
        }
        delete[] fontDescription;
        fonttypenum = 0;
        delete[] RGBValue;
}

void	PsShow::xAxisDescriptOff(void)
{
	xAxisDescript = 0;
}

void	PsShow::yAxisDescriptOff(void)
{
	yAxisDescript = 0;
}

void 	PsShow::setupsidedown(void)
{
	upsidedown = 1;
}

void 	PsShow::setMarkFontSize(double size)
{
	markFontSize = size;
}

void 	PsShow::setTextFontSize(double size)
{
	textFontSize = size;
}

void    PsShow::drawAxis(double **point, int pointNum, char *x_lab, char *y_lab)
{
        drawAxis(point, pointNum);
        double  xmiu0 = (axisXMax - axisXMin) / 2.0;
        double  xmiu1 = 35;
        double  ymiu0 = (axisYMax - axisYMin) / 2.0;
        double  ymiu1 = 20;
        if(upsidedown)	prtText(axisXMax - xmiu0, axisYMax + ymiu1, x_lab, markFontSize, 1);
	else		prtText(axisXMax - xmiu0, axisYMin - ymiu1, x_lab, markFontSize, 1);
        textAngle = 90;
        prtText(axisXMin - xmiu1, axisYMax - ymiu0, y_lab, markFontSize, 1);
        textAngle = 0;
}

void    PsShow::drawAxis(double **point, int pointNum)
{
        char str[200];
        sprintf(str, "%.1f LW %.1f G\n", lineWidth, lineGray);
        out<<str;

        double  minX, maxX, minY, maxY;
        minX = maxX = point[0][0];
        minY = maxY = point[0][1];
        int     i;
        for(i = 1; i < pointNum; i ++)  {
                minX = (minX < point[i][0])?minX:point[i][0];
                maxX = (maxX > point[i][0])?maxX:point[i][0];
                minY = (minY < point[i][1])?minY:point[i][1];
                maxY = (maxY > point[i][1])?maxY:point[i][1];
        }
        if(minX == maxX)        minX = maxX - 0.1;
        if(minY == maxY)        minY = maxY - 0.1;
	//double	edg = 0.9;
        //scaleX = (axisXMax - axisXMin) * edg / (maxX - minX);
        //scaleY = (axisYMax - axisYMin) * edg / (maxY - minY);
        //translateX = axisXMin - minX * scaleX + (axisXMax - axisXMin) * (1.0 - edg) * 0.5; 
        //translateY = axisYMin - minY * scaleY + (axisYMax - axisYMin) * (1.0 - edg) * 0.5;
                //transform information

	drawAxis(minX, maxX, minY, maxY);
}

void    PsShow::drawAxis(double minX, double maxX, double minY, double maxY)
{
	drawAxisBox();
	double	edg = 0.9;
        scaleX = (axisXMax - axisXMin) * edg / (maxX - minX);
        scaleY = (axisYMax - axisYMin) * edg / (maxY - minY);
        translateX = axisXMin - minX * scaleX + (axisXMax - axisXMin) * (1.0 - edg) * 0.5; 
        translateY = axisYMin - minY * scaleY + (axisYMax - axisYMin) * (1.0 - edg) * 0.5;
                //transform information

	int	i, j;
        double  v1, v2, v3;
        double  value;
        double  **xypoint;
        xypoint = new double* [40];
        for(i = 0; i < 40; i ++)        xypoint[i] = new double[2];
        double  **markpot;
        char    **xypointDescript;
        markpot = new double* [40];
        xypointDescript = new char* [40];
        for(i = 0; i < 40; i ++)        {
                markpot[i] = new double[2];
                xypointDescript[i] = new char[20];
        }
        int     xypointNum;
        double  beginvalue;

	//begin x axis
        v1 = (maxX - minX) / 10.0;
        digitDiv(v1, &v2, &v3);

	double	stepx;
	if(v2 > 5.0)	stepx = 10 * v3;
	else if(v2 > 2.0)	stepx = 5 * v3;
	else	stepx = 2 * v3;

	i = 0;
	if(minX >= 0)	{
		while(minX > i * stepx)	{
			i ++;
		}
		beginvalue = ((i - 1)>0?(i-1):0) * stepx;
	}
	else	{
		while(minX < -1 * i * stepx)	{
			i ++;
		}	
		beginvalue = i * stepx * -1;
	}
	//Y.Y,Mar 9,2004
	/*
	for(i = 0; i <= 10; i ++)	{
		if(minX > i * 10 * v3)	break; 	
	}
	//beginvalue = (i - 1) > 0?(i - 1):0; 
	for(i = -20; i <= 20; i ++)	{
		if(minX < i * 10 * v3)	break; 	
	}
	beginvalue = (i - 1) * 10 * v3; 
	*/

	//printf("maxX %f, minX %f, stepx %f, beginvalue %f\n", maxX, minX, stepx, beginvalue);

        j = 0;
        for(i = 0; i < 20; i ++)        {
                value = stepx * i + beginvalue;
		//printf("x-mark %d, %f\n", i, value);
                xypointDescript[j][0] = 0;
		sprintf(xypointDescript[j], "%f", value);
		remove0(xypointDescript[j]);
		//if(stepx >= 10)	sprintf(xypointDescript[j], "%.0f", value);
		//else if(stepx >= 1)	sprintf(xypointDescript[j], "%.1f", value);
                //else	sprintf(xypointDescript[j], "%.2f", value);
                value = value * scaleX + translateX;
                if(value < axisXMin)    continue;
                if(value > axisXMax)    break;
                markpot[j][0] = value;
                xypoint[2 * j][0] = value;
                xypoint[2 * j + 1][0] = value;
		if(upsidedown)	{
                	markpot[j][1] = axisYMax + 10;
	                xypoint[2 * j][1] = axisYMax + 2;
		        xypoint[2 * j + 1][1] = axisYMax;
		}
		else	{
                	markpot[j][1] = axisYMin - 10;
                	xypoint[2 * j][1] = axisYMin - 2;
	                xypoint[2 * j + 1][1] = axisYMin;
		}
                j ++;
        }
        xypointNum = j;
	
        int     ifmidalign = 1;
	if(xAxisDescript)	{
        	prtStepLine(xypoint, xypointNum * 2);
        	defaultColor();
        	prtText(markpot, xypointDescript, xypointNum, markFontSize, ifmidalign);
	}
                //x-axis
	//end of x-axis

	//begin y-axis
        v1 = (maxY - minY) / 10.0;
        digitDiv(v1, &v2, &v3);

	double	stepy;
	if(v2 > 5.0)	stepy = 10 * v3;
	else if(v2 > 2.0)	stepy = 5 * v3;
	else	stepy = 2 * v3;

	//printf("y-axis: v1 %f v2 %f, v3 %f, maxY %f, minY %f, stepy %f\n", v1, v2, v3, maxY, minY, stepy);

	i = 0;
	if(minY >= 0)	{
		while(minY > i * stepy)	{
			i ++;
		}
		beginvalue = ((i - 1)>0?(i-1):0) * stepy;
	}
	else	{
		while(minY < -1 * i * stepy)	{
			i ++;
		}	
		beginvalue = i * stepy * -1;
	}
	/*
	if(minY < 0)	{
		for(i = 0; i <= 1000; i ++)	{
			if(minY > -1 * i * 10 * v3)	break; 	
		}
		beginvalue = (-i - 1) * 10 * v3; 
	}
	else	{
		for(i = 0; i <= 1000; i ++)	{
			if(minY >= i * 10 * v3)	break; 	 //> replaced by >=, YY.7/15/03
		}
		beginvalue = ((i - 1) > 0?(i - 1):0) * 10 * v3; 
	}
	for(i = -20; i <= 20; i ++)	{
		if(minY < i * 10 * v3)	break; 	
	}
	beginvalue = (i - 1) * 10 * v3; 
	*/

        j = 0;
        for(i = 0; i < 20; i ++)        {
                value = stepy * i + beginvalue;
		//printf("y-mark %d, %f\n", i, value);
                xypointDescript[j][0] = 0;
		sprintf(xypointDescript[j], "%f", value);
		remove0(xypointDescript[j]);
                //if(stepy >= 10)	sprintf(xypointDescript[j], "%6.0f", value);
		//else if(stepy >= 1)	sprintf(xypointDescript[j], "%6.1f", value);
		//else	sprintf(xypointDescript[j], "%6.2f", value);
                //sprintf(xypointDescript[j], "%6.2f", value);
                //sprintf(xypointDescript[j], "%6.2f", value);
                value = value * scaleY + translateY;
		//printf("i %d, value %f\n", i, value);
                if(value < axisYMin)    continue;
                if(value > axisYMax)    break;
		if(upsidedown)	{
			value = axisYMax + axisYMin - value;	
		}	
                markpot[j][0] = axisXMin - 25;
                markpot[j][1] = value;
                xypoint[2 * j][0] = axisXMin;
                xypoint[2 * j][1] = value;
                xypoint[2 * j + 1][0] = axisXMin - 2;
                xypoint[2 * j + 1][1] = value;
                j ++;
        }
        xypointNum = j;
	if(yAxisDescript)	{
        	prtStepLine(xypoint, xypointNum * 2); 
		ifmidalign = 0;
        	prtText(markpot, xypointDescript, xypointNum, markFontSize, ifmidalign);
	}	//Y-axis
        for(i = 0; i < 40; i ++)        delete[] xypoint[i];
        delete[] xypoint;
        for(i = 0; i < 40; i ++)        {
                delete[] markpot[i];
                delete[] xypointDescript[i];
        }
        delete[] markpot;
        delete[] xypointDescript;

}

void PsShow::remove0(char *str)
{
	int	i, j;
	int	len = strlen(str);
	for(i = 0; i < len; i ++) {
		if(str[i] == '.')	{
			for(j = len - 1; j > i; j --)	{
				if(str[j] != '0')	break;
			}
			if(str[j] == '.')	str[j] = '\0';
			else	str[j + 1] = '\0';
			return;
		}
	}	
}

void PsShow::defaultPage(void)
{
        axisXMin = 100;
        axisXMax = 500;
        axisYMin = 150;
        axisYMax = 700; //default page, of middle size
        pointxTitle = (axisXMax - axisXMin) / 2.0 + axisXMin;
        pointyTitle = axisYMin - 60;
	axisYMaxReal = axisYMax;
	axisYMinReal = axisYMin;
}

void PsShow::PageArea(double xwide, double ywide)
{
        axisXMin = 300 - xwide / 2.0;
        axisXMax = 300 + xwide / 2.0;
        axisYMin = 400 - ywide / 2.0;
        axisYMax = 400 + ywide / 2.0;
        pointxTitle = (axisXMax - axisXMin) / 2.0 + axisXMin;
        pointyTitle = axisYMin - 40;
	axisYMaxReal = axisYMax;
	axisYMinReal = axisYMin;
}
void PsShow::PageArea(double xmin, double xwide, double ymin, double ywide)
{
        axisXMin = xmin;
        axisXMax = xmin + xwide;
        axisYMin = ymin;
        axisYMax = ymin + ywide;
        pointxTitle = (axisXMax - axisXMin) / 2.0 + axisXMin;
        pointyTitle = axisYMin - 40;
	axisYMaxReal = axisYMax;
	axisYMinReal = axisYMin;
}

void PsShow::largePage(void)
{
        axisXMin = 70;
        axisXMax = 500;
        axisYMin = 120;
        axisYMax = 800; //default page, small
        pointxTitle = (axisXMax - axisXMin) / 2.0 + axisXMin;
        pointyTitle = axisYMin - 40;
	axisYMaxReal = axisYMax;
	axisYMinReal = axisYMin;
}

void PsShow::smallPage(void)
{
        axisXMin = 120;
        axisXMax = 450;
        axisYMin = 250;
        axisYMax = 650; //default page, small
        pointxTitle = (axisXMax - axisXMin) / 2.0 + axisXMin;
        pointyTitle = axisYMin - 60;
	axisYMaxReal = axisYMax;
	axisYMinReal = axisYMin;
}

void  	PsShow::drawAxisBoxOff(void)
{
	ifDrawAxisBox = 0;
}

void PsShow::drawAxisBox(void)
{
	if(ifDrawAxisBox == 0)	return;

        char    str[200];
        str[0] = 0;

        out<<"newpath"<<endl;
        sprintf(str, "%.1f %.1f M %.1f %.1f L\n",
                axisXMin, axisYMin, axisXMax, axisYMin);
        out<<str;
        sprintf(str, "%.1f %.1f M %.1f %.1f L\n",
                axisXMax, axisYMin, axisXMax, axisYMax);
        out<<str;
        sprintf(str, "%.1f %.1f M %.1f %.1f L\n",
                axisXMax, axisYMax, axisXMin, axisYMax);
        out<<str;
        sprintf(str, "%.1f %.1f M %.1f %.1f L\n",
                axisXMin, axisYMax, axisXMin, axisYMin);
        out<<str<<"stroke"<<endl;
                //xy-box
}

void    PsShow::digitDiv(double v1, double *v2, double *v3)
{
        char    s1[50], s2[50], s3[50], s4[50];
        int     i, slen;
        s1[0] = s2[0] = s3[0] = s4[0] = 0;
        sprintf(s1, "%e", v1);
        slen = strlen(s1);
        for(i = 0; i < slen; i ++)        {
                if(s1[i] == 'e')        break;
        }
        strncpy(s2, s1, slen - i);
        s2[slen - i] = '\0';
        strcpy(s3, s1 + i);
        strcpy(s4, "1.0");
        strcat(s4, s3);
        double     v2t, v3t;
        sscanf(s2, "%lf", &v2t);
        sscanf(s4, "%lf", &v3t);
        *v2 = v2t;
        *v3 = v3t;
}

void    PsShow::codTransform(double **point, double pointNum)
{
        int     i;
        for(i = 0; i < pointNum; i ++)  {
                point[i][0] = point[i][0] * scaleX + translateX; //x-axis
                point[i][1] = point[i][1] * scaleY + translateY; //y-axis
        }
}

void    PsShow::codTransform(int ifcodtransform, double **point, double pointNum)
{
        int     i;
	if(ifcodtransform)	{
	        for(i = 0; i < pointNum; i ++)  {
       	         	point[i][0] = point[i][0] * scaleX + translateX; //x-axis
                	point[i][1] = point[i][1] * scaleY + translateY; //y-axis
		}
        }
	if(upsidedown)	{
		for(i = 0; i < pointNum; i ++)	{
			point[i][1] = axisYMaxReal + axisYMinReal - point[i][1];
		}
	}
}

void    PsShow::setTransform(void)
{
	axisYMaxReal = (axisYMax - translateY) / scaleY; 
	axisYMinReal = (axisYMin - translateY) / scaleY; 
	if(upsidedown)	{
		out<<"%upsidedown: y-value "<<axisYMaxReal + axisYMinReal<<" - real-value"<<endl;
	}
	out<<"gsave\n";
	setTranslate(translateX, translateY);
	setScale(scaleX, scaleY);
}

void    PsShow::setScale(double scalex, double scaley)
{
        char    str[200];
        str[0] = 0;
        sprintf(str, "%2.1f %2.1f scale\n", scalex, scaley);
        out<<str;
}

void    PsShow::setTranslate(double translatex, double translatey)
{
        char    str[200];
        str[0] = 0;
        sprintf(str, "%.2f %.2f translate\n", translatex, translatey);
        out<<str;
}

void    PsShow::gsave(void)
{
        out<<"gsave"<<endl;
}

void    PsShow::grestore(void)
{
        out<<"grestore"<<endl;
}

void    PsShow::creatColor(char *lab)
{
        int     i;
        strcpy(itemDescript[itemNum], lab);
        for(i = 0; i < 3; i ++)
                itemColor[itemNum][i] = RGBValue[nowColorIndex][i];
        itemNum ++;
        char    str[200];
        str[0] = 0;
        sprintf(str, "%.2f %.2f %.2f R\n", RGBValue[nowColorIndex][0],
                RGBValue[nowColorIndex][1], RGBValue[nowColorIndex][2]);
        out<<str;
        nowColorIndex ++;
        if(nowColorIndex > totalColor)    nowColorIndex = 0;
}

void    PsShow::assignColor(double  *value)
{
        char    str[200];
        str[0] = 0;
        sprintf(str, "%.2f %.2f %.2f R\n", value[0], value[1], value[2]);
        out<<str;
}

void    PsShow::defaultColor(void)
{
        char    str[200];
        str[0] = 0;
        sprintf(str, "%.2f %.2f %.2f R\n", 0.0, 0.0, 0.0);
        out<<str;
}

void    PsShow::prtContLine(double **point, int pointNum, int ifCodTransform, double linewidth)
{
	char	str[200];
	out<<"gsave"<<endl;
        sprintf(str, "%.1f LW %.1f G\n", linewidth, lineGray);
        out<<str;
	prtContLine(point, pointNum, ifCodTransform);
	out<<"grestore"<<endl;

}

void    PsShow::prtContLine(double **point, int pointNum, int ifCodTransform)
{
        int     i, j;
        char    str[200];
        double  **pointNew;
        pointNew = new double* [pointNum];
        for(i = 0; i < pointNum; i ++)  {
                pointNew[i] = new double[2];
                for(j = 0; j < 2; j ++) {
                        pointNew[i][j] = point[i][j];
                }
        }
	codTransform(ifCodTransform, pointNew, pointNum);
	if(ifPrtLine)	{
        	out<<"newpath"<<endl;
        	str[0] = 0;
        	out<<str;
        	for( i = 0; i < pointNum - 1; i ++ )    {
               	 	sprintf(str, "%.2f %.2f M %.2f %.2f L\n",
                        	pointNew[i][0], pointNew[i][1], pointNew[i + 1][0], pointNew[i + 1][1]);
                	out<<str;
        	}
        	sprintf(str, "stroke\n");
        	out<<str;
	}
        if(ifpointSymble == 1)   {
		if(!strcmp(pointSymble, "cross"))	{	
			out<<"0 G"<<endl;
			out<<"newpath"<<endl;
		}
                else	out<<"/SC {0.0 G} def"<<endl;
                for(i = 0; i < pointNum; i ++)  {
			if(!strcmp(pointSymble, "fillcircle"))
                        	sprintf(str, "%.2f %.2f %.1f SSF\n", pointNew[i][0], pointNew[i][1], symbleRad);
			else if(!strcmp(pointSymble, "redcircle"))	
                        	sprintf(str, "%.2f %.2f %.1f SSR\n", pointNew[i][0], pointNew[i][1], symbleRad);
			else if(!strcmp(pointSymble, "square"))
                        	sprintf(str, "%.2f %.2f %.1f box\n", pointNew[i][0], pointNew[i][1], symbleRad);
			else if(!strcmp(pointSymble, "circle"))	
                        	sprintf(str, "%.2f %.2f %.1f SS\n", pointNew[i][0], pointNew[i][1], symbleRad);
			else if(!strcmp(pointSymble, "cross"))	
                        	sprintf(str, "%.2f %.2f M %.2f %.2f L %.2f %.2f M %.2f %.2f L\n", 
					pointNew[i][0] - symbleRad, pointNew[i][1], 
					pointNew[i][0] + symbleRad, pointNew[i][1],
					pointNew[i][0], pointNew[i][1] - symbleRad, 
					pointNew[i][0], pointNew[i][1] + symbleRad);
                        out<<str;
                }
		if(!strcmp(pointSymble, "cross"))	out<<"stroke"<<endl;
        }
        for(i = 0; i < pointNum; i ++) delete[] pointNew[i];
        delete[] pointNew;
}

void    PsShow::prtContLine(double **point, int pointNum, int ifCodTransform, char *symble)
{
	ifpointSymble = 1;
	strcpy(pointSymble, symble);
        prtContLine(point, pointNum, ifCodTransform);
        ifpointSymble = 0;
}

void    PsShow::prtContLine(double **point, int pointNum, int ifCodTransform, char *symble, int ifline)
{
	ifPrtLine = ifline;
	prtContLine(point, pointNum, ifCodTransform, symble);
}

//expecially useful for axis drawing (rules)
void    PsShow::prtStepLine(double **point, int pointNum)
{
        int     i;
        char    str[200];
        out<<"newpath"<<endl;
        str[0] = 0;
        out<<str;
        for( i = 0; i < pointNum - 1; i += 2 )    {
                sprintf(str, "%.2f %.2f M %.2f %.2f L\n",
                        point[i][0], point[i][1], point[i + 1][0], point[i + 1][1]);
                out<<str;
        	sprintf(str, "stroke\n");
        }
        out<<str;
}

void    PsShow::prtStepLine(double **point, int pointNum, int ifCodTransform)
{
        int     i, j;
        char    str[200];
        double  **pointNew;
        pointNew = new double* [pointNum];
        for(i = 0; i < pointNum; i ++)  {
                pointNew[i] = new double[2];
                for(j = 0; j < 2; j ++) {
                        pointNew[i][j] = point[i][j];
                }
        }
        codTransform(ifCodTransform, pointNew, pointNum);
        out<<"newpath"<<endl;
        str[0] = 0;
        out<<str;
        for( i = 0; i < pointNum - 1; i += 2 )    {
                sprintf(str, "%.2f %.2f M %.2f %.2f L\n",
                        pointNew[i][0], pointNew[i][1], pointNew[i + 1][0], pointNew[i + 1][1]);
                out<<str;
        	sprintf(str, "stroke\n");
        }
        sprintf(str, "stroke\n");
        out<<str;
        for(i = 0; i < pointNum; i ++) delete[] pointNew[i];
        delete[] pointNew;
}

void    PsShow::prtGrayStepLine(double **point, int pointNum, double *gray, double *linewidth, int ifCodTransform)
{
        int     i, j;
        char    str[200];
        double  **pointNew;
        pointNew = new double* [pointNum];
        for(i = 0; i < pointNum; i ++)  {
                pointNew[i] = new double[2];
                for(j = 0; j < 2; j ++) {
                        pointNew[i][j] = point[i][j];
                }
        }
        codTransform(ifCodTransform, pointNew, pointNum);
        out<<"newpath"<<endl;
        str[0] = 0;
        out<<str;
        for( i = 0; i < pointNum - 1; i += 2 )    {
                sprintf(str, "%.1f LW %.2f G\n", linewidth[i], gray[i]);
                out<<str;
                sprintf(str, "%.2f %.2f M %.2f %.2f L\n",
                        pointNew[i][0], pointNew[i][1], pointNew[i + 1][0], pointNew[i + 1][1]);
                out<<str;
		out<<"stroke"<<endl;
        }
        sprintf(str, "stroke\n");
        out<<str;
        for(i = 0; i < pointNum; i ++) delete[] pointNew[i];
        delete[] pointNew;
}

void    PsShow::prtColorStepLine(double **point, int pointNum, double **color, int ifCodTransform)
{
        int     i, j;
        char    str[200];
        double  **pointNew;
        pointNew = new double* [pointNum];
        for(i = 0; i < pointNum; i ++)  {
                pointNew[i] = new double[2];
                for(j = 0; j < 2; j ++) {
                        pointNew[i][j] = point[i][j];
                }
        }
        codTransform(ifCodTransform, pointNew, pointNum);
        out<<"newpath"<<endl;
        str[0] = 0;
        out<<str;
        for( i = 0; i < pointNum - 1; i += 2 )    {
                sprintf(str, "%.2f %.2f %.2f R\n", color[i][0], color[i][1], color[i][2]);
                out<<str;
                sprintf(str, "%.2f %.2f M %.2f %.2f L\n",
                        pointNew[i][0], pointNew[i][1], pointNew[i + 1][0], pointNew[i + 1][1]);
                out<<str;
		out<<"stroke"<<endl;
        }
        sprintf(str, "stroke\n");
        out<<str;
        for(i = 0; i < pointNum; i ++) delete[] pointNew[i];
        delete[] pointNew;
}

void    PsShow::prtColorStepLine(double **point, int pointNum, double **color, double *linewidth, int ifCodTransform)
{
        int     i, j;
        char    str[200];
        double  **pointNew;
        pointNew = new double* [pointNum];
        for(i = 0; i < pointNum; i ++)  {
                pointNew[i] = new double[2];
                for(j = 0; j < 2; j ++) {
                        pointNew[i][j] = point[i][j];
                }
        }
        codTransform(ifCodTransform, pointNew, pointNum);
        out<<"newpath"<<endl;
        str[0] = 0;
        out<<str;
        for( i = 0; i < pointNum - 1; i += 2 )    {
                sprintf(str, "%.2f %.2f %.2f R %.1f LW\n", color[i][0], color[i][1], color[i][2], linewidth[i]);
                out<<str;
                sprintf(str, "%.2f %.2f M %.2f %.2f L\n",
                        pointNew[i][0], pointNew[i][1], pointNew[i + 1][0], pointNew[i + 1][1]);
                out<<str;
		out<<"stroke"<<endl;
        }
        sprintf(str, "stroke\n");
        out<<str;
        for(i = 0; i < pointNum; i ++) delete[] pointNew[i];
        delete[] pointNew;
}

void    PsShow::prtZhuZhuangTu(int itemnum, int linelen, double **yvalue, double *yvaluemin, double *yvaluemax, char ***xdescription, char **ydescription, int **iffill, double *gray)
{
        int     i, j, beg, end, page;
        double  yWide, xWide, ZhuXWide, ZhuYWide, xtmp, ytmp, ytmp0;
        double  **axisxcod;
        axisxcod = new double*[linelen];
        for(i = 0; i < linelen; i ++)   axisxcod[i] = new double[2];
        char    *str = new char[200];
        int     xstep = 40;
        int     ystep = 6;
        double  ysize;
        char    **text;
        text = new char*[ystep];
        for(i = 0; i < ystep; i ++)     text[i] = new char[10];

        beg = 0;
        if(linelen > xstep + 15)        {
                end = xstep;
        }
        else    {
                end = linelen;
        }
        xWide = (axisXMax - axisXMin) / double(end);
        ZhuXWide = xWide / 2.0;
        ZhuYWide = 8 * ZhuXWide;
        page = 1;
        yWide = (axisYMax - ZhuYWide - 50);
        while(beg < linelen)     {
                if(beg != 0)    {
                        sprintf(str, "%%%%Page %d\ngsave\n", page);
                        out<<str;
                }
                for(i = 0; i < itemnum; i ++)   {
                        out<<"1.0 LW 0 G"<<endl;
                        out<<"newpath"<<endl;
                        xWide = axisXMin;
                        ysize = (yvaluemax[i] - yvaluemin[i]) / double(ystep - 1);
                        for(j = 0; j < ystep; j ++)     {
                                ytmp0 = j * ysize + yvaluemin[i];
                                if(ytmp0 > yvaluemax[i])   break;
                                ytmp = ZhuYWide * (ytmp0 / yvaluemax[i]);
                                sprintf(str, "%.2f %.2f M %.2f %.2f L\n", xWide, yWide + ytmp, xWide - 2, yWide + ytmp);
                                out<<str;
                                axisxcod[j][0] = xWide - 17;
                                axisxcod[j][1] = yWide + ytmp;
                                sprintf(text[j], "%.2f", ytmp0);
                        }
                        out<<"stroke"<<endl;
                        prtText(axisxcod, text, j, 8, 0);
                        ytmp0 = double(j) * double(ysize)/2.0 + yvaluemin[i];
                        ytmp = yWide + ZhuYWide * (ytmp0 / yvaluemax[i]);
                        xtmp = xWide - 25;
                        textAngle = 90;
                        prtText(xtmp, ytmp, ydescription[i], 8, 1);
                        textAngle = 0;
                        xWide += ZhuXWide;
                        for(j = beg; j < end; j ++)     {
                                ytmp = ZhuYWide * (yvalue[i][j] / yvaluemax[i]);
                                if(iffill[i][j])        {
                                        sprintf(str, "%.2f G\n", gray[iffill[i][j] - 1]);
                                        out<<str;
                                }
                                sprintf(str, "newpath %.2f %.2f M %.2f %.2f L %.2f %.2f L %.2f %.2f L ",
                                        xWide, yWide, xWide, yWide + ytmp, xWide + ZhuXWide, yWide + ytmp, xWide + ZhuXWide, yWide);
                                out<<str;
                                if(iffill[i][j])        out<<"closepath fill"<<endl;
                                else    out<<"stroke"<<endl;
                                if(iffill[i][j])        out<<"0 G"<<endl;
                                axisxcod[j][0] = xWide + 0.5 * ZhuXWide;
                                axisxcod[j][1] = yWide - 15.0;
                                xWide += 2 * ZhuXWide;
                        }
                        sprintf(str, "%.2f %.2f M %.2f %.2f L\n",
                                axisXMin, yWide, axisXMax, yWide);
                        out<<str;
                        sprintf(str, "%.2f %.2f M %.2f %.2f L\n",
                                axisXMin, yWide, axisXMin, yWide + ZhuYWide);
                        out<<str;
                        yWide -= ZhuYWide * 1.8;
                        out<<"stroke"<<endl;
                        textAngle = 90;
                        prtText(axisxcod + beg, xdescription[i] + beg, end - beg, 8, 1);
                        textAngle = 0;
                }
                beg = end;
                end = (beg + xstep)<linelen?(beg + xstep):linelen;
                if(beg < linelen)       {
                        sprintf(str, "showpage\ngrestore\n");
                        out<<str;
                }
                page ++;
        }
        pointyTitle = yWide + 20;
        for(i = 0; i < linelen; i ++)   delete[] axisxcod[i];
        delete[] axisxcod;
        delete[] str;
        for(i = 0; i < ystep; i ++)     delete text[i];
        delete[] text;
}

void    PsShow::prtDescription(void)
{
        int     i;
        double  **point;
        point = new double* [2];
        for(i = 0; i < 2; i ++) point[i] = new double[2];
        double  yvalue = (axisYMax - axisYMin) * 0.7 + axisYMin;
        double  ystep = 10.0;
        int     iftransform = 0;
        double  **markPoint;
        markPoint = new double* [itemNum];
        for(i = 0; i < itemNum; i ++)   markPoint[i] = new double[2];
        point[0][0] = (axisXMax - axisXMin) * 0.7 + axisXMin ;
        point[1][0] = point[0][0] + 15;
        for(i = 0; i < itemNum; i ++)   {
                assignColor(itemColor[i]);
                point[0][1] = point[1][1] = yvalue - i * ystep;
                prtContLine(point, 2, iftransform);
                markPoint[i][0] = point[1][0] + 5;
                markPoint[i][1] = point[1][1];
        }
        defaultColor();
        int     ifmidalign = 0;
        prtText(markPoint, itemDescript, itemNum, textFontSize, ifmidalign);
        for(i = 0; i < itemNum; i ++)   delete[] markPoint[i];
        delete[] markPoint;

        for(i = 0; i < 2; i ++) delete[] point[i];
        delete[] point;
}

void    PsShow::prtDescription(double toMaxX, double toMaxY, char *lab)
{
        prtText(axisXMax - toMaxX, axisYMax - toMaxY, lab, 10, 1);
}

void    PsShow::prtText(double positionX, double positionY, char *text, double fontSize, int ifmidalign)
{
        char    str[200];
        out<<"gsave"<<endl;
        sprintf(str, "%.2f %.2f M\n", positionX, positionY);
        out<<str;
        if(textAngle != 0)  {
                out<<"gsave"<<endl;
                sprintf(str, "%d rotate\n", textAngle);
                out<<str;
        }
        if(ifmidalign)  {
                sprintf(str, "(%s) %.2f FR C\n", text, fontSize);
                out<<str;
        }
        sprintf(str, "(%s) FR %.2f PR\n", text, fontSize);
        out<<str;
        if(textAngle != 0)  {
                out<<"grestore"<<endl;
        }
        out<<"grestore"<<endl;
}

void    PsShow::prtText(double **position, char **text, int textNum, double fontSize, int ifmidalign)
{
        char    str[200];
        out<<"gsave"<<endl;
        for(int i = 0; i < textNum; i ++)   {
                sprintf(str, "%.2f %.2f M\n", position[i][0], position[i][1]);
                out<<str;
                if(textAngle != 0)  {
                        out<<"gsave"<<endl;
                        sprintf(str, "%d rotate\n", textAngle);
                        out<<str;
                }
                if(ifmidalign)  {
                        sprintf(str, "(%s) %.2f FR C\n", text[i], fontSize);
                        out<<str;
                }
                sprintf(str, "(%s) FR %.2f PR\n", text[i], fontSize);
                out<<str;
                if(textAngle != 0)  {
                        out<<"grestore"<<endl;
                }
        }
        out<<"grestore"<<endl;
}

//arrow types:
//uparrow, downarrow, upleftarrow, uprightarrow, downleftarrow, downrightarrow
void 	PsShow::prtArrow(double **point, int pointNum, int ifCodTransform, char *arrow)
{
        int     i;
        char    str[200];
        codTransform(ifCodTransform, point, pointNum);
        out<<"%draw arrows"<<endl;
	out<<"gsave\n0 0 0 R\n"; //black arrow
        str[0] = 0;
        out<<str;
        for( i = 0; i < pointNum; i ++)    {
                sprintf(str, "%.2f %.2f %s\n",
                        point[i][0], point[i][1], arrow);
		out<<"gsave\n"<<str<<"grestore\n";
	}
	out<<"grestore\n";
}

//allowing shape: circle, box, cross
void 	PsShow::prtPoint(double **point, int pointNum, int ifCodTransform, char *shape, double *color)
{
        int     i;
        char    str[200];
        codTransform(ifCodTransform, point, pointNum);
	sprintf(str, "gsave\n%.2f %.2f %.2f R\n", color[0], color[1], color[2]);
        out<<str;
        for( i = 0; i < pointNum; i ++)    {
                sprintf(str, "%.4f %.4f %s\n",
                        point[i][0], point[i][1], shape);
		out<<"gsave\n"<<str<<"grestore\n";
	}
	out<<"grestore\n";
}

void    PsShow::changeTitlePosition(double positionX, double positionY)
{
        pointxTitle = positionX;
        pointyTitle = positionY;
}

void    PsShow::prtTitle(char *title, double positionX, double positionY, int linelen, int ifmidalign)
{
        int     i, nowi, nowbrack;
        int     slen = strlen(title);
        nowi = nowbrack = 0;
        double  positionx = positionX;
        double  positiony = positionY;
        char    str[300];
        i = 0;
        while(i < slen) {
                for(i = nowi; i < slen && (i - nowi) < linelen; i ++)      {
                        if(title[i] == ' ')     {
                                nowbrack = i;
                        }
                }
                if(i == slen)   nowbrack = slen;
                str[0] = 0;
                strncpy(str, title + nowi, nowbrack - nowi);
                str[nowbrack - nowi] = '\0';
		out<<"gsave\n0 0 0 R"<<endl;	
                prtText(positionx, positiony, str, titleFontSize, ifmidalign);
		out<<"grestore"<<endl;
                positiony -= 15;
                nowi = nowbrack;
        }
}

void    PsShow::prtTitle(char *title)
{
        int     linelen = 80;
        int     ifmidalign = 1;
        prtTitle(title, pointxTitle, pointyTitle, linelen, ifmidalign);
}

void    PsShow::prtTitle(char *title, double positionX, double positionY)
{
        pointxTitle = positionX;
        pointyTitle = positionY;
        prtTitle(title);
}

void   PsShow::drawColorString(int lineNum, int lineLen, char ***stringLine, char **lineDescript, int **fontorient, double ***color, int **colormed, double *linewidth, int ifmidalign, double spacewide, int ifwidfit)
{
        int     i = 0;
       	int	n = 0 ;
        int     lenCut = 50;
        int     add;
        int     end;
        double  xPos;
        double  yPos;
        char    *str = new char[1000];
        char    *seq = new char[10];
        int     beg;
        int     pg, pgtotal;
        int     newpage;
        double  yold;
        double  xWide;
        double  yWide;
        int     font1, font2;
        if(lineNum > 140)       {
                xWide = 4;
                yWide = 3;
                font1 = 3;
                font2 = 3;
        }
        else if(lineNum > 80)        {
                xWide = 6;
                yWide = 5;
                font1 = 5;
                font2 = 5;
        }
        else if(lineNum > 40)        {
                xWide = 7;
                yWide = 6;
                font1 = 7;
                font2 = 7;

        }
        else    {
                xWide = 8;
                yWide = 7;
                font1 = 7;
                font2 = 7;
        }
        xWide *= spacewide;
        if(ifwidfit)    xWide = (axisXMax - axisXMin - 60) / lineLen;

        yPos = axisYMax - 50;
        newpage = 0;
        pg = 0;
        add = 0;
        while(add < lineLen)    {
                if(newpage)  {
                        yPos = axisYMax - 50;
                }
                end = (add + lenCut)<lineLen?(add + lenCut):lineLen;
                yold = yPos;
                for(n = 0; n < lineNum; n ++)   {
                        xPos = axisXMin;
                        xPos += 100;
                        for(i = add; i < end; i ++)     {
                                xPos += xWide;
                                if(xPos > axisXMax)     break;
                        }
                        yPos -= yWide * linewidth[n];
                }
                yPos -= 2.0 * yWide;
                add = i;
                if((yPos - axisYMin) < (yold - yPos))   {
                        newpage = 1;
                        pg ++;
                }
                else    {
                        newpage = 0;
                }
        }
        pgtotal = pg + 1;

        yPos = axisYMax - 50;
        beg = 1;
        newpage = 0;
        pg = 0;
        add = 0;
        end = lineLen;
        while(add < lineLen)    {
                if(beg) {
                        if(pgtotal > 1) {
                                sprintf(str, "%.1f %.1f M (Page %d) FC %d PR\n", axisXMin, axisYMax - 30, pg + 1, font1);
                                out<<str;
                        }
                        beg = 0;
                }
                else  if(newpage)  {
                        sprintf(str, "showpage\ngrestore\n\n%%%%Page %d\ngsave\n", pg);
                        out<<str;
                        sprintf(str, "%.1f %.1f M (Page %d) FC %d PR\n", axisXMin, axisYMax - 30, pg + 1, font1);
                        out<<str;
                        yPos = axisYMax - 50;
                }
                yold = yPos;
                for(n = 0; n < lineNum; n ++)   {
                        xPos = axisXMin;
                        sprintf(str, "%.1f %.1f M (%s) FC %d PR\n", xPos, yPos, lineDescript[n], font1);
                        out<<str;
                        xPos += 60;
                        for(i = add; i < end; i ++)     {
                                xPos += xWide;
                                if(xPos > axisXMax)     break;
                                seq[0] = 0;
                                sscanf(stringLine[n][i], "%s", seq);
                                if(colormed[n][i])      {
                                        sprintf(str, "%.2f %.2f %.2f R\n", color[n][i][0],
                                        color[n][i][1], color[n][i][2]);
                                        out<<str;
                                        sprintf(str, "%.2f %.2f M\n", xPos, yPos);
                                        out<<str;
                                        if(fontorient[n][i] != 0)       {
                                                sprintf(str, "gsave\n%d rotate\n", fontorient[n][i]);
                                                out<<str;
                                        }
                                        if(ifmidalign)  {
                                                sprintf(str, "(%s) %d FC C\n", seq, font2);
                                                out<<str;
                                        }
                                        sprintf(str, "(%s) FC %d PR\n", seq, font2);
                                        out<<str;
                                        if( fontorient[n][i] != 0)  {
                                                out<<"grestore"<<endl;
                                        }
                                        out<<"0 G"<<endl;
                                }   //color string itself
                                else    {
                                        sprintf(str, "(%s) %.1f %.1f %.1f %.1f %.2f %.2f %.2f color_inv\n",
                                                seq, xPos, yPos, xPos, yPos, color[n][i][0],
                                                color[n][i][1], color[n][i][2]);
                                        out<<str;
                                }   //color string box
                        }
                        yPos -= yWide * linewidth[n];

                }
                yPos -= 2.0 * yWide;
                add = i;
                if((yPos - axisYMin) < (yold - yPos))   {
                        newpage = 1;
                        pg ++;
                }
                else    {
                        newpage = 0;
                }
        }

        delete[] str;
        delete[] seq;

        pointyTitle = yPos - 40;

}

void   PsShow::drawGrayStringndZhu(
                int lineNum, int lineLen, char ***stringLine, char **lineDescript, int **fontorient, int **fontstyle, double *linewidth, int ifmidalign, double spacewide, int ifwidfit,
                int *zhuptstring, int zhuNum, double **yvalue, double *yvaluemin, double *yvaluemax, char **yDescription, int **iffill, double *gray)
{
        int     i, j, n;
        int     add, end, pg, newpage;
        double  xPos, yPos;
        char    *str = new char[1000];
        char    *str1 = new char[200];
        char    *seq = new char[10];
        int     font1, font2;
        double  yWide, xWide, ZhuXWide, ZhuYWide, ytmp, yold;

        if(lineNum > 140)       {
                xWide = 4;
                yWide = 3;
                font1 = 3;
                font2 = 3;
        }
        else if(lineNum > 80)        {
                xWide = 6;
                yWide = 5;
                font1 = 5;
                font2 = 5;
        }
        else if(lineNum > 40)        {
                xWide = 7;
                yWide = 5;
                font1 = 7;
                font2 = 7;
        }
        else    {
                xWide = 8;
                yWide = 7;
                font1 = 7;
                font2 = 7;
        }
        xWide *= spacewide;
        if(ifwidfit)    xWide = (axisXMax - axisXMin - 65) / lineLen;
        ZhuXWide = xWide / 5.0;
        ZhuYWide = 10 * ZhuXWide;

        yPos = axisYMax - 50;
        yold = yPos;
        newpage = 0;
        pg = 1;
        add = 0;
        while(add < lineLen)    {
                if(newpage)    {
                        sprintf(str, "%%%%Page %d\ngsave\n", pg);
                        out<<str;
                        yPos = axisYMax - 50;
                        yold = yPos;
                }
                end = lineLen;
                for(n = 0; n < lineNum; n ++)   {
                        xPos = axisXMin;
                        sprintf(str, "%.1f %.1f M (%s) FC %d PR\n", xPos, yPos, lineDescript[n], font1);
                        out<<str;
                        xPos += 60;
                        for(i = add; i < end; i ++)     {
                                xPos += xWide;
                                if(xPos > axisXMax)     break;
                                seq[0] = 0;
                                sscanf(stringLine[n][i], "%s", seq);
                                sprintf(str, "%.2f %.2f M\n", xPos, yPos);
                                out<<str;
                                if(fontorient[n][i] != 0)       {
                                        sprintf(str, "gsave\n%d rotate\n", fontorient[n][i]);
                                        out<<str;
                                }
                                if(ifmidalign)  {
                                        sprintf(str, "(%s) %d %s C\n", seq, font2, fontDescription[fontstyle[n][i]]);
                                        out<<str;
                                }
                                sprintf(str, "(%s) %s %d PR\n", seq, fontDescription[fontstyle[n][i]], font2);
                                out<<str;
                                if( fontorient[n][i] != 0)  {
                                        out<<"grestore"<<endl;
                                }
                        }
                        end = i;
                        yPos -= yWide * linewidth[n];
                        for(i = 0; i < zhuNum; i ++)    {
                                if(zhuptstring[i] != n)   continue;
                                out<<"1.0 LW 0 G"<<endl;
                                xPos = axisXMin + 60;
                                yPos -= ZhuYWide;
                                sprintf(str, "%.2f %.2f M %.2f %.2f L\n", xPos, yPos, axisXMax, yPos);
                                strcpy(str1, str);
                                for(j = add; j < end; j ++)     {
                                        xPos += xWide;
                                        ytmp = ZhuYWide * ((yvalue[i][j] - yvaluemin[i]) /(yvaluemax[i] - yvaluemin[i]));
                                        if(iffill[i][j])        {
                                                sprintf(str, "%.2f G\n", gray[iffill[i][j] - 1]);
                                                out<<str;
                                        }
                                        sprintf(str, "newpath %.2f %.2f M %.2f %.2f L %.2f %.2f L %.2f %.2f L ",
                                                xPos - ZhuXWide, yPos, xPos - ZhuXWide, yPos + ytmp, xPos + ZhuXWide, yPos + ytmp, xPos + ZhuXWide, yPos);
                                        out<<str;
                                        if(iffill[i][j])        out<<"closepath fill"<<endl;
                                        else    out<<"stroke"<<endl;
                                        if(iffill[i][j])        out<<"0 G"<<endl;
                                }
                                yPos -= yWide;
                                out<<"newpath\n"<<str1<<"stroke"<<endl;
                        }
                }
                yPos -= 2.0 * yWide;
                if((yPos - axisYMin) < (yold - yPos))   {
                        newpage = 1;
                        pg ++;
                }
                else    {
                        newpage = 0;
                }
                add = end;
                if(add < lineLen && newpage)       {
                        sprintf(str, "showpage\ngrestore\n");
                        out<<str;
                }
        }

        delete[] str;
        delete[] str1;
        delete[] seq;
        pointyTitle = yPos - 40;
}
