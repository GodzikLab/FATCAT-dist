#ifndef PSSHOW_H
#define PSSHOW_H
#define MAXCOLOR  100
#define MAXITEM 200

class PsShow
{
        private:
                int     totalColor;
                double  **RGBValue;
                int     nowColorIndex;
                int     itemNum;
                char    **itemDescript;
                double  **itemColor;
                int     fonttypenum;
                char    **fontDescription;
                double  lineWidth;
                double  lineGray;
                double  scaleX;
                double  scaleY;
                double  translateX;
                double  translateY;
                double  axisXMin;
                double  axisXMax;
                double  axisYMin;
                double  axisYMax;
                double  axisYMinReal;
                double  axisYMaxReal;
                double  textFontSize;
                double  titleFontSize;
                double  markFontSize;
                int     textAngle;
                int     ifpointSymble;
		double	symbleRad;
                char    pointSymble[5];
		int	ifPrtLine;
                double  pointxTitle;
                double  pointyTitle;
		int	xAxisDescript;
		int	yAxisDescript;
		int	ifDrawAxisBox;
		int	upsidedown;
                ofstream out;
        public:
                PsShow(char *filename);
                PsShow(char *filename, char *title);
                PsShow(char *filename, char *head, char *title);
                ~PsShow();
		void	Initial(char *filename, char *head, char *title);
		void	xAxisDescriptOff(void);
		void	yAxisDescriptOff(void);
                void    drawAxisBoxOff(void);
		void	drawAxis(double minx, double maxx, double miny, double maxy);
                void    drawAxis(double **point, int pointNum);
                void    drawAxis(double **point, int pointNum, char *x_lab, char *y_lab);
                void    drawAxisBox(void);
                void    defaultPage(void); //middle size page
                void    largePage(void);
                void    smallPage(void);
                void    PageArea(double xwide, double ywide);
                void    PageArea(double xmin, double xwide, double ymin, double ywide);
                void    setScale(double scalex, double scaley);
                void    setTranslate(double translatex, double translatey);
		void	setTransform(void);
		void	setupsidedown(void);
                void    gsave(void);
                void    grestore(void);
                void    creatColor(char *lab);
                void    defaultColor(void);
                void    assignColor(double *rgbvalue);
		void	setMarkFontSize(double size);
		void	setTextFontSize(double size);
                void    prtContLine(double **line, int num_line, int ifCodTranform, double linewidth);
                void    prtContLine(double **line, int num_line, int ifCodTranform);
                void    prtContLine(double **line, int num_line, int ifCodTranform, char *symbl);
                void    prtContLine(double **line, int num_line, int ifCodTranform, char *symbl, int ifline);
                void    prtStepLine(double **line, int num_line, int ifCodTranform);
                void    prtStepLine(double **line, int num_line);
		void	prtArrow(double **point, int pointNum, int ifCodTransform, char *arrow);
		void	prtPoint(double **point, int pointNum, int ifCodTransform, char *shape, double *color);
                void    prtGrayStepLine(double **line, int num_line, double *gray, double *linewidth, int ifCodTranform);
                void    prtColorStepLine(double **line, int num_line, double **color, int ifCodTranform);
		void	prtColorStepLine(double **line, int num_line, double **color, double *linewidth, int ifCodTranform);
                void    digitDiv(double v1, double *v2, double *v3);
		void	remove0(char *str);
                void    codTransform(int ifCodTransform, double **point, double pointNum);
                void    codTransform(double **point, double pointNum);
                void    prtDescription(void);
                void    prtDescription(double toMaxX, double toMaxY, char *description);
                void    prtText(double positionX, double positionY, char *text, double fontSize, int ifmidalign);
                void    prtText(double **position, char **text, int textNum, double fontSize, int ifmidalign);
                void    changeTitlePosition(double positionX, double positionY);
                void    prtTitle(char *title);
                void    prtTitle(char *title, double positionX, double positionY);
                void    prtTitle(char *title, double positionX, double positionY, int linelen, int ifmidalign);
                void    drawColorString(int lineNum, int lineLen, char ***stringLine, char **lineDescript, 
				        int **fontorient, double ***color, int **colormed, double *linwidth, 
					int ifmidalign, double spacewide, int ifwidfit);
                void    prtZhuZhuangTu(int itemnum, int linelen, double **yvalue, double *yvaluemin, 
				       double *yvaluemax, char ***xdescription, char **ydescription, 
				       int **iffill, double *gray);
                void    drawGrayStringndZhu(int lineNum, int lineLen, char ***stringLine, char **lineDescript, 
					    int **fontorient, int **fontstyle, double *linewidth, int ifmidalign, 
					    double spacewide, int ifwidfit, int *zhuptstring, int zhuNum, 
					    double **yvalue, double *yvaluemin, double *yvaluemax, 
					    char **yDescription, int **iffill, double *gray);
};

#endif

