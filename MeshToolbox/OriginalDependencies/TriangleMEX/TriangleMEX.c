/*MATLAB Interface to Triangle mesh function.*/
/*Paolo Bardella, 26/2/2016 */

#include "mex.h"

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#include <stdio.h>
#include <stdlib.h>
#include "triangle.h"

#define NUMBER_OF_FIELDS (sizeof(field_names)/sizeof(*field_names))

/*****************************************************************************/
/*                                                                           */
/*  report()   Print the input or output.                                    */
/*                                                                           */
/*****************************************************************************/

void report(struct triangulateio *io,int markers,int reporttriangles,int reportsegments,int reportedges,int reportnorms)
{
	int i, j;
	mexPrintf("Point       x           y           marker\n");
	mexPrintf("==========================================\n");
	for (i = 0; i < io->numberofpoints; i++) {
		mexPrintf("%4d    ", i);
		for (j = 0; j < 2; j++) {
			mexPrintf("%12.6f", io->pointlist[i * 2 + j]);
		}		
		if (markers) {
			mexPrintf("%8d\n", io->pointmarkerlist[i]);
		}
		else {
			mexPrintf("\n");
		}
	}
	mexPrintf("\n");

	
	mexPrintf("Segment     Pt1         Pt2         marker\n");
	mexPrintf("==========================================\n");

	if (reportsegments) {
		for (i = 0; i < io->numberofsegments; i++) {
			mexPrintf("%4d  ", i);
			mexPrintf("%8d    ", io->segmentlist[i * 2 ]);
			mexPrintf("%8d", io->segmentlist[i * 2 + 1]);
			if (markers) {
				mexPrintf("      %8d", io->segmentmarkerlist[i]);
			}
			mexPrintf("\n");
			
		}
		mexPrintf("\n");
	}
}

double mxFieldToDouble(const mxArray* mArray, const char*fieldname)
{
	mxArray* aux = mxGetField(mArray, 0, fieldname);
	if (aux == NULL)
		mexErrMsgIdAndTxt("MATLAB:triangle:missingInputField", "Field %s is missing.", fieldname);
	if (aux == NULL)
		mexErrMsgIdAndTxt("MATLAB:triangle:missingInputField", "Field %s is missing.", fieldname);
	if (!mxIsDouble(aux))
		mexErrMsgIdAndTxt("MATLAB:triangle:invalidInputFieldDouble", "%s must be double.", fieldname);
	if (mxIsComplex(aux))
		mexErrMsgIdAndTxt("MATLAB:triangle:invalidInputFieldDouble", "%s must be real.", fieldname);
	if (mxGetM(aux)*mxGetN(aux) != 1) {
		mexErrMsgIdAndTxt("MATLAB:triangle:invalidInputField",
			"%s must be a noncomplex scalar double.",fieldname);
	}
	return mxGetScalar(aux);
}


REAL* mxFieldToREALMatrix(const mxArray* mArray, const char*fieldname, const int numRows, const int numCols)
{
    double *source;
	REAL* dest;
	REAL* destIndex;
    int r, c;
	mxArray* aux = mxGetField(mArray, 0, fieldname);
	if(aux==NULL)
		mexErrMsgIdAndTxt("MATLAB:triangle:missingInputField", "Field %s is missing.", fieldname);
	if (!mxIsDouble(aux))
		mexErrMsgIdAndTxt("MATLAB:triangle:invalidInputFieldDouble","%s must contain double data.", fieldname);
	if(mxIsComplex(aux))
		mexErrMsgIdAndTxt("MATLAB:triangle:invalidInputFieldDouble", "%s must contain real data.", fieldname);
	if (mxGetM(aux)!=numRows)		
		mexErrMsgIdAndTxt("MATLAB:triangle:invalidInputField", "%s must be a matrix with %d row(s).", fieldname, numRows);
	if (mxGetN(aux) != numCols)
		mexErrMsgIdAndTxt("MATLAB:triangle:invalidInputField", "%s must be a matrix with %d column(s).", fieldname, numCols);


    source = (double*)mxGetData(aux);
	dest = (REAL *)malloc(numRows * numCols * sizeof(REAL));
	destIndex = dest;
	
	for ( r = 0; r < numRows; r++)
	{
		for ( c = 0; c < numCols; c++)
		{
			*destIndex++ = (REAL)source[c * numRows + r];
		}
	}
	return dest;
}


int* mxFieldToIntMatrix(const mxArray* mArray, const char*fieldname, const int numRows, const int numCols)
{
    int* dest;
	int*destIndex;
	double *source;
    int r,c;
	mxArray* aux = mxGetField(mArray, 0, fieldname);
	if (aux == NULL)
		mexErrMsgIdAndTxt("MATLAB:triangle:missingInputField", "Field %s is missing.", fieldname);

	if (!mxIsDouble(aux))
		mexErrMsgIdAndTxt("MATLAB:triangle:invalidInputFieldDouble", "%s must contain double data.", fieldname);
	if (mxIsComplex(aux))
		mexErrMsgIdAndTxt("MATLAB:triangle:invalidInputFieldDouble", "%s must contain real data.", fieldname);
	if (mxGetM(aux) != numRows)
		mexErrMsgIdAndTxt("MATLAB:triangle:invalidInputField", "%s must be a matrix with %d row(s).", fieldname, numRows);
	if (mxGetN(aux) != numCols)
		mexErrMsgIdAndTxt("MATLAB:triangle:invalidInputField", "%s must be a matrix with %d column(s).", fieldname, numCols);
	 dest=(int *)malloc(numRows * numCols * sizeof(int));
	destIndex = dest;
	source = (double*)mxGetData(aux);
	for ( r = 0; r < numRows; r++)
	{
		for ( c = 0; c < numCols; c++)
		{
			*destIndex++ = (int)source[c * numRows + r];
		}
	}
	return dest;
}

mxArray* toMxArrayInt(int* data, int rows, int cols)
{
    int t,c;
	mxArray* tmp = mxCreateDoubleMatrix(rows, cols, mxREAL);

	double* tmpData = mxGetPr(tmp);
	for ( t = 0; t < rows; t++)
		for ( c = 0; c<cols; c++)
			tmpData[c * rows + t] = data[t * cols + c];
	return tmp;
}

mxArray* toMxArray(double* data, int rows, int cols)
{
    int t,c;
	mxArray* tmp = mxCreateDoubleMatrix(rows, cols, mxREAL);
	double* tmpData = mxGetPr(tmp);
	for (t = 0; t < rows; t++)
		for (c = 0; c<cols; c++)
			tmpData[c * rows + t] = data[t * cols + c];
	return tmp;
}


/*****************************************************************************/
/*                                                                           */
/*  main()   Create and refine a mesh.                                       */
/*                                                                           */
/*****************************************************************************/
struct triangulateio* FreeTriangleStruct(struct triangulateio* str)
{
    if(str!=NULL)
    {
        if(str->pointlist!=NULL) 
            free(str->pointlist);
        str->pointlist = (REAL *)NULL;				

        if(str->pointattributelist!=NULL) 
            free(str->pointattributelist);
        str->pointattributelist = (REAL *)NULL;

        if(str->pointmarkerlist!=NULL) 
            free(str->pointmarkerlist);
        str->pointmarkerlist = (int *)NULL;

        if(str->trianglelist!=NULL) 
            free(str->trianglelist);
        str->trianglelist = (int *)NULL;

        if(str->triangleattributelist!=NULL) 
            free(str->triangleattributelist);
        str->triangleattributelist = (REAL *)NULL;

        if(str->segmentlist!=NULL) 
            free(str->segmentlist);
        str->segmentlist = (int *)NULL;

        if(str->segmentmarkerlist!=NULL) 
            free(str->segmentmarkerlist);
        str->segmentmarkerlist = (int *)NULL;

        if(str->edgelist!=NULL) 
            free(str->edgelist);
        str->edgelist = (int *)NULL;

        if(str->edgemarkerlist!=NULL) 
            free(str->edgemarkerlist);
        str->edgemarkerlist = (int *)NULL;

        if(str->trianglearealist!=NULL) 
            free(str->trianglearealist);
        str->trianglearealist = (REAL *)NULL;
    }

    return (struct triangulateio*)NULL;
}

struct triangulateio* AllocateTriangleStruct()
{
    struct triangulateio* str=(struct triangulateio*)malloc(sizeof(struct triangulateio));      
 /* Make necessary initializations so that Triangle can return a triangulation in `out'.*/
    str->numberoftriangles=0;
    str->numberofpoints=0;
    str->numberofsegments=0;
    str->numberofholes=0;
    str->numberofregions=0;
    str->numberofpointattributes=0;
    str->pointlist = (REAL *)NULL;				/* numberofpoints*2*/
    str->pointattributelist = (REAL *)NULL;	/* numberofpoints * numberofpointattributes*/
    str->pointmarkerlist = (int *)NULL;		/* numberofpoints*/
    str->trianglelist = (int *)NULL;			/* numberoftriangles*numberofcorners*/
    str->triangleattributelist = (REAL *)NULL; /*numberoftriangles*numberoftriangleattributes*/
    str->segmentlist = (int *)NULL;			/*numberofsegments*2*/
    str->segmentmarkerlist = (int *)NULL;		/*numberofsegments*/
    str->edgelist = (int *)NULL;				/*numberofedges*2*/
    str->edgemarkerlist = (int *)NULL;			/*numberofedges*/  
    str->trianglearealist = (double*)NULL;		/*numberoftriangles?*/
    return str;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,	const mxArray *prhs[])
{
    int k;
    int PreviousTriangles=0;
    int Refine=1;
    int VerboseInput=0;
	struct triangulateio *in, *out;
	char TriangulationSwitches[128];
	 const char *field_names[] = {"pointlist", "pointattributelist","pointmarkerlist","trianglelist",
  "triangleattributelist","segmentlist","segmentmarkerlist","edgelist","edgemarkerlist"};

    mwSize dims[2] = { 1, 1 };
	/*CHECK CONVERSION*/
	if (nrhs != 1 && nrhs != 2)
		mexErrMsgIdAndTxt("MATLAB:triangle:invalidNumInputs","One or two input required.");
	else if (!mxIsStruct(prhs[0]))
		mexErrMsgIdAndTxt("MATLAB:triangle:invalidInputType1","Input must be a structure.");
	else if (mxGetNumberOfElements(prhs[0])>1)
		mexErrMsgIdAndTxt("MATLAB:triangle:invalidInputType2", "First input must be a single structure.");
	/* get input arguments */
	else if(mxGetNumberOfFields(prhs[0])<10)
		mexErrMsgIdAndTxt("MATLAB:triangle:invalidInputType3", "The structure type is invalid.");
	Refine=(nrhs==2)*5;
	/*START INPUT CONVERSION*/

    in=AllocateTriangleStruct();
    in->numberofpoints = (int)mxFieldToDouble(prhs[0], "numberofpoints");
	
	in->pointlist = mxFieldToREALMatrix(prhs[0], "pointlist", in->numberofpoints, 2);
	in->pointmarkerlist = mxFieldToIntMatrix(prhs[0], "pointmarkerlist", in->numberofpoints,1);

	in->numberofpointattributes = (int)mxFieldToDouble(prhs[0], "numberofpointattributes");
	if (in->numberofpointattributes > 0)
		in->pointattributelist = mxFieldToREALMatrix(prhs[0], "pointattributelist", in->numberofpoints, in->numberofpointattributes);
	else
		in->pointattributelist = NULL;	
	in->numberofsegments = (int)mxFieldToDouble(prhs[0], "numberofsegments");
	if (in->numberofsegments > 0)
	{
		in->segmentlist = mxFieldToIntMatrix(prhs[0], "segmentslist", in->numberofsegments, 2);
		in->segmentmarkerlist = mxFieldToIntMatrix(prhs[0], "segmentmarkerlist", in->numberofsegments,1);
	}
	else
	{
		in->segmentlist = NULL;
		in->segmentmarkerlist = NULL;
	}

	in->numberofholes = (int)mxFieldToDouble(prhs[0], "numberofholes");
	if (in->numberofholes > 0)
		in->holelist = mxFieldToREALMatrix(prhs[0], "holelist", in->numberofholes,2);
	else
		in->holelist = NULL;

	in->numberofregions = (int)mxFieldToDouble(prhs[0], "numberofregions");
	if (in->numberofregions > 0)
		in->regionlist = mxFieldToREALMatrix(prhs[0], "regionlist", in->numberofregions,4);
	else
		in->regionlist = NULL;
	
	strcpy(TriangulationSwitches, "Aqpea");
    VerboseInput=(int)mxFieldToDouble(prhs[0], "verboseinput");
	/* END INPUT CONVERSION*/

    mexPrintf("TRIANGLE v.1.6 - MEX v.1.3\n");
	/*mexPrintf("TRIANGLE switches: %s\n", TriangulationSwitches);*/
    if (VerboseInput) 
        report(in, 1, 0, 1, 0, 0);
    /* Make necessary initializations so that Triangle can return a triangulation in `out'*/
    out=AllocateTriangleStruct();
    mexPrintf("Calculating triangulation..."); 
    triangulate(TriangulationSwitches, in, out, NULL);   
    in=FreeTriangleStruct(in);   
    mexPrintf("Done.\nTRIANGLE generated a mesh with %d points and %d triangles\n", 
          out->numberofpoints,out->numberoftriangles);
    PreviousTriangles=0;
    if (Refine>0)
        mexPrintf("Refining triangulation..."); 
    while(Refine-->0 && fabs(out->numberoftriangles-PreviousTriangles)>out->numberoftriangles*0.1)
    {  
        double* MaxAreaRef;
        mxArray* prhsTemp[3];
        mxArray* plhsTemp[1];
        double* tmpDataX;
        double* tmpDataY;
        
        mexPrintf("."); 
        PreviousTriangles=out->numberoftriangles;
        in=out;                
        prhsTemp[0]=prhs[1];
        prhsTemp[1]= mxCreateDoubleMatrix(in->numberoftriangles, 1, mxREAL);
        prhsTemp[2]= mxCreateDoubleMatrix(in->numberoftriangles, 1, mxREAL);
        tmpDataX = mxGetPr(prhsTemp[1]);
        tmpDataY = mxGetPr(prhsTemp[2]);
        for (k=0; k<in->numberoftriangles;k++)
        {                
          tmpDataX[k]=(
                  in->pointlist[2*(in->trianglelist[3*k]-1)]+
                  in->pointlist[2*(in->trianglelist[3*k+1]-1)]+
                  in->pointlist[2*(in->trianglelist[3*k+2]-1)])/3.0;
          tmpDataY[k]=(
                  in->pointlist[2*(in->trianglelist[3*k]-1)+1]+
                  in->pointlist[2*(in->trianglelist[3*k+1]-1)+1]+
                  in->pointlist[2*(in->trianglelist[3*k+2]-1)+1])/3.0;

        /*mexPrintf("T%d: %d(%f,%f) %d(%f,%f) %d(%f,%f)\n",k    ,           
                in->trianglelist[3*k],
                in->pointlist[2*(in->trianglelist[3*k]-1)],
                in->pointlist[2*(in->trianglelist[3*k]-1)+1],
                in->trianglelist[3*k+1],
                in->pointlist[2*(in->trianglelist[3*k+1]-1)],
                in->pointlist[2*(in->trianglelist[3*k+1]-1)+1],
                in->trianglelist[3*k+2],
                in->pointlist[2*(in->trianglelist[3*k+2]-1)],
                in->pointlist[2*(in->trianglelist[3*k+2]-1)+1]);*/
       }
       mexCallMATLAB(1, plhsTemp, 3,prhsTemp, "feval");
       if(in->trianglearealist!=NULL) free(in->trianglearealist);
       in->trianglearealist = (REAL *)malloc(in->numberoftriangles * sizeof(REAL));
       
       MaxAreaRef = (double*)mxGetData(plhsTemp[0]);       
       if(mxGetM(plhsTemp[0])*mxGetN(plhsTemp[0])==1)/*scalar*/
       {
           for (k=0;k<in->numberoftriangles;k++)
            {
            in->trianglearealist[k] = MaxAreaRef[0];
            }  
       }
       else if(mxGetM(plhsTemp[0])*mxGetN(plhsTemp[0])==in->numberoftriangles)/*vector*/
       {   
           for (k=0;k<in->numberoftriangles;k++)
           {
               in->trianglearealist[k] = MaxAreaRef[k];
           }
       }
       else
           mexErrMsgIdAndTxt("MATLAB:triangle:invalidAreaSize",
                   "The function describing the area must return a scalar, or a vector with %d elements.",in->numberoftriangles);
   
       mxDestroyArray(prhsTemp[1]);
       mxDestroyArray(prhsTemp[2]);   
       mxDestroyArray(plhsTemp[0]);

      /* Refine the triangulation according to the attached triangle area constraints.  */
        out=AllocateTriangleStruct();
        triangulate("rqapAe", in, out, (struct triangulateio *) NULL);
        in=FreeTriangleStruct(in);   


    }
    if (PreviousTriangles!=0)
        mexPrintf("Done.\nTRIANGLE generated a mesh with %d points and %d triangles\n", 
          out->numberofpoints,out->numberoftriangles);
    in=FreeTriangleStruct(in);
  
  
    /****************EXPORT RESULTS*******************/
 
    plhs[0] = mxCreateStructArray(2, dims, NUMBER_OF_FIELDS, field_names);

    /*mexPrintf("Exporting pointlist\n");*/
    mxSetField(plhs[0], 0, "pointlist", 
              toMxArray(out->pointlist, out->numberofpoints, 2));
  
    /*mexPrintf("Exporting pointattributelist\n");*/
    mxSetField(plhs[0], 0, "pointattributelist", 
          toMxArray(out->pointattributelist, out->numberofpoints, out->numberofpointattributes));

    /*mexPrintf("Exporting pointmarkerlist\n");*/
    mxSetField(plhs[0], 0, "pointmarkerlist", 
          toMxArrayInt(out->pointmarkerlist, out->numberofpoints,1));

    /*mexPrintf("Exporting trianglelist\n");*/
    mxSetField(plhs[0], 0, "trianglelist", 
          toMxArrayInt(out->trianglelist, out->numberoftriangles, out->numberofcorners));

    /*mexPrintf("Exporting triangleattributelist\n");*/
    mxSetField(plhs[0], 0, "triangleattributelist", 
          toMxArray(out->triangleattributelist, out->numberoftriangles, out->numberoftriangleattributes));
  
    /*mexPrintf("Exporting segmentlist\n");*/
    mxSetField(plhs[0], 0, "segmentlist",
          toMxArrayInt(out->segmentlist, out->numberofsegments,2));

    /*mexPrintf("Exporting segmentmarkerlist\n");*/
    mxSetField(plhs[0], 0, "segmentmarkerlist", 
          toMxArrayInt(out->segmentmarkerlist, out->numberofsegments,1));

    /*mexPrintf("Exporting edgemarkerlist\n");*/
    mxSetField(plhs[0], 0, "edgemarkerlist", 
          toMxArrayInt(out->edgemarkerlist, out->numberofedges,1));;

    /*mexPrintf("Exporting edgelist\n");*/
    mxSetField(plhs[0], 0, "edgelist", 
          toMxArrayInt(out->edgelist, out->numberofedges,2));
  

    /****************END EXPORT RESULTS*******************/

    /* Free all allocated arrays, including those allocated by Triangle. */
    out=FreeTriangleStruct(out); 
}
