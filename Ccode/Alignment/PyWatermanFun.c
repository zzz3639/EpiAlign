#include"Python.h"
#include"WatermanFun.h"

static PyObject* WF_SWA_Even_Naive(PyObject* self, PyObject* args)
{
    PyObject *output;
    PyObject *py_seq1,*py_seq2,*py_alpha;
    float alpha;
    if(!PyArg_ParseTuple(args, "OOO", &py_seq1,&py_seq2,&py_alpha))
	return NULL;
    if(PyFloat_Check(py_alpha)) /*check whether alpha is python float*/
	alpha = (float)PyFloat_AS_DOUBLE(py_alpha);
    else /*regard alpha as python int*/
	alpha = (float)PyInt_AS_LONG(py_alpha);
    int n1,n2;
    unsigned char *seq1,*seq2;
    n1=PyList_GET_SIZE(py_seq1);
    n2=PyList_GET_SIZE(py_seq2);
    seq1=(unsigned char *)malloc(sizeof(unsigned char)*n1);
    seq2=(unsigned char *)malloc(sizeof(unsigned char)*n2);
    int i,j;
    for(i=0;i<n1;i++)
	seq1[i] = (unsigned char)PyInt_AS_LONG(PyList_GET_ITEM(py_seq1,i));
    for(i=0;i<n2;i++)
	seq2[i] = (unsigned char)PyInt_AS_LONG(PyList_GET_ITEM(py_seq2,i));
    float **map_score;
    unsigned char ** map_trace;
    map_trace=(unsigned char **)malloc(sizeof(unsigned char*)*(n2+1));
    map_score=(float **)malloc(sizeof(float*)*(n2+1));
    for(i=0;i<n2+1;i++){
	map_trace[i]=(unsigned char *)malloc(sizeof(unsigned char)*(n1+1));
	map_score[i]=(float *)malloc(sizeof(float)*(n1+1));
    }
    SWA_Even(seq1,seq2,n1,n2,MatchScore_Naive,alpha,NULL,NULL,map_score,map_trace);
    output = PyList_New(n1+1);
    PyObject *py_temp,*py_temp2;
    for(i=0;i<n1+1;i++){
	py_temp = PyList_New(n2+1);
	for(j=0;j<n2+1;j++){
	    py_temp2 = Py_BuildValue("[fi]", map_score[j][i], (int)map_trace[j][i]);
	    PyList_SET_ITEM(py_temp,j,py_temp2);
	}
        PyList_SET_ITEM(output, i, py_temp);
    }
    Free_Map(n2,map_score,map_trace);
    return output;
}

static PyMethodDef WF_Methods[] = {
    {"SWA_Even_Naive",WF_SWA_Even_Naive,METH_VARARGS,NULL},
    {NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initPycSWF(void){
    (void) Py_InitModule("PycSWF",WF_Methods);
}


