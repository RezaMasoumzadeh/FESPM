
#ifndef __spmgriditerators_h__
#define __spmgriditerators_h__


#define NodeIdx(s,i,j) ((i) + (j)*((s)->nx))
#define CellIdx(s,i,j) ((i) + (j)*((s)->mx))

/*
#define CellGetNodeIdx(s,ci,cj,idx){\
 (idx)[0] = NodeIdx((s),(ci),(cj));\
 (idx)[1] = NodeIdx((s),(ci)+1,(cj));\
 (idx)[2] = NodeIdx((s),(ci),(cj)+1);\
 (idx)[3] = NodeIdx((s),(ci+1),(cj)+1);}
*/

#define CellGetNodeIdx(s,ci,cj,idx){\
(idx)[0] = NodeIdx((s),(ci),(cj));\
(idx)[1] = (idx)[0] + 1;\
(idx)[2] = (idx)[0] + (s)->nx;\
(idx)[3] = (idx)[2] + 1;}

#define IteratorEnd_TerminateDoubleLoop(s)\
    }\
  }\
}

#define IteratorEnd_TerminateSingleLoop(s)\
  }\
}

#define NodalIteratorBegin(s)\
{\
  int ni,nj;\
  for (nj=0; nj<(s)->ny; nj++) {\
    for (ni=0; ni<(s)->nx; ni++) {\

#define NodalIteratorEnd(s)\
    }\
  }\
}


#define NodalIteratorInteriorBegin(s)\
{\
  int ni,nj;\
  for (nj=1; nj<(s)->ny-1; nj++) {\
    for (ni=1; ni<(s)->nx-1; ni++) {\

#define NodalIteratorInteriorEnd(s)\
    }\
  }\
}


#define NodalIteratorJMINEdgeBegin(s)\
{\
  int ni,nj=0;\
  for (ni=0; ni<(s)->nx; ni++) {
#define NodalIteratorJMINEdgeEnd(s) IteratorEnd_TerminateSingleLoop(s)


#define NodalIteratorJMAXEdgeBegin(s)\
{\
  int ni,nj=(s)->ny-1;\
  for (ni=0; ni<(s)->nx; ni++) {
#define NodalIteratorJMAXEdgeEnd(s) IteratorEnd_TerminateSingleLoop(s)


#define NodalIteratorIMINEdgeBegin(s)\
{\
  int ni=0,nj;\
  for (nj=0; nj<(s)->ny; nj++) {
#define NodalIteratorIMINEdgeEnd(s) IteratorEnd_TerminateSingleLoop(s)


#define NodalIteratorIMAXEdgeBegin(s)\
{\
  int ni=(s)->nx-1,nj;\
  for (nj=0; nj<(s)->nx; nj++) {
#define NodalIteratorIMAXEdgeEnd(s) IteratorEnd_TerminateSingleLoop(s)


#define CentroidIteratorBegin(s)\
{\
  int mi,mj;\
  for (mj=0; mj<(s)->my; mj++) {\
    for (mi=0; mi<(s)->mx; mi++) {\

#define CentroidIteratorEnd(s)\
    }\
  }\
}




#endif
