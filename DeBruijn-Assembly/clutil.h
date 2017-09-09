/*Time-stamp: <2013-01-15 11:28:12 by Yu Yang>
* =======================================================================
*       Filename:  clutil.h
*        Created:  2014-09-03 11:00:35
*       Compiler:  gcc
*
*         Author:  Yu Yang
*			Email:  yy1990cn@gmail.com
* =======================================================================
*/

#ifndef _CLUTIL_H_
#define _CLUTIL_H_ 1

const char *clGetErrorString(cl_int error);
void clCheckEqWithMsg(int status, int expect, char *msg);
void clCheckNeqWithMsg(int status, int expect, char *msg);
void clCheckLtWithMsg(int status, int expect, char *msg);
void clCheckGtWithMsg(int status, int expect, char *msg);

cl_context CreateContext();

#endif /* _CLUTIL_H_ */