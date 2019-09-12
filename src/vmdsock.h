/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2008 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: vmdsock.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.12 $      $Date: 2008/03/27 19:36:52 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   socket interface layer, abstracts platform-dependent routines/APIs
 ***************************************************************************/

#if defined(VMDSOCKINTERNAL)

#if !defined(_MSC_VER)
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <sys/file.h>
#endif

typedef struct {
  struct sockaddr_in addr; /* address of socket provided by bind() */
  int addrlen;             /* size of the addr struct */
  int sd;                  /* socket file descriptor */
} vmdsocket;

#endif /* VMDSOCKINTERNAL */

#ifdef _WIN32
    #if defined(SimTK_MOLMODEL_BUILDING_SHARED_LIBRARY)
        #define SimTK_MOLMODEL_EXPORT __declspec(dllexport)
    #elif defined(SimTK_MOLMODEL_BUILDING_STATIC_LIBRARY) || defined(SimTK_USE_STATIC_LIBRARIES)
        #define SimTK_MOLMODEL_EXPORT
    #else
        #define SimTK_MOLMODEL_EXPORT __declspec(dllimport)   // i.e., a client of a shared library
    #endif
#else
    #define SimTK_MOLMODEL_EXPORT // Linux, Mac
#endif

#ifdef __cplusplus
extern "C" {
#endif

int  SimTK_MOLMODEL_EXPORT  vmdsock_init(void);
void SimTK_MOLMODEL_EXPORT *vmdsock_create(void);
int  SimTK_MOLMODEL_EXPORT  vmdsock_bind(void *, int);
int  SimTK_MOLMODEL_EXPORT  vmdsock_listen(void *);
void SimTK_MOLMODEL_EXPORT *vmdsock_accept(void *);  /* return new socket */
int  SimTK_MOLMODEL_EXPORT  vmdsock_connect(void *, const char *, int);
int  SimTK_MOLMODEL_EXPORT  vmdsock_write(void *, const void *, int);
int  SimTK_MOLMODEL_EXPORT  vmdsock_read(void *, void *, int);
int  SimTK_MOLMODEL_EXPORT  vmdsock_selread(void *, int);
int  SimTK_MOLMODEL_EXPORT  vmdsock_selwrite(void *, int);
void SimTK_MOLMODEL_EXPORT  vmdsock_shutdown(void *);
void SimTK_MOLMODEL_EXPORT  vmdsock_destroy(void *);

#ifdef __cplusplus
}
#endif

