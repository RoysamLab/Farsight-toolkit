
SET(QHULL_HDRS
 geom.h
 io.h
 libqhull.h
 mem.h
 merge.h
 poly.h
 qhull.h
 qhull_a.h
 qset.h
 random.h
 stat.h
user.h	
)


SET(QHULL_SRCS
        rboxlib.c
        user.c
        global.c
        stat.c
        io.c
        geom2.c
        poly2.c
        merge.c
        libqhull.c
        geom.c
        poly.c
        qset.c
        mem.c
        usermem.c
        userprintf.c
        random.c
        )
        
ADD_LIBRARY( QHULL
            ${QHULL_HDRS}
            ${QHULL_SRCS}
)
   






