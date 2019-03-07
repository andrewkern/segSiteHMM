#ifndef ATLAS_CACHEEDGE_H
   #define ATLAS_CACHEEDGE_H
   #define CacheEdge 917504
   #if !defined(ATL_NCPU) || ATL_NCPU <= 1
      #undef CacheEdge
      #define CacheEdge 720896
   #endif
   #ifdef SREAL
      #if !defined(ATL_NCPU) || ATL_NCPU <= 1
         #undef CacheEdge
         #define CacheEdge 458752
      #endif
   #endif
   #ifdef SCPLX
      #if !defined(ATL_NCPU) || ATL_NCPU <= 1
         #undef CacheEdge
         #define CacheEdge 524288
      #endif
   #endif
   #ifdef DCPLX
      #undef CacheEdge
      #define CacheEdge 524288
   #endif
#endif
