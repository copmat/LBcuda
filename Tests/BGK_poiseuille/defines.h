#define BGK
#define noSCP
#define noMOMBGK
#define noCG
#define noPARTICLES
#define FLUIDSPLIT
#define noHIGHGRAD
#define noDIFFDENS
#define NEWINPUT
#define noD3Q27
#define APPLYBC
#define noPBC
#define SCPD3Q19GRAD

#ifndef MYDIMESION
#define MYDIMESION 64
#endif

#ifndef TILE1
#define TILE1 8
#endif
#ifndef TILE2
#define TILE2 4
#endif
#ifndef TILE3
#define TILE3 4
#endif

#define noDEBUG_FORCE
#define noDEBUG_ROT
#define noDEBUG_MKRM
#define noDEBUG_N2P
#define noDEBUG_P2N
#define noCHECK_WRONGATOM
#define noCHECK_VOLP
#define noMPI_DEBUG

#define MINDENS   0.e-8

