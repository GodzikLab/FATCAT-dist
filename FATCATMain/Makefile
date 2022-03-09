CC = g++
CCOPT = -O2 -Wall
LIBS = -lm
INCDIRS = -I./
CFLAGS = ${CCOPT} ${INCDIRS}

#----------------all programs----------------------------
all: FATCAT FATCATDB PostFATCAT FATPAT FATCAT2Tree

#-----------------Objects--------------------------------
FUN0 = geometry.o PsShow.o matrix.o
ALN =  Align0.o Amino.o

#---------------include dependencies--------------------
geometry.o : geometry.h
PsShow.o	: PsShow.h basic.h 
Prot.o		: Prot.h basic.h tempkit.h
AFPchain.o	: AFPchain.h SigEva.h Prot.h AFP.h geometry.h tempkit.h matrix.h
FATCAT.o	: AFPchain.h 
Align0.o	: Align0.h tempkit.h
Amino.o		: Amino.h
SAlnOpt.o	: SAlnOpt.h Align0.h tempkit.h

#---------------make programs----------------------------

.C.o:
	$(CC) $(CFLAGS) -c $<

OBJ_FATCAT = FATCAT.o SAlnOpt.o SigEva.o Prot.o AFPchain.o $(ALN) $(FUN0)
FATCAT :   $(OBJ_FATCAT)
	$(CC) $(CCOPT) -o FATCAT $(OBJ_FATCAT) $(LIBS)

OBJ_FATCATDB = FATCATDBSearch.o FATCATDB.o SigEva.o $(FUN0)
FATCATDB : $(OBJ_FATCATDB)
	$(CC) $(CCOPT) -o FATCATDB $(OBJ_FATCATDB) $(LIBS)

OBJ_PostFATCAT = PostFATCAT.o PostAlign.o Prot.o $(ALN) $(FUN0)
PostFATCAT : $(OBJ_PostFATCAT)
	$(CC) $(CCOPT) -o PostFATCAT $(OBJ_PostFATCAT)

OBJ_FATPAT = FatPatMain.o FatDom.o Prot.o Amino.o $(FUN0)
FATPAT : $(OBJ_FATPAT)
	$(CC) $(CCOPT) -o FATPAT $(OBJ_FATPAT) $(LIBS)

OBJ_FATTREE = FATCAT2Tree.o SAlnOpt.o SigEva.o Prot.o AFPchain.o hiCluster.o $(ALN) $(FUN0)
FATCAT2Tree :   $(OBJ_FATTREE)
	$(CC) $(CCOPT) -o FATCAT2Tree $(OBJ_FATTREE) $(LIBS)


clean:
	rm FATCAT FATCATDB PostFATCAT FATPAT
	rm *.o
