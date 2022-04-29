#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

char *usage = (char *)"usage:%s filename increment [prefix suffix maxID]\n";

/** lock the file, read the two integers in the file, increment the first integer by the specified increment and the 2nd integer by 1, then unlock the file */
/** If prefix,suffix are specified, modify how the first integer is incremented:
   A. Repeatedly increment the first integer N by 1 then check if a non-empty file <prefix><N><suffix> exists. If not don't count the current increment.
   B. Continue until N has been increment the specified number of times, counting only increments with a non-empty file, but don't increment N past maxID
   C. If the original value of N is already >= maxID, increment N by 1 (used to signal completion of all jobs in the shell script that calls nextID
 **/


/** Returns the final value of N. Updates Group to the updated 2nd integer value */
static int nextContigId(char *filename, int increment, int &Group, int &origID, char *prefix, char* suffix, int maxID)
{
  int fd = open(filename,O_RDWR);
  if(fd < 0){
    char *err = strerror(errno);
    printf("open(%s) for read/write failed: errno=%d(%s)\n", filename, errno, err);
    exit(1);
  }
  if(increment <= 0){
    printf("increment=%d must be +ve integer\n",increment);
    exit(1);
  }

  /* lock the file for exclusive read */
#if 0 // OLD : only works if all jobs run on same server
  if(flock(fd, LOCK_EX) < 0){
    char *err = strerror(errno);
    printf("flock of %s failed: errno=%d(%s)\n", filename, errno, err);
    exit(1);
  }
#else // use fcntl which should work on NFS mounted files shared across multiple servers
  struct flock lock;
  lock.l_type = F_WRLCK;
  lock.l_whence = SEEK_SET;
  lock.l_start = 0;
  lock.l_len = 16;
  if(fcntl(fd,F_SETLKW,&lock) < 0){
    char *err = strerror(errno);
    printf("fcntl(SETLKW,F_WRLCK) of %s failed: errno=%d(%s)\n",filename,errno,err);
    exit(1);
  }
#endif

  /* read/update integer in file */
  char buf[BUFSIZ];
  ssize_t cnt = read(fd,buf,BUFSIZ);
  if(cnt < 0){
    char *err = strerror(errno);
    printf("read of %s failed: errno=%d(%s)\n", filename, errno, err);
    exit(1);
  }
  if(!cnt){
    printf("Could not read integer in file %s (file empty)\n",filename);
    exit(1);
  }
  char *pt = buf, *qt;
  int ID = strtol(pt,&qt,10);
  if(pt == qt || !(*qt==0 || isspace(*qt))){
    printf("ID value (first integer) in %s has invalid syntax: %s\n", filename, buf);
    exit(1);
  }
  if(ID < 0){
    printf("ID value in %s of %d is invalid (must be positive after adding 1)\n",filename,ID);
    exit(1);
  }
  origID = ID;

  int G = 0;

  if(prefix && suffix){    /* keep incrementing ID, until the number of non-empty files <prefix><ID><suffix> equals "increment" (or ID >= maxID) */
    pt = qt;
    G = strtol(pt,&qt,10);
    if(pt == qt || !(*qt==0 || isspace(*qt))){
      printf("Group value (second integer) in %s has invalid syntax: %s\n", filename, buf);
      exit(1);
    }
    if(G < 0){
      printf("Group value in %s of %d is invalid (must be positive after adding 1)\n",filename,G);
      exit(1);
    }
    int origG = G;
    char filename[BUFSIZ];
    int cnt = 0;
    if(origID >= maxID)
      ID++;
    else {
      G++;
      while(cnt < increment && ID < maxID){
	ID++;
	sprintf(filename,"%s%d%s",prefix,ID,suffix);
	struct stat sbuf;
	if(!stat(filename,&sbuf) && sbuf.st_size > 0)
	  cnt++;
      }
    }
    Group = G;

    fprintf(stderr,"ID=%d->%d,G=%d->%d,maxID=%d\n",origID,ID,origG,G,maxID);//HERE
    fflush(stderr);
  } else
    ID += increment;

  off_t offset = lseek(fd,0,SEEK_SET);
  if(offset < 0){
    char *err = strerror(errno);
    printf("lseek of %s to start of file failed: errno=%d(%s)\n", filename, errno, err);
    exit(1);
  }
  if(prefix && suffix)
    sprintf(buf,"%d %d\n",ID,G);
  else
    sprintf(buf,"%d\n",ID);
  cnt = write(fd,buf,strlen(buf));
  if(cnt < 0){
    char *err = strerror(errno);
    printf("write of %d to %s failed: errno=%d(%s)\n", ID, filename, errno, err);
    exit(1);
  }
  if(cnt != (ssize_t)strlen(buf)){
    printf("write to %s failed (bytes written=%ld):%s\n", filename, (long int)cnt, buf);
    exit(1);
  }

  /* unlock file */
#if 0 // OLD method, only works with all jobs on same server
  if(flock(fd,LOCK_UN) < 0){
    char *err = strerror(errno);
    printf("unlock of %s failed: errno=%d(%s)\n", filename, errno, err);
    exit(1);
  }
#else // use F_SETLKW
  // No need to unlock, since close(fd) releases all locks
#endif

  if(close(fd) < 0){
    char *err = strerror(errno);
    printf("close of file descriptor to %s failed: errno=%d(%s)\n", filename, errno, err);
    exit(1);
  }
  
  return ID;
}

int main(int argc, char **argv)
{
  if(argc < 3){
    printf(usage,argv[0]);
    exit(1);
  }

  int increment = atoi(argv[2]);

  char *prefix = 0;
  char *suffix = 0;
  int maxID = 0;
  if(argc == 6){
    prefix = argv[3];
    suffix = argv[4];
    maxID = atoi(argv[5]);
  }

  int Group, origID;
  int ID = nextContigId(argv[1], increment, Group, origID, prefix, suffix, maxID);
  if(prefix && suffix)
    printf("%d %d %d\n",Group, origID, ID);
  else
    printf("%d\n",ID);
  exit(0);
}

