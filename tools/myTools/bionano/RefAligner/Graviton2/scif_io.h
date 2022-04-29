#ifndef __SCIF_IO_H__
#define __SCIF_IO_H__

#ifdef __KERNEL__
#include <linux/types.h>
#include <linux/errno.h>
#include <linux/string.h>
#include <linux/delay.h>
//#include </usr/linux-k1om-4.7/linux-k1om/usr/include/scif.h>
//#include </home/users/vdergachev/src-bin/mpss-3.2.1/src/mpss-modules-3.2.1/include/scif.h>
#include </opt/mpss/3.8.4/sysroots/k1om-mpss-linux/usr/src/kernel/include/modules/scif.h>
#define printw(...)  printk(KERN_WARNING __VA_ARGS__)


int errno;

void *calloc(long a, long b)
{
  return vmalloc(a*b);
}

void perror(char *s)
{
  printw("Error %d\n", errno);
}

void sleep(int a)
{
  msleep(a*1000);
}

#define free(a) vfree(a)

#define SCIF_DEFAULT_REGISTER_FLAGS SCIF_MAP_KERNEL

#else
#include <string.h>
#include <errno.h>
#include <scif.h>
#define printw(...)  fprintf(stderr, __VA_ARGS__)
#define SCIF_DEFAULT_REGISTER_FLAGS 0

#endif

#ifndef RELEASE
#define RELEASE 0 /* If 1 : turns off some risky assertions. Should be 1 for customer release binaries, 0 in all other cases */
#endif

#define SCIF_IO_TYPE_NOP	0
#define SCIF_IO_TYPE_HOST_HELLO	1
#define SCIF_IO_TYPE_CLIENT_HELLO	2
#define SCIF_IO_TYPE_PING	3

#define SCIF_IO_TYPE_OPEN	16
#define SCIF_IO_TYPE_CLOSE	17
#define SCIF_IO_TYPE_WRITE	18
#define SCIF_IO_TYPE_READ	19
#define SCIF_IO_TYPE_SEEK	20
#define SCIF_IO_TYPE_UNLINK	21
#define SCIF_IO_TYPE_CHDIR	22

#define SCIF_IO_TYPE_TIFF_OPEN	36
#define SCIF_IO_TYPE_TIFF_CLOSE	37
#define SCIF_IO_TYPE_TIFF_IMG_INFO	38
#define SCIF_IO_TYPE_TIFF_IMG_READ	39
#define SCIF_IO_TYPE_TIFF_DIR_READ	40


#ifndef __KERNEL__
static int remote_port=2000;
#endif

typedef struct {
	int type;
	int fd;
	long long offset;
	long long length;
	long long aux1;
} SCIF_IO_MESSAGE;

typedef struct {
	scif_epd_t epd;
	long long transmit_window_size;
	long long transmit_window_free;
	off_t transmit_window_po;
	char *transmit_window;

	long long receive_window_size;
	long long receive_window_free;
	off_t receive_window_po;
	char *receive_window;
	
	off_t peer_po;
} SCIF_IO_CONTEXT;

#ifdef __cplusplus
extern "C" {
#endif

static inline int send_msg(SCIF_IO_CONTEXT *ctx, SCIF_IO_MESSAGE *msg)
{
  int nbytes;
// printw("Entering send_msg\n");
  while ((nbytes = scif_send(ctx->epd, msg, sizeof(*msg), SCIF_SEND_BLOCK)) >= 0) {
    if(nbytes == sizeof(*msg)) return 0;
    printw("Incomplete send message: nbytes=%d, should be %ld\n", nbytes, sizeof(*msg));
  }
  printw("Could not send message: nbytes=%d, should be %ld\n", nbytes, sizeof(*msg));
  return nbytes;
}

static inline int recv_msg(SCIF_IO_CONTEXT *ctx, SCIF_IO_MESSAGE *msg)
{
  int nbytes;
// printw("Entering recv_msg\n");
  while ((nbytes = scif_recv(ctx->epd, msg, sizeof(*msg), SCIF_RECV_BLOCK)) >= 0) {
    if(nbytes == sizeof(*msg)) return 0;
    if(nbytes<0)perror("recv_msg");
    printw( "Strange message received nbytes=%d\n", nbytes);
  }
  printw("Could not receive message: nbytes=%d, should be %ld\n", nbytes, sizeof(*msg));
  return nbytes;
}

static inline void expect_msg(SCIF_IO_CONTEXT *ctx, SCIF_IO_MESSAGE *msg, int type)
{
  int n;
  while(1) {
    if((n = recv_msg(ctx, msg)) == 0 && msg->type == type) return;
#ifdef __KERNEL__
    if(msg->type == -ERESTARTSYS) continue;
#endif
    if(msg->type != type)
      printw( "*** unexpected message type=%d (expected type=%d) fd=%d length=%lld aux1=%lld\n", msg->type, type, msg->fd, msg->length, msg->aux1);
    else
      printw( "*** Failed to receive message type=%d (expected type=%d) fd=%d length=%lld aux1=%lld\n", msg->type, type, msg->fd, msg->length, msg->aux1);
  }
}

#define SI_PAGE_SIZE 4096

static void init_buffers(SCIF_IO_CONTEXT *ctx)
{
  printw( "Initializing transmit/receive buffers\n");
  ctx->receive_window_free = 0;
 ctx->transmit_window_free = 0;
#ifndef __KERNEL__
 ctx->receive_window_size = 10*1024*1024;
 ctx->receive_window = (char *)((long long)calloc(ctx->receive_window_size+SI_PAGE_SIZE, 1) & ~(SI_PAGE_SIZE-1));

 ctx->transmit_window_size = 10*1024*1024;
 ctx->transmit_window = (char *)((long long)calloc(ctx->transmit_window_size+SI_PAGE_SIZE, 1) & ~(SI_PAGE_SIZE-1));
#else
/* Optimal transfer rate is achieved with large (10 MB) requests. However, the kernel only ever issues requests of at most 4kb size */
 ctx->receive_window_size = 16*1024;
 ctx->receive_window=kmalloc(ctx->receive_window_size, GFP_KERNEL);

 ctx->transmit_window_size = 16*1024;
 ctx->transmit_window=kmalloc(ctx->transmit_window_size, GFP_KERNEL);
#endif
 printw( "Recieve =0x%016llx Transmit =0x%016llx\n", (long long)ctx->receive_window, (long long)ctx->transmit_window);

 ctx->receive_window_po = scif_register(ctx->epd, ctx->receive_window, ctx->receive_window_size, 0, SCIF_PROT_WRITE, SCIF_DEFAULT_REGISTER_FLAGS);
 if(ctx->receive_window_po == -1) perror("receive_window");
 ctx->transmit_window_po = scif_register(ctx->epd, ctx->transmit_window, ctx->transmit_window_size, 0, SCIF_PROT_READ, SCIF_DEFAULT_REGISTER_FLAGS);
 if(ctx->transmit_window_po == -1) perror("transmit_window");
 printw( "Recieve PO=0x%016lx Transmit PO=0x%016lx\n", ctx->receive_window_po, ctx->transmit_window_po);
}

static void free_buffers(SCIF_IO_CONTEXT *ctx)
{
  scif_unregister(ctx->epd, ctx->receive_window_po, ctx->receive_window_size);
  scif_unregister(ctx->epd, ctx->transmit_window_po, ctx->transmit_window_size);
  free(ctx->receive_window);
  free(ctx->transmit_window);

  ctx->receive_window = NULL;
  ctx->receive_window_free = 0;
  ctx->receive_window_size = 0;
  ctx->transmit_window = NULL;
  ctx->transmit_window_free = 0;
  ctx->transmit_window_size = 0;
}

#ifndef __KERNEL__
static inline int get_port(void)
{
  char *s;
  int port;
  s = getenv("SCIF_PORT");
  if(s != NULL) {
    port = atoi(s);
    if(port > 0) return port;
  }
  return -1;
}

static inline int get_remote_port(void)
{
  char *s;
  int port;
  s = getenv("SCIF_DAEMON_PORT");
  if(s != NULL) {
    port = atoi(s);
    if(port > 0) return port;
  }
//port=remote_port+getuid();
  port = remote_port;
  return(port);
}

static inline long get_cookie(void)
{
  long cookie;
  char *s;
  /* Not terribly secure.. */
  s = getenv("SCIF_USERID_COOKIE");
  if(s !=NULL) {
    cookie = atol(s);
    if(cookie > 0)return cookie;
  }
  /* Hardcode for now */
  cookie = 124455 + getuid();
  return cookie;
}

#else

int get_remote_port(void);
static inline long get_cookie(void)
{
  long cookie;
  char *s;
  /* Hardcode for now */
  cookie = 124455+0;
 return cookie;
}

#endif

/*
 * Change current directory to specified path. 
 * If path is null synchronize with the current 
 * directory of the calling process 
 * 
 */
static int si_chdir(SCIF_IO_CONTEXT *ctx, const char *path)
{
  SCIF_IO_MESSAGE msg;
  msg.type = SCIF_IO_TYPE_CHDIR;
  msg.fd = 0;
  if(path != NULL) msg.length = strlen(path)+1;

#ifndef __KERNEL__
#pragma omp critical(sic)
#endif
  {
    if(path != NULL) {
      memcpy(ctx->transmit_window, path, msg.length);
    } else { /* special usage: path=NULL uses current directory */
#ifdef __KERNEL__
      ctx->transmit_window[0] = "/";
      ctx->transmit_window[1] = 0;
#else
      if(getcwd(ctx->transmit_window, ctx->transmit_window_size) == NULL) {
	perror("Failed to obtain current working directory");
	msg.fd = -ERANGE;
	msg.length = errno;
      } else {
	msg.length = strlen(ctx->transmit_window) + 1;
      }
#endif
    }
    if(msg.fd == 0) {
      scif_writeto(ctx->epd, ctx->transmit_window_po, msg.length, ctx->peer_po, SCIF_RMA_SYNC);

      send_msg(ctx, &msg);
      expect_msg(ctx, &msg, SCIF_IO_TYPE_CHDIR);
    }
  }

  if(msg.fd < 0) {
    errno = msg.length;
  }

  return(msg.fd);
}

static inline SCIF_IO_CONTEXT * si_init(void)
{
  struct scif_portID portID;
  int conn_port, err;
  SCIF_IO_CONTEXT *ctx = (SCIF_IO_CONTEXT *)calloc(1, sizeof(*ctx));
  if(ctx == NULL) return NULL;

  if ((ctx->epd = scif_open()) == SCIF_OPEN_FAILED) {
    printw("scif_open failed with error %d\n", (int)ctx->epd);
    free(ctx);
    return(NULL);
  }

  while((conn_port = scif_bind(ctx->epd, 0)) < 0) {
    printw("host: scif_bind to port %d failed with error %d\n", 0, conn_port);
    sleep(1);
  }
  printw("host: scif_bind to port %d success\n", conn_port);


#ifndef __KERNEL
  char port_name[32];
  sprintf(port_name,"%d",conn_port);

  errno = 0;
  if(setenv("SCIF_PORT",port_name,1) < 0){
    if(!RELEASE){
      int eno = errno;
      char *err = strerror(eno);
      printf("si_init(): failed to put SCIF_PORT=%s using setenv: errno=%d,%s\n",port_name,eno,err);
      fflush(stdout);
    }
  }
#endif

  portID.node = 0;
  portID.port = get_remote_port();

  int tries = 61;
  while ((err = scif_connect(ctx->epd, &portID)) < 0) {
    if (tries > 0) {
      printw("connection to node %d port %d failed : trial %d (errno=%d err=%d)\n", portID.node, portID.port, tries, errno, err);
      tries--;
#ifndef __KERNEL__
      if(!RELEASE){
	printf("si_init():connection to node %d port %d failed : trial %d (errno=%d err=%d)\n", portID.node, portID.port, tries, errno, err);
	fflush(stdout);
      }
#endif
      sleep(1);
      continue;
    }
    printw("scif_connect failed with error %d (err=%d)\n", errno, err);
#ifndef __KERNEL__
    if(!RELEASE){
      printf("scif_connect failed with error %d (err=%d)\n", errno, err);
      fflush(stdout);
    }
#endif
    scif_close(ctx->epd);
    free(ctx);
    return(NULL);
  }
  printw( "connect to node %d success\n", portID.node);

  init_buffers(ctx);

  SCIF_IO_MESSAGE msg;
  msg.type = SCIF_IO_TYPE_CLIENT_HELLO;
#ifdef __KERNEL__
  msg.fd=0;
#else
  msg.fd=getuid();
#endif

  msg.offset = get_cookie();
  msg.length = ctx->receive_window_po;

  send_msg(ctx, &msg);

  expect_msg(ctx, &msg, SCIF_IO_TYPE_HOST_HELLO);

  ctx->peer_po = msg.length;

  if(msg.fd != 0) {
    printw( "handshake failed, code=%d\n", msg.fd);
    free_buffers(ctx);
    scif_close(ctx->epd);
    free(ctx);
    return(NULL);
  }

  si_chdir(ctx, NULL);
  
  printw( "handshake successfull\n");
  return(ctx);
}

static void si_destroy_context(SCIF_IO_CONTEXT *ctx)
{
  if(ctx != NULL) {
    free_buffers(ctx);
    scif_close(ctx->epd);
    free(ctx);
  }
}
	
static int si_open(SCIF_IO_CONTEXT *ctx, const char *filename, int flags)
{
  SCIF_IO_MESSAGE msg;
  msg.type = SCIF_IO_TYPE_OPEN;
  msg.fd = flags;
  msg.length = strlen(filename)+1;

  #pragma omp critical(sic)
  {
    memcpy(ctx->transmit_window, filename, msg.length);
    scif_writeto(ctx->epd, ctx->transmit_window_po, msg.length, ctx->peer_po, SCIF_RMA_SYNC);

    send_msg(ctx, &msg);
    expect_msg(ctx, &msg, SCIF_IO_TYPE_OPEN);
  }

  if(msg.fd < 0) {
    errno = msg.length;// NOTE : scif_daemon.c replaced the msg.length with host errno value
  }

  /* Return host fd */
  return(msg.fd);
}

static int si_unlink(SCIF_IO_CONTEXT *ctx, const char *filename)
{
  SCIF_IO_MESSAGE msg;
  msg.type = SCIF_IO_TYPE_UNLINK;
  msg.fd = 0;
  msg.length = strlen(filename)+1;

  #pragma omp critical(sic)
  {
    memcpy(ctx->transmit_window, filename, msg.length);
    scif_writeto(ctx->epd, ctx->transmit_window_po, msg.length, ctx->peer_po, SCIF_RMA_SYNC);

    send_msg(ctx, &msg);
    expect_msg(ctx, &msg, SCIF_IO_TYPE_UNLINK);
  }

  if(msg.fd < 0) {
    errno = msg.length;
  }

  return(msg.fd);
}


static void si_close(SCIF_IO_CONTEXT *ctx, int fd)
{
  SCIF_IO_MESSAGE msg;
  msg.type = SCIF_IO_TYPE_CLOSE;
  msg.fd = fd;
  msg.length = 0;
  #pragma omp critical(sic)
  {
    send_msg(ctx, &msg);
  }
}

static size_t si_read(SCIF_IO_CONTEXT *ctx, int fd, void *buf, size_t count)
{
  SCIF_IO_MESSAGE msg;
  msg.type = SCIF_IO_TYPE_READ;
  msg.fd = fd;
  msg.length = count;
  if(msg.length > ctx->receive_window_size) msg.length = ctx->receive_window_size;

  #pragma omp critical(sic)
  {
    send_msg(ctx, &msg);
    // printw("send_msg complete length=%d\n", msg.length);

    expect_msg(ctx, &msg, SCIF_IO_TYPE_READ);
    // printw("expect_msg complete length=%d\n", msg.length);

    if(msg.length > 0) memcpy(buf, ctx->receive_window, msg.length);
  }
  return(msg.length);
}

static size_t si_seek(SCIF_IO_CONTEXT *ctx, int fd, off_t offset, int whence)
{
  SCIF_IO_MESSAGE msg;
  msg.type = SCIF_IO_TYPE_SEEK;
  msg.fd = fd;
  msg.length = offset;
  msg.aux1 = whence;
  #pragma omp critical(sic)
  {
    send_msg(ctx, &msg);
    // printw("seek send_msg complete length=%d\n", msg.length);

    expect_msg(ctx, &msg, SCIF_IO_TYPE_SEEK);
    // printw("seek expect_msg complete length=%d\n", msg.length);
  }
  return(msg.length);
}

static size_t si_write(SCIF_IO_CONTEXT *ctx, int fd, void *buf, size_t count)
{
  SCIF_IO_MESSAGE msg;
  msg.type = SCIF_IO_TYPE_WRITE;
  msg.fd = fd;
  msg.length = count;
  if(msg.length > ctx->transmit_window_size) msg.length = ctx->transmit_window_size;

  #pragma omp critical(sic)
  {
    memcpy(ctx->transmit_window, buf, msg.length);

    scif_writeto(ctx->epd, ctx->transmit_window_po, msg.length, ctx->peer_po, SCIF_RMA_SYNC);

    send_msg(ctx, &msg);

    expect_msg(ctx, &msg, SCIF_IO_TYPE_WRITE);
  }

  return(msg.length);
}

#ifndef __KERNEL__
FILE *sic_fopen(const char *path, const char *mode);
int sic_unlink(const char *path);
SCIF_IO_CONTEXT *get_main_context(void);
#include <fcntl.h>
#include <complex.h>

/* TIFF file I/O. R/O so far */

static int si_tiff_open(SCIF_IO_CONTEXT *ctx, const char *filename)
{
  SCIF_IO_MESSAGE msg;
  msg.type = SCIF_IO_TYPE_TIFF_OPEN;
  msg.fd = O_RDONLY;
  msg.length = strlen(filename)+1;

  #pragma omp critical(sic)
  {
    memcpy(ctx->transmit_window, filename, msg.length);
    scif_writeto(ctx->epd, ctx->transmit_window_po, msg.length, ctx->peer_po, SCIF_RMA_SYNC);

    send_msg(ctx, &msg);
    expect_msg(ctx, &msg, SCIF_IO_TYPE_TIFF_OPEN);
  }

  if(msg.fd < 0) {
    errno = msg.length;
  }

  /* Return host fd */
  return(msg.fd);
}

static void si_tiff_close(SCIF_IO_CONTEXT *ctx, int fd)
{
  SCIF_IO_MESSAGE msg;
  msg.type = SCIF_IO_TYPE_TIFF_CLOSE;
  msg.fd = fd;
  msg.length = 0;
  #pragma omp critical(sic)
  {
    send_msg(ctx, &msg);
  }
}

static int si_tiff_dir_read(SCIF_IO_CONTEXT *ctx, int fd)
{
  SCIF_IO_MESSAGE msg;
  msg.type = SCIF_IO_TYPE_TIFF_DIR_READ;
  msg.fd = fd;
  msg.length = 0;

  #pragma omp critical(sic)
  {
    send_msg(ctx, &msg);
    expect_msg(ctx, &msg, SCIF_IO_TYPE_TIFF_DIR_READ);
  }

  if(msg.fd < 0) {
    errno = msg.length;
  }

  /* Return host fd */
  return(msg.length);
}

static long long si_tiff_img_info(SCIF_IO_CONTEXT *ctx, int fd, long long attribute)
{
  SCIF_IO_MESSAGE msg;
  msg.type = SCIF_IO_TYPE_TIFF_IMG_INFO;
  msg.fd = fd;
  msg.offset = attribute;
  msg.length = 0;

  #pragma omp critical(sic)
  {
    send_msg(ctx, &msg);
    expect_msg(ctx, &msg, SCIF_IO_TYPE_TIFF_IMG_INFO);
  }

  if(msg.fd < 0) {
    errno = msg.length;
  }

  /* Return host fd */
  return(msg.length);
}

static size_t si_tiff_img_read(SCIF_IO_CONTEXT *ctx, int fd, double complex  *buf, int stride, int x0, int x1, int y0, int y1)
{
  int i,j,step;
  long long len, a;
  uint16_t *data;
  SCIF_IO_MESSAGE msg;
  msg.type = SCIF_IO_TYPE_TIFF_IMG_READ;
  msg.fd = fd;
  msg.offset = (((long long) x1) << 32) | (x0);
  msg.length = (((long long) y1) << 32) | (y0);

//   fprintf(stderr, "%llx %llx\n", msg.offset, msg.length);
  
  if(2*(y1-y0)*(x1-x0) > ctx->receive_window_size) {
	  step=(ctx->receive_window_size-1024)/(2*(x1-x0));
	  len=0;
	  for(j=y0;j<y1;j+=step) {
		  a=si_tiff_img_read(ctx, fd, &(buf[stride*j]), stride, x0, x1, j, (j+step>y1) ? y1 : (j+step));
		  if(a<0) return(a);
		  len+=a;
		}
	  return(len);
	}

  #pragma omp critical(sic)
  {
    send_msg(ctx, &msg);
    // printw("send_msg complete length=%d\n", msg.length);

    expect_msg(ctx, &msg, SCIF_IO_TYPE_TIFF_IMG_READ);
    // printw("expect_msg complete length=%d\n", msg.length);

    if(msg.length > 0) {
	    data=(uint16_t *)ctx->receive_window;
// 	    fprintf(stderr, "Copying data length=%lld\n", msg.length);
	    
	    for(j=0;j<y1-y0;j++) 
		    for(i=0;i<(x1-x0);i++) 
			    buf[stride*j+i]=data[j*(x1-x0)+i];
	    
// 	    fprintf(stderr, "Finished copying data\n");
	}
  }
  return(msg.length);
}

static size_t si_tiff_img_readF(SCIF_IO_CONTEXT *ctx, int fd, float complex  *buf, int stride, int x0, int x1, int y0, int y1)
{
  int i,j,step;
  long long len, a;
  uint16_t *data;
  SCIF_IO_MESSAGE msg;
  msg.type = SCIF_IO_TYPE_TIFF_IMG_READ;
  msg.fd = fd;
  msg.offset = (((long long) x1) << 32) | (x0);
  msg.length = (((long long) y1) << 32) | (y0);

//   fprintf(stderr, "%llx %llx\n", msg.offset, msg.length);

  if(2*(y1-y0)*(x1-x0) > ctx->receive_window_size) {
          step=(ctx->receive_window_size-1024)/(2*(x1-x0));
          len=0;
          for(j=y0;j<y1;j+=step) {
                  a=si_tiff_img_readF(ctx, fd, &(buf[stride*j]), stride, x0, x1, j, (j+step>y1) ? y1 : (j+step));
                  if(a<0) return(a);
                  len+=a;
                }
          return(len);
        }

  #pragma omp critical(sic)
  {
    send_msg(ctx, &msg);
    // printw("send_msg complete length=%d\n", msg.length);

    expect_msg(ctx, &msg, SCIF_IO_TYPE_TIFF_IMG_READ);
    // printw("expect_msg complete length=%d\n", msg.length);

    if(msg.length > 0) {
            data=(uint16_t *)ctx->receive_window;
//          fprintf(stderr, "Copying data length=%lld\n", msg.length);

            for(j=0;j<y1-y0;j++)
                    for(i=0;i<(x1-x0);i++)
                            buf[stride*j+i]=data[j*(x1-x0)+i];

//          fprintf(stderr, "Finished copying data\n");
        }
  }
  return(msg.length);
}

static size_t si_tiff_img_readFS(SCIF_IO_CONTEXT *ctx, int fd, float *buf, int stride, int x0, int x1, int y0, int y1)
{
  int i,j,step;
  long long len, a;
  uint16_t *data;
  SCIF_IO_MESSAGE msg;
  msg.type = SCIF_IO_TYPE_TIFF_IMG_READ;
  msg.fd = fd;
  msg.offset = (((long long) x1) << 32) | (x0);
  msg.length = (((long long) y1) << 32) | (y0);

//   fprintf(stderr, "%llx %llx\n", msg.offset, msg.length);

  if(2*(y1-y0)*(x1-x0) > ctx->receive_window_size) {
          step=(ctx->receive_window_size-1024)/(2*(x1-x0));
          len=0;
          for(j=y0;j<y1;j+=step) {
                  a=si_tiff_img_readFS(ctx, fd, &(buf[stride*j]), stride, x0, x1, j, (j+step>y1) ? y1 : (j+step));
                  if(a<0) return(a);
                  len+=a;
                }
          return(len);
        }

  #pragma omp critical(sic)
  {
    send_msg(ctx, &msg);
    // printw("send_msg complete length=%d\n", msg.length);

    expect_msg(ctx, &msg, SCIF_IO_TYPE_TIFF_IMG_READ);
    // printw("expect_msg complete length=%d\n", msg.length);

    if(msg.length > 0) {
            data=(uint16_t *)ctx->receive_window;
//          fprintf(stderr, "Copying data length=%lld\n", msg.length);

            for(j=0;j<y1-y0;j++)
                    for(i=0;i<(x1-x0);i++)
                            buf[stride*j+i]=data[j*(x1-x0)+i];

//          fprintf(stderr, "Finished copying data\n");
        }
  }
  return(msg.length);
}

#endif

#ifdef __cplusplus
}
#endif

#endif
