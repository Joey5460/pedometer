#include "fifo.h"
#include <math.h>
#include <string.h>

//This initializes the FIFO structure with the given buffer and size
void fifo_init(fifo_t * f, unsigned char * buffer, unsigned int size){
     f->in = f->out =0;
     f->size = size;
     f->buffer = buffer;
}

double * get_fifo_buf(fifo_t * f)
{
    return f ->buffer;
}
/**
 * fifo_read - gets some data from the FIFO, no locking version
 * @fifo: the fifo to be used.
 * @buffer: where the data must be copied.
 * @len: the size of the destination buffer.
 *
 * This function copies at most @len bytes from the FIFO into the
 * @buffer and returns the number of copied bytes.
 *
 * Note that with only one concurrent reader and one concurrent
 * writer, you don't need extra locking to use these functions.
 */
unsigned int fifo_read(fifo_t * fifo, unsigned char *buffer, unsigned int len)
{
    
        unsigned int l;

        len = (int)fmin(len, fifo->in - fifo->out);

         /*
          * Ensure that we sample the fifo->in index -before- we
          * start removing bytes from the kfifo.
          */

         //smp_rmb();

         /* first get the data from fifo->out until the end of the buffer */
         l = (int)fmin(len, fifo->size - (fifo->out & (fifo->size - 1)));
         memcpy(buffer, fifo->buffer + (fifo->out & (fifo->size - 1)), l);

         /* then get the rest (if any) from the beginning of the buffer */
         memcpy(buffer + l, fifo->buffer, (len - l));

         /*
          * Ensure that we remove the bytes from the kfifo -before-
          * we update the fifo->out index.
         */

         //smp_mb();

         fifo->out += len;

         return len;
}
unsigned int fifo_read_steps(fifo_t * fifo, unsigned char *buffer, unsigned int len,unsigned int steps)
{
    
        unsigned int l;

        len = (int)fmin(len, fifo->in - fifo->out);

         /*
          * Ensure that we sample the fifo->in index -before- we
          * start removing bytes from the kfifo.
          */

         //smp_rmb();

         /* first get the data from fifo->out until the end of the buffer */
         l = (int)fmin(len, fifo->size - (fifo->out & (fifo->size - 1)));
         memcpy(buffer, fifo->buffer + (fifo->out & (fifo->size - 1)), l);

         /* then get the rest (if any) from the beginning of the buffer */
         memcpy(buffer + l, fifo->buffer, (len - l));

         /*
          * Ensure that we remove the bytes from the kfifo -before-
          * we update the fifo->out index.
         */

         //smp_mb();

         fifo->out += steps;

         return len;
}
/**
 * fifo_write - puts some data into the FIFO, no locking version
 * @fifo: the fifo to be used.
 * @buffer: the data to be added.
 * @len: the length of the data to be added.
 *
 * This function copies at most @len bytes from the @buffer into
 * the FIFO depending on the free space, and returns the number of
 * bytes copied.
 *
 * Note that with only one concurrent reader and one concurrent
 * writer, you don't need extra locking to use these functions.
 */
unsigned int fifo_write(fifo_t * fifo, const unsigned char * buffer, unsigned int len)
{

     unsigned int l;

     len = (int)fmin(len, fifo->size - fifo->in + fifo->out);

     /*
      * Ensure that we sample the fifo->out index -before- we
      * start putting bytes into the kfifo.
      */

     //smp_mb();

     /* first put the data starting from fifo->in to buffer end */
     l = (int)fmin(len, fifo->size - (fifo->in & (fifo->size - 1)));
     memcpy(fifo->buffer + (fifo->in & (fifo->size - 1)), buffer, l);

     /* then put the rest (if any) at the beginning of the buffer */
     memcpy(fifo->buffer, buffer + l, (len - l));

     /*
      * Ensure that we add the bytes to the kfifo -before-
      * we update the fifo->in index.
      */

     //smp_wmb();

     fifo->in += len;

     return len;
}

//This initializes the FIFO structure with the given buffer and size
void int_fifo_init(int_fifo_t * f, int* buffer, unsigned int size){
     f->in = f->out =0;
     f->size = size;
     f->buffer = buffer;
}

int * get_int_fifo_buf(int_fifo_t * f)
{
    return f ->buffer;
}

/* fifo_read - gets some data from the FIFO, no locking version
 * @fifo: the fifo to be used.
 * @buffer: where the data must be copied.
 * @len: the size of the destination buffer.
 *
 * This function copies at most @len bytes from the FIFO into the
 * @buffer and returns the number of copied bytes.
 *
 * Note that with only one concurrent reader and one concurrent
 * writer, you don't need extra locking to use these functions.
 */
unsigned int int_fifo_read(int_fifo_t * fifo, int *buffer, unsigned int len)
{
        unsigned int l;

        len = (int)fmin(len, fifo->in - fifo->out);

         /*
          * Ensure that we sample the fifo->in index -before- we
          * start removing bytes from the kfifo.
          */

         //smp_rmb();

         /* first get the data from fifo->out until the end of the buffer */
         l = (int)fmin(len, fifo->size - (fifo->out & (fifo->size - 1)));
         memcpy(buffer, fifo->buffer + (fifo->out & (fifo->size - 1)), l*sizeof(int));

         /* then get the rest (if any) from the beginning of the buffer */
         memcpy(buffer + l, fifo->buffer, (len - l)*sizeof(int));

         /*
          * Ensure that we remove the bytes from the kfifo -before-
          * we update the fifo->out index.
         */

         //smp_mb();

         fifo->out += len;

         return len;
}
unsigned int int_fifo_read_steps(int_fifo_t * fifo, int *buffer, unsigned int len,unsigned int steps)
{
        unsigned int l;

        len = (int)fmin(len, fifo->in - fifo->out);

         /*
          * Ensure that we sample the fifo->in index -before- we
          * start removing bytes from the kfifo.
          */

         //smp_rmb();

         /* first get the data from fifo->out until the end of the buffer */
         l = (int)fmin(len, fifo->size - (fifo->out & (fifo->size - 1)));
         memcpy(buffer, fifo->buffer + (fifo->out & (fifo->size - 1)), l*sizeof(int));

         /* then get the rest (if any) from the beginning of the buffer */
         memcpy(buffer + l, fifo->buffer, (len - l)*sizeof(int));

         /*
          * Ensure that we remove the bytes from the kfifo -before-
          * we update the fifo->out index.
         */

         //smp_mb();

         fifo->out += steps;

         return len;
}
/*
 * fifo_write - puts some data into the FIFO, no locking version
 * @fifo: the fifo to be used.
 * @buffer: the data to be added.
 * @len: the length of the data to be added.
 *
 * This function copies at most @len bytes from the @buffer into
 * the FIFO depending on the free space, and returns the number of
 * bytes copied.
 *
 * Note that with only one concurrent reader and one concurrent
 * writer, you don't need extra locking to use these functions.
 */
unsigned int int_fifo_write(int_fifo_t * fifo, const int * buffer, unsigned int len)
{
     unsigned int l;

     len = (int)fmin(len, fifo->size - fifo->in + fifo->out);

     /*
      * Ensure that we sample the fifo->out index -before- we
      * start putting bytes into the kfifo.
      */

     //smp_mb();

     /* first put the data starting from fifo->in to buffer end */
     l = (int)fmin(len, fifo->size - (fifo->in & (fifo->size - 1)));
     memcpy(fifo->buffer + (fifo->in & (fifo->size - 1)), buffer, l*sizeof(int));

     /* then put the rest (if any) at the beginning of the buffer */
     memcpy(fifo->buffer, buffer + l, (len - l)*sizeof(int));

     /*
      * Ensure that we add the bytes to the kfifo -before-
      * we update the fifo->in index.
      */

     //smp_wmb();

     fifo->in += len;

     return len;
}
