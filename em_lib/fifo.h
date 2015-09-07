#ifndef FIFO_H
#define FIFO_H
#ifdef _cplusplus
extern "C"
{
#endif
	#ifndef  FIFO_SIZE
	#define FIFO_SIZE 8192 
	#endif
	typedef struct {
		 unsigned char * buffer; /* the buffer holding the data */
		 unsigned int in;		/* data is added at offset (in % size) */
		 unsigned int out;
		 unsigned int size; /* the size of the allocated buffer */
	} fifo_t;


	double * get_fifo_buf(fifo_t * f);

	void fifo_init(fifo_t * f, unsigned char * buffer, unsigned int size);
	unsigned int fifo_read(fifo_t * fifo, unsigned char *buffer, unsigned int len);
	unsigned int fifo_write(fifo_t * fifo, const unsigned char * buffer, unsigned int len);
	unsigned int fifo_read_steps(fifo_t * fifo, unsigned char *buffer, unsigned int len,unsigned int steps);

	/**
	 * fifo_len - returns the number of bytes available in the FIFO, no locking version
	 * @fifo: the fifo to be used.
	 */
	static __inline unsigned int fifo_len(fifo_t *fifo)
	{
		 return fifo->in - fifo->out;
	}

    typedef struct {
		 int * buffer; /* the buffer holding the data */
		 unsigned int in;		/* data is added at offset (in % size) */
		 unsigned int out;
		 unsigned int size; /* the size of the allocated buffer */
	} int_fifo_t;

	int * get_int_fifo_buf(int_fifo_t * f);

	void int_fifo_init(int_fifo_t *f, int *buffer, unsigned int size);
	unsigned int int_fifo_read(int_fifo_t *fifo, int *buffer, unsigned int len);
	unsigned int int_fifo_write(int_fifo_t *fifo, const int * buffer, unsigned int len);
	unsigned int int_fifo_read_steps(int_fifo_t *fifo, int *buffer, unsigned int len,unsigned int steps);

	/**
	 * fifo_len - returns the number of bytes available in the FIFO, no locking version
	 * @fifo: the fifo to be used.
	 */
	static __inline unsigned int int_fifo_len(int_fifo_t *fifo)
    {
		 return fifo->in - fifo->out;
	}
#ifdef _cplusplus
}
#endif
#endif
