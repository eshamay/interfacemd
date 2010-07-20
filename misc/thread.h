#ifndef THREAD_H_
#define THREAD_H_

#include <pthread.h>
#include <stdio.h>
#include <string>
#include <vector>

typedef void* (*fct_ptr)(void *);

template <class T>
class Thread {

public:
	Thread () {}
	// create a thread that will run the function and the argument given
	Thread (fct_ptr f, void * arg);

	void* Join ();
	void * Return () const { return _return; }
	pthread_t * GetThread () const { return _thread; }

	static int num_threads;
	
protected:

	// the function we're going to run
	void * (* _function) (void *);
	// and the argument we receive to run it;
	T * _arg;	

	// the thread itself
	pthread_t	* _thread;
	// attributes are usually NULL
	pthread_attr_t	* _attr;
	// errorcodes...
	void *_return;
	int _errcode;
	int _id;

};

template <class T> int Thread<T>::num_threads = 0;

template <class T>
Thread<T>::Thread (fct_ptr f, void * arg) 
:
	_function(f),
	_arg((T *)arg),
	_attr(NULL),
	_id(Thread<T>::num_threads++)
{

	printf ("making a new thread\n");

	if (_errcode = pthread_create(_thread, _attr, _function, _arg)) {
		fprintf(stderr,"%s: %s\n","Thread<T>::Thread()",strerror(_errcode));
		exit(1);
	}
}

template <class T>
void * Thread<T>::Join () {
	
	printf ("killing %d\n", _id);
	if (_errcode = pthread_join (*_thread, &_return)) {
		fprintf(stderr,"%s: %s\n", "Thread<T>::Join()",strerror(_errcode));
		exit(1);
	}
	Thread<T>::num_threads--;

return (_return);
}


#endif
