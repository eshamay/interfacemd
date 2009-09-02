#include <boost/bind.hpp>
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

#define nthreads 5
using namespace boost;

void task1(int& a) {
    printf ("task1 %d\n",a);
	a = a+3;
}

void task2() {
    printf ("task2\n");
}

int main (int argc, char ** argv) {
    using namespace boost;
	typedef boost::shared_ptr<boost::thread> Thread_Ptr_t;
	typedef std::vector < Thread_Ptr_t > Thread_Pool_t;
	Thread_Pool_t m_thread_pool;

	for (int i = 0; i < nthreads; i++) {
		Thread_Ptr_t t (new thread(bind(task1,i)));
		m_thread_pool.push_back(t);
	}
	Thread_Pool_t::iterator ti;
	for (ti = m_thread_pool.begin(); ti != m_thread_pool.end(); ti++) {
		(*ti)->join();
	}


    return 0;
}

	//int thread_count = 2; // 2-Core cpu

	// Create a pool of Apple threads

	// Create threads

/*
	for ( int index = 0; index < thread_count; ++index )
	{
		task t;
		Thread_Ptr_t new_thread_ptr (
			new boost::thread (boost::bind (&t )( index ) ) );
		m_thread_pool.push_back ( new_thread_ptr );
	}

	// Wait for all threads to finish
	for ( Thread_Pool_t::iterator pos = m_thread_pool.begin();
		pos != m_thread_pool.end();
		++pos )
	{
		(*pos)->join();
	}
*/


