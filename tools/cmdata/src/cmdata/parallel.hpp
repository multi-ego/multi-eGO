#ifndef _CMDATA_SEMAPHORE_HPP
#define _CMDATA_SEMAPHORE_HPP

#include <condition_variable>
#include <mutex>

// #define DB
#ifdef DB
#define LOG(x) std::cout << x << std::endl;
#else
#define LOG(x)
#endif

namespace cmdata::parallel
{

class Semaphore
{
private:
  std::mutex mut;
  std::condition_variable cv;
  std::uint64_t counter;

public:
  Semaphore( std::uint64_t counter = 0 ) : mut(std::mutex()), cv(std::condition_variable()), counter(counter) {}
  void acquire()
  {
    LOG("acquire out of " << counter);
    std::unique_lock<std::mutex> lock(mut);
    cv.wait(lock, [this](){ return counter > 0; });
    --counter;
    lock.unlock();
  }

  void release()
  {
    LOG("release out of " << counter);
    std::unique_lock<std::mutex> lock(mut);
    ++counter;
    lock.unlock();
    cv.notify_one();
  }

  void set_counter( std::uint64_t c )
  {
    counter = c;
  }
};
    
} // namespace cmdata::parallel

#endif // _CMDATA_SEMAPHORE_HPP