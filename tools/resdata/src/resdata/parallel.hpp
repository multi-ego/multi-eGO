#ifndef _RESDATA_SEMAPHORE_HPP
#define _RESDATA_SEMAPHORE_HPP

#include <condition_variable>
#include <mutex>

namespace resdata::parallel
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
    std::unique_lock<std::mutex> lock(mut);
    cv.wait(lock, [this](){ return counter > 0; });
    --counter;
    lock.unlock();
  }

  void release()
  {
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
    
} // namespace resdata::parallel

#endif // _RESDATA_SEMAPHORE_HPP