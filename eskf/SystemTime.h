/*
 * @Author: your name
 * @Date: 2021-05-10 09:12:21
 * @LastEditTime: 2021-07-20 09:08:06
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: \mpu\code\include\common\SystemTime.h
 */
#ifndef SYSTEM_TIME_H
#define SYSTEM_TIME_H



#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#include <unistd.h>
#endif
#include <ctime>
#include <cstdlib>
#include <chrono>


typedef unsigned char      uint8;
typedef unsigned short     uint16;
typedef unsigned int       uint32;
typedef unsigned long long uint64;
typedef float              float32;
typedef double             float64;

typedef char               int8;
typedef short              int16;
typedef int                int32;
typedef long long          int64;


inline float64 GetSystemTimeStamp()
{

#ifdef _WIN32
    auto t1 = std::chrono::system_clock::now().time_since_epoch();
    auto int_ms = std::chrono::duration_cast<std::chrono::microseconds>(t1) / 1000000.0;;
    return int_ms.count();
#elif _linux
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec / 1000000.0;
#endif
}

class TicToc
{
  public:
    TicToc()
    {
        tic();
    }

    void tic()
    {
        start = GetSystemTimeStamp();
    }

    double toc()
    {
        end = GetSystemTimeStamp();
        return end - start;
    }

  private:
    double start, end;
};

#endif


