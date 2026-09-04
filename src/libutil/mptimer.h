/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __mptimer_cuh__
#define __mptimer_cuh__

#include <ratio>
#include <chrono>

////////////////////////////////////////////////////////////////////////////
// CLASS MyMpTimer for measuring performance
//
// template<class Period = std::micro>
template<class Period = std::ratio<1>>
class MyMpTimer {
public:
    /**
     * @brief 构造 `MyMpTimer`，初始化其负责的通用工具状态。
     * @par 参数
     * 无。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    MyMpTimer() {}

    //start measuring performance
    /**
     * @brief 在通用工具中处理 `Start` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Start() {
        start_ = std::chrono::high_resolution_clock::now();
    }

    //stop the timer; measure the elapsed time
    /**
     * @brief 在通用工具中处理 `Stop` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void Stop() {
        stop_ = std::chrono::high_resolution_clock::now();
        elapsed_ = stop_- start_;
    }

    /**
     * @brief 在通用工具中读取 `GetElapsedTime` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    double GetElapsedTime() const {return elapsed_.count();}

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_, stop_;
    std::chrono::duration<double, Period> elapsed_;//time elapsed in us when using std::micro
};

#endif//__mptimer_cuh__
