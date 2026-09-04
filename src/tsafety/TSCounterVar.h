/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __TSCounterVar_h__
#define __TSCounterVar_h__

#include <mutex>
#include <condition_variable>

class TSCounterVar
{
public:
    /**
     * @brief 构造 `TSCounterVar`，初始化其负责的线程安全支持状态。
     * @par 参数
     * 无。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    TSCounterVar()
    : counter_(0)
    {}
    /**
     * @brief 在线程安全支持中重置 `reset` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void reset()
    {
        {   std::lock_guard<std::mutex> lck(mtx_);
            counter_ = 0;
        }
        cv_.notify_one();
    }
    /**
     * @brief 在线程安全支持中设置 `set` 对应的数据。
     * @param value 需要读取、写入或转换的值。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void set(int value)
    {
        {   std::lock_guard<std::mutex> lck(mtx_);
            counter_ = value;
        }
        cv_.notify_one();
    }
    /**
     * @brief 在线程安全支持中处理 `isset` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int isset() const
    {
        return (get() > 0);
    }
    /**
     * @brief 在线程安全支持中读取 `get` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int get() const
    {
        std::lock_guard<std::mutex> lck(mtx_);
        int value = counter_;
        return value;
    }
    /**
     * @brief 在线程安全支持中处理 `inc` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void inc()
    {
        std::lock_guard<std::mutex> lck(mtx_);
        counter_++;
    }
    /**
     * @brief 在线程安全支持中处理 `dec` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void dec()
    {
        {   std::lock_guard<std::mutex> lck(mtx_);
            counter_--;
        }
        cv_.notify_one();
    }
    /**
     * @brief 在线程安全支持中等待 `wait0` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void wait0()
    {
        std::unique_lock<std::mutex> lck(mtx_);
        cv_.wait(lck, [this]{return counter_ < 1;});
    }

    /**
     * @brief 在线程安全支持中读取 `get_mutex` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    std::mutex& get_mutex() {return mtx_;}

    /**
     * @brief 在线程安全支持中重置 `reset_under_lock` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void reset_under_lock() {counter_ = 0;}

    /**
     * @brief 在线程安全支持中处理 `inc_under_lock` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void inc_under_lock() {counter_++;}

    /**
     * @brief 在线程安全支持中读取 `get_under_lock` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int get_under_lock() const {return counter_;}

    /**
     * @brief 在线程安全支持中设置 `set_under_lock` 对应的数据。
     * @param value 需要读取、写入或转换的值。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void set_under_lock(int value) {counter_ = value;}

private:
    mutable std::mutex mtx_;
    std::condition_variable cv_;
    int counter_;
};

#endif//__TSCounterVar_h__
