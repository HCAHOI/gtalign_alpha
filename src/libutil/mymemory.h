/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __mymemory_h__
#define __mymemory_h__

#include <stdio.h>
#include <cstdlib>
#include "platform.h"

#ifdef OS_MS_WINDOWS
//memory allocation
/**
 * @brief 在 Windows 上分配满足指定字节对齐要求的内存。
 * @param alignment 返回地址必须满足的字节对齐量。
 * @param size 要分配的字节数。
 * @return 成功时返回对齐内存地址，失败时返回空指针。
 */
inline void* my_aligned_alloc(size_t alignment, size_t size)
{
    //NOTE: arguments order
    return _aligned_malloc(size, alignment);
}
//memory deallocation
/**
 * @brief 释放由 Windows 对齐分配器返回的内存。
 * @param memptr `my_aligned_alloc()` 返回的地址；允许为空指针。
 * @return 无返回值；`memptr` 指向的内存被释放。
 */
inline void my_aligned_free(void* memptr)
{
    return _aligned_free(memptr);
}
#else
//memory allocation
/**
 * @brief 在通用工具中处理 `my_aligned_alloc` 对应的数据。
 * @param alignment 供该函数读取或更新的 `alignment` 参数。
 * @param size 控制当前步骤范围或规模的 `size`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline void* my_aligned_alloc(size_t alignment, size_t size)
{
    return aligned_alloc(alignment, size);
}
//memory deallocation
/**
 * @brief 在通用工具中处理 `my_aligned_free` 对应的数据。
 * @param memptr 供当前步骤读取或更新的 `memptr` 缓冲区。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline void my_aligned_free(void* memptr)
{
    return free(memptr);
}
#endif

#endif//__mymemory_h__
