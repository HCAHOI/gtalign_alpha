/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __warpscan_h__
#define __warpscan_h__

#include "cutemplates.h"

// mywarpincprefixmin: inclusive prefix min perfomed in a warp
template <typename T> 
__device__ __forceinline__ 
/**
 * @brief 在CUDA 通用原语中处理 `mywarpincprefixmin` 对应的数据。
 * @param reg 供该函数读取或更新的 `reg` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
T mywarpincprefixmin( T reg )
{
    T tmp;
    int lid = threadIdx.x & 0x1f;//lane id
    tmp = __shfl_up_sync(0xffffffff, reg, 1); if (1 <= lid) reg = myhdmin(reg, tmp);
    tmp = __shfl_up_sync(0xffffffff, reg, 2); if (2 <= lid) reg = myhdmin(reg, tmp);
    tmp = __shfl_up_sync(0xffffffff, reg, 4); if (4 <= lid) reg = myhdmin(reg, tmp);
    tmp = __shfl_up_sync(0xffffffff, reg, 8); if (8 <= lid) reg = myhdmin(reg, tmp);
    tmp = __shfl_up_sync(0xffffffff, reg, 16); if (16 <= lid) reg = myhdmin(reg, tmp);
    return reg;
}

// mywarpincprefixsum: inclusive prefix sum perfomed in a warp
template <typename T> 
__device__ __forceinline__ 
/**
 * @brief 在CUDA 通用原语中处理 `mywarpincprefixsum` 对应的数据。
 * @param reg 供该函数读取或更新的 `reg` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
T mywarpincprefixsum( T reg )
{
    T tmp;
    int lid = threadIdx.x & 0x1f;//lane id
    tmp = __shfl_up_sync(0xffffffff, reg, 1); if (1 <= lid) reg += tmp;
    tmp = __shfl_up_sync(0xffffffff, reg, 2); if (2 <= lid) reg += tmp;
    tmp = __shfl_up_sync(0xffffffff, reg, 4); if (4 <= lid) reg += tmp;
    tmp = __shfl_up_sync(0xffffffff, reg, 8); if (8 <= lid) reg += tmp;
    tmp = __shfl_up_sync(0xffffffff, reg, 16); if (16 <= lid) reg += tmp;
    return reg;
}

// mywarprevincprefixsum: inclusive prefix sum accumulated in the reversed order
template <typename T> 
__device__ __forceinline__ 
/**
 * @brief 在CUDA 通用原语中处理 `mywarprevincprefixsum` 对应的数据。
 * @param reg 供该函数读取或更新的 `reg` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
T mywarprevincprefixsum( T reg )
{
    T tmp;
    int lid = threadIdx.x & 0x1f;//lane id
    tmp = __shfl_down_sync(0xffffffff, reg, 1); if (lid <= 1) reg += tmp;
    tmp = __shfl_down_sync(0xffffffff, reg, 2); if (lid <= 2) reg += tmp;
    tmp = __shfl_down_sync(0xffffffff, reg, 4); if (lid <= 4) reg += tmp;
    tmp = __shfl_down_sync(0xffffffff, reg, 8); if (lid <= 8) reg += tmp;
    tmp = __shfl_down_sync(0xffffffff, reg, 16); if (lid <= 16) reg += tmp;
    return reg;
}

// mywarpreducemax: warp reduce for the maximum value across a warp
template <typename T> 
__device__ __forceinline__ 
/**
 * @brief 在CUDA 通用原语中处理 `mywarpreducemax` 对应的数据。
 * @param value 需要读取、写入或转换的值。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
T mywarpreducemax( T value )
{
    //warp reduce
    value = myhdmax( value, __shfl_down_sync(0xffffffff, value, 16));
    value = myhdmax( value, __shfl_down_sync(0xffffffff, value, 8));
    value = myhdmax( value, __shfl_down_sync(0xffffffff, value, 4));
    value = myhdmax( value, __shfl_down_sync(0xffffffff, value, 2));
    value = myhdmax( value, __shfl_down_sync(0xffffffff, value, 1));
    return value;
}

// mywarpreducemin: warp reduce for the minimum value across a warp
template <typename T> 
__device__ __forceinline__ 
/**
 * @brief 在CUDA 通用原语中处理 `mywarpreducemin` 对应的数据。
 * @param value 需要读取、写入或转换的值。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
T mywarpreducemin( T value )
{
    //warp reduce
    value = myhdmin( value, __shfl_down_sync(0xffffffff, value, 16));
    value = myhdmin( value, __shfl_down_sync(0xffffffff, value, 8));
    value = myhdmin( value, __shfl_down_sync(0xffffffff, value, 4));
    value = myhdmin( value, __shfl_down_sync(0xffffffff, value, 2));
    value = myhdmin( value, __shfl_down_sync(0xffffffff, value, 1));
    return value;
}

// mywarpreducesum: warp reduce for the sum across a warp
template <typename T> 
__device__ __forceinline__ 
/**
 * @brief 在CUDA 通用原语中处理 `mywarpreducesum` 对应的数据。
 * @param value 需要读取、写入或转换的值。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
T mywarpreducesum( T value )
{
    //warp reduce
    value += __shfl_down_sync(0xffffffff, value, 16);
    value += __shfl_down_sync(0xffffffff, value, 8);
    value += __shfl_down_sync(0xffffffff, value, 4);
    value += __shfl_down_sync(0xffffffff, value, 2);
    value += __shfl_down_sync(0xffffffff, value, 1);
    return value;
}

#endif//__warpscan_h__
