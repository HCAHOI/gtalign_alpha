/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __cumath_cuh__
#define __cumath_cuh__

// myfastdiv3: divide dividend by 3 fast; credit: Voo (stackoverflow)
__device__ __forceinline__
/**
 * @brief 在CUDA 对齐流水线中处理 `myfastdiv3` 对应的数据。
 * @param dividend 供该函数读取或更新的 `dividend` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
uint myfastdiv3(uint dividend)
{
    enum{magic = 0x55555556};//(2^32 + 2) / 3
    return __mulhi(magic, dividend);//(dividend * (2^32 + 2) / 3) / 2^32
}

// myfastmod3: divide dividend modulo 3 fast;
__device__ __forceinline__
/**
 * @brief 在CUDA 对齐流水线中处理 `myfastmod3` 对应的数据。
 * @param dividend 供该函数读取或更新的 `dividend` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
uint myfastmod3(uint dividend)
{
    return dividend - 3 * myfastdiv3(dividend);
}

#endif