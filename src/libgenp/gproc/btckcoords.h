/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __btckcoords_h__
#define __btckcoords_h__

// invalid DP cell coordinates (positions; as a result of CombineCoords);
// or Stop marker:
#define INVALIDDPCELLCOORDS 0xffffffff
#define DPCELLSTOP 0xffffffff

// -------------------------------------------------------------------------
// CombineCoords: combine (x,y) coordinates (structure pairwise positions) 
// into one integer value;
// NOTE: the arguments x and y are supposed to contain 16-bit values!
#ifdef GPUINUSE
__host__ __device__ __forceinline__
#else
/**
 * @brief 在通用结构处理中处理 `CombineCoords` 对应的数据。
 * @param x 供该函数读取或更新的 `x` 参数。
 * @param y 供该函数读取或更新的 `y` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
#endif
unsigned int CombineCoords(unsigned int x, unsigned int y)
{
    return (x << 16) | y;//(y & 0xffff);
}

// GetCoordX and GetCoordY extract x and y coordinates from the 
// combined value
#ifdef GPUINUSE
__host__ __device__ __forceinline__
#else
/**
 * @brief 在通用结构处理中读取 `GetCoordX` 对应的数据。
 * @param xy 供该函数读取或更新的 `xy` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
#endif
unsigned int GetCoordX(unsigned int xy)
{
    return (xy >> 16) & 0xffff;
}

#ifdef GPUINUSE
__host__ __device__ __forceinline__
#else
/**
 * @brief 在通用结构处理中读取 `GetCoordY` 对应的数据。
 * @param xy 供该函数读取或更新的 `xy` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
#endif
unsigned int GetCoordY(unsigned int xy)
{
    return xy & 0xffff;
}

#endif//__btckcoords_h__
