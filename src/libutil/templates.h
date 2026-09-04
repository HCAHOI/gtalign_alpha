/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __templates_h__
#define __templates_h__

// absolute value of the difference:
template <typename T>
/**
 * @brief 在通用工具中读取 `GetAbsDiff` 对应的数据。
 * @param v1 供该函数读取或更新的 `v1` 参数。
 * @param v2 供该函数读取或更新的 `v2` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
T GetAbsDiff(T v1, T v2) {return (v1 < v2)? v2 - v1: v1 - v2;}

//simple swap template function
/**
 * @brief 在通用工具中处理 `myswap` 对应的数据。
 * @param a 供该函数读取或更新的 `a` 参数。
 * @param b 供该函数读取或更新的 `b` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
template <typename T> void myswap(T& a, T& b) { T c(a); a=b; b=c; }

//simple template functions for squares
/**
 * @brief 在通用工具中处理 `mysqrd` 对应的数据。
 * @param a 供该函数读取或更新的 `a` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
template <typename T> T mysqrd(T& a){ return a*a; }
/**
 * @brief 在通用工具中处理 `mysqrdv` 对应的数据。
 * @param a 供该函数读取或更新的 `a` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
template <typename T> T mysqrdv(T a){ return a*a; }

//min/max
/**
 * @brief 在通用工具中处理 `mymax` 对应的数据。
 * @param a 供该函数读取或更新的 `a` 参数。
 * @param b 供该函数读取或更新的 `b` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
template <typename T> T mymax(T a, T b){ return a<b?b:a; }
/**
 * @brief 在通用工具中处理 `mymin` 对应的数据。
 * @param a 供该函数读取或更新的 `a` 参数。
 * @param b 供该函数读取或更新的 `b` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
template <typename T> T mymin(T a, T b){ return a<b?a:b; }

/**
 * @brief 在通用工具中处理 `mycemax` 对应的数据。
 * @param a 供该函数读取或更新的 `a` 参数。
 * @param b 供该函数读取或更新的 `b` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
template <typename T> constexpr T mycemax(T a, T b){ return a<b?b:a; }
/**
 * @brief 在通用工具中处理 `mycemin` 对应的数据。
 * @param a 供该函数读取或更新的 `a` 参数。
 * @param b 供该函数读取或更新的 `b` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
template <typename T> constexpr T mycemin(T a, T b){ return a<b?a:b; }

//conditional assignment
template <typename T, typename T2>
/**
 * @brief 在通用工具中处理 `mymaxassgn` 对应的数据。
 * @param a 供该函数读取或更新的 `a` 参数。
 * @param b 供该函数读取或更新的 `b` 参数。
 * @param c 供该函数读取或更新的 `c` 参数。
 * @param d 供该函数读取或更新的 `d` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void mymaxassgn(T& a, T b, T2& c, T2 d){ if(a<b){a=b; c=d;} }

#endif//__templates_h__
