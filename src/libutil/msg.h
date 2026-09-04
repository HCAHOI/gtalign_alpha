/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __msg_h__
#define __msg_h__

#include <stdio.h>
#include <string>
#include <chrono>

#if 1//def __DEBUG__
#define MYMSGnonl(CSTR,lev) message(CSTR,false,lev)
#define MYMSG(CSTR,lev) message(CSTR,true,lev)
#define MYMSGBEGl(lvl) if(lvl<=VERBOSE){
#define MYMSGENDl }
#define MYMSGBEG if(1){
#define MYMSGEND }
#else
#define MYMSGnonl(ARG,lev)
#define MYMSG(ARG,lev)
#define MYMSGBEGl(lvl) if(0){
#define MYMSGENDl }
#define MYMSGBEG if(0){
#define MYMSGEND }
#endif

// typedefs
typedef int (*TPrintFunction)( void*, const char* format, ... );

// global variables used
extern const char*  PROGDIR;
extern const char*  PROGNAME;
extern const char*  PROGVERSION;
extern const char*  PROGREFERENCES[];

extern int*         __PARGC__;
extern char***      __PARGV__;

extern int          VERBOSE;
extern bool         WARNINGSRECORDED;

extern const
std::chrono::high_resolution_clock::time_point gtSTART;

// set global variables
/**
 * @brief 在通用工具中设置 `SetVerboseMode` 对应的数据。
 * @param param1 供该函数读取或更新的 `param1` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void SetVerboseMode( int );
/**
 * @brief 在通用工具中设置 `SetProgramName` 对应的数据。
 * @param name 控制当前步骤范围或规模的 `name`。
 * @param version 供该函数读取或更新的 `version` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void SetProgramName( const char* name, const char* version = NULL );
/**
 * @brief 在通用工具中设置 `SetArguments` 对应的数据。
 * @param pargc 供该函数读取或更新的 `pargc` 参数。
 * @param pargv 供该函数读取或更新的 `pargv` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void SetArguments( int* pargc, char*** pargv );

// messaging...
/**
 * @brief 在通用工具中处理 `error` 对应的数据。
 * @param param1 供该函数读取或更新的 `param1` 参数。
 * @param param2 供该函数读取或更新的 `param2` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void error( const char*, bool = true );
/**
 * @brief 在通用工具中处理 `warning` 对应的数据。
 * @param param1 供该函数读取或更新的 `param1` 参数。
 * @param param2 供该函数读取或更新的 `param2` 参数。
 * @param minlevel 供该函数读取或更新的 `minlevel` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void warning( const char*, bool = true, int minlevel = 1 );
/**
 * @brief 在通用工具中检查 `checkforwarnings` 对应的数据。
 * @par 参数
 * 无。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void checkforwarnings();
/**
 * @brief 在通用工具中处理 `message` 对应的数据。
 * @param param1 供该函数读取或更新的 `param1` 参数。
 * @param param2 供该函数读取或更新的 `param2` 参数。
 * @param minlevel 供该函数读取或更新的 `minlevel` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void message( const char*, bool = true, int minlevel = 1  );
/**
 * @brief 在通用工具中处理 `progname_and_version` 对应的数据。
 * @param FILE 供该函数读取或更新的 `FILE` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void progname_and_version( FILE* );
/**
 * @brief 在通用工具中格式化输出 `print_cmdline` 对应的数据。
 * @param outstr 接收当前步骤输出的 `outstr`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int print_cmdline(char* outstr);
/**
 * @brief 在通用工具中格式化输出 `print_dtime` 对应的数据。
 * @param minlevel 供该函数读取或更新的 `minlevel` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void print_dtime(int minlevel);
/**
 * @brief 在通用工具中读取 `getdtime` 对应的数据。
 * @param tmstr 供该函数读取或更新的 `tmstr` 参数。
 * @param sztmstr 供该函数读取或更新的 `sztmstr` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
size_t getdtime(char* tmstr, size_t sztmstr);

// string
/**
 * @brief 在通用工具中处理 `usage` 对应的数据。
 * @param progname 供该函数读取或更新的 `progname` 参数。
 * @param instructions 供该函数读取或更新的 `instructions` 参数。
 * @param version 供该函数读取或更新的 `version` 参数。
 * @param date 供该函数读取或更新的 `date` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
std::string usage( const char* progname, const char* instructions, const char* version, const char* date );

// path processing
/**
 * @brief 在通用工具中处理 `my_basename` 对应的数据。
 * @param param1 供该函数读取或更新的 `param1` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
const char* my_basename( const char* );
/**
 * @brief 在通用工具中处理 `my_dirname` 对应的数据。
 * @param param1 供该函数读取或更新的 `param1` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
const char* my_dirname( const char* );

// printing-to-stream interface
/**
 * @brief 在通用工具中处理 `file_print` 对应的数据。
 * @param param1 供该函数读取或更新的 `param1` 参数。
 * @param format 供该函数读取或更新的 `format` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int file_print( void*, const char* format, ... );
/**
 * @brief 在通用工具中处理 `string_print` 对应的数据。
 * @param param1 供该函数读取或更新的 `param1` 参数。
 * @param format 供该函数读取或更新的 `format` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int string_print( void*, const char* format, ... );

#endif//__msg_h__
