/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __localconfig_h__
#define __localconfig_h__

#include "libutil/preproc.h"
#include "libutil/platform.h"

#define CONCATSTRSEP(arg1)      CONCATSTRA(arg1, DIRSEPSTR)

#if !defined(GTALIGNSTATEDIR)
#   error   "Unable to compile: Undefined configuration directory (GTALIGNSTATEDIR)."
#endif

#define VARDIR var
#define PARAM_DIRNAME  GTALIGNSTATEDIR

#define PARAM_FILENAME gtalign.par
#define PARAM_FULLNAME CONCATSTRSEP( PARAM_DIRNAME ) TOSTR( PARAM_FILENAME )

static const char*  var_param_DIR = TOSTR( VARDIR );
static const char*  var_param_DIRNAME = TOSTR( PARAM_DIRNAME );

static const char*  var_param_FILENAME = TOSTR( PARAM_FILENAME );
static const char*  var_param_FULLPATHNAME = PARAM_FULLNAME;

/**
 * @brief 返回参数文件目录的短目录名。
 * @par 参数
 * 无。
 * @return 返回静态目录名字符串，调用方不拥有该内存。
 */
inline const char* GetParamDirectory()      {   return var_param_DIR;  }
/**
 * @brief 返回编译时确定的 GTAlign 参数文件完整目录。
 * @par 参数
 * 无。
 * @return 返回静态完整目录字符串，调用方不拥有该内存。
 */
inline const char* GetFullParamDirname()    {   return var_param_DIRNAME;  }

/**
 * @brief 返回 GTAlign 参数文件的文件名。
 * @par 参数
 * 无。
 * @return 返回静态文件名字符串，调用方不拥有该内存。
 */
inline const char* GetParamFilename()       {   return var_param_FILENAME;  }
/**
 * @brief 返回 GTAlign 参数文件的完整路径。
 * @par 参数
 * 无。
 * @return 返回静态完整路径字符串，调用方不拥有该内存。
 */
inline const char* GetFullParamFilename()   {   return var_param_FULLPATHNAME;  }

#endif//__localconfig_h__
