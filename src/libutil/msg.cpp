/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include <time.h>

#include <string>
#include <thread>
#include <chrono>

#include "platform.h"
#include "mytype.h"
#include "myexception.h"
#include "mylimits.h"
#include "msg.h"

const
std::chrono::high_resolution_clock::time_point gtSTART = 
std::chrono::high_resolution_clock::now();

const char* PROGDIR = NULL;
const char* PROGNAME = NULL;
const char* PROGVERSION = NULL;
const char* PROGREFERENCES[] = {
    "Margelevicius, M. "
    "GTalign: spatial index-driven protein structure alignment, superposition, and search. "
    "Nat Commun 15, 7305 (2024).",
    // "https://doi.org/10.1038/s41467-024-51669-z",
    NULL
};

int*        __PARGC__ = NULL;
char***     __PARGV__ = NULL;

int         VERBOSE = 0;
bool        WARNINGSRECORDED = false;
bool        ERRORSRECORDED = false;

// -------------------------------------------------------------------------
// set global variables
//
/**
 * @brief 在通用工具中设置 `SetProgramName` 对应的数据。
 * @param name 控制当前步骤范围或规模的 `name`。
 * @param version 供该函数读取或更新的 `version` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void SetProgramName( const char* name, const char* version )
{
    PROGDIR = my_dirname( name );
    PROGNAME = my_basename( name );
    if( version )
        PROGVERSION = version;
}

/**
 * @brief 在通用工具中设置 `SetArguments` 对应的数据。
 * @param pargc 供该函数读取或更新的 `pargc` 参数。
 * @param pargv 供该函数读取或更新的 `pargv` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void SetArguments( int* pargc, char*** pargv )
{
    __PARGC__ = pargc;
    __PARGV__ = pargv;
}

/**
 * @brief 在通用工具中设置 `SetVerboseMode` 对应的数据。
 * @param value 需要读取、写入或转换的值。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void SetVerboseMode( int value )
{
    VERBOSE = value;
}

// -------------------------------------------------------------------------
// messaging routines
//
/**
 * @brief 在通用工具中处理 `error` 对应的数据。
 * @param str 供该函数读取或更新的 `str` 参数。
 * @param putnl 供该函数读取或更新的 `putnl` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void error( const char* str, bool putnl )
{
    ERRORSRECORDED = true;
    if(str && PROGNAME && putnl)
        fprintf( stderr, "[%s] ERROR: %s%s%s", PROGNAME, str, NL, NL);
    else if(str && PROGNAME)
        fprintf( stderr, "[%s] ERROR: %s%s", PROGNAME, str, NL);
    else if(str && putnl)
        fprintf( stderr, "ERROR: %s%s%s", str, NL, NL);
    else if(str)
        fprintf( stderr, "ERROR: %s%s", str, NL);
}

/**
 * @brief 在通用工具中处理 `warning` 对应的数据。
 * @param str 供该函数读取或更新的 `str` 参数。
 * @param putnl 供该函数读取或更新的 `putnl` 参数。
 * @param minlevel 供该函数读取或更新的 `minlevel` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void warning( const char* str, bool putnl, int/* minlevel*/)
{
//     if( VERBOSE < minlevel )
//         return;
    WARNINGSRECORDED = true;
    if(str && PROGNAME && putnl)
        fprintf( stderr, "[%s] WARNING: %s%s%s", PROGNAME, str, NL, NL);
    else if(str && PROGNAME)
        fprintf( stderr, "[%s] WARNING: %s%s", PROGNAME, str, NL);
    else if(str && putnl)
        fprintf( stderr, "WARNING: %s%s%s", str, NL, NL);
    else if(str)
        fprintf( stderr, "WARNING: %s%s", str, NL);
}

/**
 * @brief 在通用工具中检查 `checkforwarnings` 对应的数据。
 * @par 参数
 * 无。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void checkforwarnings()
{
    if( !WARNINGSRECORDED && !ERRORSRECORDED)
        return;
    if(PROGNAME)
        fprintf( stderr, "%s[%s] There are %s.%s", NL, PROGNAME,
                ERRORSRECORDED? "ERRORS": "WARNINGS", NL);
    else 
        fprintf( stderr, "%sThere are %s.%s", NL,
                ERRORSRECORDED? "ERRORS": "WARNINGS", NL);
}

/**
 * @brief 在通用工具中处理 `message` 对应的数据。
 * @param str 供该函数读取或更新的 `str` 参数。
 * @param putnl 供该函数读取或更新的 `putnl` 参数。
 * @param minlevel 供该函数读取或更新的 `minlevel` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void message( const char* str, bool putnl, int minlevel )
{
    if( VERBOSE < minlevel )
        return;

    std::chrono::high_resolution_clock::time_point tnow = 
    std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elpsd = tnow - gtSTART;

    if(str && PROGNAME && putnl)
        fprintf( stderr, "[%s] {%.6fs} (%zu) %s%s%s", PROGNAME, elpsd.count(),
                std::hash<std::thread::id>{}(std::this_thread::get_id()), 
                str, NL, NL);
    else if(str && PROGNAME)
        fprintf( stderr, "[%s] {%.6fs} (%zu) %s%s", PROGNAME, elpsd.count(),
                std::hash<std::thread::id>{}(std::this_thread::get_id()), 
                str, NL);
    else if(str && putnl)
        fprintf( stderr, "{%.6fs} (%zu) %s%s%s", elpsd.count(),
                std::hash<std::thread::id>{}(std::this_thread::get_id()), 
                str, NL, NL);
    else if(str)
        fprintf( stderr, "{%.6fs} (%zu) %s%s", elpsd.count(),
                std::hash<std::thread::id>{}(std::this_thread::get_id()), 
                str, NL);
}

/**
 * @brief 在通用工具中处理 `progname_and_version` 对应的数据。
 * @param fp 供该函数读取或更新的 `fp` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void progname_and_version( FILE* fp )
{
    if( PROGNAME )
        fprintf( fp, "%s", PROGNAME );
    if( PROGVERSION )
        fprintf( fp, " %s", PROGVERSION );
    if( PROGNAME || PROGVERSION )
        fprintf( fp, "%s%s", NL, NL );
}

// -------------------------------------------------------------------------
// print_cmdline: print command line to string;
// NOTE: space for outstr assumed to be preallocated!
//
/**
 * @brief 在通用工具中格式化输出 `print_cmdline` 对应的数据。
 * @param outstr 接收当前步骤输出的 `outstr`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int print_cmdline(char* outstr)
{
    if(!outstr) return 0;
    int written = 0;
    if(__PARGV__ && __PARGC__  && *__PARGV__) {
        for(int n = 0; n < *__PARGC__; n++) {
            int locwrt = sprintf(outstr,"%s ",(*__PARGV__)[n]);
            written += locwrt;
            outstr += locwrt;
        }
    }
    return written;
}

// -------------------------------------------------------------------------
// print_dtime: print date and local time
//
/**
 * @brief 在通用工具中格式化输出 `print_dtime` 对应的数据。
 * @param minlevel 供该函数读取或更新的 `minlevel` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void print_dtime(int minlevel)
{
    if( VERBOSE < minlevel )
        return;
    char tmstr[BUF_MAX];
    time_t t = time(NULL);
    struct tm* ctm = localtime(&t);
    if( ctm ) {
        if( strftime(tmstr, sizeof(tmstr), "%c", ctm) != 0) {
            if(PROGNAME)
                fprintf(stderr,"[%s] %s%s%s", PROGNAME, tmstr,NL,NL);
            else
                fprintf(stderr,"%s%s%s", tmstr,NL,NL);
        }
    }
}

// -------------------------------------------------------------------------
// getdtime: thread-unsafe function for formatting and writing local 
// date and time info to string;
// NOTE: space for tmstr assumed to be preallocated!
//
/**
 * @brief 在通用工具中读取 `getdtime` 对应的数据。
 * @param tmstr 供该函数读取或更新的 `tmstr` 参数。
 * @param sztmstr 供该函数读取或更新的 `sztmstr` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
size_t getdtime(char* tmstr, size_t sztmstr)
{
    if(tmstr == NULL || sztmstr < 1) return 0;
    time_t t = time(NULL);
    struct tm* ctm = localtime(&t);
    if(ctm) return strftime(tmstr, sztmstr, "%c", ctm);
    return 0;
}

// -------------------------------------------------------------------------
// path-processing functions
//
/**
 * @brief 在通用工具中处理 `my_basename` 对应的数据。
 * @param name 控制当前步骤范围或规模的 `name`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
const char* my_basename( const char* name )
{
    const char* bp = NULL;
    if( name )
        for(bp  = name + strlen( name );
            bp != name && bp[-1] != DIRSEP;
            bp-- );
    return  bp;
}

// my_dirname: return directory name given pathname; if directory is not
//  found in the given path, return the current directory (`.')
//
/**
 * @brief 在通用工具中处理 `my_dirname` 对应的数据。
 * @param name 控制当前步骤范围或规模的 `name`。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
const char* my_dirname( const char* name )
{
    static const size_t size = KBYTE;
    static char dir[size];
    const char* bp = name;

    dir[0] = '.'; dir[1] = 0;

    if( !bp )
        return dir;

    for( bp += strlen( name ); bp != name && bp[-1] != DIRSEP; bp-- );

    if( bp != name && *bp == 0 ) {
        //the last character is dir separator
        for( ; bp != name && bp[-1] == DIRSEP; bp-- );

        if( bp != name )
            //find separator before the last one
            for( ; bp != name && bp[-1] != DIRSEP; bp-- );
    }

    if( bp == name )
        //no directory separator has been found
        return dir;

    for( ; bp != name && bp[-1] == DIRSEP; bp-- );

    if( bp == name ) {
        bp++;

        if( *bp == DIRSEP && bp[1] != DIRSEP )
            //there are exactly two separators at the beginning; save them
            bp++;
    }

    if( bp - name < (ssize_t)size ) {
        //if enough memory allocated
        memcpy( dir, name, bp - name );
        dir[ bp - name ] = 0;
    }

    return dir;
}

// -------------------------------------------------------------------------
// usage: returns <instructions> translated by appropriately inserting
//     program name, version, and date information
//
/**
 * @brief 在通用工具中处理 `usage` 对应的数据。
 * @param path 供该函数读取或更新的 `path` 参数。
 * @param instructions 供该函数读取或更新的 `instructions` 参数。
 * @param version 供该函数读取或更新的 `version` 参数。
 * @param date 供该函数读取或更新的 `date` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
std::string usage( const char* path, const char* instructions, const char* version, const char* date )
{
    size_t pos;
    std::string instr = instructions;
    std::string prog = path;
    std::string full;
    bool first = true;

    if(!path || !instructions )
        return instr;

    if((pos = prog.rfind( DIRSEP )) != std::string::npos){
        prog = prog.substr( pos+1 );
    }

    full = prog + std::string(" ");

    if(version && strlen( version ))
        full += version;

    if(date && strlen( date ))
        full += std::string(" (" ) + date + std::string( ")");

    while((pos = instr.find("<>")) != std::string::npos){
        instr.erase( pos, 2 );
        instr.insert( pos, first? full: prog );
        if( first ) 
            first = false;
    }

    while((pos = instr.find("-*")) != std::string::npos){
        instr.erase( pos, 2 );
        for( size_t len = full.length(); len; len-- )
            instr.insert( pos++, "-");
    }

    if((pos = instr.find("[]")) != std::string::npos){
        instr.erase(pos, 2);
#ifdef GPUINUSE
        instr.insert(pos, "(compiled with GPU support)");
#else
        instr.insert(pos, "(CPU/multiprocessing compilation without support for GPUs)");
#endif
    }

//     if((pos = instr.find("{}")) != std::string::npos){
//         instr.erase(pos, 2);
// #ifdef GPUINUSE
// #else
// #endif
//     }

    return instr;
}

// -------------------------------------------------------------------------
// file_print: redirect formatted output string to file; file
//     pointer (vpn) must be preintialized and must be valid
//
/**
 * @brief 在通用工具中处理 `file_print` 对应的数据。
 * @param vpn 供该函数读取或更新的 `vpn` 参数。
 * @param format 供该函数读取或更新的 `format` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int file_print( void* vpn, const char* format, ... )
{
    if( !vpn )
        return -1;

    FILE*   fp = ( FILE* )vpn;
    va_list ap;
    int     ret;

    va_start( ap, format );

    ret = vfprintf( fp, format, ap );

    va_end( ap );

    if( ret < 0 )
        throw myruntime_error("file_print: Formatted print to file failed.", __EXCPOINT__ );

    return ret;
}

// -------------------------------------------------------------------------
// string_print: same as file_print except that redirection is to chr. 
//     string, which is assumed to have enough space to contain information;
//     string pointer (vpn) must be preallocated and must be valid
//
/**
 * @brief 在通用工具中处理 `string_print` 对应的数据。
 * @param vpn 供该函数读取或更新的 `vpn` 参数。
 * @param format 供该函数读取或更新的 `format` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int string_print( void* vpn, const char* format, ... )
{
    if( !vpn )
        return -1;

    char*   sp = ( char* )vpn;
    va_list ap;
    int     ret;

    va_start( ap, format );

    ret = vsprintf( sp + strlen( sp ), format, ap );

    va_end( ap );

    if( ret < 0 )
        throw myruntime_error("string_print: Formatted print to string failed.", __EXCPOINT__ );

    return ret;
}

