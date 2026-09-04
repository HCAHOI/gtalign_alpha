/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

// #include <ctype.h>
// #include <stdlib.h>
// #include <string.h>

#include <string>
#include <algorithm>

#include "mylimits.h"
#include "myexception.h"
#include "mygetopt.h"

// file globals:

const size_t ggonMaxopts = BUF_MAX;
const char* ggosTer = "--";

// -------------------------------------------------------------------------
// CLASS MyGetopt
//
// constructor
//
/**
 * @brief 构造 `MyGetopt`，初始化其负责的通用工具状态。
 * @param opts 供该函数读取或更新的 `opts` 参数。
 * @param argv 命令行参数字符串数组。
 * @param argc 命令行参数个数。
 * @return 无返回值；完成对象构造与初始状态设置。
 */
MyGetopt::MyGetopt( const myoption* opts, const char *argv[], int argc )
:   ncnt_( 0 ),
    stopped_( false ),
    argv_( NULL ),
    argc_( 0 )
{
    srtoptions_.reserve(50);
    Init( opts, argv, argc );
}

// destructor
//
/**
 * @brief 销毁 `MyGetopt`，释放其管理的通用工具资源。
 * @par 参数
 * 无。
 * @return 无返回值；对象持有的资源在返回前完成释放。
 */
MyGetopt::~MyGetopt()
{
    argv_ = NULL;
}

// -------------------------------------------------------------------------
// Init: initialize object
//
/**
 * @brief 在通用工具中初始化 `MyGetopt::Init` 对应的数据。
 * @param options 供该函数读取或更新的 `options` 参数。
 * @param argv 命令行参数字符串数组。
 * @param argc 命令行参数个数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void MyGetopt::Init( const myoption* options, const char *argv[], int argc )
{
    size_t  n;
    if( options == NULL || argv == NULL )
        throw myruntime_error("MyGetopt::Init: Memory access error.", __EXCPOINT__ );
    argv_ = argv;
    argc_ = argc;
    for( n = 0; ; n++) {
        if( ggonMaxopts <= n )
            throw myruntime_error("MyGetopt::Init: Too long list of options.", __EXCPOINT__ );
        if( options[n].sLongname_ == NULL )
            break;
        //make a sorted list of option names
        srtoptions_.insert(
            std::lower_bound(srtoptions_.begin(), srtoptions_.end(), options[n]),
            options[n]);
    }
}

// -------------------------------------------------------------------------
// GetNextOption: get the next option from processing the command line;
//      return value:
//          option's value if the option has been found;
//          '?' if option has not been found;
//          -1 if the processing of the command line has to be 
//      stopped; the processing is stopped on the end of the command line or 
//      options terminator `--'.
//
/**
 * @brief 在通用工具中读取 `MyGetopt::GetNextOption` 对应的数据。
 * @param argument 供该函数读取或更新的 `argument` 参数。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
int MyGetopt::GetNextOption( std::string* argument )
{
    if( stopped_ )
        return my_gorv_term;
    if( ncnt_ < 0 || ncnt_ > argc_ )
        return my_gorv_term;
    const char* carg = argv_[ncnt_];
    ncnt_++;
    if( carg == NULL || *carg == 0 )
        return my_gorv_term;
    if( strcmp( carg, ggosTer ) == 0 ) {
        stopped_ = true;
        return my_gorv_term;
    }
    if( *carg != '-' ) {
        if( argument )
            *argument = carg;
        return my_gorv_value;
    }
    carg++;
    if( *carg == 0 )
        return my_gorv_illformed;
    if( *carg == '-' ) {
        ++carg;
        if( *carg == 0 || *carg == '-' )
            return my_gorv_illformed;
    }

    myoption privopt = {NULL,my_n_targflags,0};//option read from the command line
    std::string arg = carg;
    std::string key = arg;
    size_t p = arg.find ('=');

    if( p != std::string::npos ) {
        //command line option is provided with '='
        key = arg.substr( 0, p );
        if( argument )
            *argument = arg.substr( p+1 );
    }

    privopt.sLongname_ = key.c_str();

    auto it_option = std::lower_bound(srtoptions_.begin(), srtoptions_.end(), privopt);

    if( it_option == srtoptions_.end() || !(*it_option == privopt))
        //option not found
        return my_gorv_notfound;

    while(( it_option->eFlag_ == my_required_argument || 
            it_option->eFlag_ == my_optional_argument ) && 
            p == std::string::npos ) {
        //command line option found; next, read argument (value) if required
        if( ncnt_ > argc_ ) {
            if( it_option->eFlag_ == my_optional_argument ) {
                if( argument )
                    argument->erase();
                break;
            }
            return my_gorv_noarg;
        }
        carg = argv_[ncnt_];
        ncnt_++;
        if( carg == NULL || *carg == 0 || *carg == '-') {
            if( it_option->eFlag_ == my_optional_argument ) {
                if( argument )
                    argument->erase();
                ncnt_--;
                break;
            }
            return my_gorv_noarg;
        }
        if( argument )
            *argument = carg;
        break;
    }

    return it_option->nRet_;
}
