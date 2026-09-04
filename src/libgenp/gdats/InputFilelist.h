/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __InputFilelist_h__
#define __InputFilelist_h__

#include <stdio.h>
#include <string.h>

#include <string>
#include <vector>
#include <algorithm>
#include <functional>

#include <fstream>

#include "libutil/mybase.h"
#include "libutil/CLOptions.h"

// -------------------------------------------------------------------------
// _________________________________________________________________________
// Class InputFilelist
//
// Controller for conducting reading of data
//
class InputFilelist
{
public:
    enum {
        DEFNUMCHAINS = 16,//default number of chains per structure file
        CUDBREADER_TMPBUFFSIZE = 4096,
        MAXCOMPLETESTRUCTSIZE = 128 * ONEM
    };
    //extensions
    enum TKnownExts{
            FDRKETar,//.tar file
            FDRKEGz,//.gz file
            FDRKEPDBGz,//.pdb.gz file
            FDRKEENTGz,//.ent.gz file
            FDRKECIFGz,//.cif.gz file
            FDRKEPDB,//.pdb file
            FDRKEENT,//.ent file
            FDRKECIF,//.cif file
            FDRKEN
    };
    //data file type...
    enum TDataFile{
            FDRFlTar,//file is a tar archive
            FDRFlZip,//file is a compressed structure file
            FDRFlPDBZip,//file is a compressed structure PDB file
            FDRFlPDBxmmCIFZip,//file is a compressed structure PDBx/mmCIF file
            FDRFlStructure,//file is a structure file either in PDB or PDBx/mmCIF format
            FDRFlStructurePDB,//file is a structure file in PDB format
            FDRFlStructurePDBxmmCIF,//file is a structure file in PDBx/mmCIF format
            FDRFlN
    };

    /**
     * @brief 构造 `InputFilelist`，初始化其负责的结构数据读取与布局状态。
     * @param dnamelst 供该函数读取或更新的 `dnamelst` 参数。
     * @param sfxlst 供该函数读取或更新的 `sfxlst` 参数。
     * @param clustering 供该函数读取或更新的 `clustering` 参数。
     * @param construct 供该函数读取或更新的 `construct` 参数。
     * @return 无返回值；完成对象构造与初始状态设置。
     */
    InputFilelist(const std::vector<std::string>& dnamelst,
                const std::vector<std::string>& sfxlst,
                const bool clustering = false,
                const bool construct = true);
    /**
     * @brief 销毁 `InputFilelist`，释放其管理的结构数据读取与布局资源。
     * @par 参数
     * 无。
     * @return 无返回值；对象持有的资源在返回前完成释放。
     */
    ~InputFilelist();

    /**
     * @brief 在结构数据读取与布局中处理 `ConstructFileList` 对应的数据。
     * @par 参数
     * 无。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void ConstructFileList();
    /**
     * @brief 在结构数据读取与布局中读取 `GetFileTypeFromFilename` 对应的数据。
     * @param string 供该函数读取或更新的 `string` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    int GetFileTypeFromFilename(const std::string&);

    /**
     * @brief 在结构数据读取与布局中读取 `GetStrFilelist` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const std::vector<std::string>& GetStrFilelist() const {return strfilelist_;}
    /**
     * @brief 在结构数据读取与布局中读取 `GetPntFilelist` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const std::vector<std::string>& GetPntFilelist() const {return pntfilelist_;}
    /**
     * @brief 在结构数据读取与布局中读取 `GetStrFilePositionlist` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const std::vector<size_t>& GetStrFilePositionlist() const {return strfilepositionlist_;}
    /**
     * @brief 在结构数据读取与布局中读取 `GetStrFilesizelist` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const std::vector<size_t>& GetStrFilesizelist() const {return strfilesizelist_;}
    /**
     * @brief 在结构数据读取与布局中读取 `GetStrParenttypelist` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const std::vector<int>& GetStrParenttypelist() const {return strparenttypelist_;}
    /**
     * @brief 在结构数据读取与布局中读取 `GetStrFiletypelist` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const std::vector<int>& GetStrFiletypelist() const {return strfiletypelist_;}
    /**
     * @brief 在结构数据读取与布局中读取 `GetFilendxlist` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const std::vector<size_t>& GetFilendxlist() const {return filendxlist_;}

    /**
     * @brief 在结构数据读取与布局中读取 `GetGlobalIds` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    const std::vector<std::vector<int>>& GetGlobalIds() const {return globalids_;}
    /**
     * @brief 在结构数据读取与布局中读取 `GetGlobalIds` 对应的数据。
     * @par 参数
     * 无。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    std::vector<std::vector<int>>& GetGlobalIds() {return globalids_;}

protected:
    /**
     * @brief 在结构数据读取与布局中处理 `ProcessEntry` 对应的数据。
     * @param string 供该函数读取或更新的 `string` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    template<int LEVEL> void ProcessEntry(const std::string&);
    /**
     * @brief 在结构数据读取与布局中处理 `AddFile` 对应的数据。
     * @param string 供该函数读取或更新的 `string` 参数。
     * @param param2 供该函数读取或更新的 `param2` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void AddFile(const std::string&, bool);
    /**
     * @brief 在结构数据读取与布局中处理 `AddFilesFromTAR` 对应的数据。
     * @param string 供该函数读取或更新的 `string` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void AddFilesFromTAR(const std::string&);

    /**
     * @brief 在结构数据读取与布局中处理 `SuffixFound` 对应的数据。
     * @param entryname 供该函数读取或更新的 `entryname` 参数。
     * @return 返回该步骤计算、查询或状态判断的结果。
     */
    bool SuffixFound(const std::string& entryname) const
    {
        //verify whether the file suffix is among the specified ones
        auto sit = std::find_if(sfxlst_.begin(), sfxlst_.end(),
            [&entryname](const std::string& sfx){
                return//file has this extension/suffix if true:
                    sfx.size() <= entryname.size() && 
                    entryname.compare(entryname.size()-sfx.size(), 
                        sfx.size(), sfx) == 0;
            });
        if(sit == sfxlst_.end())
            //file does not have a valid extension
            return false;
        return true;
    }

    /**
     * @brief 在结构数据读取与布局中处理 `AddEntry` 对应的数据。
     * @param entryname 供该函数读取或更新的 `entryname` 参数。
     * @param parentname 供该函数读取或更新的 `parentname` 参数。
     * @param position 供该函数读取或更新的 `position` 参数。
     * @param filesize 控制当前步骤范围或规模的 `filesize`。
     * @param parenttype 供该函数读取或更新的 `parenttype` 参数。
     * @param filetype 供该函数读取或更新的 `filetype` 参数。
     * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
     */
    void AddEntry(
        const std::string& entryname, const std::string& parentname,
        size_t position, size_t filesize, int parenttype, int filetype)
    {
        strfilelist_.push_back(entryname);
        pntfilelist_.push_back(parentname);
        strfilepositionlist_.push_back(position);
        strfilesizelist_.push_back(filesize);
        strparenttypelist_.push_back(parenttype);
        strfiletypelist_.push_back(filetype);
        filendxlist_.push_back(filendxlist_.size());
        if(clustering_) {
            globalids_.push_back(std::vector<int>());
            if(!globalids_.empty()) globalids_.back().reserve(DEFNUMCHAINS);
        }
    }

private:
    const std::vector<std::string>& dnamelst_;//input
    const std::vector<std::string>& sfxlst_;//input
    const bool clustering_;//file list for clustering

    std::vector<std::string> strfilelist_;//structure file list
    std::vector<std::string> pntfilelist_;//parent file list
    std::vector<size_t> strfilepositionlist_;//list of structure file positions within an archive
    std::vector<size_t> strfilesizelist_;//list of structure file sizes
    std::vector<int> strparenttypelist_;//parent file type list of structure files (e.g, tar)
    std::vector<int> strfiletypelist_;//list of structure file types
    std::vector<size_t> filendxlist_;//list of the indices of files sorted by filesize

    std::vector<std::vector<int>> globalids_;//global ids for structures across all files

    char tmpbuff_[CUDBREADER_TMPBUFFSIZE];//buffer for temporary data

public:
    static const std::string knownexts_[FDRKEN];//known extension
};


// /////////////////////////////////////////////////////////////////////////
// INLINES
//
// -------------------------------------------------------------------------
// GetFileTypeFromFilename: get file type from the filename
//
/**
 * @brief 在结构数据读取与布局中读取 `InputFilelist::GetFileTypeFromFilename` 对应的数据。
 * @param filename 输入或输出文件路径。
 * @return 返回该步骤计算、查询或状态判断的结果。
 */
inline
int InputFilelist::GetFileTypeFromFilename(const std::string& filename)
{
    //verify file extension for match
    std::function<bool(int)> lfExtMatched = [&filename](int extcode) {
        if(knownexts_[extcode].size() < filename.size() && 
           filename.compare(filename.size()-knownexts_[extcode].size(),
                knownexts_[extcode].size(), knownexts_[extcode]) == 0)
            return true;
        return false;
    };

    if(lfExtMatched(FDRKEGz)) {
        if(CLOptions::GetI_INFMT() == CLOptions::iifPDB) return FDRFlPDBZip;
        if(CLOptions::GetI_INFMT() == CLOptions::iifmmCIF) return FDRFlPDBxmmCIFZip;
        if(lfExtMatched(FDRKEPDBGz) || lfExtMatched(FDRKEENTGz)) return FDRFlPDBZip;
        if(lfExtMatched(FDRKECIFGz)) return FDRFlPDBxmmCIFZip;
        return FDRFlZip;
    }
    if(lfExtMatched(FDRKETar)) return FDRFlTar;
    if(CLOptions::GetI_INFMT() == CLOptions::iifPDB) return FDRFlStructurePDB;
    if(CLOptions::GetI_INFMT() == CLOptions::iifmmCIF) return FDRFlStructurePDBxmmCIF;
    if(lfExtMatched(FDRKEPDB)) return FDRFlStructurePDB;
    if(lfExtMatched(FDRKEENT)) return FDRFlStructurePDB;
    if(lfExtMatched(FDRKECIF)) return FDRFlStructurePDBxmmCIF;
    //if no extension is recognized, assume it is a structure file whose 
    //type is to be determined
    return FDRFlStructure;
}

// -------------------------------------------------------------------------
// AddFile: add file to a list of structure files
//
/**
 * @brief 在结构数据读取与布局中处理 `InputFilelist::AddFile` 对应的数据。
 * @param entryname 供该函数读取或更新的 `entryname` 参数。
 * @param sfxcheck 供该函数读取或更新的 `sfxcheck` 参数。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
inline
void InputFilelist::AddFile(const std::string& entryname, bool sfxcheck)
{
    if(sfxcheck && sfxlst_.size())
        if(!SuffixFound(entryname))
            return;

    size_t filesize = 0;

    if(file_size(entryname.c_str(), &filesize) != 0) {
        warning(("Unable to get file information; Ignored: " + entryname).c_str());
        return;
    }

    if(filesize < 10)
        return;

    int filetype = GetFileTypeFromFilename(entryname);
    int parenttype = filetype;
    size_t position = 0;

    if(filetype == FDRFlTar) {
        AddFilesFromTAR(entryname);
        return;
    }

    AddEntry(entryname, entryname/*parentname*/,
        position, filesize, parenttype, filetype);
}

#endif//__InputFilelist_h__
