/***************************************************************************
 *   Copyright (C) 2021-2023 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <omp.h>

#include <math.h>
#include <string>

#include "libutil/cnsts.h"
#include "libutil/macros.h"
#include "libutil/templates.h"
#include "libutil/CLOptions.h"
#include "libgenp/gproc/gproc.h"
#include "libgenp/gdats/PM2DVectorFields.h"
#include "libgenp/gdats/PMBatchStrData.h"
#include "libmymp/mpproc/mpprocconf.h"
#include "libmymp/mplayout/MpGlobalMemory.h"
#include "MpReform.h"

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// SelectAndReformatKernel: reformat (with a store and load) a reference
// database chunk to include candidates proceeding to stages of more
// detailed superposition search and refinement;
//
/**
 * @brief 在CPU 对齐流水线中并行计算 `MpReform::SelectAndReformatKernel` 对应的数据。
 * @param ndbCstrs2 控制当前步骤范围或规模的 `ndbCstrs2`。
 * @param maxndbCposs 描述参考结构的 `maxndbCposs`。
 * @param filterdata 供该函数读取或更新的 `filterdata` 参数。
 * @param querypmbeg 查询结构打包字段的起始指针数组。
 * @param bdbCpmbeg 参考结构打包字段的起始指针数组。
 * @param tfmmem 保存最终刚体变换矩阵的缓冲区。
 * @param wrkmemaux 保存分数、收敛标记等辅助状态的工作缓冲区。
 * @param tmpdpdiagbuffers 供当前步骤读取或更新的 `tmpdpdiagbuffers` 缓冲区。
 * @return 无返回值；结果写入传入缓冲区、输出参数或对象状态。
 */
void MpReform::SelectAndReformatKernel(
    const int ndbCstrs2,
    const int maxndbCposs,
    const uint* const __RESTRICT__ filterdata,
    const char* const * const __RESTRICT__ querypmbeg,
    const char* const * const __RESTRICT__ bdbCpmbeg,
    float* const __RESTRICT__ tfmmem,
    float* const __RESTRICT__ wrkmemaux,
    float* const __RESTRICT__ tmpdpdiagbuffers)
{
    enum {
        lXNDX = fdNewReferenceIndex,//index for reference structure new indices
        lXADD = fdNewReferenceAddress//index for reference structure new addresses
    };

    MYMSG("MpReform::SelectAndReformatKernel", 4);
    // static const std::string preamb = "MpReform::SelectAndReformatKernel: ";
    static const int nthreads = CLOptions::GetCPU_THREADS();
    constexpr int memalignment = CuMemoryBase::GetMinMemAlignment();

    const int nblocks_x_1 = ndbCstrs_;
    const int nblocks_x_2 = ndbCstrs2;
    const int nblocks_y = nqystrs_;

    size_t chunksize_helper1 = ((size_t)nblocks_y * (size_t)nblocks_x_1 + (size_t)nthreads - 1) / nthreads;
    size_t chunksize_helper2 = ((size_t)nblocks_y * (size_t)nblocks_x_2 + (size_t)nthreads - 1) / nthreads;
    const int chunksize1 = (int)mymin(chunksize_helper1, (size_t)MPFL_STORECANDIDATEDATA_CHSIZE);
    const int chunksize2 = (int)mymin(chunksize_helper2, (size_t)MPFL_LOADCANDIDATEDATA_CHSIZE);

    #pragma omp parallel num_threads(nthreads) default(shared)
    {
        //store data
        #pragma omp for collapse(2) schedule(dynamic, chunksize1)
        for(int qi = 0; qi < nblocks_y; qi++)
            for(int ri = 0; ri < nblocks_x_1; ri++)
            {
                int offset = (qi/*qryndx*/) * maxndbCposs;
                int newndx = filterdata[lXNDX * ndbCstrs_ + ri/*dbstrndx*/];
                if(newndx == 0) continue;
                //adjust index appropriately (NOTE: newaddr is already valid):
                newndx--;

                #pragma omp simd aligned(wrkmemaux,tmpdpdiagbuffers:memalignment)
                for(int fi = 0; fi < nTAuxWorkingMemoryVars; fi++) {
                    //NOTE: uncoalesced read, coalesced write:
                    int mloc = ((qi * maxnsteps_ + 0) * nTAuxWorkingMemoryVars + fi) * ndbCstrs_;
                    //NOTE: nTAuxWorkingMemoryVars * max(ndbCstrs) < max(ndbCposs) by definition (memory layout)
                    int tloc = offset + nTAuxWorkingMemoryVars * newndx;
                    tmpdpdiagbuffers[tloc + fi] = wrkmemaux[mloc + ri/*dbstrndx*/];
                }

                offset = (nqystrs_ + qi/*qryndx*/) * maxndbCposs;
                #pragma omp simd aligned(tfmmem,tmpdpdiagbuffers:memalignment)
                for(int fi = 0; fi < nTTranformMatrix; fi++) {
                    //NOTE: coalesced read, coalesced write:
                    int mloc = (qi * ndbCstrs_ + ri/*dbstrndx*/) * nTTranformMatrix;
                    //NOTE: nTTranformMatrixv * max(ndbCstrs) < max(ndbCposs) by definition (memory layout)
                    int tloc = offset + nTTranformMatrix * newndx;
                    tmpdpdiagbuffers[tloc + fi] = tfmmem[mloc + fi];
                }
            }//omp for
        //implicit barrier here

        //load data
        #pragma omp for collapse(2) schedule(dynamic, chunksize2)
        for(int qi = 0; qi < nblocks_y; qi++)
            for(int ri = 0; ri < nblocks_x_2; ri++)
            {
                const int qrydst = PMBatchStrData::GetAddressAt(querypmbeg, qi);
                const int typexq = PMBatchStrData::GetFieldAt<INTYPE,pmv2D_Ins_Ch_Ord>(querypmbeg, qrydst);
                const int dbstrdst = PMBatchStrData::GetAddressAt(bdbCpmbeg, ri);
                const int typexr = PMBatchStrData::GetFieldAt<INTYPE,pmv2D_Ins_Ch_Ord>(bdbCpmbeg, dbstrdst);
                const int typeqry = GetMoleculeType(typexq);
                const int typerfn = GetMoleculeType(typexr);

                int offset = (qi/*qryndx*/) * maxndbCposs;
                #pragma omp simd aligned(wrkmemaux,tmpdpdiagbuffers:memalignment)
                for(int fi = 0; fi < nTAuxWorkingMemoryVars; fi++) {
                    //NOTE: coalesced read, uncoalesced write:
                    int mloc = ((qi * maxnsteps_ + 0) * nTAuxWorkingMemoryVars + fi) * ndbCstrs2;
                    int tloc = offset + nTAuxWorkingMemoryVars * ri/*dbstrndx*/;
                    wrkmemaux[mloc + ri/*dbstrndx*/] = tmpdpdiagbuffers[tloc + fi];
                }

                {   const int mloc0 = ((qi * maxnsteps_ + 0) * nTAuxWorkingMemoryVars + tawmvConverged) * ndbCstrs2;
                    wrkmemaux[mloc0 + ri] = (typeqry == typerfn)? 0.0f: (float)CONVERGED_LOWTMSC_bitval;
                }

                offset = (nqystrs_ + qi/*qryndx*/) * maxndbCposs;
                #pragma omp simd aligned(tfmmem,tmpdpdiagbuffers:memalignment)
                for(int fi = 0; fi < nTTranformMatrix; fi++) {
                    //NOTE: coalesced read, coalesced write:
                    int mloc = (qi * ndbCstrs2 + ri/*dbstrndx*/) * nTTranformMatrix;
                    int tloc = offset + nTTranformMatrix * ri/*dbstrndx*/;
                    tfmmem[mloc + fi] = tmpdpdiagbuffers[tloc + fi];
                }
            }
    }//omp parallel
}
