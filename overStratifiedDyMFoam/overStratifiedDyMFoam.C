/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    overPimpleDyMFoam

Group
    grpIncompressibleSolvers grpMovingMeshSolvers

Description
    Transient solver for incompressible flow of Newtonian fluids
    on a moving mesh using the PIMPLE (merged PISO-SIMPLE) algorithm.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"   //到此处与PimpleFoam完全一致（少了一个CorrectPhi文件）
#include "radiationModel.H"  //引入辐射传热模型

#include "cellCellStencilObject.H" // 构建重叠网格插值模版（决定供体、受体）
#include "zeroGradientFvPatchFields.H" // 零梯度边界条件
#include "localMin.H"  //并行计算中的最小值操作
#include "interpolationCellPoint.H" //插值方法相关
#include "transform.H" //坐标变化
#include "fvMeshSubset.H" //创建子网格
#include "oversetAdjustPhi.H" //(重要）修正质量通量，使其在重叠网格上守恒

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " on a moving mesh."
    );

    #include "postProcess.H" // 后处理相关支持

    #include "setRootCaseLists.H" //解析case路径
    #include "createTime.H" //创建时间对象（读取controdict
    #include "createDynamicFvMesh.H"//创建动态网格对象（重叠网格）
    #include "initContinuityErrs.H"  //初始化连续性误差

    pimpleControl pimple(mesh); // 创建Pimple算法控制对象

    #include "createFields.H"  //*（重要！）读取所有物理场，该文件在本求解器文件夹中，读取0文件夹
    #include "createUf.H" //创建面速度Uf
    //#include "createMRF.H" //创建旋转坐标系
    //#include "createFvOptions.H" //创建fvOption对象
    #include "createControls.H" //读取控制参数*
    #include "CourantNo.H" //计算库朗数*
    #include "setInitialDeltaT.H"//设置初始时间步

    turbulence->validate(); //湍流模型合法性检查

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())  //最外层的时间循环
    {
        #include "readControls.H" //读取控制参数*
        #include "CourantNo.H" //计算库朗数*

        #include "setDeltaT.H" //自动调整时间步

        ++runTime; //时间步推进

        Info<< "Time = " << runTime.timeName() << nl << endl;

        bool changed = mesh.update(); //（核心）判断网格是否发生变化

        if (changed) //如果网格发生变化
        {
            #include "setCellMask.H" //设置洞单元、计算单元
            #include "setInterpolatedCells.H" //设置插值单元

            surfaceScalarField faceMaskOld
            (
                localMin<scalar>(mesh).interpolate(cellMask.oldTime())
            );  //将旧时间步的cellMask插值到面（关键步骤）

            // Zero Uf on old faceMask (H-I)
            Uf *= faceMaskOld; //旧hole face上的速度被清除
            // Update Uf and phi on new C-I faces
            Uf += (1-faceMaskOld)*fvc::interpolate(U); //新face上重新计算速度
            phi = mesh.Sf() & Uf; //重新计算通量

            // Zero phi on current H-I
            surfaceScalarField faceMask
            (
                localMin<scalar>(mesh).interpolate(cellMask)
            ); //将当前cellMask插值到face
            phi *= faceMask; //清除当前holeface上的通量
        }


        if (mesh.changing() && correctPhi)
        {
            // Calculate absolute flux from the mapped surface velocity
            #include "correctPhi.H" //修正通量*
        }

        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi, U); //转化为相对通量（动网格）

        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H" //计算网格库朗数
        }

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop()) //pimple外循环
        {
            #include "UEqn.H"//(核心）求解动量方程，文件在本文件夹，求解预测速度而不是最终速度
            #include "TEqn.H"
            #include "SEqn.H"
            rhok = 1.0 - betaT*(T - TRef) + betaS * (S - SRef);
            
            // --- Pressure corrector loop
            while (pimple.correct()) //压力修正循环（内循环）
            {
                #include "pEqn.H"//求解压力泊松方程，强制满足连续性方程
            }

            if (pimple.turbCorr()) //湍流修正，判断是否更新湍流模型
            {
                laminarTransport.correct();
                turbulence->correct();  // 更新湍流相关物理量
            }
        }

        runTime.write(); // 按照writeInterval写入结果

        runTime.printExecutionTime(Info); //打印实际计算时间
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
