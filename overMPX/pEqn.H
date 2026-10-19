// Option 1: interpolate rAU, do not block out rAU on blocked cells
volScalarField rAU("rAU", 1.0/UEqn.A()); //计算动量方程逆对角系数
mesh.interpolate(rAU);  //重叠网格插值，数据传递
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU)); //插值生成面心场
volVectorField H("H", UEqn.H());
volVectorField HbyA("HbyA", U); //无压力梯度的预估速度
HbyA = constrainHbyA(rAU*H, U, p_rgh);

surfaceScalarField phig(-rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf()); //(独有)浮力产生的通量项

if (massFluxInterpolation) //重叠网格交界面高精度通量插值
{
    #include "interpolatedFaces.H"
}

if (runTime.outputTime()) //调试输出
{
    H.write();
    rAU.write();
    HbyA.write();
}

if (pimple.nCorrPISO() <= 1) //释放内存
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA //（区别）预估通量=对流扩散通量+时间项通量+浮力通量
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + phig
    );

if (ddtCorr) //时间导数通量校正 faceMaskOld强制被遮挡网格不参与通量计算
{
    surfaceScalarField faceMaskOld
    (
        localMin<scalar>(mesh).interpolate(cellMask.oldTime())
    );
    phiHbyA += rAUf*faceMaskOld*fvc::ddtCorr(U, Uf);
}

MRF.makeRelative(phiHbyA); //旋转坐标系修正

constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);


if (adjustFringe) //重叠网格通量修正
{
    fvc::makeRelative(phiHbyA, U);
    oversetAdjustPhi(phiHbyA, U);
    fvc::makeAbsolute(phiHbyA, U);
}

while (pimple.correctNonOrthogonal()) //求解压力泊松方程
{
    fvScalarMatrix p_rghEqn
    (
        fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
    );

    p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

    p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
       phi = phiHbyA - p_rghEqn.flux();
       
       p_rgh.relax(); 

       U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
       
       U.correctBoundaryConditions();
       
       fvOptions.correct(U);
    }
}

#include "continuityErrs.H" //连续性误差检查

p = p_rgh + rhok*gh;

if (p_rgh.needReference()) //压力参考值修正
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rhok*gh;
    }

{
    Uf = fvc::interpolate(U);
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    Uf += n*(phi/mesh.magSf() - (n & Uf));
} //面心速度强修正

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

surfaceScalarField faceMask
(
    localMin<scalar>(mesh).interpolate(cellMask)
);
phi *= faceMask;
