{
    kappaS = turbulence->nut()/Sct;
    kappaS.correctBoundaryConditions();

    volScalarField Kst("Kst", turbulence->nu()/Sc + kappaS);

    fvScalarMatrix SEqn
    (
        fvm::ddt(S)
      + fvm::div(phi, S)
      - fvm::laplacian(Kst, S)      //ysj,待添加盐分子扩散系数
     ==
       fvOptions(S)                             //ysj,待添加盐度辐射源项
    );

    SEqn.relax();

    fvOptions.constrain(SEqn);

    SEqn.solve();

    radiation->correct();

    fvOptions.correct(S);

    
}
