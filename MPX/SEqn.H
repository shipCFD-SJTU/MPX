{
    kappaS = turbulence->nut()/Sct;
    kappaS.correctBoundaryConditions();

    volScalarField Kst("Kst", turbulence->nu()/Sc + kappaS);

    fvScalarMatrix SEqn
    (
        fvm::ddt(S)
      + fvm::div(phi, S)
      - fvm::laplacian(Kst, S)   
     ==
       fvOptions(S)                         
    );

    SEqn.relax();

    fvOptions.constrain(SEqn);

    SEqn.solve();

    radiation->correct();

    fvOptions.correct(S);

    
}
