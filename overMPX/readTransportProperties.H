singlePhaseTransportModel laminarTransport(U, phi);

// Thermal expansion coefficient [1/K]
/*
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport
);
*/
// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);

// Reference sailty [-]
dimensionedScalar SRef("SRef", dimless, laminarTransport);

// Laminar Schmidt number
dimensionedScalar Sc("Sc", dimless, laminarTransport);

// Turbulent Schmidt number
dimensionedScalar Sct("Sct", dimless, laminarTransport);

// Coefficient of Thermal Expansion
dimensionedScalar betaT("betaT", dimless/dimTemperature, laminarTransport);

// Coefficient of Salinity Contraction
dimensionedScalar betaS("betaS", dimless, laminarTransport);

//Reference density
//dimensionedScalar rhoRef("rhoRef", dimless, laminarTransport);