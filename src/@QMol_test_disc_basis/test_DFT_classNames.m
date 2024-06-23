function test_DFT_classNames(obj)
%test_DFT_classNames test class names for DFT simulations

    % Component names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obj.showSection('DFT component names');

    obj.showResult('DFT_classDensity',          strcmpi(QMol_disc_basis.DFT_classDensity,          'QMol_DFT_density'));
    obj.showResult('DFT_classOrbital',          strcmpi(QMol_disc_basis.DFT_classOrbital,          'QMol_DFT_orbital_basis'));
    obj.showResult('DFT_classPotential',        strcmpi(QMol_disc_basis.DFT_classPotential,        'QMol_DFT_Vks_basis'));
    obj.showResult('DFT_classPotentialGradient',strcmpi(QMol_disc_basis.DFT_classPotentialGradient,'N/A'));

end